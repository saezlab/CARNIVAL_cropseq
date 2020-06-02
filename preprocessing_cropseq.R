#!/usr/bin/env Rscript  
# coding: utf-8
# Copyright (C) 2019 Olga Ivanova
#
# Contact: olga.ivanova@bioquant.uni-heidelberg.de
#
# ====================
# GNU-GLPv3:
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation.
#
# This program is distributed in the hope that it will be useful, but
# WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
# General Public License for more details.
#
# A full copy of the GNU General Public License can be found on
# http://www.gnu.org/licenses/.
#
# ====================
# DESCRIPTION:
# Preprocessing and running VIPER/DoRothEA for CropSeq single cell data.
# ====================

if ( !require("here") ) { 
  tryCatch( {
    install.packages("here")  
  }, error = function(e) {
    print("Unable to install package 'here'. Continuing without it.")
  })
}

########################################################################################
### ------ SETTING UP DEFAULT PARAMETERS (if running as a standalone script) ------- ###
######################################################################################## 

if ( (!exists("source_folder") || source_folder == "") && "here" %in% (.packages()) ) {
  source_folder = here()
} else if ( !exists("source_folder") ) {
  source_folder = ""  
}

if ( !exists("settings_run") ) {
  
  if ( !require("yaml") ) { 
    install.packages("yaml")
  }
  
  if ("here" %in% (.packages()) ) {
    source_folder = here()   
  } else {
    source_folder = ""
  }
  
  settings_file = file.path( source_folder, "settings.yml" )
  
  if ( file.exists(settings_file) ) { 
    settings_run = yaml.load_file( settings_file )$local
  } else {
    stop("Can't continue, please specify correct path to settings file")
  }
  
  if ( dir.exists(source_folder) ) {
    source( file.path(source_folder, "setting_up_pipeline.R") )
  } else {
    stop("processing_cropseq.R: can't continue, please specify correct source folder")
  }
}

########################################################################################
### ------------ INSTALLING/LOADING NECESSARY PACKAGES ----------------------------- ###
########################################################################################

if ( dir.exists(source_folder) ) {
  source( file.path(source_folder, "packages_utils.R") )  
} else {
  stop("run_carnival_cropseq.R: can't continue, please specify correct source folder")
}
  
cran_list_packages = c("dplyr","readr", "stringr", "purrr", "tibble", "tidyr",
                       "Seurat", "logging", "BiocManager", "devtools")
bioc_list_packages = c("viper", "biomaRt", "UniProt.ws", "rhdf5", "hdf5r")

CheckAndLoadLibraries( cran_list_packages, bioc_list_packages )

basicConfig(level = "DEBUG")
addHandler(writeToFile, logger = "preprocessing_run", file = logfile)
loginfo("Preprocessing script started", logger = "preprocessing_run.module")

source( file.path(source_folder, "utils_CARNIVAL.R") )
source( file.path(source_folder, "utils_cropseq.R") )

########################################################################################
### ------------ READING THE DATA -------------------------------------------------- ###
########################################################################################
loginfo( paste( "Reading the data files:", annotation_file, 
                ";", raw_file ), 
         logger = "preprocessing_run.module" )

tcr_cell_annotation = read.delim( annotation_file )
trc_crispr_data     = Read10X_h5( raw_file )

loginfo( "Running standard seurat pipeline", logger = "preprocessing_run.module" ) 

# Following standard seurat pipeline from 
# https://satijalab.org/seurat/v3.1/pbmc3k_tutorial.html
tcr_cropseq_seurat = CreateSeuratObject(counts = trc_crispr_data, 
                                        project = "tcr_crispr_data", 
                                        min.cells = 3, 
                                        min.features = 200)

tcr_cropseq_seurat[["percent.mt"]] = PercentageFeatureSet( tcr_cropseq_seurat, pattern = "^MT-" )

# Filtering cells with low feature counts, too high fecture counts (), and high percentage of mitochondrial genes 
tcr_cropseq_seurat_subset = subset(tcr_cropseq_seurat, subset = nFeature_RNA > 200 & 
                                     nFeature_RNA < 2500 & 
                                     percent.mt < 5)
tcr_cropseq_seurat_subset = NormalizeData( tcr_cropseq_seurat_subset )

# Keeping only annotated cells
tcr_cropseq_seurat_annotated = subset(tcr_cropseq_seurat_subset, 
                                      cells = tcr_cell_annotation$Cell.barcode)

#########################################################################################
### ------------ COLLECTING PERTURBED BY CRISPR-CAS9 GENE NAMES --------------------= ###
#########################################################################################
loginfo( "Collecting gene names perturbed by CRISPR-CAS9", logger = "preprocessing_run.module" )

# Cleans the gene names from library identifiers in CRISPR-CAS9 experiment
CleanGeneNames = function( genes, pattern_lib_identifier, pattern_to_clean ) {
  
  gene_names_cleaned =  genes %>% keep( str_detect(., pattern = pattern_lib_identifier) ) %>% 
    set_names(., .) %>% 
    modify(., function(x) { 
      gsub(pattern_to_clean, x, replacement = "", perl = TRUE)
    })
  return( gene_names_cleaned )
}

# All features ending with "-gene" are CRISPR-CAS9 perturbed genes. 
# 3 types of gRNA are used: controls non-target, essential genes and TCR library genes
edited_genes = rownames(tcr_cropseq_seurat_annotated) %>% keep( str_detect(., pattern = "-gene$") )

controls_edited_genes    = CleanGeneNames(edited_genes,  "^CTRL", "-gene")
essential_edited_genes   = CleanGeneNames(edited_genes, "^Essential", "Essential-library-|-\\d|-gene|")
tcr_library_edited_genes = CleanGeneNames(edited_genes, "^Tcrlibrary-", "Tcrlibrary-|-\\d|-gene")

# With these selection, we exclude Cas9-blast-gene edits (do not need them for analysis)
combined_edited_genes = c(controls_edited_genes, essential_edited_genes, tcr_library_edited_genes)

##################################################################################################################
### ------------ SELECTING CELLS WHERE ONLY ONE GENE IS PERTURBED BY CRISPR/CAS9 ----------------------------- ###
##################################################################################################################
loginfo( "Selecting cells with only one perturbed gene", logger = "preprocessing_run.module" )

# From provided dataset, returns the barcodes for cells where only one edit of CRISPR/CAS9 happened
SelectSingleEdits = function( dataset, cell_annotations, edited_genes, one_gene ) {
  
  barcodes = colnames( dataset ) 
  # select only the cells with non-zero RNA counts for a perturbed gene, with only this specified gene perturbed
  edited_cells = FetchData( dataset, edited_genes ) %>% mutate( sum = rowSums(.) ) %>%
                                                        rownames_to_column( "barcode" ) %>%
                                                        filter_at( vars(one_gene), all_vars(. != 0) ) %>%
                                                        dplyr::filter( get(one_gene) == sum ) %>% 
                                                        dplyr::select( barcode, one_gene, sum )
  
  if ( dim(edited_cells)[1] ) { 
    assigned_barcodes = barcodes[as.numeric(edited_cells$barcode)]
    annotations = cell_annotations %>% dplyr::filter(Cell.barcode %in% assigned_barcodes)
    return( list(annotations) )
  }
  
  return( c() )
}

seurat_edited_genes = subset( tcr_cropseq_seurat_annotated, features = edited_genes )

# Subset the cells, with only one gene edited, for all genes
cells_with_one_edits_sep = list()
for ( i in names(combined_edited_genes) ) {
  cells_with_one_edits_sep[i] = SelectSingleEdits( tcr_cropseq_seurat_annotated, 
                                                   tcr_cell_annotation, edited_genes, 
                                                   one_gene = i )
}

# Combine edits with different gRNA but the same genes (for essential and TCR library). 
list_names_to_reassign = combined_edited_genes %>% keep( names(.) %in% names(cells_with_one_edits_sep) )
names( cells_with_one_edits_sep ) = list_names_to_reassign

cells_with_one_edits_merged = list()

for ( i in unique(list_names_to_reassign) ) { 
  new_elem = cells_with_one_edits_sep %>% keep( names(.) == i )
  new_elem_combined = list()
  for ( elem in new_elem ) {
    new_elem_combined = rbind( new_elem_combined, elem )
  }
  
  cells_with_one_edits_merged[[i]] = new_elem_combined
}

#############################################################################################################
### ------------ SELECTING CELLS WHERE BOTH (TCR) NAIVE AND STIMULATED SET ARE PRESENT------------------- ###
#############################################################################################################
loginfo( "Selecting cells where both (TCR) naive and stimulated set are present",
         logger = "preprocessing_run.module" ) 

counts_conditions    = sapply( cells_with_one_edits_merged, function(x) table( x$Sample ) )
conditions           = rownames( counts_conditions )
counts_conditions    = as_tibble( counts_conditions ) %>% select_if( ~all(.!= 0) )
edited_genes_to_keep = colnames( counts_conditions )

##############################################################################################################
### ------------ SELECTING DATA FOR 3 SEPARATE SETS: CONTROLS, ESSENTIAL GENES, TCR LIBRARY -------------- ###
##############################################################################################################
loginfo("Selecting data for 3 separate sets: controls, essential genes, TCR library")

# Select barcodes for postprocessed gene list: CONTROLS
control_genes_to_keep = edited_genes_to_keep %>% keep( . %in% controls_edited_genes )
control_genes_to_keep = cells_with_one_edits_merged[control_genes_to_keep]                              
control_genes_to_keep = do.call(rbind, control_genes_to_keep)

# Select barcodes for postprocessed gene list : ESSENTIAL
essential_genes_to_keep = edited_genes_to_keep %>% keep( . %in% essential_edited_genes )
essential_genes_to_keep = cells_with_one_edits_merged[essential_genes_to_keep]

# Select barcodes for postprocessed gene list: TCR signaling
names_tcr_genes_to_keep = edited_genes_to_keep %>% keep( . %in% tcr_library_edited_genes )
tcr_genes_to_keep = cells_with_one_edits_merged[names_tcr_genes_to_keep] 

# -- Add extra column to each list with gene names, for convenience in further processing
tcr_genes_to_keep = lapply( c(1:length( tcr_genes_to_keep) ), function(x) {
  gene_name = names( tcr_genes_to_keep[x] )
  gene_name_rep = tcr_genes_to_keep %>% pluck(x) %>% dim() %>% pluck(1)
  tcr_genes_to_keep %>% pluck(x) %>% bind_cols(gene = rep(gene_name, gene_name_rep))
})
names(tcr_genes_to_keep) = names_tcr_genes_to_keep


###############################################################################################################
### -------- RUNNING VIPER/DOROTHEA ----------------------------------------------------------------------- ###
###############################################################################################################
loginfo( "Running VIPER/DoRothEA", logger = "preprocessing_run.module" )

# Collect the averaged signatures for different conditions on the basis of provided barcodes
CollectAverageSignatures = function( dataset, preselected_annotations, gene_name, conditions ) {
  
  CollectAverageSignature = function( dataset, preselected_annotations, gene_name, condition ) {
    #Choose barcodes for a specific condition, specific CRISPR-CAS9 edit
    preselected_annotations = preselected_annotations %>% keep( names(.) == gene_name ) %>% 
      pluck( gene_name ) %>%
      dplyr::filter( Sample == condition )
    
    genes_data_one_gene = subset( dataset, cells = preselected_annotations$Cell.barcode )
    averaged_expression_one_gene_condition = AverageExpression( genes_data_one_gene, return.seurat = TRUE )
    
    return( averaged_expression_one_gene_condition[["RNA"]]@data )
  }
  
  signatures_conditions = lapply(conditions, function(x) { CollectAverageSignature(dataset, preselected_annotations,
                                                                                   gene_name, x) 
  })
  
  return( signatures_conditions )
}

### -------- Collecting averaged signatures for 2 conditions------------------------------------------ ###
averaged_gene_expression_naive = c()
averaged_gene_expression_stim = c()

for ( i in names(tcr_genes_to_keep ) ) { 
  tmp_res = invisible( CollectAverageSignatures( tcr_cropseq_seurat_annotated, tcr_genes_to_keep, i, conditions ) ) 
  naive_res = tmp_res %>% pluck(1) %>% as.matrix()
  stim_res = tmp_res %>% pluck(2) %>% as.matrix()
  
  averaged_gene_expression_naive = cbind( averaged_gene_expression_naive, naive_res ) 
  averaged_gene_expression_stim = cbind( averaged_gene_expression_stim, stim_res )
} 

colnames( averaged_gene_expression_naive ) = unique( names(tcr_genes_to_keep) )
colnames( averaged_gene_expression_stim ) = unique( names(tcr_genes_to_keep) )


### -------- LOADING DOROTHEA ----------------------------------------------------------------------- ###
loginfo("Loading DoRothEA and starting VIPER", logger = "preprocessing_run.module")
dorothea_regulon_human = read_csv( dorothea_path )
regulon = dorothea_regulon_human %>% dplyr::filter( confidence %in% c("A","B","C") ) %>% df2regulon()

### -------- RUNNING_VIPER ----------------------------------------------------------------------- ###
tcr_genes_viper_naive      = viper(eset = averaged_gene_expression_naive, regulon = regulon, method = "scale", 
                                   minsize = 4, 
                                   eset.filter = FALSE, 
                                   cores = 1, 
                                   verbose = FALSE)

tcr_genes_viper_stimulated = viper(eset = averaged_gene_expression_stim,  regulon = regulon, method = "scale", 
                                   minsize = 4, 
                                   eset.filter = FALSE, 
                                   cores = 1, 
                                   verbose = FALSE)

loginfo( "Translating IDs", logger = "preprocessing_run.module" )
uniprot_ids = TranslateIds( names(tcr_genes_to_keep) )

loginfo("Saving the results of preprocessing", logger = "preprocessing_run.module")

# Saving intermediate results to RData file
save(trc_crispr_data,
     tcr_cropseq_seurat_subset,
     tcr_cropseq_seurat_annotated,
     combined_edited_genes,
     counts_conditions,
     control_genes_to_keep,
     essential_genes_to_keep,
     tcr_genes_to_keep,
     averaged_gene_expression_naive,
     averaged_gene_expression_stim,
     tcr_genes_viper_naive,
     tcr_genes_viper_stimulated,
     uniprot_ids,
     file = Rdata_file)

loginfo("PREPROCESSING DONE", logger = "preprocessing_run.module")
