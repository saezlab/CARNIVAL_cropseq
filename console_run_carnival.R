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
# VIPER/DoRothEA + CARNIVAL Analysis for CropSeq single cell data.
# ====================
source("packages_utils.R")
source("paths.R")

########################################################################################
### ------------ INSTALLING/LOADING NECESSARY PACKAGES ----------------------------- ###
########################################################################################

cran_list_packages = c("dplyr","readr", "stringr", "purrr", "tibble", "tidyr",
                       "Seurat", "logging", "BiocManager", "devtools")
bioc_list_packages = c("viper", "biomaRt", "UniProt.ws", "rhdf5", "OmnipathR")
github_packages    = c("CARNIVAL" = CARNIVAL_installation_path)

CheckAndLoadLibraries( cran_list_packages, bioc_list_packages, github_packages )

basicConfig(level = "DEBUG")
addHandler(writeToFile, logger = "CARNIVAL_run", file = logfile)
loginfo("CARNIVAL script started", logger = "CARNIVAL_run.module")

source("utils_CARNIVAL.R")
source("utils_cropseq.R")

########################################################################################
### ------------ READING THE DATA -------------------------------------------------- ###
########################################################################################
loginfo( paste( "Reading the data files:", annotated_filename, 
                ";", raw_crispr_hdf5 ), 
         logger = "CARNIVAL_run.module" )

tcr_cell_annotation = read.delim( annotated_filename )
trc_crispr_data     = Read10X_h5( raw_crispr_hdf5 )

loginfo( "Running standard seurat pipeline", logger = "CARNIVAL_run.module" ) 

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
loginfo( "Collecing gene names perturbed by CRISPR-CAS9", logger = "CARNIVAL_run.module" )

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
loginfo( "Selecting cells with only one perturbed gene", logger = "CARNIVAL_run.module" )

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
    annotations = cell_annotations %>% filter(Cell.barcode %in% assigned_barcodes)
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
        logger = "CARNIVAL_run.module" ) 

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
loginfo( "Running VIPER/DoRothEA", logger = "CARNIVAL_run.module" )

# Collect the averaged signatures for different conditions on the basis of provided barcodes
CollectAverageSignatures = function( dataset, preselected_annotations, gene_name, conditions ) {
  
  CollectAverageSignature = function( dataset, preselected_annotations, gene_name, condition ) {
    #Choose barcodes for a specific condition, specific CRISPR-CAS9 edit
    preselected_annotations = preselected_annotations %>% keep( names(.) == gene_name ) %>% 
                                                          pluck( gene_name ) %>%
                                                          filter( Sample == condition )
    
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
dorothea_regulon_human = read_csv( dorothea_path )
regulon = dorothea_regulon_human %>% filter( confidence %in% c("A","B","C") ) %>% df2regulon()


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


#########################################################################################################
###  ---------- RUNNING CARNIVAL ON AVERAGE SIGNATURES FOR EACH GENE PERTURBATION ------------------- ###
#########################################################################################################

RunCarnivalOneTime = function( edited_gene_name, prior_knowledge_network, viper_scores, 
                               perturbations, save_outfile = FALSE, 
                               output_filename = "out_carnival.csv",
                               output_dir = "", dot_figures = FALSE ) { 
  
  res_carnival = runCARNIVAL(solverPath = cplex_solver_path, 
                             netObj = prior_knowledge_network, 
                             measObj = viper_scores, 
                             inputObj = perturbations,
                             solver = "cplex", 
                             dir_name = output_dir,
                             DOTfig = FALSE)
  if ( save_outfile ) { 
    results_carnival = res_carnival$weightedSIF %>% select("Node1", "Node2", "Sign")
    write.csv(results_carnival, file = paste0(output_dir, output_filename), 
              quote = FALSE, row.names = FALSE) 
  }
  
  return( res_carnival )
}

loginfo( "Reading/requesting prior knowledge network", logger = "CARNIVAL_run.module" )
prior_knowledge_network = LoadPKNForCarnival( omnipath_filename )

loginfo( "Translating IDs", logger = "CARNIVAL_run.module" )
uniprot_ids = TranslateIds( names(tcr_genes_to_keep) )

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
#TODO may be will be used on later stages
#load(Rdata_file)

for ( i in uniprot_ids$GENES[1] ) { 
  loginfo( paste("Running CARNIVAL for naive (TCR) data, gene: ", i), 
           logger = "CARNIVAL_run.module" )
  uniprot_id = uniprot_ids %>% filter( GENES == i ) 
  perturbations = data.frame( "1" )
  names(perturbations) = uniprot_id$UNIPROTKB
  
  viper_scores = tcr_genes_viper_naive %>% as_tibble( rownames = "id" ) %>% 
                                            dplyr::select( id, all_of(i) ) %>% 
                                            drop_na() %>% 
                                            spread(id, i) %>% 
                                            as.data.frame()
  
  tf_names = colnames(viper_scores)
  colnames(viper_scores) = ReadDorotheaMapping( tf_names, dorothea_tf_mapping_filename )$UNIPROT
  
  res_carn = RunCarnivalOneTime(i, prior_knowledge_network, viper_scores, perturbations,
                                output_dir = output_directory_carnival,
                                save_outfile = TRUE, 
                                output_filename = paste0("out_carnival_naive", i, ".csv"),
                                dot_figures = TRUE)
}


for ( i in uniprot_ids$GENES[1] ) { 
  loginfo( paste("Running CARNIVAL for perturbed (with TCR) data, gene,", i), 
           logger = "CARNIVAL_run.module" )
  uniprot_id = uniprot_ids %>% filter( GENES == i ) 
  perturbations = data.frame( "1" )
  names(perturbations) = uniprot_id$UNIPROTKB
  
  viper_scores = tcr_genes_viper_stimulated %>% as_tibble( rownames = "id" ) %>% 
                                            dplyr::select( id, all_of(i) ) %>% 
                                            drop_na() %>% 
                                            spread(id, i) %>% 
                                            as.data.frame()
  
  tf_names = colnames(viper_scores)
  colnames(viper_scores) = ReadDorotheaMapping( tf_names, dorothea_tf_mapping_filename )$UNIPROT
  
  res_carn = RunCarnivalOneTime(i, prior_knowledge_network, viper_scores, perturbations,
                                output_dir = output_directory_carnival,
                                save_outfile = TRUE, 
                                output_filename = paste0("out_carnival_perturbed", i, ".csv"),
                                dot_figures = TRUE)
}

