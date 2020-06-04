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
# Load preprocessed data and run CARNIVAL analysis for CropSeq single cell data.
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
    # After reading yml file, this will do preliminary set ups and checks for folders/files. 
    source( file.path(source_folder, "setting_up_pipeline.R") )
    source( file.path(source_folder, "utils_cropseq.R") )
  } else {
    stop("run_carnival_cropseq.R: can't continue, please specify correct source folder")
  }
}

########################################################################################
### ------------ INSTALLING/LOADING NECESSARY PACKAGES ----------------------------- ###
########################################################################################
setwd(working_dir)
print( paste("Current working directory:", getwd()) ) 

if ( dir.exists(source_folder) ) {
  source( file.path(source_folder, "packages_utils.R") )  
} else {
  stop("run_carnival_cropseq.R: can't continue, please specify correct source folder")
}

cran_list_packages = c("dplyr", "logging", "tidyr", "yaml")
bioc_list_packages = c("OmnipathR")
# TODO this will be removed once there will be a setup of directory for temp files in Bioconductor version of CARNIVAL
github_packages    = c("CARNIVAL" = CARNIVAL_installation_path)
CheckAndLoadLibraries(  cran_list_packages, bioc_list_packages, github_packages )

basicConfig(level = "DEBUG")
addHandler(writeToFile, logger = "CARNIVAL_run", file = logfile)
loginfo("CARNIVAL script started", logger = "CARNIVAL_run.module")
loginfo( paste0("Running CARNIVAL script with a setup: source path: ", source_folder, ";", 
               " Is it a test run: ", test_run, ";",
               " N threads:", carnival_threads, ";",
               " Start id: ", start_id, ";",
               " End id: ", end_id), 
         logger = "CARNIVAL_run.module" )

########################################################################################
### ------------ READING PREPROCESSED DATA OR RUNNING PREPROCESSING IF NEEDED ------ ###
########################################################################################
if ( file.exists(Rdata_file) ) {
    loginfo( "Loading preprocessed Rdata file", logger = "CARNIVAL_run.module" )
    load(Rdata_file)
} else { 
    loginfo( "Cannot load Rdata file with preprocessed data, trying to run preprocessing first",
             logger = "CARNIVAL_run.module" )
    source( file.path(source_folder, "preprocessing_cropseq.R") )
}

#########################################################################################################
###  ---------- RUNNING CARNIVAL ON AVERAGE SIGNATURES FOR EACH GENE PERTURBATION ------------------- ###
#########################################################################################################

RunCarnivalOneTime = function( edited_gene_name, prior_knowledge_network, viper_scores, 
                               perturbations, output_filename = "out_carnival.csv",
                               output_dir = "", produce_dot_figure = FALSE, 
                               threads = 0, save_outfile = TRUE,
                               save_image = TRUE) { 
  
  res_carnival = runCARNIVAL(solverPath = solver_path, 
                             netObj = prior_knowledge_network, 
                             measObj = viper_scores, 
                             inputObj = perturbations,
                             solver = "cplex", 
                             dir_name = output_dir,
                             DOTfig = produce_dot_figure,
                             timelimit = carnival_timelimit,
                             threads = threads)
  if ( save_outfile ) { 
    results_carnival = res_carnival$weightedSIF %>% as_tibble() 
    results_carnival_nodes_attributes = results_nodes_attributes$nodesAttributes %>% as_tibble()
    
    write.csv( results_carnival, file = file.path( output_dir, output_filename ), 
               quote = FALSE, row.names = FALSE ) 
    
    attribute_file_name = paste0( strsplit(output_filename, "\\.")[[1]][0], "nodes_attributes" ) #get the filename without extension
    attribute_file_extension = strsplit(output_filename, "\\.")[[1]][1]
    write.csv( results_carnival_nodes_attributes, file = file.path( output_dir, attribute_file_name, attribute_file_extension), 
               quote = FALSE, row.names = FALSE ) 
  } 

  if ( save_image ) {
    save.image( file = paste0( strsplit(output_filename, "\\.")[[1]][0], ".RData") ) 
  }

  return( res_carnival )
}

CollectPerturbationsData = function(crispr_gene) {
  perturbations_file = list.files( path = perturbations_folder, pattern = crispr_gene, 
                                   full.names = T)[1]
  perturbations_data = read.csv2( perturbations_file )
  genes = colnames( perturbations_data ) 
  
  translated_genes = lapply( genes, function(x) {
    if ( x %in% uniprot_ids$GENES ) {
      translated_id = ( uniprot_ids %>% filter(GENES == x) )$UNIPROTKB
      return( translated_id )
    } else {
      return( TranslateIds( x ) ) 
    }
  })
  
  uniprot_ids %>% bind_rows()
  
  colnames( perturbations_data ) = unlist( translated_genes )
  return( perturbations_data )
}

RunCarnivalOnListGenes = function( uniprot_ids, perturbations, 
                                   tcr_genes_viper, prior_knowledge_network, 
                                   threads = 0, 
                                   naive = FALSE) {
  res_carnivals_genes = list() 
  
  if ( naive ) {
    file_prefix = "naive"
  } else { 
    file_prefix = "stimulated"
  }
  
  for ( i in uniprot_ids$GENES ) { 
    loginfo( paste("Running CARNIVAL for (TCR) data, gene: ", i), 
             logger = "CARNIVAL_run.module" )
    tryCatch({
      viper_scores = tcr_genes_viper %>% as_tibble( rownames = "id" ) %>% 
                                         dplyr::select( id, all_of(i) ) %>% 
                                         drop_na() %>% 
                                         spread(id, i) %>% 
                                         as.data.frame()
    
      tf_names = colnames(viper_scores)
      colnames(viper_scores) = ReadDorotheaMapping( tf_names, dorothea_mapping_file )$UNIPROT
  
      res_carn = RunCarnivalOneTime( i, prior_knowledge_network, viper_scores, perturbations,
                                    output_dir = file.path( output_folder, i ),
                                    output_filename = paste0("out_carnival_", file_prefix, i, ".csv"),
                                    produce_dot_figure = TRUE,
                                    threads = threads )
        
      res_carnivals_genes[[i]] = res_carn
    }, error = function( e ) {
        loginfo( paste("Cannot process data for", i, ":", e ), 
               logger = "CARNIVAL_run.module" )
    })
  }
  return( res_carnivals_genes )
}

loginfo( "Reading/requesting prior knowledge network", logger = "CARNIVAL_run.module" )
prior_knowledge_network = LoadPKNForCarnival( PKN_file, filter_by_references = 1 )
loginfo( paste("Prior knowledge network contains", dim(prior_knowledge_network)[1], "interactions" ), 
         logger = "CARNIVAL_run.module" )

if ( test_run ) {
  loginfo( "Test run of CARNIVAL is initiated", logger = "CARNIVAL_run.module" )
  # Test CARNIVAL with running on LAT gene (finishes fast)
  RunCarnivalOnListGenes( uniprot_ids[28, ], tcr_genes_viper_naive, prior_knowledge_network, carnival_threads)
  loginfo( "Test run of CARNIVAL is finished", logger = "CARNIVAL_run.module" )
} else {
  
  if ( run_naive ) {
    RunCarnivalOnListGenes( uniprot_ids[ c(start_id:end_id), ], tcr_genes_viper_naive, 
                            prior_knowledge_network, 
                            carnival_threads, naive = TRUE )  
  }
  
  if ( run_stimulated ) { 
    RunCarnivalOnListGenes( uniprot_ids[ c(start_id:end_id), ], tcr_genes_viper_stimulated, 
                            prior_knowledge_network, 
                            carnival_threads )  
  }
  
}

loginfo( "ALL CARNIVAL RUNS DONE", logger = "CARNIVAL_run.module" )
