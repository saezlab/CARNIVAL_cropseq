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

# Reassigne source path if the script is being executed from run_pipeline_cropseq.R
if ( exists( "opt" ) && opt["source-path"] != "" )  {
  source_path = opt['source-path']
} else if ( "here" %in% (.packages()) ) {
  source_path = here()
  run_carnival = TRUE
} else {
  source_path = ""  
  run_carnival = TRUE
}

if ( exists( "opt" ) ) {
  carnival_threads =  opt["carnival-threads"]
  test_run = opt["test"]
} else {
  carnival_threads = 0
  test_run = FALSE
}


########################################################################################
### ------------ INSTALLING/LOADING NECESSARY PACKAGES ----------------------------- ###
########################################################################################
source( file.path(source_path, "packages_utils.R") )
source( file.path(source_path, "paths.R") )

cran_list_packages = c("dplyr", "logging")
bioc_list_packages = c("biomaRt", "UniProt.ws", "OmnipathR")
github_packages    = c("CARNIVAL" = CARNIVAL_installation_path)
CheckAndLoadLibraries(  cran_list_packages, bioc_list_packages, github_packages )

basicConfig(level = "DEBUG")
addHandler(writeToFile, logger = "CARNIVAL_run", file = logfile)
loginfo("CARNIVAL script started", logger = "CARNIVAL_run.module")
loginfo( paste0("Running with the setup: source path: ", source_path, ";", 
               " Is it a test run: ", test_run, ";",
               " N threads:", carnival_threads), 
         logger = "CARNIVAL_run.module" )

source( file.path(source_path, "utils_cropseq.R") )

########################################################################################
### ------------ RUNNING PREPROCESSING IF NEEDED ----------------------------------- ###
########################################################################################

if ( exists("opt") && unlist( opt["run-preprocessing"] ) ) { 
  source( file.path(source_path, "preprocessing_cropseq.R") )
} else {
  
  if( file.exists(Rdata_file) ) {
    loginfo( "Loading preprocessed Rdata file...", logger = "CARNIVAL_run.module" )
    load(Rdata_file)
  } else { 
    loginfo( "Cannot load RData file with preprocessed data, trying to run preprocessing first...",
             logger = "CARNIVAL_run.module" )
    source( file.path(source_path, "preprocessing_cropseq.R") )
  }
}

#########################################################################################################
###  ---------- RUNNING CARNIVAL ON AVERAGE SIGNATURES FOR EACH GENE PERTURBATION ------------------- ###
#########################################################################################################

RunCarnivalOneTime = function( edited_gene_name, prior_knowledge_network, viper_scores, 
                               perturbations, save_outfile = FALSE, 
                               output_filename = "out_carnival.csv",
                               output_dir = "", dot_figures = FALSE, 
                               threads = 0) { 
  
  res_carnival = runCARNIVAL(solverPath = cplex_solver_path, 
                             netObj = prior_knowledge_network, 
                             measObj = viper_scores, 
                             inputObj = perturbations,
                             solver = "cplex", 
                             dir_name = output_dir,
                             DOTfig = FALSE,
                             threads = threads)
  if ( save_outfile ) { 
    results_carnival = res_carnival$weightedSIF %>% select("Node1", "Node2", "Sign")
    write.csv(results_carnival, file = paste0(output_dir, output_filename), 
              quote = FALSE, row.names = FALSE) 
  }
  
  return( res_carnival )
}

RunCarnivalOnListGenes = function( uniprot_ids, tcr_genes_viper, prior_knowledge_network, 
                                   threads = 0) {
  res_carnivals_genes = list() 
  for ( i in uniprot_ids$GENES ) { 
    loginfo( paste("Running CARNIVAL for naive (TCR) data, gene: ", i), 
             logger = "CARNIVAL_run.module" )
    uniprot_id = uniprot_ids %>% filter( GENES == i ) 
    perturbations = data.frame( "1" )
    names(perturbations) = uniprot_id$UNIPROTKB
    
    viper_scores = tcr_genes_viper %>% as_tibble( rownames = "id" ) %>% 
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
                                  threads = threads,
                                  dot_figures = TRUE)
    
    res_carnivals_genes[[i]] = res_carn
  }
  return( res_carnivals_genes )
}

loginfo( "Reading/requesting prior knowledge network", logger = "CARNIVAL_run.module" )
prior_knowledge_network = LoadPKNForCarnival( omnipath_filename )

loginfo( "Translating IDs", logger = "CARNIVAL_run.module" )
uniprot_ids = TranslateIds( names(tcr_genes_to_keep) )

if ( test_run ) {
  loginfo( "Test run of CARNIVAL is initiated", logger = "CARNIVAL_run.module" )
  RunCarnivalOnListGenes( uniprot_ids[1], tcr_genes_viper_naive, prior_knowledge_network, carnival_threads)
  loginfo( "Test run of CARNIVAL is finished", logger = "CARNIVAL_run.module" )
} else {
  RunCarnivalOnListGenes( uniprot_ids, tcr_genes_viper_naive, prior_knowledge_network, carnival_threads )
  RunCarnivalOnListGenes( uniprot_ids, tcr_genes_viper_stimulated, prior_knowledge_network, carnival_threads )  
}

loginfo( "ALL CARNIVAL RUNS DONE", logger = "CARNIVAL_run.module" )
