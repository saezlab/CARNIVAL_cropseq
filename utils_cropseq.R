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
# DESCRIPTION: Utility functions for CropSeq analysis with CARNIVAL.
# ====================

if ( !require("here") ) { 
  install.packages("here")
}

# Reassign source path if the script is being executed from run_pipeline_cropseq.R
# if ( exists( "opt" ) && opt["source-path"] != "" )  {
#   source_path = opt['source-path']
# } else if ( "here" %in% (.packages()) ) {
#   source_path = here()
# } else {
#   source_path = ""  
# }

source( file.path(source_folder, "packages_utils.R") )

cran_list_packages = c("dplyr", "stringr", "purrr")
bioc_list_packages = c("biomaRt", "UniProt.ws", "OmnipathR")
CheckAndLoadLibraries( cran_list_packages, bioc_list_packages )


# Reads Omnipath as a prior knowledge network for CARNIVAL run
LoadPKNForCarnival = function( path_file = "", filter_by_references = 0 ) {

  if ( !file.exists(path_file) ) {
    omnipath = import_Omnipath_Interactions(from_cache_file = NULL,
                                            filter_databases = get_interaction_databases(),
                                            select_organism = 9606)

    write.csv(omnipath, file = paste0( "omnipath", "_", format(Sys.time(), "%d_%m_%Y_%H_%M"), ".txt" ),
                        quote = FALSE, row.names = FALSE)
  } else {
    omnipath = read.csv(file = path_file)
  }

  #TODO check with Denes/Enio that I am using the right columns here
  omnipath_pkn = omnipath %>% dplyr::filter( is_directed == 1 ) %>%
                              dplyr::mutate_at( c("nrefs"), 
                                                ~as.numeric(as.character(.))
                                                ) %>%
                              dplyr::filter( nrefs >= filter_by_references ) %>%
                              dplyr::select( "source", "target",
                                            "consensus_stimulation",
                                            "consensus_inhibition" ) %>%
                              dplyr::mutate_at( c("consensus_stimulation",
                                                 "consensus_inhibition"),
                                                ~as.numeric(as.character(.))
                                              )

  omnipath_pkn$consensus_stimulation = as.numeric(as.character(omnipath_pkn$consensus_stimulation))
  omnipath_pkn$consensus_inhibition = as.numeric(as.character(omnipath_pkn$consensus_inhibition))

  # Keep only interactions with -1 or 1
  omnipath_pkn = omnipath_pkn %>% mutate( consensus_inhibition_upd = ifelse(consensus_inhibition == 1, -1, 0) ) %>%
                                  mutate( effect = consensus_stimulation + consensus_inhibition_upd ) %>%
                                  dplyr::select( "source", "effect", "target" ) %>%
                                  dplyr::filter( effect != 0 )

  omnipath_pkn = omnipath_pkn %>% filter_at( vars(source, target),
                                             all_vars(!str_detect(., "COMPLEX:")) )
  

  return( omnipath_pkn )
}


# Translates gene names to uniprot IDs. It will select only curated IDs (that have a corresponding Swissprot id). 
# If multiple IDs are available, it will return the one with the highest annotation score.
# Because each request to the database takes a long time, it is better to request all gene_names at the same time
TranslateIds = function( gene_names, only_first = FALSE){

  ensembl = useMart( 'ensembl', dataset = "hsapiens_gene_ensembl" )
  up = UniProt.ws(taxId = 9606)

  # Returns the uniprot IDs of the corresponding proteins listed in VIPER/DoRothEA results
  # Because BioMart doesn't allow to request more than 3 fields, there are two request
  requested_mart = getBM(attributes = c("uniprot_gn_id",
                                        "uniprot_gn_symbol",
                                        "uniprotswissprot"),
                             values = gene_names,
                             filters = "uniprot_gn_symbol",
                             mart = ensembl)

  # select only curated ids
  requested_mart = requested_mart %>% filter(uniprotswissprot != '')

  # collect annotation scores from uniprot database
  # more information on annotation:
  # https://www.ebi.ac.uk/training/online/course/uniprot-exploring-protein-sequence-and-functional/exploring-uniprotkb-entry/annotation-score
  ids = select(up, columns = c("UNIPROTKB", "GENES", "SCORE"), keys = requested_mart$uniprot_gn_id)
  ids = ids %>% mutate( SCORE_upd = str_remove(SCORE, " out of \\d") ) %>%
                mutate_at( vars(SCORE_upd), as.numeric )

  # If mapping on gene names is not unique, keep only the first gene from the list
  # The first in the list should be the one originally requested.
  ids$GENES = sapply( ids$GENES, function(x) strsplit(x, " ") %>% pluck(1, 1) )

  uniprot_id = ids %>% group_by(GENES) %>% filter( SCORE_upd == max(SCORE_upd) )

  return( uniprot_id )
}

# Reads the curated mapping between dorothea regulon and uniprots, 
# selects only the first value from this list if several mapping are available
#TODO update: it should accept all values, duplicating the vipers scores for CARNIVAl run then
ReadDorotheaMapping = function( ids, dorothea_tf_mapping_filename ) {

  dorothea_tf_mapping = read.csv( dorothea_tf_mapping_filename )
  dorothea_unique_mappings = dorothea_tf_mapping[unique(dorothea_tf_mapping$SYMBOL), ]
  dorothea_unique_mappings = dorothea_unique_mappings %>% filter(SYMBOL %in% ids)

  return( dorothea_unique_mappings )
}
