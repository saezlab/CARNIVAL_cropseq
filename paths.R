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
# DESCRIPTION
# File contains the paths needed for CARNIVAL analysis for CropSeq data. 
# To download the last version of the data use the link below: 
# https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE137554 
# ====================

if ( !require("here") ) { 
  install.packages("here")
}

if ( !exists( "opt" ) & "here" %in% (.packages()) ) {
  base_folder = here()
} else if ( !exists( "opt" ) ) {
  base_folder = ""
}

if ( exists( "opt" ) ) { 
  run_preprocessing = unlist( opt["run-preprocessing"] )
  run_carnival = unlist( opt["run-carnival"] )
  test_run = unlist( opt["test"] )
} 

directories_to_check = c()
files_to_check = c()

########################################################################################
### ------------ Update if running either preprocessing/carnival part or both------- ###
########################################################################################

if ( exists("opt") && unlist(opt["default-input-output"]) ) {
  input_folder = file.path( getwd(), "input" )
  output_folder = file.path( getwd(), "output" )
} else { 
  
  if ( exists( "opt" ) && unlist(opt["input-folder"]) != "" )  {
    input_folder = unlist( opt["input-folder"] )
  } else {
    input_folder  = file.path( base_folder, "input" )  
  }
  
  if ( exists( "opt" ) && unlist(opt["output-folder"]) != "" )  {
    output_folder = unlist( opt["output-folder"] )
  } else {
    output_folder  = file.path( base_folder, "output" )  
  }
  
}

if ( exists( "opt" ) )  {
  Rdata_file = opt["rdata-file"]
} 

intermediate_results_folder = file.path( output_folder, "intermediate_results/" )
logfile = file.path( output_folder, paste0( "carnival_run_", format(Sys.time(), "%d_%m_%Y_%H_%M"), ".log" ) )

directories_to_check = c( input_folder, output_folder, intermediate_results_folder )

########################################################################################
### ------------ Update only if running preprocessing part ------------------------- ###
########################################################################################
# Update the path to your own path with downloaded data
data_path = "/Users/olgaivanova/_WORK/_PhD_Heidlbrg_data/CROPseq_data_updated"

raw_crispr_hdf5 = file.path( data_path, "GSE137554_raw_gene_bc_matrices_h5.h5" )
annotated_filename = file.path( data_path, "GSE137554_CellAnnotation.tsv" )

# DoRothEA link to download human regulon
dorothea_path = "https://raw.githubusercontent.com/saezlab/ConservedFootprints/master/data/dorothea_benchmark/regulons/dorothea_regulon_human_v1.csv"

if ( run_preprocessing ) {
  Rdata_file = file.path( intermediate_results_folder, 
                          paste0("carnival_run_", format(Sys.time(), "%d_%m_%Y_%H_%M"), ".Rdata") )
  files_to_check = c( raw_crispr_hdf5,
                      annotated_filename )  
}

########################################################################################
### ------------ Update if running carnival part ----------------------------------  ###
########################################################################################
# PKN file to run CARNIVAL
omnipath_filename = file.path( input_folder, "PKNs/omnipath_10_02_2020_13_46.txt" )
dorothea_tf_mapping_filename = file.path( input_folder, "dorothea_TF_mapping.csv" )

# CARNIVAL will be installed from github using the link below. 
# TODO Once CARNIVAL is accepted on Bioconductor, it can be loaded/installed as other packages
CARNIVAL_installation_path = "https://github.com/saezlab/CARNIVAL-Bioconductor-Dev"
cplex_solver_path = "/Applications/CPLEX_Studio1210/cplex/bin/x86-64_osx/cplex"
output_directory_carnival = file.path( output_folder, "Results_CARNIVAL_massive/" )

if ( !run_preprocessing ) {
  if ( Rdata_file == "" ){
    print( "No preprocessing requested and no Rdata file provided." )
    print( paste0("Trying to find R data file in the output folder: ", output_folder) )
    all_Rdata_files = list.files( path = output_folder, pattern = "\\.Rdata", 
                                  recursive = TRUE, full.names = TRUE )
    Rdata_file  = all_Rdata_files[ length( all_Rdata_files ) ]
  }
  files_to_check = c( files_to_check, Rdata_file )
}

print(Rdata_file)

if ( run_carnival ) {
  directories_to_check = c( directories_to_check, 
                            output_directory_carnival )
  
  files_to_check = c( files_to_check, cplex_solver_path, omnipath_filename, 
                      dorothea_tf_mapping_filename )  
}

########################################################################################
### ------------ CHECKING PROVIDED DIRECTORIES AND FILES  -------------------------- ###
########################################################################################

CreateDirectoriesIfDontExist = function( sub_dir, main_dir = "", needed_main_dir = F ) {
  if ( main_dir == "" & needed_main_dir ) {
    main_dir = getwd()  
  }
  dir.create( file.path(main_dir, sub_dir), showWarnings = FALSE )
}

CheckFile = function( filename, path = "" ) {
  full_filename = file.path( path, filename )
  if( !file.exists(full_filename) ) {
    stop( paste( "Cannot find the speficied file! Check provided path:", full_filename ) )
  }
}

tryCatch( 
  expr = { 
      invisible( sapply(directories_to_check, CreateDirectoriesIfDontExist) )
  }, 
  error = function(e) { 
      message( "Cannot create directories. Please specify correct paths." )
      print(e)
      stop()
  }
)

invisible( sapply(files_to_check, CheckFile) )

print( paste0("All directories for a run:", paste0(directories_to_check, collapse = "\n, ")) )
print( paste0("All specified input files for a run:", paste0(files_to_check, collapse = "\n, ")) )