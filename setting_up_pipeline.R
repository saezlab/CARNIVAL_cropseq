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
# File contains the set up for the pipeline of DoRothEA/VIPER-CARNIVAL.  
# ====================
if ( !exists("working_dir") ) {
  # if working_dir is not specified from the console arguments, read from the .yml file
  working_dir = settings_run$working_dir  
}

base_path = settings_run$base_path
input_folder = file.path( base_path, settings_run$input_folder )

perturbations_folder = settings_run$perturbations_folder
raw_file = file.path( input_folder, settings_run$raw_file )
annotation_file = file.path( input_folder, settings_run$annotation_file )

dorothea_path = settings_run$dorothea_path
dorothea_mapping_file = file.path( input_folder, settings_run$dorothea_mapping_file )
PKN_file = file.path( input_folder, settings_run$PKN_file )

output_folder = file.path( base_path, settings_run$output_folder )
intermediate_results_folder = file.path( output_folder, "intermediate_results")

preprocessing = settings_run$preprocessing 
test_run = settings_run$test_run
carnival_run = settings_run$carnival_run


CARNIVAL_installation_path =  settings_run$CARNIVAL_installation_path
carnival_threads = settings_run$carnival_threads
carnival_timelimit = settings_run$carnival_timelimit
solver_path = settings_run$solver_path
Rdata_file = file.path ( output_folder, settings_run$Rdata_file )

PKN_filter_references = settings_run$PKN_filter_references

if ( !exists("start_id") || !exists("end_id") ) {
  # if start_id and end_id are not specified from the console,arguments read from the .yml file
  start_id = settings_run$start_id
  end_id = settings_run$end_id
  
}

run_naive = settings_run$run_naive
run_stimulated = settings_run$run_stimulated

########################################################################################
### ------------ FURTHER SETTING UP THE PARAMETERS --------------------------------- ###
########################################################################################
logfile = file.path( output_folder, paste0( "carnival_run_", format(Sys.time(), "%d_%m_%Y_%T"), ".log" ) )

if ( start_id == -1  ) {
  start_id = 1
}

if ( end_id == -1) {
  end_id = 30
}

directories_to_check = c( base_path, input_folder, perturbations_folder, 
                          output_folder, intermediate_results_folder )
directories_to_check = directories_to_check[ directories_to_check != '']

files_to_check = c()
if ( preprocessing ) {
  files_to_check = c( raw_file, annotation_file )
  Rdata_file = file.path( intermediate_results_folder, 
                          paste0("carnival_run_", format(Sys.time(), "%d_%m_%Y_%H_%M_%S"), ".Rdata") )
} else {
  files_to_check = c( Rdata_file )
}

files_to_check = c( files_to_check, dorothea_mapping_file, PKN_file, solver_path )

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
    stop( paste( "Cannot find the speficied file! Check provided files in settings:", full_filename ) )
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

cat("All directories for a run:", " \n ", paste0(directories_to_check, collapse = " \n "), " \n ")
cat("All specified input files for a run:", " \n ", paste0(files_to_check, collapse = " \n "), " \n ")