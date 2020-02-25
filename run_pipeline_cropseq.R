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
# DESCRIPTION: Reading console arguments for CARNIVAL runs.
# ====================

if ( !require("here") ) { 
  tryCatch( {
    install.packages("here")  
  }, error = function( e ) {
    print("Unable to install package 'here'. Continuing without it.")
  })
}

if ( !require("optparse") ) { 
  install.packages("optparse")
}

if ( !require("yaml") ) { 
  install.packages("yaml")
}

option_list = list(
  make_option( c("-f", "--settings-file"), type = "character", default = "settings.yml", 
               help = "Specify settings file location", 
               metavar = "character" ), 
  
  make_option( c("-t", "--test"), type = "logical", default = FALSE, action = "store_true",
               help = "Add the option if running on the server if running in a test mode",
               metavar = "logical" ),
  
  make_option( c("-l", "--local"), type = "logical", default = FALSE, action = "store_true",
               help = "Add the option if running on the local machine",
               metavar = "logical" ),
  
  make_option( c("-s", "--server"), type = "logical", default = FALSE, action = "store_true",
                help = "Add the option if running on the server",
                metavar = "logical" )
)

opt_parser = OptionParser(option_list = option_list);
opt = parse_args(opt_parser)

print( "Initiating the script with parameters: " )
print( opt )

server_run = unlist( opt["server"] )
test_run = unlist( opt["test"] )
local_run = unlist( opt["local"] )
settings_file = unlist( opt["settings-file"])

if ( !( server_run || test_run || local_run ) )   {
  print_help( opt_parser )
  stop("Nothing to run. At least one argument must be supplied: -r, -t, or -l", call. = FALSE)
}

if ( server_run & test_run & local_run ) {
  print_help( opt_parser )
  stop("Conflicting arguments provided. Please specify only one option from these: -r, -t, or -l", call. = FALSE)
}

if ( file.exists(settings_file) ) { 
  settings = yaml.load_file( settings_file )
} else { 
  print_help( opt_parser )
  stop("Please specify a correct location for settings file using -f option.", call. = FALSE)
}

if ( test_run ) { 
  settings_run = settings$test
} else if ( local_run ) {
  settings_run = settings$local
} else if ( server_run ) {
  settings_run = settings$server
}

print( cat("Specified options for a run:", paste0(settings_run, collapse = " \n "))  )

source_folder = settings_run$source_folder

if ( source_folder == "" && "here" %in% (.packages()) ) {
  source_folder = here()
  print( paste0("Source folder was not provided. Using the path provided by here() library:", source_folder) )
} else if (source_folder == "" && !"here" %in% (.packages()) ) {
  print( paste0("Source folder was not provided. Using: ", source_folder) )
}

if ( dir.exists(source_folder) ) { 
  source( file.path(source_folder, "setting_up_pipeline.R") )
  if ( settings_run$preprocessing ) {
    source( file.path(source_folder, "preprocessing_cropseq.R") )
  } 
  source( file.path(source_folder, "run_carnival_cropseq.R") )
} else {
  stop("Please specify a correct location for source R files in settings.yml", call. = FALSE)
}

