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
  }, error = function(e) {
    print("Unable to install package 'here'. Continuing without it.")
  })
}

if ( !require("optparse") ) { 
  install.packages("optparse")
}

option_list = list(
  make_option( c("-p", "--run-preprocessing"), type = "logical", default = FALSE, action = "store_true",
                help = "Specify if the data should be preprocessed first (Seurat + VIPER/DoRothEA)", 
                metavar = "logical" ), 
  make_option( c("-t", "--test"), type = "logical", default = FALSE, action = "store_true",
               help = "Run CARNIVAL for one gene. Should be run together with -p option or RData file (in paths.R) 
                       with preprocessed data should be prodived", metavar = "logical" ),
  make_option( c("-c", "--run-carnival"), type = "logical", default = FALSE, action = "store_true",
               help = "Provide source path for the run", metavar = "logical" ), 
  make_option( c("-s", "--source-path"), type = "character", default = "", 
                help = "Provide source path for the run", metavar = "character" ), 
  make_option( c("-r", "--carnival-threads"), type = "numeric", default = 0, 
                help = "Provide number of threads to run CARNIVAL", metavar = "numeric" )
)

opt_parser = OptionParser(option_list = option_list);
opt = parse_args(opt_parser)

print( "Initiating the script with parameters: " )
print( opt )

if ( !( unlist(opt["run-preprocessing"]) |
        unlist(opt["test"]) | 
        unlist(opt["run-carnival"]) )    ) {
  print_help( opt_parser )
  stop("Nothing to run. At least one argument must be supplied: -p, -t or -c", call. = FALSE)
}

if ( opt["source-path"] != "") {
  source_path = opt["source-path"]
} else if ( "here" %in% (.packages()) ) {
  source_path = here()
} else {
  source_path = ""  
}

print(source_path)
if ( unlist( opt["run-preprocessing"] ) ) {
  source( file.path(source_path, "preprocessing_cropseq.R") )
}

if ( unlist( opt["run-carnival"] ) | 
     unlist( opt["test"] )       ) {
  source( file.path(source_path, "run_carnival_cropseq.R") )
}