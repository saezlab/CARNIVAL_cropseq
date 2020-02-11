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

if ( !require("optparse") ) { 
  install.packages("optparse")
}

library("optparse")

option_list = list(
  make_option( c("-t", "--test"), type = "character", default = FALSE, action = "store_true",
                help = "run CARNIVAL for one gene", metavar = "character" ),
  make_option( c("-p", "--preprocess"), type = "character", default = FALSE, action = "store_true",
                help = "Specify if the data should be preprocessed first (Seurat + VIPER/DoRothEA", 
                metavar = "character"),
  make_option( c("-s", "--source_path"), type = "character", default = "", 
                help = "Provide source path for the run", metavar = "character" ), 
  make_option( c("-r", "--carnival_threads"), type = "numeric", default = 0, 
                help = "Provide number of threads to run CARNIVAL", metavar = "numeric" )
)

opt_parser = OptionParser(option_list = option_list);
opt = parse_args(opt_parser)

print(opt)