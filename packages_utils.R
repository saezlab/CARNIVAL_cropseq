#!/usr/bin/env Rscript  
# coding: utf-8
# Copyright (C) 2020 Olga Ivanova
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
# DESCRIPTION: Utility functions for handling libraries in R.
# ====================
# NB setting a cran mirror is necessary for running on the cluster
cran_mirrors = "https://cran.uni-muenster.de/"

CheckAndLoadLibraries = function(cran_list_packages = "", bioc_list_packages = "", github_packages = "") {
  all_packages = c( cran_list_packages, bioc_list_packages, names(github_packages) )
  all_packages = all_packages[all_packages != ""]
  new_packages = all_packages[!(all_packages %in% installed.packages()[,"Package"])]
  
  if ( length(new_packages) > 0 ) {
    print( paste0("New packages to install:", paste0(new_packages, collapse = ", ")) )  
  }

  cran_packages_to_install = cran_list_packages[cran_list_packages %in% new_packages]
  bioc_packages_to_install = bioc_list_packages[bioc_list_packages %in% new_packages]
  github_packages_to_install = github_packages[names(github_packages) %in% new_packages]
  
  if ( length(cran_packages_to_install) ) {
    
    utils::install.packages( cran_packages_to_install, repos = cran_mirrors ) 
  } else if ( length(bioc_packages_to_install) ) {
    
    BiocManager::install(bioc_packages_to_install)
  } else if ( length(github_packages_to_install) ) {
    
    for (package in github_packages_to_install) {
      devtools::install_github(package, build_vignettes = FALSE)
    }
  }
    
  invisible( sapply( all_packages, function(x) {
    if ( !x %in% github_packages & x != "") {
      library( x, character.only = TRUE )
    }
  }))
  
  invisible( sapply( names(github_packages), function(x) {
    if (x != "") {
      library( x, character.only = TRUE )
    }
  }))
  
}

cran_list_packages = c("BiocManager", "devtools")
CheckAndLoadLibraries( cran_list_packages = c("BiocManager", "devtools") )
