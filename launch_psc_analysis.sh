#!/bin/zsh
# launch an instance of R.app and run analysis
open -n -a R --args -e '
  rm(list = ls(all = TRUE));
  graphics.off();

  load_required_packages <- function(packages) {
    new.packages <- packages[!(packages %in% installed.packages()[, "Package"])]
    if (length(new.packages)) install.packages(new.packages)
    invisible(lapply(packages, library, character.only = TRUE))
  }

  required.packages <- c("dbscan", "minpack.lm", "Rcpp", "robustbase",
    "shiny", "signal", "readABF", "readxl", "tcltk", "tkrplot", "openxlsx")
  load_required_packages(required.packages)

  UserName <- Sys.getenv("USER")
  path_repository <- "/Documents/Repositories/Rfits"
  file_path <- paste0("/Users/", UserName, path_repository)
  source(paste0(file_path, "/nNLS functions.R"))

  analyseABFtk()
'