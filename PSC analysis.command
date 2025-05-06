#!/bin/zsh
# launch PSC Analysis via Rscript so the tcltk GUI stays alive

# explicit path to the Rscript binary on macOS
RSCRIPT="/Library/Frameworks/R.framework/Resources/bin/Rscript"

"$RSCRIPT" --vanilla -e "
  # load/install packages
  load_required_packages <- function(pkgs) {
    new.pkgs <- setdiff(pkgs, rownames(installed.packages()))
    if (length(new.pkgs)) install.packages(new.pkgs)
    invisible(lapply(pkgs, library, character.only=TRUE))
  }
  load_required_packages(c(
    'dbscan','minpack.lm','Rcpp','robustbase',
    'shiny','signal','readABF','readxl',
    'tcltk','tkrplot','openxlsx'
  ))

  # source your GUI code
  source('~/Documents/Repositories/Rfits/nNLS functions.R')

  # launch the GUI (this blocks until you close the window)
  analysePSCtk()
"

