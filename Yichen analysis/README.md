# <center>`R` Code guide</center>


## Table of Contents
- [Step-by-Step Guide](#Step-by-Step-Guide)
   - [analysing `2P_PD 90-100d.xlsx`](#analysing-2p_pd-90-100dxlsx)
   - [analysing `2P_PD 90-100d.xlsx` (an alternate approach)](#analysing-2p_pd-90-100dxlsx-an-alternate-approach)
- [Contact](#Contact)

## Step-by-Step Guide

The following code will analyse the datasets contained in the EXCEL spreadsheet `2P_PD 90-100d.xlsx`

#### analysing `2P_PD 90-100d.xlsx`

**Run this code to analyse `2P_PD 90-100d.xlsx`**

   ```R
   # Remove all objects from the environment
   rm(list = ls(all = TRUE))

   # remove any graphs
   graphics.off()

   # Load and install necessary packages
   load_required_packages <- function(packages) {
    new.packages <- packages[!(packages %in% installed.packages()[, 'Package'])]
    if (length(new.packages)) install.packages(new.packages)
    invisible(lapply(packages, library, character.only = TRUE))
   }

   required.packages <- c('robustbase', 'minpack.lm', 'Rcpp', 'signal', 'writexl')
   load_required_packages(required.packages)

   # User defined values


   # Load and install necessary packages
   load_required_packages <- function(packages) {
    new.packages <- packages[!(packages %in% installed.packages()[, 'Package'])]
    if (length(new.packages)) install.packages(new.packages)
    invisible(lapply(packages, library, character.only = TRUE))
   }

   required.packages <- c('robustbase', 'minpack.lm', 'Rcpp', 'signal', 'writexl')
   load_required_packages(required.packages)

   # User defined values
   UserName <- 'YourUserName' # substitute your UserName here
   root_dir <- paste0('/Users/', UserName, '/Documents/Repositories/Rfits') # substitute your directory structure here
   
   # load custom functions which must be contained in the root directory root_dir
   path <- file.path(root_dir, 'nNLS functions.R')
   source(path)
   
   # create and change working directory wd
   wd <- paste0('/Users/', UserName, '/Documents/Repositories/Rfits/Yichen analysis') # substitute your directory structure here
   setwd(wd)

   # name of EXCEL xlsx dataset
   name <- '2P_PD 90-100d'

   # use load_data2 to load xlsx sheet file with multiple sheets
   example_data <- load_data2(wd=wd, name=name)

   # using names(example_data) will show the names of the sheets; in this case only 'Sheet1'
   # if multiple sheets then remember to analyse all; in this eg analysing first one

   # create data_sheet to for the first sheet
   data_sheet <- example_data$Sheet1

   # if required can easily plot using next 3 lines of code:
   # y <-  data_sheet[,1]; dt <- 1/40 * 1000
   # x <- seq(0, length(y)-1) * dt
   # plot(x, y, type='l')

   # metadata/settings for analysis
   stimulation_time <- 20000                 # time of stimulation
   baseline <- 15000                         # required baseline period
   fps <- 9.757                              # collection rate in frames per second
   dt <- 1/fps * 1e3                         # sample rate in ms
   func <- product1N                         # function used to model the response: product1N models 1 product eqn (with N reponses in a train)
   latency.limit <- NULL                     # default is NULL; will constrain the latency of the delay from the stimulation time
   method <- 'BF.LM'                         # default fitting algorithm
   ylab <- ''                                # chose appopriate y axis label eg ylab <- 'PSC amplitude (pA)'
   n <- 30                                   # number of performed fits; the function will chose the best fit by comparing an appropriate 'goodness of fit' 

   # use this code to produce simple plot of the data
   # y <-  data_sheet[,4]
   # x <- seq(0, length(y)-1) * dt
   # plot(x, y, type='l')

   # use this code to determine fit.limits for each (if necessary)
   # analyse_PSC(response=data_sheet[,1], dt=dt, n=n, func=func, ylab=ylab, stimulation_time=stimulation_time, baseline=baseline, method=method, return.output=FALSE) 

   # analyse all
   out1  <- analyse_PSC(response=data_sheet[,1], dt=dt, n=n, func=func, ylab=ylab, stimulation_time=stimulation_time, baseline=baseline, method=method, fit.limits=6e5) 

   out2  <- analyse_PSC(response=data_sheet[,2], dt=dt, n=n, func=func, ylab=ylab, stimulation_time=stimulation_time, baseline=baseline, method=method, fit.limits=1.5e5)

   out3  <- analyse_PSC(response=data_sheet[,3], dt=dt, n=n, func=func, ylab=ylab, stimulation_time=stimulation_time, baseline=baseline, method=method, fit.limits=4e4) 

   out4  <- analyse_PSC(response=data_sheet[,4], dt=dt, n=n, func=func, ylab=ylab, stimulation_time=stimulation_time, baseline=baseline, method=method, fit.limits=4e4) 

   out5  <- analyse_PSC(response=data_sheet[,5], dt=dt, n=n, func=func, ylab=ylab, stimulation_time=stimulation_time, baseline=baseline, method=method, fit.limits=5.5e4) 

   out6  <- analyse_PSC(response=data_sheet[,6], dt=dt, n=n, func=func, ylab=ylab, stimulation_time=stimulation_time, baseline=baseline, method=method, fit.limits=6e4) 
   ```

   The R output will look like this:

   ```
      out1  <- analyse_PSC(response=data_sheet[,1], dt=dt, n=n, func=func, ylab=ylab, stimulation_time=stimulation_time, baseline=baseline, method=method, fit.limits=6e5) 
        A1    τrise   τdecay    tpeak   r10_90   d90_10 delay    area1
   1 0.236 35905.15 357171.6 91704.94 48305.57 791210.7     0 108988.2 
     
   out2  <- analyse_PSC(response=data_sheet[,2], dt=dt, n=n, func=func, ylab=ylab, stimulation_time=stimulation_time, baseline=baseline, method=method, fit.limits=1.5e5)
     A1  τrise   τdecay   tpeak  r10_90   d90_10   delay    area1
   1 0.081 71.836 59540.49 483.322 153.958 130823.8 284.223 4888.096
     
   out3  <- analyse_PSC(response=data_sheet[,3], dt=dt, n=n, func=func, ylab=ylab, stimulation_time=stimulation_time, baseline=baseline, method=method, fit.limits=4e4) 
     A1   τrise   τdecay    tpeak  r10_90   d90_10   delay    area1
   1 0.458 444.553 8843.184 1399.726 697.667 19454.63 478.431 4749.758
     
      out4  <- analyse_PSC(response=data_sheet[,4], dt=dt, n=n, func=func, ylab=ylab, stimulation_time=stimulation_time, baseline=baseline, method=method, fit.limits=4e4) 
        A1   τrise   τdecay    tpeak  r10_90   d90_10   delay   area1
   1 0.59 612.806 8532.573 1738.762 893.226 18814.52 229.674 6174.95
     
   out5  <- analyse_PSC(response=data_sheet[,5], dt=dt, n=n, func=func, ylab=ylab, stimulation_time=stimulation_time, baseline=baseline, method=method, fit.limits=5.5e4) 
     A1   τrise   τdecay    tpeak  r10_90   d90_10   delay    area1
   1 0.697 321.125 14049.55 1241.757 572.083 30871.34 380.769 10694.72
    
   out6  <- analyse_PSC(response=data_sheet[,6], dt=dt, n=n, func=func, ylab=ylab, stimulation_time=stimulation_time, baseline=baseline, method=method, fit.limits=6e4) 
     A1   τrise   τdecay   tpeak r10_90   d90_10   delay    area1
   1 0.875 216.931 8700.299 821.285 381.95 19117.83 279.996 8368.254
   ```

   <div style="display: flex; justify-content: space-between;">
       <img src="../images/Yichen%20analysis/trace1.svg" alt="trace 1" style="width: 30%;"/>
       <img src="../images/Yichen%20analysis/trace2.svg" alt="trace 2" style="width: 30%;"/>
       <img src="../images/Yichen%20analysis/trace3.svg" alt="trace 3" style="width: 30%;"/>
   </div>

  <div style="display: flex; justify-content: space-between;">
       <img src="../images/Yichen%20analysis/trace4.svg" alt="trace 4" style="width: 30%;"/>
       <img src="../images/Yichen%20analysis/trace5.svg" alt="trace 5" style="width: 30%;"/>
       <img src="../images/Yichen%20analysis/trace6.svg" alt="trace 6" style="width: 30%;"/>
   </div>


   ```R   
   # plot smooth plots with upsampled fits on same x and y axes; upsampling can be increased usinf upsample.fit=c(upsample=TRUE, factor=100) 
   # simply increase factor; higher upsampling may slow graphing process
   
   # set axes limits and labels for the graph
   ylim <- c(-0.1, 1)
   xlab <- 'time (ms)'; ylab=''

   save <- FALSE # if FALSE graphs will be shown; if TRUE graphs are NOT shown but are saved as `svg` in the current working directory

   smooth.plots(y=data_sheet[,1], fits=out1$fits, dt=dt,  stimulation_time=stimulation_time, baseline=baseline, func=func, 
     upsample.fit=c(upsample=TRUE, factor=100), xlab=xlab, ylab=ylab, xlim=c(0, 600000), ylim=ylim, filename='smoothtrace1.svg', save=save)

   smooth.plots(y=data_sheet[,2], fits=out2$fits, dt=dt,  stimulation_time=stimulation_time, baseline=baseline, func=func, 
     upsample.fit=c(upsample=TRUE, factor=100), xlab=xlab, ylab=ylab, xlim=c(0, 250000), ylim=ylim, filename='smoothtrace2.svg', save=save)

   smooth.plots(y=data_sheet[,3], fits=out3$fits, dt=dt,  stimulation_time=stimulation_time, baseline=baseline, func=func, 
     upsample.fit=c(upsample=TRUE, factor=100), xlab=xlab, ylab=ylab, xlim=c(0, 150000), ylim=ylim, filename='smoothtrace3.svg', save=save)
   
   smooth.plots(y=data_sheet[,4], fits=out4$fits, dt=dt,  stimulation_time=stimulation_time, baseline=baseline, func=func, 
     upsample.fit=c(upsample=TRUE, factor=100), xlab=xlab, ylab=ylab, xlim=c(0, 150000), ylim=ylim, filename='smoothtrace4.svg', save=save)

   smooth.plots(y=data_sheet[,5], fits=out5$fits, dt=dt,  stimulation_time=stimulation_time, baseline=baseline, func=func, 
     upsample.fit=c(upsample=TRUE, factor=100), xlab=xlab, ylab=ylab, xlim=c(0, 150000), ylim=ylim, filename='smoothtrace5.svg', save=save)

   smooth.plots(y=data_sheet[,6], fits=out6$fits, dt=dt,  stimulation_time=stimulation_time, baseline=baseline, func=func,
     upsample.fit=c(upsample=TRUE, factor=100), xlab=xlab, ylab=ylab, xlim=c(0, 120000), ylim=ylim, filename='smoothtrace6.svg', save=save)


   ```
   
   If using smooth.plots and save=TRUE no graphical output will appear but graphs will be saved in svg format in the current working directory (i.e. wd as defined above):

   <div style="display: flex; justify-content: space-between;">
       <img src="../images/Yichen%20analysis/smoothtrace1.svg" alt="smooth trace 1" style="width: 30%;"/>
       <img src="../images/Yichen%20analysis/smoothtrace2.svg" alt="smooth trace 2" style="width: 30%;"/>
       <img src="../images/Yichen%20analysis/smoothtrace3.svg" alt="smooth trace 3" style="width: 30%;"/>
   </div>

   <div style="display: flex; justify-content: space-between;">
       <img src="../images/Yichen%20analysis/smoothtrace4.svg" alt="smooth trace 4" style="width: 30%;"/>
       <img src="../images/Yichen%20analysis/smoothtrace5.svg" alt="smooth trace 5" style="width: 30%;"/>
       <img src="../images/Yichen%20analysis/smoothtrace6.svg" alt="smooth trace 6" style="width: 30%;"/>
   </div>

   The fits have been upsampled 100-fold and appear smooth in comparison to the original fits (which were plotted at the same sample rate as the original response).


#### analysing `2P_PD 90-100d.xlsx` (an alternate approach)

**Run this code to analyse `2P_PD 90-100d.xlsx`**

   In this approach, rather than store the data from each fit as out1, out2 etc, simply store all the fits in a single `list`:

   ```R
   # Remove all objects from the environment
   rm(list = ls(all = TRUE))

   # remove any graphs
   graphics.off()

   # Load and install necessary packages
   load_required_packages <- function(packages) {
    new.packages <- packages[!(packages %in% installed.packages()[, 'Package'])]
    if (length(new.packages)) install.packages(new.packages)
    invisible(lapply(packages, library, character.only = TRUE))
   }

   required.packages <- c('robustbase', 'minpack.lm', 'Rcpp', 'signal', 'writexl')
   load_required_packages(required.packages)

   # User defined values

   # Load and install necessary packages
   load_required_packages <- function(packages) {
    new.packages <- packages[!(packages %in% installed.packages()[, 'Package'])]
    if (length(new.packages)) install.packages(new.packages)
    invisible(lapply(packages, library, character.only = TRUE))
   }

   required.packages <- c('robustbase', 'minpack.lm', 'Rcpp', 'signal', 'writexl')
   load_required_packages(required.packages)

   # User defined values
   UserName <- 'YourUserName' # substitute your UserName here
   root_dir <- paste0('/Users/', UserName, '/Documents/Repositories/Rfits') # substitute your directory structure here
   
   # load custom functions which must be contained in the root directory root_dir
   path <- file.path(root_dir, 'nNLS functions.R')
   source(path)
   
   # create and change working directory wd
   wd <- paste0('/Users/', UserName, '/Documents/Repositories/Rfits/Yichen analysis') # substitute your directory structure here
   setwd(wd)

   # name of EXCEL xlsx dataset
   name <- '2P_PD 90-100d'

   # use load_data2 to load xlsx sheet file with multiple sheets
   example_data <- load_data2(wd=wd, name=name)

   # using names(example_data) will show the names of the sheets; in this case only 'Sheet1'
   # if multiple sheets then remember to analyse all; in this eg analysing first one

   # create data_sheet to for the first sheet
   data_sheet <- example_data$Sheet1

   # if required can easily plot using next 3 lines of code:
   # y <-  data_sheet[,1]; dt <- 1/40 * 1000
   # x <- seq(0, length(y)-1) * dt
   # plot(x, y, type='l')

   # metadata/settings for analysis
   stimulation_time <- 20000                 # time of stimulation
   baseline <- 15000                         # required baseline period
   fps <- 9.757                              # collection rate in frames per second
   dt <- 1/fps * 1e3                         # sample rate in ms
   func <- product1N                         # function used to model the response: product1N models 1 product eqn (with N reponses in a train)
   latency.limit <- NULL                     # default is NULL; will constrain the latency of the delay from the stimulation time
   method <- 'BF.LM'                         # default fitting algorithm
   ylab <- ''                                # chose appopriate y axis label eg ylab <- 'PSC amplitude (pA)'
   n <- 30                                   # number of performed fits; the function will chose the best fit by comparing an appropriate 'goodness of fit' 

   # use this code to produce simple plot of the data
   # y <-  data_sheet[,4]
   # x <- seq(0, length(y)-1) * dt
   # plot(x, y, type='l')

   # use this code to determine fit.limits for each (if necessary)
   # analyse_PSC(response=data_sheet[,1], dt=dt, n=n, func=func, ylab=ylab, stimulation_time=stimulation_time, baseline=baseline, method=method, return.output=FALSE) 

   time.limits <- c(6e5, 1.5e5, 4e4, 4e4, 5.5e4, 6e4)

   # Create an empty list to store results
   out_list <- list()

   # Loop over the columns and store the results in the list
   for (ii in 1:ncol(data_sheet)) {
     out_list[[ii]] <- analyse_PSC(response=data_sheet[,ii], dt=dt, n=n, func=func, ylab=ylab, 
        stimulation_time=stimulation_time, baseline=baseline, method=method, fit.limits=time.limits[ii])  
   }

   # Create a summary table of the results
   # Loop over the list 'out_list', extracting the 'output' element from each result
   # 't()' is used to transpose the extracted results into rows 
   # 'as.vector' flattens the matrix to ensure all values are in a single vector per row
   summary <- t(sapply(1:length(out_list), function(ii){
     X <- out_list[[ii]]$output
     as.vector(t(X))
   }))

   # Create new column names by appending 1 and 2 to the original names
   # The 'rep()' function duplicates each column name from the first output twice
   # This is because the summary now has two rows for each original column
   new_colnames <- rep(colnames(out_list[[1]]$output), 1)

   # Assign the new column names to the summary table
   colnames(summary) <- new_colnames

   # Set row names as the index from 1 to the length of out_list, representing the row numbers (or the index of the analysis results)
   rownames(summary) <- 1:length(out_list)

   # 'summary' now holds the flattened results from all fits contained in 'out_list'
   summary
   ```

   `summary` now contains the main output in a single matrix.

   ```
        A1     τrise     τdecay     tpeak    r10_90    d90_10   delay      area1
   1 0.236 35905.150 357171.636 91704.944 48305.572 791210.75   0.000 108988.239
   2 0.081    71.836  59540.488   483.322   153.958 130823.82 284.223   4888.096
   3 0.458   444.553   8843.184  1399.726   697.667  19454.63 478.431   4749.758
   4 0.590   612.806   8532.573  1738.762   893.226  18814.52 229.674   6174.950
   5 0.697   321.125  14049.549  1241.757   572.083  30871.34 380.769  10694.720
   6 0.875   216.931   8700.299   821.285   381.950  19117.83 279.996   8368.254
   ```

   **Saving the the `summary` output to a `XLSX` file**

   ```R
   # Change the working directory to 'wd' before saving (this ensures saved environment is within the working directory)
   setwd(wd)
   # Save the summary output to XLSX file (nb as.data.frame(summary) conversion to data frame is necessary)
   write.xlsx(as.data.frame(summary), 'summary.xlsx')
   ```
   The working directory (use `getwd()` to see the current working directory) now contains an excel spreadsheet named `summary.xlsx` 

   **Saving the objects and variables within an image file**

   ```R
   # Change the working directory to 'wd' before saving (this ensures saved environment is within the working directory)
   setwd(wd)
   # Save the entire R environment including all objects and variables
   # This will save all objects created during this session and any other data or functions in the environment
   save.image(file = 'example.RData')  
   ```


   **Starting a new session from an image file**

   Having saved the current environment to an image file, `RData` the session can be resumed in a new instance of `R`.

   Note that any required packages will still need to be loaded into the opened environemnt.

   ```R
   # Remove all objects from the environment
   rm(list = ls(all = TRUE))

   UserName <- 'YourUserName' # substitute your UserName here
   
   # create and change working directory wd
   wd <- paste0('/Users/', UserName, '/Documents/Repositories/Rfits/Yichen analysis') # substitute your directory structure here
   setwd(wd)

   # path to data
   path <- file.path(wd, 'example.RData')

   # load the data
   load(path)

   # This should load the entire environment which will include the custom functions contained in 'nNLS functions.R'
   # However, any required packages must still be loaded
   # To load and install necessary packages, run the following code as before:
   
   load_required_packages <- function(packages) {
    new.packages <- packages[!(packages %in% installed.packages()[, 'Package'])]
    if (length(new.packages)) install.packages(new.packages)
    invisible(lapply(packages, library, character.only = TRUE))
   }

   required.packages <- c('robustbase', 'minpack.lm', 'Rcpp', 'signal', 'writexl')
   load_required_packages(required.packages)   
   ```

   This session is now restored to the saved point and can be used for further analysis


## Contact

If any bug fixes are necessary (most likely related to providing help on other operating systems), it will be provided as an update on the parent [`GitHub` page](https://github.com/vernonclarke/Rfits).

For queries related to this repository, please [open an issue](https://github.com/vernonclarke/Rfits/issues) or [email](mailto:WOPR2@proton.me) directly 


