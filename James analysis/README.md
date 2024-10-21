# <center>`R` Code guide</center>


## Table of Contents
- [Step-by-Step Guide](#Step-by-Step-Guide)
   - [analysing `KineticsDataset_Single.xlsx`](#analysing-KineticsDataset_Singlexlsx)
   - [analysing `KineticsDataset_5P(5Hz).xlsx`](#analysing-KineticsDataset_5P(5Hz)xlsx)
- [Contact](#Contact)


## Step-by-Step Guide

The following code will analyse the datasets contained in the EXCEL spreadsheets `KineticsDataset_Single.xlsx` and `KineticsDataset_5P(5Hz).xlsx`

#### analysing KineticsDataset_Single.xlsx

**Run this code to analyse `KineticsDataset_Single.xlsx`**

   ```R
   # Remove all objects from the environment
   rm(list = ls(all = TRUE))

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
   wd <- paste0('/Users/', UserName, '/Documents/Repositories/Rfits/James analysis') # substitute your directory structure here
   setwd(wd)

   # name of EXCEL xlsx dataset
   name <- 'KineticsDataset_Single'

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
   stimulation_time <- 2950                  # time of stimulation
   baseline <- 500                           # required baseline period
   dt <- 1/40 * 1000                         # sampling interval in ms (only)
   func <- product1N                         # function used to model the response: product1N models 1 product eqn (with N reponses in a train)
   latency.limit <- NULL                     # default is NULL; will constrain the latency of the delay from the stimulation time
   method <- 'BF.LM'                         # default fitting algorithm
   ylab <- ''                                # chose appopriate y axis label eg ylab <- 'PSC amplitude (pA)'
   n <- 30                                   # number of performed fits; the function will chose the best fit by comparing an appropriate 'goodness of fit' 


   # use this code to determine fit.limits for each (if necessary)
   # analyse_PSC(response=data_sheet[,1], dt=dt, n=n, func=func, ylab=ylab, stimulation_time=stimulation_time, baseline=baseline, method=method, return.output=FALSE) 

   # analyse all
   out1  <- analyse_PSC(response=data_sheet[,1], dt=dt, n=n, func=func, ylab=ylab, stimulation_time=stimulation_time, baseline=baseline, method=method, fit.limits=1000) 

   out2  <- analyse_PSC(response=data_sheet[,2], dt=dt, n=n, func=func, ylab=ylab, stimulation_time=stimulation_time, baseline=baseline, method=method, fit.limits=2000) 

   out3  <- analyse_PSC(response=data_sheet[,3], dt=dt, n=n, func=func, ylab=ylab, stimulation_time=stimulation_time, baseline=baseline, method=method, fit.limits=1250) 

   out4  <- analyse_PSC(response=data_sheet[,4], dt=dt, n=n, func=func, ylab=ylab, stimulation_time=stimulation_time, baseline=baseline, method=method, fit.limits=1500)

   out5  <- analyse_PSC(response=data_sheet[,5], dt=dt, n=n, func=func, ylab=ylab, stimulation_time=stimulation_time, baseline=baseline, method=method, fit.limits=1500) 
 


   ```

   The R output will look like this:

   ```
   out1  <- analyse_PSC(response=data_sheet[,1], dt=dt, n=n, func=func, ylab=ylab, stimulation_time=stimulation_time, baseline=baseline, method=method, fit.limits=1000) 
        A1  τrise  τdecay  tpeak r10_90  d90_10  delay   area1
   1 0.666 22.551 171.618 52.691 28.246 382.665 23.132 155.311
 
   out2  <- analyse_PSC(response=data_sheet[,2], dt=dt, n=n, func=func, ylab=ylab, stimulation_time=stimulation_time, baseline=baseline, method=method, fit.limits=2000) 
        A1  τrise  τdecay   tpeak r10_90  d90_10  delay   area1
   1 1.972 57.459 271.347 113.157  62.24 619.036 19.586 812.116
   out3  <- analyse_PSC(response=data_sheet[,3], dt=dt, n=n, func=func, ylab=ylab, stimulation_time=stimulation_time, baseline=baseline, method=method, fit.limits=1250) 
        A1  τrise  τdecay  tpeak r10_90  d90_10  delay   area1
   1 1.142 35.251 219.562 76.812 41.658 493.207 23.263 355.866

   out4  <- analyse_PSC(response=data_sheet[,4], dt=dt, n=n, func=func, ylab=ylab, stimulation_time=stimulation_time, baseline=baseline, method=method, fit.limits=1500)
        A1  τrise  τdecay   tpeak r10_90  d90_10  delay   area1
   1 2.088 58.931 274.842 115.511 63.572 627.564 20.717 873.671
 
   out5  <- analyse_PSC(response=data_sheet[,5], dt=dt, n=n, func=func, ylab=ylab, stimulation_time=stimulation_time, baseline=baseline, method=method, fit.limits=1500) 
        A1  τrise  τdecay  tpeak r10_90  d90_10 delay   area1
   1 1.818 31.051 184.339 66.509 36.165 414.993 20.08 480.834
   ```

   <div style="display: flex; justify-content: space-between;">
       <img src="../images/James%20analysis/trace1.svg" alt="trace 1" style="width: 30%;"/>
       <img src="../images/James%20analysis/trace2.svg" alt="trace 2" style="width: 30%;"/>
       <img src="../images/James%20analysis/trace3.svg" alt="trace 3" style="width: 30%;"/>
   </div>

  <div style="display: flex; justify-content: space-between;">
       <img src="../images/James%20analysis/trace4.svg" alt="trace 4" style="width: 30%;"/>
       <img src="../images/James%20analysis/trace5.svg" alt="trace 5" style="width: 30%;"/>
   </div>


   ```R   
   # plot smooth plots with upsampled fits on same x and y axes; upsampling can be increased usinf upsample.fit=c(upsample=TRUE, factor=100) 
   # simply increase factor; higher upsampling may slow graphing process
   
   # set axes limits and labels for the graph
   xlim=c(0, 2500); ylim=c(-0.05, 2)
   xlab='time (ms)'; ylab=''

   smooth.plots(y=data_sheet[,1], fits=out1$fits, dt=dt,  stimulation_time=stimulation_time, baseline=baseline, func=func, 
     upsample.fit=c(upsample=TRUE, factor=100), xlab=xlab, ylab=ylab, xlim=xlim, ylim=ylim, filename='smoothtrace1.svg', save=FALSE)

   smooth.plots(y=data_sheet[,2], fits=out2$fits, dt=dt,  stimulation_time=stimulation_time, baseline=baseline, func=func, 
     upsample.fit=c(upsample=TRUE, factor=100), xlab=xlab, ylab=ylab, xlim=xlim, ylim=ylim, filename='smoothtrace2.svg', save=FALSE)

   smooth.plots(y=data_sheet[,3], fits=out3$fits, dt=dt,  stimulation_time=stimulation_time, baseline=baseline, func=func, 
     upsample.fit=c(upsample=TRUE, factor=100), xlab=xlab, ylab=ylab, xlim=xlim, ylim=ylim, filename='smoothtrace3.svg', save=FALSE)
   
   smooth.plots(y=data_sheet[,4], fits=out4$fits, dt=dt,  stimulation_time=stimulation_time, baseline=baseline, func=func, 
     upsample.fit=c(upsample=TRUE, factor=100), xlab=xlab, ylab=ylab, xlim=xlim, ylim=ylim, filename='smoothtrace4.svg', save=FALSE)

   smooth.plots(y=data_sheet[,5], fits=out5$fits, dt=dt,  stimulation_time=stimulation_time, baseline=baseline, func=func, 
     upsample.fit=c(upsample=TRUE, factor=100), xlab=xlab, ylab=ylab, xlim=xlim, ylim=ylim, filename='smoothtrace5.svg', save=FALSE)


   ```
   
   If using smooth.plots and save=TRUE no graphical output will appear but graphs will be saved in svg format in the current working directory (i.e. wd as defined above):

   <div style="display: flex; justify-content: space-between;">
       <img src="../images/James%20analysis/smoothtrace1.svg" alt="smooth trace 1" style="width: 30%;"/>
       <img src="../images/James%20analysis/smoothtrace2.svg" alt="smooth trace 2" style="width: 30%;"/>
       <img src="../images/James%20analysis/smoothtrace3.svg" alt="smooth trace 3" style="width: 30%;"/>
   </div>

   <div style="display: flex; justify-content: space-between;">
       <img src="../images/James%20analysis/smoothtrace4.svg" alt="smooth trace 4" style="width: 30%;"/>
       <img src="../images/James%20analysis/smoothtrace5.svg" alt="smooth trace 5" style="width: 30%;"/>
   </div>

   The fits have been upsampled 100-fold and appear smooth in comparison to the original fits (which were plotted at the same sample  rate as the original response).



**Run this code to analyse `KineticsDataset_Single.xlsx`; an alternate approach**

   In this approach, rather than store the data from each fit as out1, out2 etc, simply store all the fits in a single `list`:

   ```R
   # Remove all objects from the environment
   rm(list = ls(all = TRUE))

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
   wd <- paste0('/Users/', UserName, '/Documents/Repositories/Rfits/James analysis') # substitute your directory structure here
   setwd(wd)

   # name of EXCEL xlsx dataset
   name <- 'KineticsDataset_Single'

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
   stimulation_time <- 2950                  # time of stimulation
   baseline <- 500                           # required baseline period
   dt <- 1/40 * 1000                         # sampling interval in ms (only)
   func <- product1N                         # function used to model the response: product1N models 1 product eqn (with N reponses in a train)
   latency.limit <- NULL                     # default is NULL; will constrain the latency of the delay from the stimulation time
   method <- 'BF.LM'                         # default fitting algorithm
   ylab <- ''                                # chose appopriate y axis label eg ylab <- 'PSC amplitude (pA)'
   n <- 30                                   # number of performed fits; the function will chose the best fit by comparing an appropriate 'goodness of fit' 


   # use this code to determine fit.limits for each (if necessary)
   # analyse_PSC(response=data_sheet[,1], dt=dt, n=n, func=func, ylab=ylab, stimulation_time=stimulation_time, baseline=baseline, method=method, return.output=FALSE) 

   time.limits <- c(1000, 2000, 1250, 1500, 1500)

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
        A1  τrise  τdecay   tpeak r10_90  d90_10  delay   area1
   1 0.666 22.551 171.618  52.691 28.246 382.665 23.132 155.311
   2 1.972 57.459 271.347 113.157 62.240 619.036 19.586 812.116
   3 1.142 35.251 219.562  76.812 41.658 493.207 23.263 355.866
   4 2.088 58.931 274.842 115.511 63.572 627.564 20.717 873.671
   5 1.818 31.051 184.339  66.509 36.165 414.993 20.080 480.834
   ```

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
   wd <- paste0('/Users/', UserName, '/Documents/Repositories/Rfits/James analysis') # substitute your directory structure here
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

   This session is now restored to the saved point.

   
   #### analysing `KineticsDataset_5P(5Hz).xlsx`

   **Run this code to analyse `KineticsDataset_5P(5Hz).xlsx`**

   ```R
   # Remove all objects from the environment
   rm(list = ls(all = TRUE))

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
   wd <- paste0('/Users/', UserName, '/Documents/Repositories/Rfits/James analysis') # substitute your directory structure here
   setwd(wd)

   # name of EXCEL xlsx dataset
   name <- 'KineticsDataset_5P(5Hz)'

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
   stimulation_time <- 2950                  # time of stimulation
   baseline <- 500                           # required baseline period
   dt <- 1/40 * 1000                         # sampling interval in ms (only)
   func <- product1N                         # function used to model the response: product1N models 1 product eqn (with N reponses in a train)
   latency.limit <- NULL                     # default is NULL; will constrain the latency of the delay from the stimulation time
   method <- 'BF.LM'                         # default fitting algorithm
   ylab <- ''                                # chose appopriate y axis label eg ylab <- 'PSC amplitude (pA)'
   n <- 30                                   # number of performed fits; the function will chose the best fit by comparing an appropriate 'goodness of fit' 

   N <- 5
   IEI <- 200

   # use this code to determine fit.limits for each (if necessary)
   # analyse_PSC(response=data_sheet[,1], dt=dt, n=n, func=func, N=N, IEI=IEI, ylab=ylab, stimulation_time=stimulation_time, baseline=baseline, method=method, return.output=FALSE) 

   # use same time.limits for all i.e. 2000

   # Create an empty list to store results
   out_list <- list()

   # Loop over the columns and store the results in the list
   for (ii in 1:ncol(data_sheet)) {
     out_list[[ii]] <- analyse_PSC(response=data_sheet[,ii], dt=dt, n=n, N=N, IEI=IEI, func=func, ylab=ylab, 
      stimulation_time=stimulation_time, baseline=baseline, method=method, fit.limits=2000) 
   }

   # Change the working directory to 'wd' before saving (this ensures saved environment is within the working directory)
   setwd(wd)
   # Save the entire R environment including all objects and variables
   # This will save all objects created during this session and any other data or functions in the environment
   save.image(file = 'example1.RData')  
   ```   


   **Starting a new session from an image file**

   Having saved the current environment to an image file, `RData` the session can be resumed in a new instance of `R`.

   Note that any required packages will still need to be loaded into the opened environemnt.

   ```R
   # Remove all objects from the environment
   rm(list = ls(all = TRUE))

   UserName <- 'YourUserName' # substitute your UserName here
   
   # create and change working directory wd
   wd <- paste0('/Users/', UserName, '/Documents/Repositories/Rfits/James analysis') # substitute your directory structure here
   setwd(wd)

   # path to data
   path <- file.path(wd, 'example1.RData')

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


   ```R
   # plot smooth plots with upsampled fits on same x and y axes
   # set axes limits and labels for the graph
   xlim=c(0, 2500); ylim=c(-0.05, 2)
   xlab='time (ms)'; ylab=''

   # use a loop to plot all as smooth plots
   for (ii in 1:length(out_list)) {
      smooth.plots(y=data_sheet[,ii], fits=out_list[[ii]]$fits, N=N, IEI=IEI, dt=dt,  stimulation_time=stimulation_time, baseline=baseline, func=func, 
     upsample.fit=c(upsample=TRUE, factor=100), xlab=xlab, ylab=ylab, xlim=xlim, ylim=ylim, filename=paste0('smoothtrain', ii,  '.svg'), save=FALSE)
   }
    ```

  If using smooth.plots and save=TRUE no graphical output will appear but graphs will be saved in svg format in the current working directory (i.e. wd as defined above):

   <div style="display: flex; justify-content: space-between;">
       <img src="../images/James%20analysis/smoothtrain1.svg" alt="smooth train 1" style="width: 30%;"/>
       <img src="../images/James%20analysis/smoothtrain2.svg" alt="smooth train 2" style="width: 30%;"/>
       <img src="../images/James%20analysis/smoothtrain3.svg" alt="smooth train 3" style="width: 30%;"/>
   </div>

   <div style="display: flex; justify-content: space-between;">
       <img src="../images/James%20analysis/smoothtrain4.svg" alt="smooth train 4" style="width: 30%;"/>
       <img src="../images/James%20analysis/smoothtrain5.svg" alt="smooth train 5" style="width: 30%;"/>
   </div>


   ```R
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
   
   Entering `summary` gives:

   ```
        A1    A2    A3    A4    A5  τrise  τdecay  tpeak r10_90  d90_10  delay   area1   area2   area3   area4   area5
   1 0.447 0.174 0.139 0.160 0.172 12.488 195.044 36.670 18.663 429.664 23.004 105.229  40.886  32.714  37.591  40.542
   2 1.332 0.315 0.248 0.320 0.340 29.780 213.776 68.199 36.688 477.570 20.617 391.679  92.649  72.877  94.202 100.121
   3 0.810 0.384 0.257 0.232 0.240 24.825 205.589 59.689 31.828 457.304 21.729 222.622 105.606  70.730  63.706  65.900
   4 1.522 0.330 0.241 0.291 0.247 37.524 221.928 80.265 43.654 499.705 22.288 485.015 105.206  76.906  92.612  78.815
   5 1.888 0.830 0.670 0.687 0.684 31.750 220.724 71.908 38.756 493.650 18.079 577.106 253.821 204.838 210.061 209.027
   ```

   Box plot of data sets

   ```R
   # Calculate the number of rows in the 'summary' table (n), then create a new data frame 'A'
   # 's' represents a subjects of 1:n repeated twice (for the two categories in 'x')
   # 'x' represents two groups: the first n rows are labeled 1 for the fast component 
   # amplitude, and the second n rows are labeled 2 for the slow component amplitude
   # 'y' combines the first column of 'summary' and the 9th column (without names) 
   # into a single vector i.e the fast amd slow amplitudes
   n <- dim(summary)[1]
   A <- data.frame(
      's' = rep(1:n,5),
      'x' = c(rep(1,n), rep(2,n), rep(3,n), rep(4,n), rep(5,n)),
      'y' = c(summary[,1], summary[,2], summary[,3], summary[,4], summary[,5])
      )
   
   # show structure of A is standard 'longitudinal' format
   A

   BoxPlot(A, xrange=c(0.75, 5.25), yrange=c(0, 2), ylab='amplitude (units)', wid=0.3, cap=0.1, y_tick_interval=0.2, 
   width=8, height=5, tick_length=0.2, lwd=1.25, amount=0.1, p.cex=1)

   # use subset(A, x %in% c(1, 5)) to isolate the first and last amplitudes in the train

   ScatterPlot(subset(A, x %in% c(1, 2)), xlim=c(0,2), ylim=c(-2,2), x_tick_interval=0.5, y_tick_interval=1, 
         xlab=expression(A[fast] * ' ' * (pA)), ylab=expression(A[slow] * ' ' * (pA)), 
         lwd=1.25, p.cex=0.6, reg=TRUE, plot.CI=TRUE, width=5, height=5)

   ScatterPlot(subset(A, x %in% c(1, 2)), xlim=c(0,2), ylim=c(-2,2), x_tick_interval=0.5, y_tick_interval=1, 
         xlab=expression(A[fast] * ' ' * (pA)), ylab=expression(A[slow] * ' ' * (pA)), 
         lwd=1.25, p.cex=0.6, reg=TRUE, reg.method='Theil-Sen', plot.CI=TRUE, width=5, height=5)

   # STATISTICS

   # non-parametric equivalent of a one-way repeated measures ANOVA is the Friedman test
   result <- friedman.test(y ~ x | s, data = A)

   # Perform pairwise Wilcoxon signed-rank tests with Holm correction, specifying paired data
   pairwise_result <- pairwise.wilcox.test(A$y, A$x, paired = TRUE, p.adjust.method = 'holm')

   # Alternative method
   # Get all unique pairs of conditions
   condition_pairs <- combn(unique(A$x), 2, simplify = FALSE)

   # Perform pairwise Wilcoxon signed-rank tests with Holm correction (alternative method)
   p_values <- sapply(condition_pairs, function(pair) {
     # Subset data for the two conditions
     data1 <- A$y[A$x == pair[1]]
     data2 <- A$y[A$x == pair[2]]
     
     # Perform the paired Wilcoxon signed-rank test
     test_result <- wilcox.test(data1, data2, paired = TRUE)
     
     # Return the p-value
     test_result$p.value
   })

   # Apply Holm correction
   adjusted_p_values <- p.adjust(p_values, method = 'holm')

   # Combine results into a data frame
   pairwise_results2 <- data.frame(
     contrast = sapply(condition_pairs, function(pair) paste(pair, collapse = " vs ")),
     p = p_values,
     'p adj' = adjusted_p_values
   )

   ```

   ### `friedman.test`
   - **s** represents the subject (or block)
   - **x** represents the conditions (repeated measures)
   - **y** is the dependent variable
   - **y ~ x | s** the formula indicates that y is measured across different levels of x for each subject s

   ### `pairwise.wilcox.test` 
   paired Wilcoxon signed-rank test with a multiple comparisons correction
   - **A$y** represents the dependent variable, 
   - **A$x** represents the groups being compared.
   - **paired = TRUE** argument indicates that the data is paired (repeated measures for the same subjects). 
   - **method** can be set to 'holm', 'hochberg', 'hommel', 'bonferroni', 'BH', 'BY', 'fdr' or 'none'

   ### alternative method for pairwise Wilcoxon signed-rank tests with Holm correction

   - Create all unique pairs of conditions using `combn`.
   - For each pair, extract the corresponding values for the two conditions and perform a paired Wilcoxon signed-rank test.
   - Apply the Bonferroni correction to the resulting p-values.
   - Compiled results into a data frame showing the original and adjusted p-values for each pairwise comparison.

   This ensures that the tests account for the paired nature of the data based on the subject (s). 




## Contact

If any bug fixes are necessary (most likely related to providing help on other operating systems), it will be provided as an update on the parent [`GitHub` page](https://github.com/vernonclarke/Rfits).

For queries related to this repository, please [open an issue](https://github.com/vernonclarke/Rfits/issues) or [email](mailto:WOPR2@proton.me) directly 


