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

# set axes limits and labels for the graph
xlim=c(0, 2500); ylim=c(-0.05, 2)
xlab='time (ms)'; ylab=''

smooth.plots(y=data_sheet[,1], fits=out1$fits, dt=dt,  stimulation_time=stimulation_time, baseline=baseline, func=func, 
 upsample.fit=c(upsample=TRUE, factor=100), xlab=xlab, ylab=ylab, xlim=xlim, ylim=ylim, filename='smoothtrace1.svg', save=FALSE)

smooth.plots(y=data_sheet[,2], fits=out2$fits, dt=dt,  stimulation_time=stimulation_time, baseline=baseline, func=func, 
 upsample.fit=c(upsample=TRUE, factor=100), xlab=xlab, ylab=ylab, xlim=xlim, ylim=ylim, filename='smoothtrace2.svg', save=FALSE)

smooth.plots(y=data_sheet[,3], fits=out3$fits, dt=dt,  stimulation_time=stimulation_time, baseline=baseline, func=func, 
 upsample.fit=c(upsample=TRUE, factor=100), xlab=xlab, ylab=ylab, xlim=xlim, ylim=ylim, filename='smoothtrace3.svg', save=FALSE)
 
 # to save (and not show) change save to TRUE


### TRAINS

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
UserName <- 'euo9382' # substitute your UserName here
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


# analyse all
out1  <- analyse_PSC(response=data_sheet[,1], dt=dt, n=n, N=N, IEI=IEI, func=func, ylab=ylab, stimulation_time=stimulation_time, baseline=baseline, method=method, fit.limits=2000)

out2  <- analyse_PSC(response=data_sheet[,2], dt=dt, n=n, N=N, IEI=IEI, func=func, ylab=ylab, stimulation_time=stimulation_time, baseline=baseline, method=method, fit.limits=2000) 

out3  <- analyse_PSC(response=data_sheet[,3], dt=dt, n=n, N=N, IEI=IEI, func=func, ylab=ylab, stimulation_time=stimulation_time, baseline=baseline, method=method, fit.limits=2000) 

out4  <- analyse_PSC(response=data_sheet[,4], dt=dt, n=n, N=N, IEI=IEI, func=func, ylab=ylab, stimulation_time=stimulation_time, baseline=baseline, method=method, fit.limits=2000) 

out5  <- analyse_PSC(response=data_sheet[,5], dt=dt, n=n, N=N, IEI=IEI, func=func, ylab=ylab, stimulation_time=stimulation_time, baseline=baseline, method=method, fit.limits=2000)


# plot smooth plots with upsampled fits on same x and y axes
# set axes limits and labels for the graph
xlim=c(0, 2500); ylim=c(-0.05, 2)
xlab='time (ms)'; ylab=''

smooth.plots(y=data_sheet[,1], fits=out1$fits, N=N, IEI=IEI, dt=dt,  stimulation_time=stimulation_time, baseline=baseline, func=func, 
 upsample.fit=c(upsample=TRUE, factor=100), xlab=xlab, ylab=ylab, xlim=xlim, ylim=ylim, filename='smoothtrain1.svg', save=FALSE)

smooth.plots(y=data_sheet[,2], fits=out2$fits, N=N, IEI=IEI, dt=dt,  stimulation_time=stimulation_time, baseline=baseline, func=func, 
 upsample.fit=c(upsample=TRUE, factor=100), xlab=xlab, ylab=ylab, xlim=xlim, ylim=ylim, filename='smoothtrain2.svg', save=FALSE)

smooth.plots(y=data_sheet[,3], fits=out3$fits, N=N, IEI=IEI, dt=dt,  stimulation_time=stimulation_time, baseline=baseline, func=func, 
 upsample.fit=c(upsample=TRUE, factor=100), xlab=xlab, ylab=ylab, xlim=xlim, ylim=ylim, filename='smoothtrain3.svg', save=FALSE)

smooth.plots(y=data_sheet[,4], fits=out4$fits, N=N, IEI=IEI, dt=dt,  stimulation_time=stimulation_time, baseline=baseline, func=func, 
 upsample.fit=c(upsample=TRUE, factor=100), xlab=xlab, ylab=ylab, xlim=xlim, ylim=ylim, filename='smoothtrain4.svg', save=FALSE)

smooth.plots(y=data_sheet[,5], fits=out5$fits, N=N, IEI=IEI, dt=dt,  stimulation_time=stimulation_time, baseline=baseline, func=func, 
 upsample.fit=c(upsample=TRUE, factor=100), xlab=xlab, ylab=ylab, xlim=xlim, ylim=ylim, filename='smoothtrain5.svg', save=FALSE)

# to save (and not show) change save to TRUE



