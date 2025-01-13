# Remove all objects from the environment
rm(list = ls(all = TRUE))

# Load and install necessary packages
load_required_packages <- function(packages) {
  new.packages <- packages[!(packages %in% installed.packages()[, 'Package'])]
  if (length(new.packages)) install.packages(new.packages)
  invisible(lapply(packages, library, character.only = TRUE))
}
required.packages <- c('openxlsx', 'robustbase', 'minpack.lm', 'Rcpp', 'signal')
load_required_packages(required.packages)

source('/Users/euo9382/Documents/Repositories/Rfits//nNLS functions.R')

wd <- '/Users/euo9382/Documents/Repositories/Rfits/Yichen analysis'
name <- '2P_PD 90-100d'

setwd(wd)


# use load_data2 to load xlsx sheet file with multiple sheets
example_data <- load_data2(wd=wd, name=name)

# if multiple sheets then remember to analyse all

# in this eg analysing first one

data_sheet <- example_data$Sheet1

fps <- 9.757 
dt <- 1/fps * 1e3

y <-  data_sheet[,4]
x <- seq(0, length(y)-1) * dt
plot(x, y, type='l')


# metadata/settings for analysis
stimulation_time <- 20000
baseline <- 15000


func <- product1N
latency.limit <- NULL
method <- 'BF.LM'
ylab=''
n <- 30


# use this code to determine fit.limits for each (if necessary)
# analyse_PSC(response=data_sheet[,1], dt=dt, n=n, func=func, ylab=ylab, stimulation_time=stimulation_time, baseline=baseline, method=method, return.output=FALSE) 

# analyse all
out1  <- analyse_PSC(response=data_sheet[,1], dt=dt, n=n, func=func, ylab=ylab, stimulation_time=stimulation_time, baseline=baseline, method=method, fit.limits=6e5) 

out2  <- analyse_PSC(response=data_sheet[,2], dt=dt, n=n, func=func, ylab=ylab, stimulation_time=stimulation_time, baseline=baseline, method=method, fit.limits=6e5)

out3  <- analyse_PSC(response=data_sheet[,3], dt=dt, n=n, func=func, ylab=ylab, stimulation_time=stimulation_time, baseline=baseline, method=method, fit.limits=4e4) 

out4  <- analyse_PSC(response=data_sheet[,4], dt=dt, n=n, func=func, ylab=ylab, stimulation_time=stimulation_time, baseline=baseline, method=method, fit.limits=4e4) 

out5  <- analyse_PSC(response=data_sheet[,5], dt=dt, n=n, func=func, ylab=ylab, stimulation_time=stimulation_time, baseline=baseline, method=method, fit.limits=5e4) 

   

# plot smooth plots with upsampled fits on same x and y axes

xlim=c(0, 2500); ylim=c(-0.05, 2)
xlab='time (ms)'; ylab=''

smooth.plots(y=data_sheet[,1], fits=out1$fits, dt=dt,  stimulation_time=stimulation_time, baseline=baseline, func=func, 
  upsample.fit=c(upsample=TRUE, factor=100), xlab=xlab, ylab=ylab, xlim=xlim, ylim=ylim, filename='trace.svg1', save=FALSE)

smooth.plots(y=data_sheet[,2], fits=out2$fits, dt=dt,  stimulation_time=stimulation_time, baseline=baseline, func=func, 
  upsample.fit=c(upsample=TRUE, factor=100), xlab=xlab, ylab=ylab, xlim=xlim, ylim=ylim, filename='trace.svg2', save=FALSE)

smooth.plots(y=data_sheet[,3], fits=out3$fits, dt=dt,  stimulation_time=stimulation_time, baseline=baseline, func=func, 
  upsample.fit=c(upsample=TRUE, factor=100), xlab=xlab, ylab=ylab, xlim=xlim, ylim=ylim, filename='trace.svg3', save=FALSE)

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
required.packages <- c('openxlsx', 'robustbase', 'minpack.lm', 'Rcpp', 'signal')
load_required_packages(required.packages)

source('/Users/euo9382/Documents/Repositories/Rfits//nNLS functions.R')

wd <- '/Users/euo9382/Documents/Repositories/Rfits/James analysis'
name <- 'KineticsDataset_5P(5Hz)_VC'

setwd(wd)


# use load_data2 to load xlsx sheet file with multiple sheets
example_data <- load_data2(wd=wd, name=name)

# load 1st data sheet (in theory you could put data into mutiple excel spreadsheets

data_sheet <- example_data$Sheet1

# y <-  data_sheet[,1]; dt <- 1/40 * 1000 # ms
# x <- seq(0, length(y)-1) * dt
# plot(x, y, type='l')


# metadata/settings for analysis
stimulation_time <- 2950
baseline <- 500
dt <- 1/40 * 1000 # ms
func <- product1N
latency.limit <- NULL
method <- 'LM'
ylab=''
n <- 30

N <- 5
IEI <- 200

# use this code to determine fit.limits for each (if necessary)
# analyse_PSC(response=data_sheet[,1], dt=dt, n=n, func=func, ylab=ylab, stimulation_time=stimulation_time, baseline=baseline, method=method, return.output=FALSE) 

# analyse all
out1  <- analyse_PSC(response=data_sheet[,1], dt=dt, n=n, N=N, IEI=IEI, func=func, ylab=ylab, stimulation_time=stimulation_time, baseline=baseline, method=method, fit.limits=2000) 

out2  <- analyse_PSC(response=data_sheet[,2], dt=dt, n=n, N=N, IEI=IEI, func=func, ylab=ylab, stimulation_time=stimulation_time, baseline=baseline, method=method, fit.limits=2000) 

out3  <- analyse_PSC(response=data_sheet[,3], dt=dt, n=n, N=N, IEI=IEI, func=func, ylab=ylab, stimulation_time=stimulation_time, baseline=baseline, method=method, fit.limits=2000) 

out4  <- analyse_PSC(response=data_sheet[,4], dt=dt, n=n, N=N, IEI=IEI, func=func, ylab=ylab, stimulation_time=stimulation_time, baseline=baseline, method=method, fit.limits=2000) 

out5  <- analyse_PSC(response=data_sheet[,5], dt=dt, n=n, N=N, IEI=IEI, func=func, ylab=ylab, stimulation_time=stimulation_time, baseline=baseline, method=method, fit.limits=2000) 


# plot smooth plots with upsampled fits on same x and y axes

xlim=c(0, 2500); ylim=c(-0.05, 2)
xlab='time (ms)'; ylab=''

smooth.plots(y=data_sheet[,1], fits=out1$fits, N=N, IEI=IEI, dt=dt,  stimulation_time=stimulation_time, baseline=baseline, func=func, 
  upsample.fit=c(upsample=TRUE, factor=100), xlab=xlab, ylab=ylab, xlim=xlim, ylim=ylim, filename='trace.svg1', save=FALSE)

smooth.plots(y=data_sheet[,2], fits=out2$fits, N=N, IEI=IEI, dt=dt,  stimulation_time=stimulation_time, baseline=baseline, func=func, 
  upsample.fit=c(upsample=TRUE, factor=100), xlab=xlab, ylab=ylab, xlim=xlim, ylim=ylim, filename='trace.svg2', save=FALSE)

smooth.plots(y=data_sheet[,3], fits=out3$fits, N=N, IEI=IEI, dt=dt,  stimulation_time=stimulation_time, baseline=baseline, func=func, 
  upsample.fit=c(upsample=TRUE, factor=100), xlab=xlab, ylab=ylab, xlim=xlim, ylim=ylim, filename='trace.svg3', save=FALSE)

smooth.plots(y=data_sheet[,4], fits=out4$fits, N=N, IEI=IEI, dt=dt,  stimulation_time=stimulation_time, baseline=baseline, func=func, 
  upsample.fit=c(upsample=TRUE, factor=100), xlab=xlab, ylab=ylab, xlim=xlim, ylim=ylim, filename='trace.svg3', save=FALSE)

smooth.plots(y=data_sheet[,5], fits=out5$fits, N=N, IEI=IEI, dt=dt,  stimulation_time=stimulation_time, baseline=baseline, func=func, 
  upsample.fit=c(upsample=TRUE, factor=100), xlab=xlab, ylab=ylab, xlim=xlim, ylim=ylim, filename='trace.svg3', save=FALSE)

# to save (and not show) change save to TRUE



