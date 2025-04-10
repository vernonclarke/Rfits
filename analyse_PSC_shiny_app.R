# These steps 
# By manually creating the /tmp/.X11-unix directory with proper permissions, 
# ensuring no conflicting X server processes are running, restarting XQuartz, 
# and setting the DISPLAY variable, you’ve set up the correct environment 
# for your Tcl/Tk applications in R.

# #  1. Create the /tmp/.X11-unix directory manually with correct permissions
# sudo mkdir /tmp/.X11-unix
# sudo chmod 1777 /tmp/.X11-unix

# # 2. Verify that no other X server is running
# ps aux | grep X11
# # If you see any running X11-related processes, terminate them
# sudo killall Xquartz

# # 3. Restart XQuartz
# open -a XQuartz

# # 4. Set the DISPLAY environment variable
# export DISPLAY=:0

# # 5. Allow connections from localhost
# xhost +localhost

# # 6. Check XQuartz Security Settings (manually):
# # Go to XQuartz > Preferences > Security, and ensure "Allow connections from network clients" is checked

# # 7. Check the Console for XQuartz-related errors:
# # - Open the Console app (found via Spotlight search).
# # - Filter for "XQuartz" and look for any error messages that could provide additional clues.

# # 8. Test if the XQuartz configuration is working by running a simple application
# xterm


######################################################################################
# Final version
# for xquartz to work properly in (some) systems open R from terminal:
# open -n -a R

rm(list=ls(all=TRUE))
graphics.off()

# load and install necessary packages 
load_required_packages <- function(packages) {
  new.packages <- packages[!(packages %in% installed.packages()[, 'Package'])]
  if (length(new.packages)) install.packages(new.packages)
  invisible(lapply(packages, library, character.only=TRUE))
}

required.packages <- c('robustbase', 'minpack.lm', 'Rcpp', 'signal',
                       'dbscan', 'tkrplot', 'tcltk', 'readxl', 'shiny', 'readxl')
load_required_packages(required.packages)

# insert your username and repository path
username <- 'euo9382'
path_repository <- '/Documents/Repositories/Rfits'
file_path1 <- paste0('/Users/', username, path_repository)
source(paste0(file_path1, '/nNLS functions.R'))

PSC_analysis_tk <- function() {
  tt <- tktoplevel()
  tkwm.title(tt, 'PSC Analysis')
  
  # divide window into sidebar and main panels
  sidebarFrame <- tkframe(tt)
  mainFrame <- tkframe(tt)
  tkgrid(sidebarFrame, row=0, column=0, sticky='ns')
  tkgrid(mainFrame, row=0, column=1, sticky='nsew')
  tkgrid.rowconfigure(tt, 0, weight=0)
  tkgrid.columnconfigure(tt, 1, weight=1)
  
  # sidebar controls
  fileLabel <- tklabel(sidebarFrame, text='Upload CSV or XLSX:')
  tkgrid(fileLabel, row=0, column=0, sticky='w')
  filePathVar <- tclVar('')
  fileEntry <- tkentry(sidebarFrame, textvariable=filePathVar, width=30)
  tkgrid(fileEntry, row=0, column=1, sticky='w')
  browseButton <- tkbutton(sidebarFrame, text='Browse', command=function() {
    filePath <- tclvalue(tkgetOpenFile(filetypes='{{CSV Files} {.csv}} {{Excel Files} {.xlsx .xls}}'))
    if (nchar(filePath) > 0) {
      tclvalue(filePathVar) <- filePath
      ext <- tools::file_ext(filePath)
      if (tolower(ext) == 'csv') {
        uploaded_data <<- read.csv(filePath)
      } else {
        uploaded_data <<- readxl::read_excel(filePath)
      }
      columns <<- names(uploaded_data)
      tkconfigure(columnCombo, values=columns)
    }
  })
  tkgrid(browseButton, row=0, column=2, padx=5)
  
  colLabel <- tklabel(sidebarFrame, text='Select column:')
  tkgrid(colLabel, row=1, column=0, sticky='w')
  columnVar <- tclVar('')
  columnCombo <- ttkcombobox(sidebarFrame, textvariable=columnVar, values='', width=20)
  tkgrid(columnCombo, row=1, column=1, columnspan=2, sticky='w')
  
  # Notebook for option tabs
  nb <- ttknotebook(sidebarFrame)
  tkgrid(nb, row=2, column=0, columnspan=3, pady=5, sticky='nsew')
  
  mainOptionsFrame   <- tkframe(nb)
  fitOptionsFrame    <- tkframe(nb)
  mleSettingsFrame   <- tkframe(nb)
  advancedFrame      <- tkframe(nb)
  graphSettingsFrame <- tkframe(nb)
  
  tkadd(nb, mainOptionsFrame, text='Main Options')
  tkadd(nb, fitOptionsFrame, text='Fit Options')
  tkadd(nb, mleSettingsFrame, text='MLE Settings')
  tkadd(nb, advancedFrame, text='Advanced')
  tkadd(nb, graphSettingsFrame, text='Plot Settings')
  
  # Main Options Tab
  dtVar <- tclVar('0.1')
  stimTimeVar <- tclVar('100')
  baselineVar <- tclVar('50')
  nVar <- tclVar('30')
  yAblineVar <- tclVar('0.1')
  funcVar <- tclVar('product1N')
  tkgrid(tklabel(mainOptionsFrame, text='dt (ms):'), row=0, column=0, sticky='w')
  tkgrid(tkentry(mainOptionsFrame, textvariable=dtVar, width=10), row=0, column=1)
  tkgrid(tklabel(mainOptionsFrame, text='Stimulation Time:'), row=1, column=0, sticky='w')
  tkgrid(tkentry(mainOptionsFrame, textvariable=stimTimeVar, width=10), row=1, column=1)
  tkgrid(tklabel(mainOptionsFrame, text='Baseline:'), row=2, column=0, sticky='w')
  tkgrid(tkentry(mainOptionsFrame, textvariable=baselineVar, width=10), row=2, column=1)
  tkgrid(tklabel(mainOptionsFrame, text='n:'), row=3, column=0, sticky='w')
  tkgrid(tkentry(mainOptionsFrame, textvariable=nVar, width=10), row=3, column=1)
  tkgrid(tklabel(mainOptionsFrame, text='Fit cutoff:'), row=4, column=0, sticky='w')
  tkgrid(tkentry(mainOptionsFrame, textvariable=yAblineVar, width=10), row=4, column=1)
  tkgrid(tklabel(mainOptionsFrame, text='Function:'), row=5, column=0, sticky='w')
  funcChoices <- c('product1N', 'product2N', 'product3N')
  funcCombo <- ttkcombobox(mainOptionsFrame, textvariable=funcVar, values=funcChoices, width=10)
  tkgrid(funcCombo, row=5, column=1)
  
  dsVar <- tclVar('1')
  tkgrid(tklabel(mainOptionsFrame, text='Downsample Factor:'), row=6, column=0, sticky='w')
  tkgrid(tkentry(mainOptionsFrame, textvariable=dsVar, width=10), row=6, column=1)
  
  # Fit Options Tab
  NVar <- tclVar('1')
  IEIVar <- tclVar('50')
  smoothVar <- tclVar('5')
  methodVar <- tclVar('BF.LM')
  weightMethodVar <- tclVar('none')
  sequentialFitVar <- tclVar('0')
  intervalMinVar <- tclVar('0.1')
  intervalMaxVar <- tclVar('0.9')
  lowerVar <- tclVar('')
  upperVar <- tclVar('')
  latencyLimitVar <- tclVar('')
  tkgrid(tklabel(fitOptionsFrame, text='N:'), row=0, column=0, sticky='w')
  tkgrid(tkentry(fitOptionsFrame, textvariable=NVar, width=10), row=0, column=1)
  tkgrid(tklabel(fitOptionsFrame, text='IEI:'), row=1, column=0, sticky='w')
  tkgrid(tkentry(fitOptionsFrame, textvariable=IEIVar, width=10), row=1, column=1)
  tkgrid(tklabel(fitOptionsFrame, text='Smooth:'), row=2, column=0, sticky='w')
  tkgrid(tkentry(fitOptionsFrame, textvariable=smoothVar, width=10), row=2, column=1)
  tkgrid(tklabel(fitOptionsFrame, text='Method:'), row=3, column=0, sticky='w')
  methodChoices <- c('BF.LM', 'LM', 'GN', 'port', 'robust', 'MLE')
  methodCombo <- ttkcombobox(fitOptionsFrame, textvariable=methodVar, values=methodChoices, width=10)
  tkgrid(methodCombo, row=3, column=1)
  tkgrid(tklabel(fitOptionsFrame, text='Weighting:'), row=4, column=0, sticky='w')
  weightChoices <- c('none', '~y_sqrt', '~y')
  weightCombo <- ttkcombobox(fitOptionsFrame, textvariable=weightMethodVar, values=weightChoices, width=10)
  tkgrid(weightCombo, row=4, column=1)
  tkgrid(tklabel(fitOptionsFrame, text='Sequential Fit:'), row=5, column=0, sticky='w')
  sequentialFitCheck <- tkcheckbutton(fitOptionsFrame, variable=sequentialFitVar)
  tkgrid(sequentialFitCheck, row=5, column=1)
  tkgrid(tklabel(fitOptionsFrame, text='Min interval:'), row=6, column=0, sticky='w')
  tkgrid(tkentry(fitOptionsFrame, textvariable=intervalMinVar, width=10), row=6, column=1)
  tkgrid(tklabel(fitOptionsFrame, text='Max interval:'), row=7, column=0, sticky='w')
  tkgrid(tkentry(fitOptionsFrame, textvariable=intervalMaxVar, width=10), row=7, column=1)
  tkgrid(tklabel(fitOptionsFrame, text='Lower bounds:'), row=8, column=0, sticky='w')
  tkgrid(tkentry(fitOptionsFrame, textvariable=lowerVar, width=10), row=8, column=1)
  tkgrid(tklabel(fitOptionsFrame, text='Upper bounds:'), row=9, column=0, sticky='w')
  tkgrid(tkentry(fitOptionsFrame, textvariable=upperVar, width=10), row=9, column=1)
  tkgrid(tklabel(fitOptionsFrame, text='Latency limit:'), row=10, column=0, sticky='w')
  tkgrid(tkentry(fitOptionsFrame, textvariable=latencyLimitVar, width=10), row=10, column=1)
  
  # MLE Settings Tab
  iterVar <- tclVar('1000')
  metropolisScaleVar <- tclVar('1.5')
  fitAttemptsVar <- tclVar('10')
  RWmVar <- tclVar('0')
  tkgrid(tklabel(mleSettingsFrame, text='MLE Iterations:'), row=0, column=0, sticky='w')
  tkgrid(tkentry(mleSettingsFrame, textvariable=iterVar, width=10), row=0, column=1)
  tkgrid(tklabel(mleSettingsFrame, text='Metropolis Scale:'), row=1, column=0, sticky='w')
  tkgrid(tkentry(mleSettingsFrame, textvariable=metropolisScaleVar, width=10), row=1, column=1)
  tkgrid(tklabel(mleSettingsFrame, text='Fit Attempts:'), row=2, column=0, sticky='w')
  tkgrid(tkentry(mleSettingsFrame, textvariable=fitAttemptsVar, width=10), row=2, column=1)
  tkgrid(tklabel(mleSettingsFrame, text='Random Walk Metropolis:'), row=3, column=0, sticky='w')
  RWmCheck <- tkcheckbutton(mleSettingsFrame, variable=RWmVar)
  tkgrid(RWmCheck, row=3, column=1)
  
  # Advanced Tab
  filterVar <- tclVar('0')
  fcVar <- tclVar('1000')
  # relDecayFitLimitVar <- tclVar('0.1')
  halfWidthFitLimitVar <- tclVar('500')
  seedVar <- tclVar('42')
  dpVar <- tclVar('3')
  fastConstraintVar <- tclVar('0')
  fastConstraintMethodVar <- tclVar('rise')
  fastDecayLimitVar <- tclVar('')
  firstDelayConstraintVar <- tclVar('0')
  tkgrid(tklabel(advancedFrame, text='Filter:'), row=0, column=0, sticky='w')
  filterCheck <- tkcheckbutton(advancedFrame, variable=filterVar)
  tkgrid(filterCheck, row=0, column=1)
  tkgrid(tklabel(advancedFrame, text='Filter cutoff (Hz):'), row=1, column=0, sticky='w')
  tkgrid(tkentry(advancedFrame, textvariable=fcVar, width=10), row=1, column=1)
  tkgrid(tklabel(advancedFrame, text='Half-width fit limit:'), row=3, column=0, sticky='w')
  tkgrid(tkentry(advancedFrame, textvariable=halfWidthFitLimitVar, width=10), row=3, column=1)
  tkgrid(tklabel(advancedFrame, text='Seed:'), row=4, column=0, sticky='w')
  tkgrid(tkentry(advancedFrame, textvariable=seedVar, width=10), row=4, column=1)
  tkgrid(tklabel(advancedFrame, text='Decimal points:'), row=5, column=0, sticky='w')
  tkgrid(tkentry(advancedFrame, textvariable=dpVar, width=10), row=5, column=1)
  tkgrid(tklabel(advancedFrame, text='Fast constraint:'), row=6, column=0, sticky='w')
  fastConstraintCheck <- tkcheckbutton(advancedFrame, variable=fastConstraintVar)
  tkgrid(fastConstraintCheck, row=6, column=1)
  tkgrid(tklabel(advancedFrame, text='Fast constraint method:'), row=7, column=0, sticky='w')
  fastConstraintChoices <- c('rise', 'peak')
  fastConstraintCombo <- ttkcombobox(advancedFrame, textvariable=fastConstraintMethodVar, values=fastConstraintChoices, width=10)
  tkgrid(fastConstraintCombo, row=7, column=1)
  tkgrid(tklabel(advancedFrame, text='Fast decay limit(s):'), row=8, column=0, sticky='w')
  tkgrid(tkentry(advancedFrame, textvariable=fastDecayLimitVar, width=10), row=8, column=1)
  tkgrid(tklabel(advancedFrame, text='First delay constraint:'), row=9, column=0, sticky='w')
  firstDelayCheck <- tkcheckbutton(advancedFrame, variable=firstDelayConstraintVar)
  tkgrid(firstDelayCheck, row=9, column=1)
  
  # Graph Settings Tab
  lwdVar <- tclVar('1.2')
  xbarVar <- tclVar('50')
  ybarVar <- tclVar('50')
  xbarLabVar <- tclVar('ms')
  ybarLabVar <- tclVar('pA')
  tkgrid(tklabel(graphSettingsFrame, text='Line width:'), row=0, column=0, sticky='w')
  tkgrid(tkentry(graphSettingsFrame, textvariable=lwdVar, width=10), row=0, column=1)
  tkgrid(tklabel(graphSettingsFrame, text='x-bar length:'), row=1, column=0, sticky='w')
  tkgrid(tkentry(graphSettingsFrame, textvariable=xbarVar, width=10), row=1, column=1)
  tkgrid(tklabel(graphSettingsFrame, text='x-bar units:'), row=2, column=0, sticky='w')
  tkgrid(tkentry(graphSettingsFrame, textvariable=xbarLabVar, width=10), row=2, column=1)
  tkgrid(tklabel(graphSettingsFrame, text='y-bar length:'), row=3, column=0, sticky='w')
  tkgrid(tkentry(graphSettingsFrame, textvariable=ybarVar, width=10), row=3, column=1)
  tkgrid(tklabel(graphSettingsFrame, text='y-bar units:'), row=4, column=0, sticky='w')
  tkgrid(tkentry(graphSettingsFrame, textvariable=ybarLabVar, width=10), row=4, column=1)
  
  # Additional sidebar controls
  userTmaxVar <- tclVar('')
  tkgrid(tklabel(sidebarFrame, text='User maximum time for fit:'), row=3, column=0, sticky='w', pady=5)
  tkgrid(tkentry(sidebarFrame, textvariable=userTmaxVar, width=10), row=3, column=1, pady=5)
  
  repeatConstraintVar <- tclVar('0')
  tkgrid(tklabel(sidebarFrame, text='Add fast constraint:'), row=4, column=0, sticky='w')
  repeatConstraintCheck <- tkcheckbutton(sidebarFrame, variable=repeatConstraintVar)
  tkgrid(repeatConstraintCheck, row=4, column=1)
  
  # Analysis action buttons
  runAnalysisButton <- tkbutton(sidebarFrame, text='Run Initial Analysis', command=function() {
    filePath <- tclvalue(filePathVar)
    if (nchar(filePath) == 0) {
      tkmessageBox(message='Please select a file first')
      return()
    }
    if (nchar(tclvalue(columnVar)) == 0) {
      tkmessageBox(message='Please select a column')
      return()
    }
    ext <- tools::file_ext(filePath)
    if (tolower(ext) == 'csv') {
      uploaded_data <<- read.csv(filePath)
    } else {
      uploaded_data <<- readxl::read_excel(filePath)
    }
    response_data <<- uploaded_data[[tclvalue(columnVar)]]
    ds <- as.numeric(tclvalue(dsVar))
    if (ds > 1) {
      response_data <<- response_data[seq(1, length(response_data), by=ds)]
    }
    tkrreplot(plotWidget, fun=drawPlot1)
  })
  tkgrid(runAnalysisButton, row=5, column=0, columnspan=3, pady=5)
    
  runMainAnalysisButton <- tkbutton(sidebarFrame, text='Run Main Analysis', command=function() {
    fast.constraint        <- as.logical(as.numeric(tclvalue(repeatConstraintVar)))
    ds                     <- as.numeric(tclvalue(dsVar))
    dt                     <- as.numeric(tclvalue(dtVar)) * ds
    stimulation_time       <- as.numeric(tclvalue(stimTimeVar))
    baseline               <- as.numeric(tclvalue(baselineVar))
    smooth                 <- as.numeric(tclvalue(smoothVar))
    n                      <- as.numeric(tclvalue(nVar))
    N                      <- as.numeric(tclvalue(NVar))
    IEI                    <- as.numeric(tclvalue(IEIVar))
    func                   <- get(tclvalue(funcVar))
    method                 <- tclvalue(methodVar)
    weight_method          <- tclvalue(weightMethodVar)
    sequential.fit         <- as.logical(as.numeric(tclvalue(sequentialFitVar)))
    fit.limits             <- as.numeric(tclvalue(userTmaxVar))
    rel.decay.fit.limit    <- as.numeric(tclvalue(yAblineVar))
    lwd                    <- as.numeric(tclvalue(lwdVar))
    fc                     <- as.numeric(tclvalue(fcVar))
    interval               <- c(as.numeric(tclvalue(intervalMinVar)), as.numeric(tclvalue(intervalMaxVar)))
    lower                  <- if (nchar(tclvalue(lowerVar)) > 0) as.numeric(unlist(strsplit(tclvalue(lowerVar), ','))) else NULL
    upper                  <- if (nchar(tclvalue(upperVar)) > 0) as.numeric(unlist(strsplit(tclvalue(upperVar), ','))) else NULL
    iter                   <- as.numeric(tclvalue(iterVar))
    metropolis.scale       <- as.numeric(tclvalue(metropolisScaleVar))
    fit.attempts           <- as.numeric(tclvalue(fitAttemptsVar))
    RWm                    <- as.logical(as.numeric(tclvalue(RWmVar)))
    fast.decay.limit       <- if (nchar(tclvalue(fastDecayLimitVar)) > 0) as.numeric(unlist(strsplit(tclvalue(fastDecayLimitVar), ','))) else NULL
    fast.constraint.method <- tclvalue(fastConstraintMethodVar)
    first.delay.constraint <- as.logical(as.numeric(tclvalue(firstDelayConstraintVar)))
    dp                     <- as.numeric(tclvalue(dpVar))
    seed                   <- as.numeric(tclvalue(seedVar))
    filter                 <- as.logical(as.numeric(tclvalue(filterVar)))
    
    y <- response_data
    if (all(is.na(y[(which(!is.na(y))[length(which(!is.na(y)))] + 1):length(y)]))) {
      y <- y[!is.na(y)]
    }
    x <- seq(0, (length(y) - 1) * dt, by=dt)
    
    if (!sequential.fit) {
      tmax <- fit.limits
      x_limit <- determine_tmax2(y=y, N=N, stimulation_time=stimulation_time, baseline=baseline, lwd=lwd, 
                                smooth=smooth, tmax=tmax, y_abline=rel.decay.fit.limit, xbar=as.numeric(tclvalue(xbarVar)),
                                ybar=as.numeric(tclvalue(ybarVar)), xbar_lab=tclvalue(xbarLabVar), ybar_lab=tclvalue(ybarLabVar))
      adjusted_response <- y[x < x_limit]
      
      out <- nFIT(response=adjusted_response, n=n, N=N, IEI=IEI, dt=dt, func=func, method=method,
                  weight_method=weight_method, MLEsettings=list(iter=iter, metropolis.scale=metropolis.scale, 
                  fit.attempts=fit.attempts, RWm=RWm), stimulation_time=stimulation_time, baseline=baseline,
                  filter=filter, fc=fc, interval=interval, fast.decay.limit=fast.decay.limit, 
                  fast.constraint=fast.constraint, fast.constraint.method=fast.constraint.method, 
                  first.delay.constraint=first.delay.constraint, lower=lower, upper=upper,
                  latency.limit=if (nchar(tclvalue(latencyLimitVar)) > 0) as.numeric(unlist(strsplit(tclvalue(latencyLimitVar), ','))) else NULL,
                  return.output=TRUE, show.plot=FALSE, half_width_fit_limit=as.numeric(tclvalue(halfWidthFitLimitVar)),
                  dp=dp, height=5, width=5,seed=seed)
      
      out$traces <- traces_fun2(y=y, fits=out$fits, dt=dt, N=N, IEI=IEI, stimulation_time=stimulation_time,
                                baseline=baseline, func=func, filter=filter, fc=fc)
      tkrreplot(plotWidget, fun=function() {
        drawPlot2(traces=out$traces, func=func, lwd=lwd, cex=0.6, filter=filter, xbar=as.numeric(tclvalue(xbarVar)),
                  ybar=as.numeric(tclvalue(ybarVar)), xbar_lab=tclvalue(xbarLabVar), ybar_lab=tclvalue(ybarLabVar))
      })
    } else {
      out <- nFIT_sequential(response=y, n=n, dt=dt, func=func, method=method, weight_method=weight_method,
                  stimulation_time=stimulation_time, baseline=baseline, fit.limits=fit.limits, fast.decay.limit=fast.decay.limit,
                  fast.constraint=as.logical(as.numeric(tclvalue(fastConstraintVar))), fast.constraint.method=fast.constraint.method,
                  first.delay.constraint=first.delay.constraint, latency.limit=if (nchar(tclvalue(latencyLimitVar)) > 0) as.numeric(unlist(strsplit(tclvalue(latencyLimitVar), ','))) else NULL,
                  lower=lower, upper=upper, filter=filter, fc=fc, interval=interval,
                  MLEsettings=list(iter=iter, metropolis.scale=metropolis.scale, fit.attempts=fit.attempts, RWm=RWm),
                  MLE.method=method, half_width_fit_limit=as.numeric(tclvalue(halfWidthFitLimitVar)),
                  dp=dp, lwd=lwd, xlab='', ylab='', width=5, height=5, return.output=TRUE, show.output=TRUE,
                  show.plot=TRUE, seed=seed)
              }
    
    analysis_output <<- out
    df_out <- out$output
    if(sum(grepl('^A\\d+$', names(df_out))) == 1) names(df_out)[which(grepl('^A\\d+$', names(df_out)))] <- 'A'
    if(sum(grepl('^area\\d+$', names(df_out)))==1) names(df_out)[which(grepl('^area\\d+$', names(df_out)))] <- 'area'
    names(df_out) <- gsub("^r(\\d+)[_-](\\d+)$", "r\\1-\\2", names(df_out))
    names(df_out) <- gsub("^d(\\d+)[_-](\\d+)$", "d\\1-\\2", names(df_out))
    names(df_out)[names(df_out) == 'half_width'] <- 'half width'
    tkdelete(consoleText, '1.0', 'end')
    tkinsert(consoleText, 'end', 'Analysis complete.')
    tkdelete(fitOutputText, '1.0', 'end')
    tkinsert(fitOutputText, 'end', paste(capture.output(print(df_out)), collapse='\n'))
  })
  tkgrid(runMainAnalysisButton, row=7, column=0, columnspan=3, pady=5)
  
  downloadOutputButton <- tkbutton(sidebarFrame, text='Download Output', 
    command=function() {
      if (!exists('analysis_output') || is.null(analysis_output)) {
        tkmessageBox(message='No analysis output available!')
        return()
      }
      saveFile <- tclvalue(tkgetSaveFile(filetypes='{{Rdata Files} {.Rdata}} {{All Files} *}'))
      if (nchar(saveFile) > 0) {
        # Build metadata list using the tk variable values.
        metadata <- list(
          dt=as.numeric(tclvalue(dtVar)),
          stimTime=as.numeric(tclvalue(stimTimeVar)),
          baseline=as.numeric(tclvalue(baselineVar)),
          n=as.numeric(tclvalue(nVar)),
          y_abline=as.numeric(tclvalue(yAblineVar)),
          func=tclvalue(funcVar),
          ds=as.numeric(tclvalue(dsVar)),
          N=as.numeric(tclvalue(NVar)),
          IEI=as.numeric(tclvalue(IEIVar)),
          smooth=as.numeric(tclvalue(smoothVar)),
          method=tclvalue(methodVar),
          weight_method=tclvalue(weightMethodVar),
          sequential_fit=as.logical(as.numeric(tclvalue(sequentialFitVar))),
          interval=c(as.numeric(tclvalue(intervalMinVar)), as.numeric(tclvalue(intervalMaxVar))),
          lower=if (nchar(tclvalue(lowerVar)) > 0)
                    as.numeric(unlist(strsplit(tclvalue(lowerVar), ",")))
                  else NULL,
          upper=if (nchar(tclvalue(upperVar)) > 0)
                    as.numeric(unlist(strsplit(tclvalue(upperVar), ",")))
                  else NULL,
          latency_limit=if (nchar(tclvalue(latencyLimitVar)) > 0)
                            as.numeric(unlist(strsplit(tclvalue(latencyLimitVar), ",")))
                          else NULL,
          iter=as.numeric(tclvalue(iterVar)),
          metropolis_scale=as.numeric(tclvalue(metropolisScaleVar)),
          fit_attempts=as.numeric(tclvalue(fitAttemptsVar)),
          RWm=as.logical(as.numeric(tclvalue(RWmVar))),
          fast_decay_limit=if (nchar(tclvalue(fastDecayLimitVar)) > 0)
                               as.numeric(unlist(strsplit(tclvalue(fastDecayLimitVar), ",")))
                             else NULL,
          fast_constraint=as.logical(as.numeric(tclvalue(fastConstraintVar))),
          fast_constraint_method=tclvalue(fastConstraintMethodVar),
          first_delay_constraint=as.logical(as.numeric(tclvalue(firstDelayConstraintVar))),
          dp=as.numeric(tclvalue(dpVar)),
          seed=as.numeric(tclvalue(seedVar)),
          filter=as.logical(as.numeric(tclvalue(filterVar))),
          fc=as.numeric(tclvalue(fcVar)),
          userTmax=as.numeric(tclvalue(userTmaxVar)),
          data_col=tclvalue(columnVar)
        )
        # Combine analysis output and metadata into one list.
        results <- list(
          analysis=analysis_output,
          metadata=metadata
        )
        save(results, file=saveFile)
        tkmessageBox(message='Output saved successfully.')
      }
    }
  )

  tkgrid(downloadOutputButton, row=8, column=0, columnspan=3, pady=5)
  
  clearOutputButton <- tkbutton(sidebarFrame, text='Clear Output', command=function() {
    analysis_output <<- NULL
    tkdelete(consoleText, '1.0', 'end')
    tkdelete(fitOutputText, '1.0', 'end')
    tkrreplot(plotWidget, fun=drawPlot1)
  })
  tkgrid(clearOutputButton, row=9, column=0, columnspan=3, pady=5)
  
  # Main Panel plot
  drawPlot1 <- function() {
    
    ds <- as.numeric(tclvalue(dsVar))
    dt <- as.numeric(tclvalue(dtVar)) * ds
    lwd <- as.numeric(tclvalue(lwdVar))
    stimTime <- as.numeric(tclvalue(stimTimeVar))
    baseline <- as.numeric(tclvalue(baselineVar))
    smooth <- as.numeric(tclvalue(smoothVar))
    y_abline <- as.numeric(tclvalue(yAblineVar))
    y_val <- if (exists('response_data') && !is.null(response_data)) response_data else rnorm(10000, 0.1)
    
    determine_tmax2(y=y_val, N=1, dt=dt, stimulation_time=stimTime, baseline=baseline, smooth=smooth, lwd=lwd,
      tmax=NULL, y_abline=y_abline, xbar=as.numeric(tclvalue(xbarVar)), ybar=as.numeric(tclvalue(ybarVar)),
      xbar_lab=tclvalue(xbarLabVar), ybar_lab=tclvalue(ybarLabVar))
  }
  
  plotWidget <- tkrplot(tt, fun=drawPlot1)
  tkgrid(plotWidget, row=0, column=1, sticky='nsew')
  
  consoleText <- tktext(mainFrame, width=80, height=4)
  tkgrid(consoleText, row=1, column=1, sticky='nsew')
  
  fitOutputLabel <- tklabel(sidebarFrame, text='Fit Output:')
  tkgrid(fitOutputLabel, row=10, column=0, columnspan=3, sticky='w', pady=c(10,2), padx=20)
  
  fitOutputText <- tktext(sidebarFrame, width=90, height=5)
  tkgrid(fitOutputText, row=11, column=0, columnspan=3, sticky='w', padx=20)
  
  tkfocus(tt)
}

# Launch the PSC Analysis interface
PSC_analysis_tk()




###########################
# Shiny version of analyse_PSC function

rm(list=ls(all=TRUE))
graphics.off()

# load and install necessary packages # 
load_required_packages <- function(packages) {
  new.packages <- packages[!(packages %in% installed.packages()[, 'Package'])]
  if (length(new.packages)) install.packages(new.packages)
  invisible(lapply(packages, library, character.only=TRUE))
}

required.packages <- c('robustbase', 'minpack.lm', 'Rcpp', 'signal',
                       'dbscan', 'tkrplot', 'tcltk', 'readxl', 'shiny', 'readxl')
load_required_packages(required.packages)

# Define your username and repository path, and source your custom functions.
username <- 'euo9382'
path_repository <- '/Documents/Repositories/Rfits'
file_path1 <- paste0('/Users/', username, path_repository)
source(paste0(file_path1, '/nNLS functions.R'))


# Shiny User Interface
ui <- fluidPage(
  titlePanel('PSC Analysis'),
  sidebarLayout(
    sidebarPanel(
      fileInput('file', 'Upload CSV or XLSX', accept=c('.csv', '.xlsx')),
      uiOutput('column_selector'),
      
      tabsetPanel(
        tabPanel('Main Options',
                 numericInput('dt', 'dt (ms):', 0.1),
                 numericInput('stimulation_time', 'Stimulation Time:', 100),
                 numericInput('baseline', 'Baseline:', 50),
                 numericInput('n', 'n:', 30),
                 numericInput('y_abline', 'Fit Cutoff:', 0.1),
                 selectInput('func', 'Function:', choices=c('product1N', 'product2N', 'product3N')),
                 checkboxInput('fast_constraint', 'Fast Constraint', FALSE),
                 numericInput('ds', 'Downsample Factor:', 1, min=1)
        ),
        tabPanel('Fit Options',
                 numericInput('N', 'N:', 1),
                 numericInput('IEI', 'IEI:', 50),
                 numericInput('smooth', 'Smooth:', 5),
                 selectInput('method', 'Method:', choices=c('BF.LM', 'LM', 'GN', 'port', 'robust', 'MLE')),
                 selectInput('weight_method', 'Weighting:', choices=c('none', '~y_sqrt', '~y')),
                 checkboxInput('sequential_fit', 'Sequential Fit', FALSE),
                 numericInput('interval_min', 'Min Interval:', 0.1),
                 numericInput('interval_max', 'Max Interval:', 0.9),
                 textInput('lower', 'Lower Bounds (comma-separated):', ''),
                 textInput('upper', 'Upper Bounds (comma-separated):', ''),
                 textInput('latency_limit', 'Latency Limit:', '')
        ),
        tabPanel('MLE Settings',
                 numericInput('iter', 'MLE Iterations:', 1000),
                 numericInput('metropolis_scale', 'Metropolis Scale:', 1.5),
                 numericInput('fit_attempts', 'Fit Attempts:', 10),
                 checkboxInput('RWm', 'Random Walk Metropolis', FALSE)
        ),
        tabPanel('Advanced',
                 checkboxInput('filter', 'Filter', FALSE),
                 numericInput('fc', 'Filter Cutoff (Hz):', 1000),
                 numericInput('half_width_fit_limit', 'Half-width Fit Limit:', 500),
                 numericInput('seed', 'Seed:', 42),
                 numericInput('dp', 'Decimal Points:', 3),
                 checkboxInput('fast_constraint', 'Fast Constraint', FALSE),
                 selectInput('fast_constraint_method', 'Fast Constraint Method:', choices=c('rise', 'peak')),
                 textInput('fast_decay_limit', 'Fast Decay Limit(s) (comma-separated):', ''),
                 checkboxInput('first_delay_constraint', 'First Delay Constraint', FALSE)
        ),
        tabPanel('Plot Settings',
                 numericInput('lwd', 'Line Width:', 1.2),
                 numericInput('xbar', 'x-bar Length:', 50),
                 numericInput('ybar', 'y-bar Length:', 50),
                 textInput('xbar_lab', 'x-axis Units:', 'ms'),
                 textInput('ybar_lab', 'y-axis Units:', 'pA')
        )
      ),
      
      numericInput('userTmax', 'User Maximum Time for Fit:', NA),
      actionButton('run_initial', 'Run Initial Analysis'),
      actionButton('run_main', 'Run Main Analysis'),
      actionButton('clear_output', 'Clear Output'),
      downloadButton('download_output', 'Download Output')
    ),
    mainPanel(
      plotOutput('plot', height='500px'),
      verbatimTextOutput('console')
    )
  )
)

server <- function(input, output, session) {
  
  # Reactive values to store data and analysis results.
  state <- reactiveValues(
    response=NULL,
    analysis=NULL
  )
  
  # upload file
  uploaded_data <- reactive({
    req(input$file)
    ext <- tools::file_ext(input$file$name)
    if (tolower(ext) == 'csv') {
      read.csv(input$file$datapath)
    } else {
      readxl::read_excel(input$file$datapath)
    }
  })
  
  # update the column selector
  output$column_selector <- renderUI({
    req(uploaded_data())
    selectInput('data_col', 'Select Column to Analyse', choices=names(uploaded_data()))
  })
  
  # run initial analysis
  observeEvent(input$run_initial, {
    req(uploaded_data(), input$data_col)
    # Clear any previous response and analysis.
    state$response <- NULL
    state$analysis <- NULL
    # Extract the column from the uploaded data.
    data_col <- uploaded_data()[[input$data_col]]
    ds <- as.numeric(input$ds)
    if (ds > 1) {
      data_col <- data_col[seq(1, length(data_col), by=ds)]
    }
    state$response <- data_col
  })
  
  # update when downsampled
  observeEvent(input$ds, {
    req(uploaded_data(), input$data_col)
    # Only proceed if a response is already loaded.
    if (!is.null(state$response)) {
      data_col <- uploaded_data()[[input$data_col]]
      ds <- as.numeric(input$ds)
      if (ds > 1) {
        data_col <- data_col[seq(1, length(data_col), by=ds)]
      }
      state$response <- data_col
      # Also clear any analysis result to force re-running the main analysis.
      state$analysis <- NULL
      cat('Downsample factor changed: Updated response with length =', length(data_col), '\n')
    }
  }, ignoreInit=TRUE)
  
  # plot output
  output$plot <- renderPlot({
    req(state$response)
    # Compute effective dt using the current ds.
    dt <- as.numeric(input$dt) * as.numeric(input$ds)
    lwd <- as.numeric(input$lwd)
    stim_time <- as.numeric(input$stimulation_time)
    baseline <- as.numeric(input$baseline)
    smooth <- as.numeric(input$smooth)
    y_abline <- as.numeric(input$y_abline)
    xbar <- as.numeric(input$xbar)
    ybar <- as.numeric(input$ybar)
    xbar_lab <- input$xbar_lab
    ybar_lab <- input$ybar_lab
    
    if (is.null(state$analysis)) {
      determine_tmax2(y=state$response, N=as.numeric(input$N), dt=dt, stimulation_time=stim_time, baseline=baseline, smooth=smooth,
                      lwd=lwd, cex=1, tmax=NULL, y_abline=y_abline, xbar=xbar, ybar=ybar, xbar_lab=xbar_lab, ybar_lab=ybar_lab)
    } else {
      req(state$analysis$traces)
      func <- switch(input$func,
                     'product1N'=product1N,
                     'product2N'=product2N,
                     'product3N'=product3N,
                     product1N)
      drawPlot2(traces=state$analysis$traces, func=func, lwd=lwd,
                filter=input$filter, xbar=xbar, ybar=ybar,
                xbar_lab=xbar_lab, ybar_lab=ybar_lab)
    }
  })
  
  # run main analysis
  observeEvent(input$run_main, {
    req(state$response)
    
    dt <- as.numeric(input$dt) * as.numeric(input$ds)
    stim_time <- as.numeric(input$stimulation_time)
    baseline <- as.numeric(input$baseline)
    smooth <- as.numeric(input$smooth)
    n <- as.numeric(input$n)
    N <- as.numeric(input$N)
    IEI <- as.numeric(input$IEI)
    method <- input$method
    weight_method <- input$weight_method
    sequential_fit <- input$sequential_fit
    interval <- c(as.numeric(input$interval_min), as.numeric(input$interval_max))
    lower <- if (nchar(input$lower) > 0) as.numeric(unlist(strsplit(input$lower, ','))) else NULL
    upper <- if (nchar(input$upper) > 0) as.numeric(unlist(strsplit(input$upper, ','))) else NULL
    latency_limit <- if (nchar(input$latency_limit) > 0) as.numeric(unlist(strsplit(input$latency_limit, ','))) else NULL
    iter <- as.numeric(input$iter)
    metropolis_scale <- as.numeric(input$metropolis_scale)
    fit_attempts <- as.numeric(input$fit_attempts)
    RWm <- input$RWm
    fast_decay_limit <- if (nchar(input$fast_decay_limit) > 0) as.numeric(unlist(strsplit(input$fast_decay_limit, ','))) else NULL
    fast_constraint <- input$fast_constraint
    fast_constraint_method <- input$fast_constraint_method
    first_delay_constraint <- input$first_delay_constraint
    dp <- as.numeric(input$dp)
    seed <- as.numeric(input$seed)
    filter_flag <- input$filter
    fc <- as.numeric(input$fc)
    
    y <- state$response
    if (any(is.na(y))) y <- y[!is.na(y)]
    x <- seq(0, (length(y) - 1) * dt, by=dt)
    
    tmax_value <- if (is.na(as.numeric(input$userTmax))) {
      determine_tmax2(y=y, N=N, dt=dt, stimulation_time=stim_time, baseline=baseline,
                      smooth=smooth, tmax=NULL, y_abline=as.numeric(input$y_abline),
                      xbar=as.numeric(input$xbar), ybar=as.numeric(input$ybar),
                      xbar_lab=input$xbar_lab, ybar_lab=input$ybar_lab)
    } else {
      as.numeric(input$userTmax)
    }
    x_limit <- tmax_value
    
    adjusted_response <- y[x < x_limit]

        func <- switch(input$func,
                   'product1N'=product1N,
                   'product2N'=product2N,
                   'product3N'=product3N,
                   product1N)
    
    if (!sequential_fit) {
      result <- nFIT(response=adjusted_response, n=n, N=N, IEI=IEI, dt=dt, func=func,
                     method=method, weight_method=weight_method,
                     MLEsettings=list(iter=iter, metropolis.scale=metropolis_scale, fit.attempts=fit_attempts, RWm=RWm),
                     stimulation_time=stim_time, baseline=baseline, filter=filter_flag, fc=fc,
                     interval=interval, fast.decay.limit=fast_decay_limit, fast.constraint=fast_constraint,
                     fast.constraint.method=fast_constraint_method, first.delay.constraint=first_delay_constraint,
                     lower=lower, upper=upper, latency.limit=latency_limit,
                     return.output=TRUE, show.plot=FALSE, half_width_fit_limit=as.numeric(input$half_width_fit_limit),
                     dp=dp, height=5, width=5, seed=seed)
      result$traces <- traces_fun2(y=y, fits=result$fits, dt=dt, N=N, IEI=IEI,
                                   stimulation_time=stim_time, baseline=baseline, func=func,
                                   filter=filter_flag, fc=fc)
    } else {
      result <- nFIT_sequential(response=y, n=n, dt=dt, func=func, method=method, weight_method=weight_method,
                                stimulation_time=stim_time, baseline=baseline, fit.limits=as.numeric(input$userTmax),
                                fast.decay.limit=fast_decay_limit, fast.constraint=fast_constraint,
                                fast.constraint.method=fast_constraint_method, first.delay.constraint=first_delay_constraint,
                                latency.limit=latency_limit, lower=lower, upper=upper, filter=filter_flag, fc=fc, interval=interval,
                                MLEsettings=list(iter=iter, metropolis.scale=metropolis_scale, fit.attempts=fit_attempts, RWm=RWm),
                                MLE.method=method, half_width_fit_limit=as.numeric(input$half_width_fit_limit),
                                dp=dp, lwd=as.numeric(input$lwd), xlab='', ylab='', width=5, height=5,
                                return.output=TRUE, show.output=TRUE, show.plot=TRUE, seed=seed)
    }
    
    state$analysis <- result
  })
  
  # clear output
  observeEvent(input$clear_output, {
    state$analysis <- NULL
  })
  
  # output to console
  output$console <- renderPrint({
    if (!is.null(state$analysis)) {

      df_out <- state$analysis$output
      if(sum(grepl('^A\\d+$', names(df_out))) == 1) names(df_out)[which(grepl('^A\\d+$', names(df_out)))] <- 'A'
      if(sum(grepl('^area\\d+$', names(df_out)))==1) names(df_out)[which(grepl('^area\\d+$', names(df_out)))] <- 'area'
      names(df_out) <- gsub("^r(\\d+)[_-](\\d+)$", "r\\1-\\2", names(df_out))
      names(df_out) <- gsub("^d(\\d+)[_-](\\d+)$", "d\\1-\\2", names(df_out))
      names(df_out)[names(df_out) == 'half_width'] <- 'half width'

      print(df_out)
    } else {
      cat('No analysis output performed')
    }
  })
  
  # download output
  output$download_output <- downloadHandler(
    filename=function() {
      req(input$file)
      paste0(
        tools::file_path_sans_ext(basename(input$file$name)),
        "_", input$data_col,
        "_PSC_analysis.RData"
      )
    },
    content=function(file) {
      # all settings
      metadata <- list(
        dt=as.numeric(input$dt),
        ds=as.numeric(input$ds),
        stimulation_time=as.numeric(input$stimulation_time),
        baseline=as.numeric(input$baseline),
        n=as.numeric(input$n),
        y_abline=as.numeric(input$y_abline),
        func=input$func,
        N=as.numeric(input$N),
        IEI=as.numeric(input$IEI),
        smooth=as.numeric(input$smooth),
        method=input$method,
        weight_method=input$weight_method,
        sequential_fit=input$sequential_fit,
        interval=c(as.numeric(input$interval_min),
                     as.numeric(input$interval_max)),
        lower=if(nchar(input$lower) > 0)
                  as.numeric(unlist(strsplit(input$lower, ",")))
                else NULL,
        upper=if(nchar(input$upper) > 0)
                  as.numeric(unlist(strsplit(input$upper, ",")))
                else NULL,
        latency_limit=if(nchar(input$latency_limit) > 0)
                          as.numeric(unlist(strsplit(input$latency_limit, ",")))
                        else NULL,
        iter=as.numeric(input$iter),
        metropolis_scale=as.numeric(input$metropolis_scale),
        fit_attempts=as.numeric(input$fit_attempts),
        RWm=input$RWm,
        fast_decay_limit=if(nchar(input$fast_decay_limit) > 0) as.numeric(unlist(strsplit(input$fast_decay_limit, ","))) else NULL,
        fast_constraint=input$fast_constraint,
        fast_constraint_method=input$fast_constraint_method,
        first_delay_constraint=input$first_delay_constraint,
        dp=as.numeric(input$dp),
        seed=as.numeric(input$seed),
        filter=input$filter,
        fc=as.numeric(input$fc),
        userTmax=as.numeric(input$userTmax),
        data_col=input$data_col
      )
      
      # save analysis and metadata
      results <- list(
        analysis=state$analysis,
        metadata=metadata
      )
      
      save(results, file=file)
    }
  )
  
}
  
# Run the Shiny App
shinyApp(ui=ui, server=server)
