# have to run R gui for tcl packages to execute
# open -n -a R

# this version combines versions A and B below:
# Remove all objects from the environment
rm(list = ls(all = TRUE))

load_required_packages <- function(packages) {
  new.packages <- packages[!(packages %in% installed.packages()[, 'Package'])]
  if (length(new.packages)) install.packages(new.packages)
  invisible(lapply(packages, library, character.only = TRUE))
}
required.packages <- c('robustbase', 'minpack.lm', 'Rcpp', 'shiny', 'signal', 'readABF', 'tcltk', 'tkrplot')
load_required_packages(required.packages)



analyseABFtk <- function() {

  # helper Functions
  extract_metadata <- function(abf_dataset) {
    list(
      path                  = abf_dataset$path,
      formatVersion         = abf_dataset$formatVersion,
      channelNames          = abf_dataset$channelNames,
      channelUnits          = abf_dataset$channelUnits,
      samplingIntervalInSec = abf_dataset$samplingIntervalInSec,
      header                = abf_dataset$header,
      tags                  = abf_dataset$tags,
      sections              = abf_dataset$sections
    )
  }

  choose_data_column <- function(channelUnits, experiment) {
    if (experiment == 'voltage clamp') {
      idx <- grep('A', channelUnits, ignore.case = TRUE)
    } else if (experiment == 'current clamp') {
      idx <- grep('V', channelUnits, ignore.case = TRUE)
    } else {
      idx <- integer(0)
    }
    if (length(idx) > 0) return(idx[1])
    else return(NA)
  }

  check_consistency <- function(metadata) {
    dt_values <- sapply(metadata, function(meta) meta$samplingIntervalInSec * 1000)
    traces_values <- sapply(metadata, function(meta) meta$header$lActualEpisodes)
    expType <- tclvalue(experimentVar)
    unit_values <- sapply(metadata, function(meta) {
      col_idx <- choose_data_column(meta$channelUnits, expType)
      if (!is.na(col_idx)) meta$channelUnits[col_idx] else NA_character_
    })
    dt_good <- (length(unique(dt_values)) == 1)
    traces_good <- (length(unique(traces_values)) == 1)
    unit_good <- (length(unique(unit_values)) == 1)
    if (dt_good && traces_good && unit_good) {
      return('Data is consistent')
    } else {
      error_msgs <- c()
      if (!dt_good) error_msgs <- c(error_msgs, paste('Inconsistent dt values:', paste(dt_values, collapse = ', ')))
      if (!unit_good) error_msgs <- c(error_msgs, paste('Inconsistent Units:', paste(unit_values, collapse = ', ')))
      if (!traces_good) error_msgs <- c(error_msgs, paste('Inconsistent Traces:', paste(traces_values, collapse = ', ')))
      return(paste(error_msgs, collapse = '; '))
    }
  }

  egs_plot <- function(x, y, sign = -1, xlim = NULL, ylim = NULL, lwd = 1, show_text = FALSE, 
                       xbar = 100, ybar = 50, color = '#4C77BB', show_bar = FALSE, cex = 0.6) {
    if (is.null(ylim))
      ylim <- if (sign == 1) c(0, max(y)) else c(-max(-y), 0)
    if (is.null(xlim))
      xlim <- c(min(x), max(x))
    idx1 <- which.min(abs(x - xlim[1]))
    idx2 <- which.min(abs(x - xlim[2]))
    plot(x[idx1:idx2], y[idx1:idx2], type = 'l', col = color, xlim = xlim, ylim = ylim, bty = 'n', 
      lwd = lwd, lty = 1, axes = FALSE, frame = FALSE, xlab = '', ylab = '')
    if (show_bar) {
      ybar_start <- min(ylim) + (max(ylim) - min(ylim)) / 20
      x_start <- max(xlim) - xbar - 50
      y_start <- ybar_start
      x_end <- x_start + xbar
      y_end <- y_start + ybar
      segments(x_start, y_start, x_end, y_start, lwd = lwd, col = 'black')
      segments(x_start, y_start, x_start, y_end, lwd = lwd, col = 'black')
      if (show_text) {
        text(x = (x_start + x_end) / 2, y = y_start - ybar / 20, labels = paste(xbar, 'ms'), 
             adj = c(0.5, 1), cex = cex)
        text(x = x_start - xbar / 4, y = (y_start + y_end) / 2, labels = paste(ybar, 'pA'), 
             adj = c(0.5, 0.5), srt = 90, cex = cex)
      }
    }
  }

  load_abf_data <- function(abf_files = NULL, abf_path = NULL) {
    abf_path <- if (is.null(abf_path)) getwd() else abf_path
    setwd(abf_path)
    N <- length(abf_files)
    datasets <- lapply(seq_len(N), function(ii) readABF(abf_files[ii]))
    names(datasets) <- abf_files
    metadata <- lapply(datasets, extract_metadata)
    return(list(datasets = datasets, metadata = metadata))
  }

  download_data <- function() {
    # check for averaged data
    if (is.null(averaged_data) || length(averaged_data) == 0) {
      tkmessageBox(message = "No averaged data available.")
      return()
    }

    df <- as.data.frame(do.call(cbind, averaged_data))
    colnames(df) <- as.character(seq_len(length(averaged_data)))

    download_folder <- tclvalue(folderPathVar)
    if (nchar(download_folder) == 0) {
      tkmessageBox(message = "No folder selected for download.")
      return()
    }
    file_path <- file.path(download_folder, 'averaged_data.csv')

    # write CSV (no time column, just the numbered averages)
    write.csv(df, file = file_path, row.names = FALSE)
    tkmessageBox(message = paste('Averaged data saved to', file_path))
  }

  abf_averages <- function(datasets, baseline = 100, stimulation_time = 350, traces2average = NULL, dataCol = 1, ylim = NULL, xlim = NULL, 
    color = 'darkgray', xbar = 100, ybar = 50, width = 5.25, height = 2.75, save = FALSE, plotIt = TRUE) {
    
    N <- length(datasets)
    sampling_intervals <- sapply(datasets, function(ds) ds$samplingIntervalInSec * 1000)
    responses <- lapply(seq_len(N), function(iii) {
      sapply(seq_along(datasets[[iii]]$data), function(ii) {
        datasets[[iii]]$data[[ii]][, dataCol]
      })
    })
    names(responses) <- names(datasets)
    baseline2zero <- function(y, dt, stim, baseline) {
      idx_baseline <- round(baseline / dt)
      idx_start    <- round((stim - baseline) / dt) + 1
      y0 <- y - mean(y[1:idx_baseline])
      y0[idx_start:length(y0)]
    }
    responses0 <- lapply(seq_len(N), function(iii) {
      sapply(seq_len(ncol(responses[[iii]])), function(jj) {
        baseline2zero(responses[[iii]][, jj],
                      dt = sampling_intervals[iii],
                      stim = stimulation_time,
                      baseline = baseline)
      })
    })
    names(responses0) <- names(responses)
    responses0_mean <- if(is.null(traces2average)) {
      lapply(seq_len(N), function(iii) apply(responses0[[iii]], 1, mean))
    } else {
      lapply(seq_len(N), function(iii)
        apply(responses0[[iii]][, traces2average[[iii]], drop = FALSE], 1, mean))
    }
    time <- lapply(seq_len(N), function(iii) {
      dt_val <- sampling_intervals[iii]
      stim_time <- stimulation_time
      base_val <- baseline
      seq(from = stim_time - base_val, by = dt_val, length.out = length(responses0_mean[[iii]]))
    })
    if(plotIt){
      par(mfrow = c(1, N))
      show_bar <- rep(FALSE, N)
      if (N > 0) show_bar[N] <- TRUE
      for(ii in seq_len(N)) {
        for(ii in seq_len(N)) {
          egs_plot(x = time[[ii]], y = responses0_mean[[ii]], color = 'darkgray',
                   show_bar = FALSE, show_text = FALSE)
        }
      }
    }
    return(list(raw_data = responses,
                baseline_corrected_data = responses0,
                baseline_corrected_mean_data = responses0_mean,
                datasets = datasets))
  }

  combine_abf_data <- function(result) {
    master_abf <- list()
    master_abf$data <- list()
    master_abf$source_files <- c()
    master_abf$samplingIntervalInSec <- result$datasets[[1]]$samplingIntervalInSec
    for(i in seq_along(result$datasets)) {
      ds <- result$datasets[[i]]
      n_traces <- length(ds$data)
      master_abf$data <- c(master_abf$data, ds$data)
      master_abf$source_files <- c(master_abf$source_files, rep(names(result$datasets)[i], n_traces))
    }
    return(master_abf)
  }

  smart_axis_limits <- function(vec, n_steps = 5) {
    rng <- range(vec)
    spread <- diff(rng)
    
    # Pick a base that's a nice round number (1, 2, 5, 10, 20, 50, 100, etc.)
    raw_step <- spread / n_steps
    base <- 10^floor(log10(raw_step))
    
    # Refine to nicer step (1, 2, or 5 × 10^n)
    nice_steps <- c(1, 2, 5, 10)
    best_step <- base * nice_steps[which.min(abs(nice_steps * base - raw_step))]
    
    lower <- floor(rng[1] / best_step) * best_step
    upper <- ceiling(rng[2] / best_step) * best_step
    c(lower, upper)
  }

  # global variables
  master_abf <<- NULL      # Will hold either a concatenated master object or the original structure.
  averaged_data <<- NULL   # Will hold the averaged (baseline_corrected_mean) data.
  traces2average <<- list()  # Used in separate mode.
  # For concatenated (master) mode:
  current_trace <<- 1      
  total_traces <<- 0       
  current_group_selected <<- integer(0)  
  groups_list <<- list()  
  # For separate (non-concatenated) mode:
  current_dataset <<- 1  

  # review functions
  # for concatenated
  review_master_recordings <- function() {
    if (is.null(master_abf)) {
      tkmessageBox(message = "No master ABF data available. Please load data first.")
      return()
    }
    total_traces <<- length(master_abf$data)
    current_trace <<- 1L
    current_group_selected <<- integer(0)
    groups_list <<- list()

    # clear out old widgets
    children <- as.character(tkwinfo('children', plotPanel))
    if (length(children) > 0) {
      sapply(children, function(ch) tcl("destroy", ch))
    }

    reviewFrame <<- tkframe(plotPanel)
    tkgrid(reviewFrame, row = 0, column = 0, sticky = 'nsew')

    infoLabel <<- tklabel(
      reviewFrame,
      text = paste('Trace', current_trace, 'of', total_traces)
    )
    tkgrid(infoLabel, row = 0, column = 0, columnspan = 3)

    reviewPlot <<- tkrplot(
      reviewFrame,
      fun = function() {
        mat   <- master_abf$data[[current_trace]]
        dt    <- master_abf$samplingIntervalInSec * 1000
        time  <- seq(0, by = dt, length.out = nrow(mat))
        dc    <- as.numeric(tclvalue(dataColVar))
        if (is.na(dc) || dc < 1 || dc > ncol(mat)) dc <- 1
        trace <- mat[, dc]
        par(cex.lab = 0.6, cex.axis = 0.6, cex.main = 0.6)
        plot(
          time, trace, type = 'l', col = 'darkgray', xlab = 'time (ms)',
          xlim = smart_axis_limits(time),
          ylim = smart_axis_limits(trace),
          ylab = tclvalue(unitVar),
          axes = FALSE, bty = 'l',
          main = paste('Trace', current_trace, 'of', total_traces)
        )
        axis(1); axis(2, las = 1)
      },
      hscale = 1, vscale = 1
    )
    tkgrid(reviewPlot, row = 1, column = 0, columnspan = 3)

    # helper to advance and redraw
    move_next_master <- function() {
      tkconfigure(acceptButton, state = 'normal')
      tkconfigure(rejectButton, state = 'normal')
      if (current_trace < total_traces) {
        current_trace <<- current_trace + 1L
        tkconfigure(infoLabel, text = paste('Trace', current_trace, 'of', total_traces))
        tkrreplot(reviewPlot)
      } else {
        tkmessageBox(message = "Review complete for all traces.")
      }
    }

    acceptButton <<- tkbutton(
      reviewFrame, text = 'Accept',
      command = function() {
        current_group_selected <<- c(current_group_selected, current_trace)
        tkconfigure(acceptButton, state = 'disabled')
        tkconfigure(rejectButton, state = 'normal')
        move_next_master()
      }
    )
    tkgrid(acceptButton, row = 2, column = 0)

    rejectButton <<- tkbutton(
      reviewFrame, text = 'Reject',
      command = function() {
        tkconfigure(rejectButton, state = 'disabled')
        tkconfigure(acceptButton, state = 'normal')
        move_next_master()
      }
    )
    tkgrid(rejectButton, row = 2, column = 1)

    nextTraceButton <<- tkbutton(
      reviewFrame, text = 'Next Trace',
      command = move_next_master
    )
    tkgrid(nextTraceButton, row = 2, column = 2)

    averageGroupButton <<- tkbutton(
      reviewFrame, text = 'Add Selected Group',
      command = function() {
        if (length(current_group_selected) == 0) {
          tkinsert(consoleText, 'end', 'No traces selected in current group.\n')
        } else {
          groups_list[[length(groups_list) + 1]] <<- current_group_selected
          msg <- paste(
            'Group', length(groups_list),
            'selected with traces:', paste(current_group_selected, collapse = ', ')
          )
          tkinsert(consoleText, 'end', paste0(msg, '\n'))
          current_group_selected <<- integer(0)
        }
      }
    )
    tkgrid(averageGroupButton, row = 3, column = 0, columnspan = 3)

    selectionCompleteButton <<- tkbutton(
      reviewFrame, text = 'Selection Complete',
      command = function() {
        tkinsert(consoleText, 'end', 'Review complete: Approved traces stored.\n')
      }
    )
    tkgrid(selectionCompleteButton, row = 4, column = 0, columnspan = 3)
  }

  # for non-concatenated
  review_recordings <- function() {
    children <- as.character(tkwinfo('children', plotPanel))
    if (length(children)>0) sapply(children, function(ch) tcl("destroy", ch))

    if (!exists('abf_analysis_result', envir=.GlobalEnv)) {
      tkmessageBox(message="No analysis result available for review.")
      return()
    }
    result <- get('abf_analysis_result', envir=.GlobalEnv)
    datasets <- result$datasets
    traces2average <<- vector('list', length(datasets))
    for (i in seq_along(traces2average)) traces2average[[i]] <<- integer(0)
    current_dataset <<- 1L
    current_trace   <<- 1L

    reviewFrame <<- tkframe(plotPanel)
    tkgrid(reviewFrame, row=0, column=0, sticky='nsew')

    infoLabel <<- tklabel(reviewFrame, text = paste(names(datasets)[1],'trace 1'))
    tkgrid(infoLabel, row=0, column=0, columnspan=2)

    reviewPlot <<- tkrplot(reviewFrame, fun = function() {
      ds    <- datasets[[current_dataset]]
      fname <- names(datasets)[current_dataset]
      tkconfigure(infoLabel, text=paste(fname,'trace',current_trace))
      if (current_trace>length(ds$data)) {
        plot.new(); text(0.5,0.5,paste('No more recordings in',fname))
      } else {
        mat <- ds$data[[current_trace]]
        dc  <- as.numeric(tclvalue(dataColVar))
        if (is.na(dc)||dc<1||dc>ncol(mat)) dc <- 1
        dt  <- ds$samplingIntervalInSec*1000
        time<- seq(0,by=dt,length.out=nrow(mat))
        trace <- mat[,dc]
        par(cex.lab=0.6, cex.axis=0.6, cex.main=0.6)
        plot(time, trace, type='l', col='darkgray', xlab='time (ms)',
             xlim=smart_axis_limits(time), ylim=smart_axis_limits(trace),
             ylab=tclvalue(unitVar), axes=FALSE,
             main=paste(fname,'trace',current_trace), bty='l')
        axis(1); axis(2, las=1)
      }
    }, hscale=1, vscale=1)
    tkgrid(reviewPlot, row=1, column=0, columnspan=2)

    move_next <- function() {
      tkconfigure(acceptButton, state='normal', relief='raised')
      tkconfigure(rejectButton, state='normal', relief='raised')
      ds    <- datasets[[current_dataset]]
      fname <- names(datasets)[current_dataset]
      if (current_trace<length(ds$data)) {
        current_trace <<- current_trace + 1L
      } else {
        tkinsert(consoleText,'end',paste0(fname,' complete\n'))
        if (current_dataset<length(datasets)) {
          current_dataset <<- current_dataset+1L
          current_trace   <<- 1L
        } else {
          tkinsert(consoleText,'end','Review complete: Approved recordings stored.\n')
          return()
        }
      }
      tkrreplot(reviewPlot)
    }

    acceptButton <- tkbutton(reviewFrame, text='Accept', command=function() {
      traces2average[[current_dataset]] <<- c(traces2average[[current_dataset]], current_trace)
      move_next()
    })
    tkgrid(acceptButton, row=2, column=0)

    rejectButton <- tkbutton(reviewFrame, text='Reject', command=move_next)
    tkgrid(rejectButton, row=2, column=1)

    nextTraceButton <<- tkbutton(reviewFrame, text='Next Recording', command=move_next)
    tkgrid(nextTraceButton, row=3, column=0, columnspan=2)
  }

  # averaging Functions
  # function to average selected groups for concatenated mode.
  average_selected_groups <- function() {
    if (length(groups_list) == 0) {
      tkmessageBox(message = "No groups available for averaging. Please select groups first.")
      return()
    }

    dt_val <- master_abf$samplingIntervalInSec * 1000
    stim_time <- as.numeric(tclvalue(stimTimeVar))
    base_val <- as.numeric(tclvalue(baselineVar))
    data_column <- as.numeric(tclvalue(dataColVar))
    if (is.na(data_column) || data_column < 1) data_column <- 1

    baseline2zero <- function(y, dt, stim, baseline) {
      idx_baseline <- round(baseline / dt)
      idx_start    <- round((stim - baseline) / dt) + 1
      y0 <- y - mean(y[1:idx_baseline])
      y0[idx_start:length(y0)]
    }

    group_corrected_mean <- lapply(groups_list, function(group_indices) {
      traces_corrected <- lapply(group_indices, function(i) {
        trace <- master_abf$data[[i]][, data_column]
        baseline2zero(trace, dt = dt_val, stim = stim_time, baseline = base_val)
      })
      trace_mat <- do.call(cbind, traces_corrected)
      rowMeans(trace_mat)
    })

    averaged_data <<- group_corrected_mean

    # Ensure that this update is reflected properly when cycling
    current_avg_index <<- 1

    # Code for clearing existing UI elements
    children <- as.character(tkwinfo('children', plotPanel))
    for (child in children) {
      tryCatch({ tkdestroy(.Tk.ID[[child]]) }, error = function(e) {}, silent = TRUE)
    }

    # Create a new frame for displaying averages
    avgFrame <<- tkframe(plotPanel)
    tkgrid(avgFrame, row = 0, column = 0, sticky = 'nsew')

    # Function to draw the current average
    drawAvgPlot <- function() {
      num_groups <- length(group_corrected_mean)
      if (num_groups < 1) {
        plot.new()
        text(0.5, 0.5, 'No averaged data available')
        return()
      }

      dt_val <- master_abf$samplingIntervalInSec * 1000
      stim_time <- as.numeric(tclvalue(stimTimeVar))
      all_y <- unlist(group_corrected_mean)
      shared_ylim <- range(all_y)
      max_time <- max(sapply(group_corrected_mean, function(avg) dt_val * (length(avg) - 1)))
      shared_xlim <- c(0, max_time)

      par(mfrow = c(1, 1))
      cex <- 0.6
      par(cex.lab = cex, cex.axis = cex, cex.main = cex)

      avg_trace <- group_corrected_mean[[current_avg_index]]
      time <- seq(from = stim_time - base_val, by = dt_val, length.out = length(avg_trace))
      stim_y <- avg_trace[which.min(abs(time - stim_time))]

      egs_plot(x = time, y = avg_trace, color = 'darkgray', show_bar = TRUE, show_text = TRUE, 
        xbar = as.numeric(tclvalue(xbarVar)), ybar = as.numeric(tclvalue(ybarVar)), 
        xlim = shared_xlim, ylim = shared_ylim, cex = cex)
      points(stim_time, stim_y, pch = 8, col = 'black')
      text(stim_time, stim_y, labels = 'stim', pos = 3, cex = 0.6)
    }

    avgPlot <<- tkrplot(avgFrame, fun = drawAvgPlot, hscale = 1, vscale = 1)
    tkgrid(avgPlot, row = 0, column = 0)

    navFrame <- tkframe(avgFrame)
    tkgrid(navFrame, row = 1, column = 0)

    tkgrid(tklabel(navFrame, text = 'Average:'), row = 0, column = 0, padx = 5)
    avgLabel <- tklabel(navFrame, text = paste(current_avg_index, 'of', length(averaged_data)))
    tkgrid(avgLabel, row = 0, column = 1, padx = 5)

    nextAvgButton <- tkbutton(navFrame, text = 'Next', command = function(){
      if (current_avg_index < length(averaged_data)) {
        current_avg_index <<- current_avg_index + 1
      } else {
        current_avg_index <<- 1
      }
      tkconfigure(avgLabel, text = paste(current_avg_index, 'of', length(averaged_data)))
      tkrreplot(avgPlot)
    })
    tkgrid(nextAvgButton, row = 0, column = 2, padx = 5)

    tkdelete(consoleText, '1.0', 'end')
    tkinsert(consoleText, 'end', 'Averaging complete. Check the updated plot.')
  }

  # for separate mode (non-concatenated):
  averageApprovedTraces_sep <- function() {
    if(length(traces2average) == 0 || all(sapply(traces2average, length) == 0)){
      tkmessageBox(message = "No approved traces available. Please review recordings first.")
      return()
    }
    folderPath <- tclvalue(folderPathVar)
    if(nchar(folderPath) == 0){
      tkmessageBox(message = "Please select an ABF folder first.")
      return()
    }
    selIndices <- as.integer(tkcurselection(abfListBox))
    allFiles <- as.character(tkget(abfListBox, 0, 'end'))
    abf_files <- if(length(selIndices)==0) allFiles else allFiles[selIndices+1]
    if(length(abf_files)==0){
      tkmessageBox(message = "No ABF files selected.")
      return()
    }
    baseline <- as.numeric(tclvalue(baselineVar))
    stimTime <- as.numeric(tclvalue(stimTimeVar))
    xbar <- as.numeric(tclvalue(xbarVar))
    ybar <- as.numeric(tclvalue(ybarVar))
    
    result <- tryCatch({
      abf_out <<- abf_averages(
        datasets = abf_analysis_result$datasets,
        traces2average = traces2average,
        baseline = baseline,
        stimulation_time = stimTime,
        dataCol = as.numeric(tclvalue(dataColVar)),
        xlim = NULL, ylim = NULL,
        color = 'darkgray',
        xbar = xbar, ybar = ybar,
        width = 5.25, height = 2.75,
        plotIt = FALSE
      )
      abf_out
    }, error = function(e){
      tkmessageBox(message = paste('Error during averaging of approved traces:', e$message))
      NULL
    })
    if(!is.null(result)){
      tkdelete(consoleText, '1.0', 'end')
      msg <- sprintf("Averaging on approved traces complete.\nProcessed %d file(s).", length(abf_files))
      tkinsert(consoleText, 'end', paste0(msg, '\n'))
      abf_analysis_result <<- result
      
      averaged_data <<- result$baseline_corrected_mean_data
      datasets <- result$datasets
      current_avg_index <<- 1
      
      children <- as.character(tkwinfo('children', plotPanel))
      for(child in children){
        tryCatch({ tkdestroy(.Tk.ID[[child]]) }, error = function(e){}, silent = TRUE)
      }
      
      avgFrame <- tkframe(plotPanel)
      tkgrid(avgFrame, row = 0, column = 0, sticky = 'nsew')
      
      drawSingleAvg <- function(){
        if(length(averaged_data) == 0){
          plot.new()
          text(0.5, 0.5, 'No averaged data')
          return()
        }
        cex <- 0.6
        par(cex.lab = cex, cex.axis = cex, cex.main = cex)
        dt_val <- datasets[[current_avg_index]]$samplingIntervalInSec * 1000
        time <- seq(from = stimTime - baseline, by = dt_val, length.out = length(averaged_data[[current_avg_index]]))
        all_y <- unlist(averaged_data)
        ylim <- range(all_y)
        xlim <- c(0, max(sapply(averaged_data, function(y) length(y) * dt_val)))
        stim_time <- as.numeric(tclvalue(stimTimeVar))
        trace_y <- averaged_data[[current_avg_index]]
        stim_y <- trace_y[which.min(abs(time - stim_time))]

        egs_plot(x = time, y = trace_y, color = 'darkgray',
                 show_bar = TRUE, show_text = TRUE,
                 xbar = xbar, ybar = ybar,
                 xlim = xlim, ylim = ylim, cex = cex)

        points(stim_time, stim_y, pch = 8, col = 'black')
        text(stim_time, stim_y, labels = 'stim', pos = 3, cex = 0.6)
      }
      
      reviewPlot <<- tkrplot(avgFrame, fun = drawSingleAvg, hscale = 1, vscale = 1)
      tkgrid(reviewPlot, row = 0, column = 0, columnspan = 3)
      
      tkgrid(tklabel(avgFrame, text = 'Average:'), row = 1, column = 0)
      avgLabel <- tklabel(avgFrame, text = paste(current_avg_index, 'of', length(averaged_data)))
      tkgrid(avgLabel, row = 1, column = 1)
      
      nextAvgButton <- tkbutton(avgFrame, text = 'Next', command = function(){
        if(current_avg_index < length(averaged_data)){
          current_avg_index <<- current_avg_index + 1
        } else {
          current_avg_index <<- 1
        }
        tkconfigure(avgLabel, text = paste(current_avg_index, 'of', length(averaged_data)))
        tkrreplot(reviewPlot)
      })
      tkgrid(nextAvgButton, row = 1, column = 2)
    }
  }

  # UI Setup
  ABF_analysis_tk <- function() {
    tt <- tktoplevel()
    tkwm.title(tt, 'ABF Analysis')
    sidebarFrame <- tkframe(tt)
    mainFrame   <- tkframe(tt)
    tkgrid(sidebarFrame, row = 0, column = 0, sticky = 'ns')
    tkgrid(mainFrame,   row = 0, column = 1, sticky = 'nsew')
    tkgrid.rowconfigure(tt, 0, weight = 1)
    tkgrid.columnconfigure(tt, 1, weight = 1)

    # save mainFrame as the global plot panel
    plotPanel <<- mainFrame

    folderLabel <- tklabel(sidebarFrame, text = 'Select ABF Folder:')
    tkgrid(folderLabel, row = 0, column = 0, sticky = 'w')
    folderPathVar <<- tclVar('')
    folderEntry <- tkentry(sidebarFrame, textvariable = folderPathVar, width = 30)
    tkgrid(folderEntry, row = 0, column = 1, sticky = 'w')

    browseFolderButton <- tkbutton(sidebarFrame, text = 'Browse', command = function(){
      folderPath <- tclvalue(tkchooseDirectory())
      if(nchar(folderPath) > 0){
        tclvalue(folderPathVar) <<- folderPath
        abf_list <- list.files(path = folderPath, pattern = '\\.abf$', ignore.case = TRUE)
        if(length(abf_list) == 0){
          tkmessageBox(message = "No ABF files found in the selected folder.")
        } else {
          tkdelete(abfListBox, 0, 'end')
          for(f in abf_list){ tkinsert(abfListBox, 'end', f) }
          firstFilePath <- file.path(folderPath, abf_list[1])
          ds <- readABF(firstFilePath)
          dummy_result <- list(metadata = list(extract_metadata(ds)))
          updateAdditionalParams(dummy_result)
        }
      }
    })
    tkgrid(browseFolderButton, row = 0, column = 2, padx = 5)

    abfListLabel <- tklabel(sidebarFrame, text = 'ABF Files:')
    tkgrid(abfListLabel, row = 1, column = 0, sticky = 'w', pady = 5)
    
    tkgrid.columnconfigure(sidebarFrame, 0, weight = 1)
    tkgrid.columnconfigure(sidebarFrame, 1, weight = 1)
    tkgrid.columnconfigure(sidebarFrame, 2, weight = 1)
    abfListBox <<- tklistbox(sidebarFrame, height = 5, selectmode = 'multiple')
    tkgrid(abfListBox, row = 2, column = 0, columnspan = 3, sticky = 'we', padx = 10, pady = 5)

    paramFrame <- tkframe(sidebarFrame)
    tkgrid(paramFrame, row = 3, column = 0, columnspan = 3, sticky = 'w')
    experimentVar <<- tclVar('voltage clamp')
    unitVar       <<- tclVar('')
    dataColVar    <<- tclVar('')
    dtVar         <<- tclVar('')
    ntracesVar    <<- tclVar('')

    tkgrid(tklabel(paramFrame, text = 'Experiment:'), row = 0, column = 0, sticky = 'w')
    experimentCombo <- ttkcombobox(paramFrame, textvariable = experimentVar,
                                   values = c('voltage clamp', 'current clamp'), width = 15)
    tkgrid(experimentCombo, row = 0, column = 1, sticky = 'w')

    tkgrid(tklabel(paramFrame, text = 'Units:'), row = 1, column = 0, sticky = 'w')
    unitEntry <- tkentry(paramFrame, textvariable = unitVar, width = 10)
    tkgrid(unitEntry, row = 1, column = 1, sticky = 'w')

    tkgrid(tklabel(paramFrame, text = 'Data Column:'), row = 2, column = 0, sticky = 'w')
    dataColEntry <- tkentry(paramFrame, textvariable = dataColVar, width = 10)
    tkgrid(dataColEntry, row = 2, column = 1, sticky = 'w')
    tkbind(dataColEntry, '<FocusOut>', function() {
      dc <- as.numeric(tclvalue(dataColVar))
      if (!exists('abf_analysis_result', envir = .GlobalEnv)) return()
      cu <- abf_analysis_result$datasets[[1]]$channelUnits
      if (!is.na(dc) && dc >= 1 && dc <= length(cu)) {
        tclvalue(unitVar) <<- cu[dc]
      }
    })

    tkgrid(tklabel(paramFrame, text = 'dt (ms):'), row = 3, column = 0, sticky = 'w')
    dtEntry <- tkentry(paramFrame, textvariable = dtVar, width = 10)
    tkgrid(dtEntry, row = 3, column = 1, sticky = 'w')

    tkgrid(tklabel(paramFrame, text = '# traces:'), row = 4, column = 0, sticky = 'w')
    ntracesEntry <- tkentry(paramFrame, textvariable = ntracesVar, width = 10)
    tkgrid(ntracesEntry, row = 4, column = 1, sticky = 'w')

    baselineVar <<- tclVar('100')
    stimTimeVar <<- tclVar('150')
    xbarVar     <<- tclVar('100')
    ybarVar     <<- tclVar('50')
    concatMode  <<- tclVar('0')

    tkgrid(tklabel(sidebarFrame, text = 'Baseline:'), row = 4, column = 0, sticky = 'w')
    tkgrid(tkentry(sidebarFrame, textvariable = baselineVar, width = 10), row = 4, column = 1, sticky = 'w')
    tkgrid(tklabel(sidebarFrame, text = 'Stimulation Time:'), row = 5, column = 0, sticky = 'w')
    tkgrid(tkentry(sidebarFrame, textvariable = stimTimeVar, width = 10), row = 5, column = 1, sticky = 'w')
    tkgrid(tklabel(sidebarFrame, text = 'x-bar length:'), row = 6, column = 0, sticky = 'w')
    tkgrid(tkentry(sidebarFrame, textvariable = xbarVar, width = 10), row = 6, column = 1, sticky = 'w')
    tkgrid(tklabel(sidebarFrame, text = 'y-bar length:'), row = 7, column = 0, sticky = 'w')
    tkgrid(tkentry(sidebarFrame, textvariable = ybarVar, width = 10), row = 7, column = 1, sticky = 'w')

    concatButton <- tkcheckbutton(sidebarFrame, variable = concatMode,
                                  text = 'Concatenate Imported ABFs')
    tkgrid(concatButton, row = 8, column = 0, columnspan = 3)

    tkgrid.columnconfigure(sidebarFrame, 0, weight = 1)
    tkgrid.columnconfigure(sidebarFrame, 1, weight = 1)
    tkgrid.columnconfigure(sidebarFrame, 2, weight = 1)

    consoleText <<- tktext(sidebarFrame, height = 5)
    tkgrid(consoleText, row = 9, column = 0, columnspan = 3, sticky = 'we', padx = 10, pady = 5)

    updateAdditionalParams <<- function(result) {
      if (!is.null(result) && length(result$metadata) >= 1) {
        meta1 <- result$metadata[[1]]
        tclvalue(dtVar) <<- as.character(meta1$samplingIntervalInSec * 1000)
        if (as.character(tclvalue(concatMode)) != '1') {
          if (!is.null(meta1$header$lActualEpisodes))
            tclvalue(ntracesVar) <<- as.character(meta1$header$lActualEpisodes)
          else
            tclvalue(ntracesVar) <<- 'N/A'
        }
        expType <- tclvalue(experimentVar)
        col_idx <- choose_data_column(meta1$channelUnits, expType)
        if (!is.na(col_idx)) {
          tclvalue(unitVar) <<- meta1$channelUnits[col_idx]
          tclvalue(dataColVar) <<- as.character(col_idx)
        } else {
          tclvalue(unitVar) <<- 'N/A'
          tclvalue(dataColVar) <<- 'N/A'
        }
      }
    }
    tkbind(experimentCombo, '<<ComboboxSelected>>', function() {
      if (exists('abf_analysis_result', envir = .GlobalEnv)) {
        updateAdditionalParams(get('abf_analysis_result', envir = .GlobalEnv))
      }
    })

    runAnalysis <<- function() {
        # clear previous right‑panel widgets (graphs, metadata, table)
        children <- as.character(tkwinfo('children', plotPanel))
        if (length(children) > 0) {
          sapply(children, function(ch) tcl("destroy", ch))
        }

        # ... now your existing load_abf_data logic follows ...
        folderPath <- tclvalue(folderPathVar)
        if (nchar(folderPath) == 0) {
          tkmessageBox(message = "Please select an ABF folder first.")
          return()
        }

      folderPath <- tclvalue(folderPathVar)
      if (nchar(folderPath) == 0) {
        tkmessageBox(message = "Please select an ABF folder first.")
        return()
      }
      selIndices <- as.integer(tkcurselection(abfListBox))
      allFiles    <- as.character(tkget(abfListBox, 0, 'end'))
      abf_files   <- if (length(selIndices) == 0) allFiles else allFiles[selIndices + 1]
      if (length(abf_files) == 0) {
        tkmessageBox(message = "No ABF files selected.")
        return()
      }
      result <- tryCatch({
        load_abf_data(abf_files = abf_files, abf_path = folderPath)
      }, error = function(e) {
        tkmessageBox(message = paste('Error during data loading:', e$message))
        NULL
      })
      if (!is.null(result)) {
        tkdelete(consoleText, '1.0', 'end')
        tkinsert(consoleText, 'end', paste0('Data loaded. Processed ', length(abf_files), ' file(s).\n'))
        assign('abf_analysis_result', result, envir = .GlobalEnv)
        updateAdditionalParams(result)

        # display metadata
        meta1 <- result$metadata[[1]]
        first <- result$datasets[[1]]$data[[1]]
        length_sweep <- nrow(first)
        metaText <- paste(
          paste0("Format version: ", meta1$formatVersion),
          paste0("Sampling interval: ", meta1$samplingIntervalInSec, " s"),
          paste0("Channel names: ", paste(meta1$channelNames, collapse = " ")),
          paste0("Channel units: ", paste(meta1$channelUnits, collapse = " ")),
          paste0("Number of sweeps: ", meta1$header$lActualEpisodes),
          paste0("Length of first sweep: ", length_sweep),
          paste0("Path: ", meta1$path),
          sep = "\n"
        )
        kids <- as.character(tkwinfo('children', plotPanel))
        for (k in kids) tryCatch(tkdestroy(.Tk.ID[[k]]), error = function(e) {}, silent = TRUE)
        metaFrame  <- tkframe(plotPanel)
        tkgrid(metaFrame, row = 0, column = 0, sticky = 'w', pady = 2)
        metaLabel  <- tklabel(metaFrame, text = metaText, justify = 'left')
        tkgrid(metaLabel)

        # display first 10 rows of first trace
        out <- first[1:10, ]
        colnames(out) <- meta1$channelUnits
        rownames(out) <- seq(nrow(out))
        tableFrame  <- tkframe(plotPanel)
        tkgrid(tableFrame, row = 1, column = 0, sticky = 'nsew')
        textWidget  <<- tktext(tableFrame, width = 50, height = 11, wrap = 'none')
        tkgrid(textWidget, row = 0, column = 0)
        for (line in capture.output(print(out))) {
          tkinsert(textWidget, 'end', paste0(line, '\n'))
        }

        cons_msg <- check_consistency(result$metadata)
        if (cons_msg == 'Data is consistent') {
          tkinsert(consoleText, 'end', paste0(cons_msg, '\n'))
        } else {
          tkinsert(consoleText, 'end', paste0('ERROR: ', cons_msg, '\n'))
        }

        if (as.character(tclvalue(concatMode)) == '1') {
          master_abf <<- combine_abf_data(result)
          tclvalue(ntracesVar) <<- as.character(length(master_abf$data))
        } else {
          master_abf <<- result
        }
        tkconfigure(runAnalysisButton, text = 'Load Data')
      }
    }

    runAnalysisButton        <<- tkbutton(sidebarFrame, text = 'Load Data',               command = runAnalysis)
    reviewButton             <<- tkbutton(sidebarFrame, text = 'Review Recordings',        command = function() {
                                  if (as.character(tclvalue(concatMode)) == '1') review_master_recordings()
                                  else review_recordings()
                                })
    avgApprovedTracesButton  <<- tkbutton(sidebarFrame, text = 'Average Approved Traces', command = function() {
                                  if (as.character(tclvalue(concatMode)) == '1') average_selected_groups()
                                  else averageApprovedTraces_sep()
                                })
    tkDownloadBtn            <<- tkbutton(sidebarFrame, text = 'Download Data',            command = download_data)

    tkgrid(runAnalysisButton,       row = 10, column = 0, columnspan = 3, pady = 5)
    tkgrid(reviewButton,            row = 11, column = 0, columnspan = 3, pady = 5)
    tkgrid(avgApprovedTracesButton, row = 12, column = 0, columnspan = 3, pady = 5)
    tkgrid(tkDownloadBtn,           row = 13, column = 0, columnspan = 3, pady = 5)

    tkfocus(tt)
  }

  # launch UI
  ABF_analysis_tk()
}

analyseABFtk()

############################################################################################
analyseABFshiny <- function() {

  extract_metadata <- function(abf_dataset) {
    list(
      path                  = abf_dataset$path,
      formatVersion         = abf_dataset$formatVersion,
      channelNames          = abf_dataset$channelNames,
      channelUnits          = abf_dataset$channelUnits,
      samplingIntervalInSec = abf_dataset$samplingIntervalInSec,
      header                = abf_dataset$header
    )
  }

  choose_data_column <- function(channelUnits, experiment) {
    if (experiment == 'voltage clamp') {
      idx <- grep('A', channelUnits, ignore.case = TRUE)
    } else if (experiment == 'current clamp') {
      idx <- grep('V', channelUnits, ignore.case = TRUE)
    } else {
      idx <- integer(0)
    }
    if (length(idx) > 0) return(idx[1])
    NA_integer_
  }

  combine_abf_data <- function(result) {
    master_abf <- list(data = list(), samplingIntervalInSec = result$datasets[[1]]$samplingIntervalInSec)
    for (i in seq_along(result$datasets)) {
      ds <- result$datasets[[i]]
      master_abf$data <- c(master_abf$data, ds$data)
    }
    master_abf
  }

  egs_plot <- function(x, y, sign = -1, xlim = NULL, ylim = NULL, lwd = 1,
                       show_text = FALSE, xbar = 100, ybar = 50,
                       color = 'darkgray', show_bar = FALSE, cex = 0.6) {
    if (is.null(ylim))
      ylim <- if (sign==1) c(0,max(y)) else c(-max(-y),0)
    if (is.null(xlim))
      xlim <- c(min(x), max(x))
    idx1 <- which.min(abs(x - xlim[1]))
    idx2 <- which.min(abs(x - xlim[2]))
    plot(x[idx1:idx2], y[idx1:idx2], type='l', col=color, xlim=xlim, ylim=ylim,
         axes=FALSE, xlab='time (ms)', ylab='')
    if (show_bar) {
      segments(max(x)-xbar, min(ylim), max(x), min(ylim))
      segments(max(x), min(ylim), max(x), min(ylim)+ybar)
      if (show_text) {
        text(max(x)-xbar/2, min(ylim)-0.05*diff(ylim), paste(xbar,'ms'))
        text(max(x)+0.02*diff(xlim), min(ylim)+ybar/2, paste(ybar,'pA'), srt=90)
      }
    }
  }

  ui <- fluidPage(
    titlePanel("ABF Analysis"),
    sidebarLayout(
      sidebarPanel(
        tabsetPanel(
          id = "sideTabs",
          tabPanel("Main",
            fileInput("abfFiles","Upload ABF Files", multiple=TRUE, accept=".abf"),
            checkboxInput("concatenate","Concatenate ABFs",FALSE),
            actionButton("load","Load Data"), br(), br(),
            actionButton("review","Review Recordings"), br(), br(),
            actionButton("accept","Accept"), actionButton("reject","Reject"),
            actionButton("nextReview","Next"), br(), br(),
            actionButton("addGroup","Add Selected Group"), br(), br(),
            actionButton("completeSel","Selection Complete"), br(), br(),
            actionButton("average","Average Approved Traces"), br(), br(),
            actionButton("nextAvg","Next Average"), br(), br(),
            downloadButton("downloadData","Download Averaged CSV")
          ),
          tabPanel("Settings",
            selectInput("experiment","Experiment", c("voltage clamp","current clamp")),
            uiOutput("columnUI"),
            verbatimTextOutput("unitsText"),
            numericInput("dt","dt (ms)", NA),
            numericInput("ntraces","# traces", NA),
            numericInput("baseline","Baseline (ms)", 100),
            numericInput("stimTime","Stimulation time (ms)", 150),
            numericInput("xbar","x-bar length (ms)", 100),
            numericInput("ybar","y-bar length (pA)", 50)
          )
        )
      ),
      mainPanel(
        tabsetPanel(
          id = "mainTabs",
          tabPanel("Metadata",
            verbatimTextOutput("metaText"),
            tableOutput("firstTable")
          ),
          tabPanel("Review",
            fluidRow(
              column(8,   # ~66% for the plot
                plotOutput("reviewPlot", height="600px", width="100%")
              ),
              column(4,   # ~33% for the console box
                wellPanel(
                  style="height:600px; overflow:auto; padding:10px;",
                  verbatimTextOutput("console")
                )
              )
            )
          ),

          tabPanel("Average",
            fluidRow(
              column(8,
                plotOutput("avgPlot", height="600px", width="100%")
              ),
              column(4,
                wellPanel(
                  style="height:600px; overflow:auto; padding:10px;",
                  verbatimTextOutput("avgInfo")
                )
              )
            )
          )
        )
      )
    )
  )

  server <- function(input, output, session) {
    vals <- reactiveValues(
      datasets=NULL, metadata=NULL,
      master=NULL, traces2avg=NULL,
      mode=NULL, ct=1, total=0,
      curGroup=NULL, groups=list(),
      avg=NULL, ca=1,
      log="", avgLog=""
    )

    baseline2zero <- function(y, dt, stim, baseline) {
      idx_baseline <- round(baseline / dt)
      idx_start    <- round((stim - baseline) / dt) + 1
      y0 <- y - mean(y[1:idx_baseline])
      y0[idx_start:length(y0)]
    }

    observeEvent(input$load, {
      req(input$abfFiles)
      paths <- input$abfFiles$datapath
      names(paths) <- input$abfFiles$name
      vals$datasets <- lapply(paths, readABF)
      vals$metadata <- lapply(vals$datasets, extract_metadata)
      m1 <- vals$metadata[[1]]
      updateNumericInput(session, "dt", value = m1$samplingIntervalInSec * 1000)
      updateNumericInput(session, "ntraces", value = m1$header$lActualEpisodes)
      sel <- choose_data_column(m1$channelUnits, input$experiment)
      updateSelectInput(session, "column", choices = seq_along(m1$channelUnits), selected = sel)
      vals$log <- ""
      updateTabsetPanel(session, "sideTabs", selected = "Settings")
    })

    observe({
      req(vals$datasets)    # only run once data are loaded

      total <- if (input$concatenate) {
        # sum traces across every ABF
        sum(vapply(vals$datasets, function(ds) length(ds$data), integer(1)))
      } else {
        # original sweep count in the first ABF
        vals$metadata[[1]]$header$lActualEpisodes
      }

      updateNumericInput(session, "ntraces", value = total)
    })

    output$columnUI <- renderUI({
      req(vals$metadata)
      selectInput("column","Data Column", seq_along(vals$metadata[[1]]$channelUnits))
    })
    output$unitsText <- renderText({
      req(input$column)
      paste0("Units: ", vals$metadata[[1]]$channelUnits[as.integer(input$column)])
    })
    output$metaText <- renderText({
      req(vals$metadata)
      m <- vals$metadata[[1]]; first <- vals$datasets[[1]]$data[[1]]
      c(
        paste0("Format version: ", m$formatVersion),
        paste0("Sampling interval: ", m$samplingIntervalInSec, " s"),
        paste0("Channel names: ", paste(m$channelNames, collapse=" ")),
        paste0("Channel units: ", paste(m$channelUnits, collapse=" ")),
        paste0("Number of sweeps: ", m$header$lActualEpisodes),
        paste0("Length of first sweep: ", nrow(first)),
        paste0("Path: ", m$path)
      )
    })
    output$firstTable <- renderTable({
      req(vals$datasets, input$column)
      out <- vals$datasets[[1]]$data[[1]][1:10, ]
      colnames(out) <- vals$metadata[[1]]$channelUnits
      head(out, 10)
    })

    # Review
    observeEvent(input$review, {
      req(vals$datasets, input$column)
      vals$mode <- if (input$concatenate) "concat" else "sep"
      if (vals$mode == "concat") {
        vals$master <- combine_abf_data(list(datasets = vals$datasets))
        vals$total  <- length(vals$master$data)
        vals$ct     <- 1
        vals$curGroup <- integer(0)
        vals$groups   <- list()
      } else {
        vals$curFile   <- 1; vals$ct <- 1
        vals$traces2avg <- vector("list", length(vals$datasets))
        for (i in seq_along(vals$traces2avg)) vals$traces2avg[[i]] <- integer(0)
      }
      vals$log <- ""
      updateTabsetPanel(session, "mainTabs", selected = "Review")
    })

    output$reviewPlot <- renderPlot({
      req(vals$mode, vals$datasets, input$column)
      par(mar = c(6, 4, 2, 2)) 
      colIdx <- as.integer(input$column)
      if (vals$mode == "concat") {
        mat <- vals$master$data[[vals$ct]]
        dt  <- vals$master$samplingIntervalInSec * 1000
        time <- seq(0, by = dt, length.out = nrow(mat))
        trace <- mat[, colIdx]
      } else {
        ds   <- vals$datasets[[vals$curFile]]
        mat  <- ds$data[[vals$ct]]
        dt   <- ds$samplingIntervalInSec * 1000
        time <- seq(0, by = dt, length.out = nrow(mat))
        trace<- mat[, colIdx]
      }
      egs_plot(time, trace, show_bar = FALSE)
      axis(1); axis(2)
      # usr <- par("usr")
      # y0  <- usr[3] + 0.05*(usr[4]-usr[3])
      # text(input$stimTime, y0, "*", col = "black", cex = 2.5)
    })

    observeEvent(input$accept, {
      if (is.null(vals$mode)) {
        showNotification(
          "Error: No review in progress. Click 'Review Recordings' first.",
          type = "error",
          duration = 5
        )
        return()
      }

      if (vals$mode == "concat") {
        # concatenated mode
        vals$curGroup <- union(vals$curGroup, vals$ct)
        vals$log <- paste0(vals$log, "Accepted trace ", vals$ct, "\n")

        isolate({ input$nextReview })
        if (vals$ct < vals$total) {
          vals$ct <- vals$ct + 1
        } else {
          vals$log <- paste0(vals$log, "Review complete for all traces.\n")
        }

      } else {
        # separate mode
        fidx <- vals$curFile
        ds   <- vals$datasets[[fidx]]
        fname <- names(vals$datasets)[fidx]

        # log the accept
        vals$traces2avg[[fidx]] <- union(vals$traces2avg[[fidx]], vals$ct)
        vals$log <- paste0(vals$log, "Accepted ", fname, " trace ", vals$ct, "\n")

        isolate({ input$nextReview })
        if (vals$ct < length(ds$data)) {
          vals$ct <- vals$ct + 1
        } else if (vals$curFile < length(vals$datasets)) {
          # finished this file (but not the last) → log and advance
          vals$log <- paste0(vals$log, fname, " complete\n")
          vals$curFile <- vals$curFile + 1
          vals$ct      <- 1
        } else {
          # last file → log complete and final message
          vals$log <- paste0(vals$log, fname, " complete\n")
          vals$log <- paste0(vals$log, "Review complete: Approved recordings stored.\n")
        }
      }
    })

    observeEvent(input$reject, {
      if (is.null(vals$mode)) {
        showNotification(
          "Error: No review in progress. Click 'Review Recordings' first.",
          type = "error",
          duration = 5
        )
        return()
      }

      if (vals$mode == "concat") {
        # concatenated mode
        vals$log <- paste0(vals$log, "Rejected trace ", vals$ct, "\n")

        isolate({ input$nextReview })
        if (vals$ct < vals$total) {
          vals$ct <- vals$ct + 1
        } else {
          vals$log <- paste0(vals$log, "Review complete for all traces.\n")
        }

      } else {
        # separate mode
        fidx <- vals$curFile
        ds   <- vals$datasets[[fidx]]
        fname <- names(vals$datasets)[fidx]

        # log the reject
        vals$log <- paste0(vals$log, "Rejected ", fname, " trace ", vals$ct, "\n")

        isolate({ input$nextReview })
        if (vals$ct < length(ds$data)) {
          vals$ct <- vals$ct + 1
        } else if (vals$curFile < length(vals$datasets)) {
          # finished this file → log and advance
          vals$log <- paste0(vals$log, fname, " complete\n")
          vals$curFile <- vals$curFile + 1
          vals$ct      <- 1
        } else {
          # last file → log complete and final message
          vals$log <- paste0(vals$log, fname, " complete\n")
          vals$log <- paste0(vals$log, "Review complete: Approved recordings stored.\n")
        }
      }
    })

    output$console <- renderText(vals$log)

    observeEvent(input$addGroup, {
      if (length(vals$curGroup) == 0) {
        vals$log <- paste0(vals$log, "No traces selected in current group.\n")
      } else {
        vals$groups[[length(vals$groups) + 1]] <- vals$curGroup
        vals$log <- paste0(vals$log,
                           "Group ", length(vals$groups),
                           " selected: ",
                           paste(vals$curGroup, collapse = ","),
                           "\n")
        vals$curGroup <- integer(0)
      }
    })
    observeEvent(input$completeSel, {
      vals$log <- paste0(vals$log, "Review complete: Approved traces stored.\n")
    })

    # average traces
    observeEvent(input$average, {
      req(vals$mode)
      updateTabsetPanel(session, "mainTabs", selected = "Average")

      # common dt in ms
      dt <- if (vals$mode == "concat") {
        vals$master$samplingIntervalInSec * 1000
      } else {
        vals$datasets[[1]]$samplingIntervalInSec * 1000
      }

      if (vals$mode == "concat") {
        if (length(vals$groups) == 0) {
          vals$avgLog <- "No groups to average.\n"
          return()
        }
        vals$avg <- lapply(vals$groups, function(gr) {
          # combine accepted traces, compute mean
          y_full <- rowMeans(
            do.call(cbind, lapply(gr, function(i)
              vals$master$data[[i]][, as.integer(input$column)]))
          )
          baseline2zero(y_full, dt, input$stimTime, input$baseline)
        })

      } else {
        # separate files
        bc2 <- mapply(function(ds, idxs) {
          if (length(idxs) == 0) return(NULL)
          y_full <- rowMeans(
            do.call(cbind, lapply(idxs, function(i)
              ds$data[[i]][, as.integer(input$column)]))
          )
          baseline2zero(y_full, dt, input$stimTime, input$baseline)
        }, vals$datasets, vals$traces2avg,
        SIMPLIFY = FALSE)

        vals$avg <- bc2[!sapply(bc2, is.null)]
      }

      vals$ca     <- 1
      vals$avgLog <- "Averaging complete.\n"
    })

    observeEvent(input$nextAvg, {
      req(vals$avg)
      n <- length(vals$avg)
      vals$ca <- if (vals$ca < n) vals$ca + 1 else 1
    })

    output$avgPlot <- renderPlot({
      req(vals$avg)
      par(mar = c(2, 2, 1, 1))

      # data
      y   <- vals$avg[[vals$ca]]
      dt  <- vals$metadata[[1]]$samplingIntervalInSec * 1000
      time <- seq(0, by = dt, length.out = length(y))

      # plot trace
      egs_plot(time, y, show_bar  = FALSE, show_text = FALSE, color = 'darkgray')

      usr    <- par("usr")
      x_min  <- usr[1]; x_max <- usr[2]
      y_min  <- usr[3]; y_max <- usr[4]
      x_span <- x_max - x_min
      y_span <- y_max - y_min

      margin_x <- 0.05 * x_span
      margin_y <- 0.05 * y_span

      # bottom-right origin for bars
      x0 <- x_max - input$xbar - margin_x
      y0 <- y_min + margin_y

      # draw bars
      segments(x0, y0, x0 + input$xbar, y0, lwd = 1)           # xbar
      segments(x0, y0, x0, y0 + input$ybar, lwd = 1)           # ybar

      text(x0 + input$xbar/2,
           y0 - 0.03 * y_span,
           paste0(input$xbar, " ms"),
           adj = c(0.5, 1),
           cex = 1.2)

      text(
        x = x0 - margin_x/2,
        y = y0 + input$ybar/2,
        labels = paste0(input$ybar, " pA"),
        srt    = 90,
        adj    = c(0.5, 0.5),  # center in both directions
        cex    = 1.2)

      # stimulation marker at baseline
      text(input$baseline, 0, "*", col = "black", cex = 2.5)
      text(input$baseline, 0, labels = "stim", pos = 3, cex = 1)
    })

    output$avgInfo <- renderText({
      req(vals$avg)
      paste0("Average ", vals$ca, " of ", length(vals$avg), "\n", vals$avgLog)
    })

    # download as csv
    output$downloadData <- downloadHandler(
      filename = function() "averaged_data.csv",
      content = function(file) {
        # combine only the averaged traces (no time column)
        df <- as.data.frame(do.call(cbind, lapply(vals$avg, as.vector)))
        colnames(df) <- as.character(seq_along(vals$avg))
        write.csv(df, file, row.names = FALSE)
      }
    )
  }

  shinyApp(ui, server)

}

analyseABFshiny()

