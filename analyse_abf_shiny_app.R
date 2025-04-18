# this version combines versions A and B below:
# Remove all objects from the environment
rm(list = ls(all = TRUE))

load_required_packages <- function(packages) {
  new.packages <- packages[!(packages %in% installed.packages()[, 'Package'])]
  if (length(new.packages)) install.packages(new.packages)
  invisible(lapply(packages, library, character.only = TRUE))
}
required.packages <- c('robustbase', 'minpack.lm', 'Rcpp', 'signal', 'readABF', 'tcltk', 'tkrplot')
load_required_packages(required.packages)

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
  plot(x[idx1:idx2], y[idx1:idx2], type = 'l', col = color,
       xlim = xlim, ylim = ylim, bty = 'n', lwd = lwd, lty = 1,
       axes = FALSE, frame = FALSE, xlab = '', ylab = '')
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
  # do.call(cbind, averaged_data)

  if (is.null(averaged_data) || length(averaged_data) == 0) {
    tkmessageBox(message = "No averaged data available.")
    return()
  }
  # Determine source of dt_val
  dt_val <- if (!is.null(master_abf)) {
    master_abf$samplingIntervalInSec * 1000
  } else if (exists('abf_analysis_result', envir = .GlobalEnv)) {
    result <- get('abf_analysis_result', envir = .GlobalEnv)
    result$datasets[[1]]$samplingIntervalInSec * 1000
  } else {
    stop("No valid data source found for dt")
  }
  download_folder <- tclvalue(folderPathVar)
  if (nchar(download_folder) == 0) {
    tkmessageBox(message = "No folder selected for download.")
    return()
  }
  file_path <- file.path(download_folder, 'averaged_data.csv')
  write.csv(df, file = file_path, row.names = FALSE)
  tkmessageBox(message = paste('Averaged data saved to', file_path))
}

abf_averages <- function(datasets, baseline = 100, stimulation_time = 350, traces2average = NULL, dataCol = 1, ylim = NULL, xlim = NULL, 
  color = '#4C77BB', xbar = 100, ybar = 50, width = 5.25, height = 2.75, save = FALSE, plotIt = TRUE) {
  
  N <- length(datasets)
  sampling_intervals <- sapply(datasets, function(ds) ds$samplingIntervalInSec * 1000)
  responses <- lapply(seq_len(N), function(iii) {
    sapply(seq_along(datasets[[iii]]$data), function(ii) {
      datasets[[iii]]$data[[ii]][, dataCol]
    })
  })
  names(responses) <- names(datasets)
  baseline2zero <- function(y, dt, stimulation_time, baseline) {
    idx1 <- (stimulation_time - baseline) / dt
    idx2 <- baseline / dt
    y1 <- y[idx1:length(y)]
    y1 <- y1 - mean(y1[1:idx2])
    y1 - mean(y1[1:idx2])
  }
  responses0 <- lapply(seq_len(N), function(iii) {
    sapply(seq_len(ncol(responses[[iii]])), function(jj) {
      baseline2zero(responses[[iii]][, jj],
                    dt = sampling_intervals[iii],
                    stimulation_time = stimulation_time,
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
  if(is.null(master_abf)) {
    tkmessageBox(message = "No master ABF data available. Please load data first.")
    return()
  }
  total_traces <<- length(master_abf$data)
  current_trace <<- 1
  current_group_selected <<- integer(0)
  groups_list <<- list()
  
  children <- as.character(tkwinfo('children', plotPanel))
  for(child in children) {
    tryCatch({ tkdestroy(.Tk.ID[[child]]) }, error = function(e) {}, silent = TRUE)
  }
  
  reviewFrame <<- tkframe(plotPanel)
  tkgrid(reviewFrame, row = 0, column = 0, sticky = 'nsew')
  
  infoLabel <<- tklabel(reviewFrame, text = paste('Trace', current_trace, 'of', total_traces))
  tkgrid(infoLabel, row = 0, column = 0, columnspan = 3)
  
  reviewPlot <<- tkrplot(reviewFrame, fun = function() {
    if(current_trace > total_traces){
      plot.new()
      text(0.5, 0.5, 'No more traces to review.')
    } else {
      cex <- 0.6
      par(cex.lab = cex, cex.axis = cex, cex.main = cex)  # Force text sizes before plot()
      trace_matrix <- master_abf$data[[current_trace]]
      dt_val <- master_abf$samplingIntervalInSec * 1000
      time <- seq(0, by = dt_val, length.out = nrow(trace_matrix))
      data_column <- as.numeric(tclvalue(dataColVar))
      if(is.na(data_column) || data_column < 1 || data_column > ncol(trace_matrix))
        data_column <- 1
      trace <- trace_matrix[, data_column]
      plot(time, trace, col = 'darkgray', xlab = 'time (ms)',
            xlim = smart_axis_limits(time), ylim = smart_axis_limits(trace),
            ylab = tclvalue(unitVar), type = 'l', bty = 'l',
           axes = FALSE, main = paste('trace', current_trace))
      axis(1, tcl = -0.2)
      axis(2, las = 1, tcl = -0.2)
    }
  }, hscale = 1, vscale = 1)
  tkgrid(reviewPlot, row = 1, column = 0, columnspan = 3)
  
  acceptButton <<- tkbutton(reviewFrame, text = 'Accept', command = function() {
    current_group_selected <<- c(current_group_selected, current_trace)
    tkconfigure(acceptButton, state = 'disabled')
    tkconfigure(rejectButton, state = 'normal')
  })
  tkgrid(acceptButton, row = 2, column = 0)
  
  rejectButton <<- tkbutton(reviewFrame, text = 'Reject', command = function() {
    tkconfigure(rejectButton, state = 'disabled')
    tkconfigure(acceptButton, state = 'normal')
  })
  tkgrid(rejectButton, row = 2, column = 1)
  
  nextTraceButton <<- tkbutton(reviewFrame, text = 'Next Trace', command = function() {
    tkconfigure(acceptButton, state = 'normal')
    tkconfigure(rejectButton, state = 'normal')
    if(current_trace < total_traces){
      current_trace <<- current_trace + 1
      tkconfigure(infoLabel, text = paste('Trace', current_trace, 'of', total_traces))
      tkrreplot(reviewPlot)
    } else {
      tkmessageBox(message = "Review complete for all traces.")
    }
  })
  tkgrid(nextTraceButton, row = 2, column = 2)
  
  averageGroupButton <<- tkbutton(reviewFrame, text = 'Add Selected Group', 
                                   command = function() {
    if(length(current_group_selected) == 0) {
      tkmessageBox(message = "No traces selected in current group.")
    } else {
      groups_list[[length(groups_list) + 1]] <<- current_group_selected
      tkmessageBox(message = paste('Group', length(groups_list), 'selected with traces:',
                                   paste(current_group_selected, collapse = ', ')))
      current_group_selected <<- integer(0)
    }
  })
  tkgrid(averageGroupButton, row = 3, column = 0, columnspan = 3)
  
  selectionCompleteButton <<- tkbutton(reviewFrame, text = 'Selection Complete', 
                                        command = function() {
    tkdelete(consoleText, '1.0', 'end')
    tkinsert(consoleText, 'end', 'Review complete: Approved traces stored.')
  })
  tkgrid(selectionCompleteButton, row = 4, column = 0, columnspan = 3)
}

# for non-concatenated
review_recordings <- function() {
  if(!exists('abf_analysis_result', envir = .GlobalEnv)){
    tkmessageBox(message = "No analysis result available for review.")
    return()
  }
  result <- get('abf_analysis_result', envir = .GlobalEnv)
  datasets <- result$datasets
  traces2average <<- vector('list', length = length(datasets))
  for(i in seq_along(datasets)){
    traces2average[[i]] <<- integer(0)
  }
  current_dataset <<- 1
  current_trace <<- 1
  children <- as.character(tkwinfo('children', plotPanel))
  for(child in children){
    tryCatch({ tkdestroy(.Tk.ID[[child]]) }, error = function(e){}, silent = TRUE)
  }
  reviewFrame <<- tkframe(plotPanel)
  tkgrid(reviewFrame, row = 0, column = 0, sticky = 'nsew')
  current_filename <- names(datasets)[current_dataset]
  infoLabel <<- tklabel(reviewFrame, text = paste(current_filename, 'trace', current_trace))
  tkgrid(infoLabel, row = 0, column = 0, columnspan = 2)
  reviewPlot <<- tkrplot(reviewFrame, fun = function() {
    ds <- datasets[[current_dataset]]
    if(current_trace > length(ds$data)){
      plot.new()
      text(0.5, 0.5, paste('No more recordings in', current_filename))
    } else {
      trace_matrix <- ds$data[[current_trace]]
      data_column <- as.numeric(tclvalue(dataColVar))
      if(is.na(data_column) || data_column < 1 || data_column > ncol(trace_matrix))
        data_column <- 1
      dt_val <- ds$samplingIntervalInSec * 1000
      time <- seq(0, by = dt_val, length.out = nrow(trace_matrix))
      trace <- trace_matrix[, data_column]
      
      cex <- 0.6
      par(cex.lab = cex, cex.axis = cex, cex.main = cex)
      plot(time, trace, col = 'darkgray', xlab = 'time (ms)',
           xlim = smart_axis_limits(time),
           ylim = smart_axis_limits(trace),
           ylab = tclvalue(unitVar), type = 'l', bty = 'l',
           axes = FALSE, main = paste(current_filename, 'trace', current_trace))
      axis(1, tcl = -0.2)
      axis(2, las = 1, tcl = -0.2)

    }
  }, hscale = 1, vscale = 1)
  tkgrid(reviewPlot, row = 1, column = 0, columnspan = 2)
  acceptButton <- tkbutton(reviewFrame, text = 'Accept', command = function() {
    traces2average[[current_dataset]] <<- c(traces2average[[current_dataset]], current_trace)
    tkconfigure(acceptButton, state = 'disabled', relief = 'sunken')
    tkconfigure(rejectButton, state = 'normal', relief = 'raised')
  })
  tkgrid(acceptButton, row = 2, column = 0)
  rejectButton <- tkbutton(reviewFrame, text = 'Reject', command = function() {
    tkconfigure(rejectButton, state = 'disabled', relief = 'sunken')
    tkconfigure(acceptButton, state = 'normal', relief = 'raised')
  })
  tkgrid(rejectButton, row = 2, column = 1)
  nextTraceButton <- tkbutton(reviewFrame, text = 'Next Recording', command = function() {
    tkconfigure(acceptButton, state = 'normal', relief = 'raised')
    tkconfigure(rejectButton, state = 'normal', relief = 'raised')
    ds <- datasets[[current_dataset]]
    numRecordings <- length(ds$data)
    if(current_trace < numRecordings){
      current_trace <<- current_trace + 1
    } else {
      tkmessageBox(message = paste('Finished reviewing', current_filename))
      if(current_dataset < length(datasets)){
        current_dataset <<- current_dataset + 1
        current_trace <<- 1
      } else {
        tkmessageBox(message = "Review complete. Approved recordings are in 'traces2average'.")
        return()
      }
    }
    current_filename <- names(datasets)[current_dataset]
    tkconfigure(infoLabel, text = paste(current_filename, 'trace', current_trace))
    tkrreplot(reviewPlot)
  })
  tkgrid(nextTraceButton, row = 3, column = 0, columnspan = 2)
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

  baseline2zero <- function(y, dt, stimulation_time, baseline) {
    idx1 <- (stimulation_time - baseline) / dt
    idx2 <- baseline / dt
    y1 <- y[idx1:length(y)]
    y1 <- y1 - mean(y1[1:idx2])
    y1 - mean(y1[1:idx2])
  }

  group_corrected_mean <- lapply(groups_list, function(group_indices) {
    traces_corrected <- lapply(group_indices, function(i) {
      trace <- master_abf$data[[i]][, data_column]
      baseline2zero(trace, dt = dt_val, stimulation_time = stim_time, baseline = base_val)
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

    egs_plot(x = time, y = avg_trace, color = 'darkgray',
             show_bar = TRUE, show_text = TRUE,
             xbar = as.numeric(tclvalue(xbarVar)),
             ybar = as.numeric(tclvalue(ybarVar)),
             xlim = shared_xlim, ylim = shared_ylim,
             cex = cex)
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
    if(current_avg_index < length(averaged_data)){
      current_avg_index <<- current_avg_index + 1
    } else {
      current_avg_index <<- 1
    }
    tkconfigure(avgLabel, text = paste(current_avg_index, 'of', length(averaged_data)))
    tkrreplot(reviewPlot)
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
  abfListBox <<- tklistbox(sidebarFrame, height = 5, width = 40, selectmode = 'multiple')
  tkgrid(abfListBox, row = 2, column = 0, columnspan = 3, pady = 5, sticky = 'we')

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

  consoleText <<- tktext(sidebarFrame, width = 40, height = 4)
  tkgrid(consoleText, row = 9, column = 0, columnspan = 3, pady = 5)

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
      tkinsert(consoleText, 'end', paste('Data loaded. Processed',
                                        length(abf_files), 'file(s).'))
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
      textWidget  <<- tktext(tableFrame, width = 50, height = 10, wrap = 'none')
      tkgrid(textWidget, row = 0, column = 0)
      for (line in capture.output(print(out))) {
        tkinsert(textWidget, 'end', paste0(line, '\n'))
      }

      cons_msg <- check_consistency(result$metadata)
      if (cons_msg == 'Data is consistent') {
        tkmessageBox(message = cons_msg)
        if (as.character(tclvalue(concatMode)) == '1') {
          master_abf <<- combine_abf_data(result)
          tclvalue(ntracesVar) <<- as.character(length(master_abf$data))
        } else {
          master_abf <<- result
        }
      } else {
        tkmessageBox(message = paste('ERROR:', cons_msg))
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
  downloadButton           <<- tkbutton(sidebarFrame, text = 'Download Data',            command = download_data)

  tkgrid(runAnalysisButton,       row = 10, column = 0, columnspan = 3, pady = 5)
  tkgrid(reviewButton,            row = 11, column = 0, columnspan = 3, pady = 5)
  tkgrid(avgApprovedTracesButton, row = 12, column = 0, columnspan = 3, pady = 5)
  tkgrid(downloadButton,          row = 13, column = 0, columnspan = 3, pady = 5)

  tkfocus(tt)
}

# launch UI
ABF_analysis_tk()


# # this version combines versions A and B below:
# # Remove all objects from the environment
# rm(list = ls(all = TRUE))

# load_required_packages <- function(packages) {
#   new.packages <- packages[!(packages %in% installed.packages()[, 'Package'])]
#   if (length(new.packages)) install.packages(new.packages)
#   invisible(lapply(packages, library, character.only = TRUE))
# }
# required.packages <- c('robustbase', 'minpack.lm', 'Rcpp', 'signal', 'readABF', 'tcltk', 'tkrplot')
# load_required_packages(required.packages)

# # helper Functions
# extract_metadata <- function(abf_dataset) {
#   list(
#     path                  = abf_dataset$path,
#     formatVersion         = abf_dataset$formatVersion,
#     channelNames          = abf_dataset$channelNames,
#     channelUnits          = abf_dataset$channelUnits,
#     samplingIntervalInSec = abf_dataset$samplingIntervalInSec,
#     header                = abf_dataset$header,
#     tags                  = abf_dataset$tags,
#     sections              = abf_dataset$sections
#   )
# }

# choose_data_column <- function(channelUnits, experiment) {
#   if (experiment == 'voltage clamp') {
#     idx <- grep('A', channelUnits, ignore.case = TRUE)
#   } else if (experiment == 'current clamp') {
#     idx <- grep('V', channelUnits, ignore.case = TRUE)
#   } else {
#     idx <- integer(0)
#   }
#   if (length(idx) > 0) return(idx[1])
#   else return(NA)
# }

# check_consistency <- function(metadata) {
#   dt_values <- sapply(metadata, function(meta) meta$samplingIntervalInSec * 1000)
#   traces_values <- sapply(metadata, function(meta) meta$header$lActualEpisodes)
#   expType <- tclvalue(experimentVar)
#   unit_values <- sapply(metadata, function(meta) {
#     col_idx <- choose_data_column(meta$channelUnits, expType)
#     if (!is.na(col_idx)) meta$channelUnits[col_idx] else NA_character_
#   })
#   dt_good <- (length(unique(dt_values)) == 1)
#   traces_good <- (length(unique(traces_values)) == 1)
#   unit_good <- (length(unique(unit_values)) == 1)
#   if (dt_good && traces_good && unit_good) {
#     return('Data is consistent')
#   } else {
#     error_msgs <- c()
#     if (!dt_good) error_msgs <- c(error_msgs, paste('Inconsistent dt values:', paste(dt_values, collapse = ', ')))
#     if (!unit_good) error_msgs <- c(error_msgs, paste('Inconsistent Units:', paste(unit_values, collapse = ', ')))
#     if (!traces_good) error_msgs <- c(error_msgs, paste('Inconsistent Traces:', paste(traces_values, collapse = ', ')))
#     return(paste(error_msgs, collapse = '; '))
#   }
# }

# egs_plot <- function(x, y, sign = -1, xlim = NULL, ylim = NULL, lwd = 1, show_text = FALSE, 
#                      xbar = 100, ybar = 50, color = '#4C77BB', show_bar = FALSE, cex = 0.6) {
#   if (is.null(ylim))
#     ylim <- if (sign == 1) c(0, max(y)) else c(-max(-y), 0)
#   if (is.null(xlim))
#     xlim <- c(min(x), max(x))
#   idx1 <- which.min(abs(x - xlim[1]))
#   idx2 <- which.min(abs(x - xlim[2]))
#   plot(x[idx1:idx2], y[idx1:idx2], type = 'l', col = color,
#        xlim = xlim, ylim = ylim, bty = 'n', lwd = lwd, lty = 1,
#        axes = FALSE, frame = FALSE, xlab = '', ylab = '')
#   if (show_bar) {
#     ybar_start <- min(ylim) + (max(ylim) - min(ylim)) / 20
#     x_start <- max(xlim) - xbar - 50
#     y_start <- ybar_start
#     x_end <- x_start + xbar
#     y_end <- y_start + ybar
#     segments(x_start, y_start, x_end, y_start, lwd = lwd, col = 'black')
#     segments(x_start, y_start, x_start, y_end, lwd = lwd, col = 'black')
#     if (show_text) {
#       text(x = (x_start + x_end) / 2, y = y_start - ybar / 20, labels = paste(xbar, 'ms'), 
#            adj = c(0.5, 1), cex = cex)
#       text(x = x_start - xbar / 4, y = (y_start + y_end) / 2, labels = paste(ybar, 'pA'), 
#            adj = c(0.5, 0.5), srt = 90, cex = cex)
#     }
#   }
# }

# load_abf_data <- function(abf_files = NULL, abf_path = NULL) {
#   abf_path <- if (is.null(abf_path)) getwd() else abf_path
#   setwd(abf_path)
#   N <- length(abf_files)
#   datasets <- lapply(seq_len(N), function(ii) readABF(abf_files[ii]))
#   names(datasets) <- abf_files
#   metadata <- lapply(datasets, extract_metadata)
#   return(list(datasets = datasets, metadata = metadata))
# }

# download_data <- function() {
#   # do.call(cbind, averaged_data)

#   if (is.null(averaged_data) || length(averaged_data) == 0) {
#     tkmessageBox(message = "No averaged data available.")
#     return()
#   }
#   # Determine source of dt_val
#   dt_val <- if (!is.null(master_abf)) {
#     master_abf$samplingIntervalInSec * 1000
#   } else if (exists('abf_analysis_result', envir = .GlobalEnv)) {
#     result <- get('abf_analysis_result', envir = .GlobalEnv)
#     result$datasets[[1]]$samplingIntervalInSec * 1000
#   } else {
#     stop("No valid data source found for dt")
#   }
#   download_folder <- tclvalue(folderPathVar)
#   if (nchar(download_folder) == 0) {
#     tkmessageBox(message = "No folder selected for download.")
#     return()
#   }
#   file_path <- file.path(download_folder, 'averaged_data.csv')
#   write.csv(df, file = file_path, row.names = FALSE)
#   tkmessageBox(message = paste('Averaged data saved to', file_path))
# }

# abf_averages <- function(datasets, baseline = 100, stimulation_time = 350, traces2average = NULL, dataCol = 1, ylim = NULL, xlim = NULL, 
#   color = '#4C77BB', xbar = 100, ybar = 50, width = 5.25, height = 2.75, save = FALSE, plotIt = TRUE) {
  
#   N <- length(datasets)
#   sampling_intervals <- sapply(datasets, function(ds) ds$samplingIntervalInSec * 1000)
#   responses <- lapply(seq_len(N), function(iii) {
#     sapply(seq_along(datasets[[iii]]$data), function(ii) {
#       datasets[[iii]]$data[[ii]][, dataCol]
#     })
#   })
#   names(responses) <- names(datasets)
#   baseline2zero <- function(y, dt, stimulation_time, baseline) {
#     idx1 <- (stimulation_time - baseline) / dt
#     idx2 <- baseline / dt
#     y1 <- y[idx1:length(y)]
#     y1 <- y1 - mean(y1[1:idx2])
#     y1 - mean(y1[1:idx2])
#   }
#   responses0 <- lapply(seq_len(N), function(iii) {
#     sapply(seq_len(ncol(responses[[iii]])), function(jj) {
#       baseline2zero(responses[[iii]][, jj],
#                     dt = sampling_intervals[iii],
#                     stimulation_time = stimulation_time,
#                     baseline = baseline)
#     })
#   })
#   names(responses0) <- names(responses)
#   responses0_mean <- if(is.null(traces2average)) {
#     lapply(seq_len(N), function(iii) apply(responses0[[iii]], 1, mean))
#   } else {
#     lapply(seq_len(N), function(iii)
#       apply(responses0[[iii]][, traces2average[[iii]], drop = FALSE], 1, mean))
#   }
#   time <- lapply(seq_len(N), function(iii) {
#     dt_val <- sampling_intervals[iii]
#     stim_time <- stimulation_time
#     base_val <- baseline
#     seq(from = stim_time - base_val, by = dt_val, length.out = length(responses0_mean[[iii]]))
#   })
#   if(plotIt){
#     par(mfrow = c(1, N))
#     show_bar <- rep(FALSE, N)
#     if (N > 0) show_bar[N] <- TRUE
#     for(ii in seq_len(N)) {
#       for(ii in seq_len(N)) {
#         egs_plot(x = time[[ii]], y = responses0_mean[[ii]], color = 'darkgray',
#                  show_bar = FALSE, show_text = FALSE)
#       }
#     }
#   }
#   return(list(raw_data = responses,
#               baseline_corrected_data = responses0,
#               baseline_corrected_mean_data = responses0_mean,
#               datasets = datasets))
# }

# combine_abf_data <- function(result) {
#   master_abf <- list()
#   master_abf$data <- list()
#   master_abf$source_files <- c()
#   master_abf$samplingIntervalInSec <- result$datasets[[1]]$samplingIntervalInSec
#   for(i in seq_along(result$datasets)) {
#     ds <- result$datasets[[i]]
#     n_traces <- length(ds$data)
#     master_abf$data <- c(master_abf$data, ds$data)
#     master_abf$source_files <- c(master_abf$source_files, rep(names(result$datasets)[i], n_traces))
#   }
#   return(master_abf)
# }

# smart_axis_limits <- function(vec, n_steps = 5) {
#   rng <- range(vec)
#   spread <- diff(rng)
  
#   # Pick a base that's a nice round number (1, 2, 5, 10, 20, 50, 100, etc.)
#   raw_step <- spread / n_steps
#   base <- 10^floor(log10(raw_step))
  
#   # Refine to nicer step (1, 2, or 5 × 10^n)
#   nice_steps <- c(1, 2, 5, 10)
#   best_step <- base * nice_steps[which.min(abs(nice_steps * base - raw_step))]
  
#   lower <- floor(rng[1] / best_step) * best_step
#   upper <- ceiling(rng[2] / best_step) * best_step
#   c(lower, upper)
# }

# # global variables
# master_abf <<- NULL      # Will hold either a concatenated master object or the original structure.
# averaged_data <<- NULL   # Will hold the averaged (baseline_corrected_mean) data.
# traces2average <<- list()  # Used in separate mode.
# # For concatenated (master) mode:
# current_trace <<- 1      
# total_traces <<- 0       
# current_group_selected <<- integer(0)  
# groups_list <<- list()  
# # For separate (non-concatenated) mode:
# current_dataset <<- 1  

# # review functions
# # for concatenated
# review_master_recordings <- function() {
#   if(is.null(master_abf)) {
#     tkmessageBox(message = "No master ABF data available. Please load data first.")
#     return()
#   }
#   total_traces <<- length(master_abf$data)
#   current_trace <<- 1
#   current_group_selected <<- integer(0)
#   groups_list <<- list()
  
#   children <- as.character(tkwinfo('children', plotPanel))
#   for(child in children) {
#     tryCatch({ tkdestroy(.Tk.ID[[child]]) }, error = function(e) {}, silent = TRUE)
#   }
  
#   reviewFrame <<- tkframe(plotPanel)
#   tkgrid(reviewFrame, row = 0, column = 0, sticky = 'nsew')
  
#   infoLabel <<- tklabel(reviewFrame, text = paste('Trace', current_trace, 'of', total_traces))
#   tkgrid(infoLabel, row = 0, column = 0, columnspan = 3)
  
#   reviewPlot <<- tkrplot(reviewFrame, fun = function() {
#     if(current_trace > total_traces){
#       plot.new()
#       text(0.5, 0.5, 'No more traces to review.')
#     } else {
#       cex <- 0.6
#       par(cex.lab = cex, cex.axis = cex, cex.main = cex)  # Force text sizes before plot()
#       trace_matrix <- master_abf$data[[current_trace]]
#       dt_val <- master_abf$samplingIntervalInSec * 1000
#       time <- seq(0, by = dt_val, length.out = nrow(trace_matrix))
#       data_column <- as.numeric(tclvalue(dataColVar))
#       if(is.na(data_column) || data_column < 1 || data_column > ncol(trace_matrix))
#         data_column <- 1
#       trace <- trace_matrix[, data_column]
#       plot(time, trace, col = 'darkgray', xlab = 'time (ms)',
#             xlim = smart_axis_limits(time), ylim = smart_axis_limits(trace),
#             ylab = tclvalue(unitVar), type = 'l', bty = 'l',
#            axes = FALSE, main = paste('trace', current_trace))
#       axis(1, tcl = -0.2)
#       axis(2, las = 1, tcl = -0.2)
#     }
#   }, hscale = 1, vscale = 1)
#   tkgrid(reviewPlot, row = 1, column = 0, columnspan = 3)
  
#   acceptButton <<- tkbutton(reviewFrame, text = 'Accept', command = function() {
#     current_group_selected <<- c(current_group_selected, current_trace)
#     tkconfigure(acceptButton, state = 'disabled')
#     tkconfigure(rejectButton, state = 'normal')
#   })
#   tkgrid(acceptButton, row = 2, column = 0)
  
#   rejectButton <<- tkbutton(reviewFrame, text = 'Reject', command = function() {
#     tkconfigure(rejectButton, state = 'disabled')
#     tkconfigure(acceptButton, state = 'normal')
#   })
#   tkgrid(rejectButton, row = 2, column = 1)
  
#   nextTraceButton <<- tkbutton(reviewFrame, text = 'Next Trace', command = function() {
#     tkconfigure(acceptButton, state = 'normal')
#     tkconfigure(rejectButton, state = 'normal')
#     if(current_trace < total_traces){
#       current_trace <<- current_trace + 1
#       tkconfigure(infoLabel, text = paste('Trace', current_trace, 'of', total_traces))
#       tkrreplot(reviewPlot)
#     } else {
#       tkmessageBox(message = "Review complete for all traces.")
#     }
#   })
#   tkgrid(nextTraceButton, row = 2, column = 2)
  
#   averageGroupButton <<- tkbutton(reviewFrame, text = 'Add Selected Group', 
#                                    command = function() {
#     if(length(current_group_selected) == 0) {
#       tkmessageBox(message = "No traces selected in current group.")
#     } else {
#       groups_list[[length(groups_list) + 1]] <<- current_group_selected
#       tkmessageBox(message = paste('Group', length(groups_list), 'selected with traces:',
#                                    paste(current_group_selected, collapse = ', ')))
#       current_group_selected <<- integer(0)
#     }
#   })
#   tkgrid(averageGroupButton, row = 3, column = 0, columnspan = 3)
  
#   selectionCompleteButton <<- tkbutton(reviewFrame, text = 'Selection Complete', 
#                                         command = function() {
#     tkdelete(consoleText, '1.0', 'end')
#     tkinsert(consoleText, 'end', 'Review complete: Approved traces stored.')
#   })
#   tkgrid(selectionCompleteButton, row = 4, column = 0, columnspan = 3)
# }

# # for non-concatenated
# review_recordings <- function() {
#   if(!exists('abf_analysis_result', envir = .GlobalEnv)){
#     tkmessageBox(message = "No analysis result available for review.")
#     return()
#   }
#   result <- get('abf_analysis_result', envir = .GlobalEnv)
#   datasets <- result$datasets
#   traces2average <<- vector('list', length = length(datasets))
#   for(i in seq_along(datasets)){
#     traces2average[[i]] <<- integer(0)
#   }
#   current_dataset <<- 1
#   current_trace <<- 1
#   children <- as.character(tkwinfo('children', plotPanel))
#   for(child in children){
#     tryCatch({ tkdestroy(.Tk.ID[[child]]) }, error = function(e){}, silent = TRUE)
#   }
#   reviewFrame <<- tkframe(plotPanel)
#   tkgrid(reviewFrame, row = 0, column = 0, sticky = 'nsew')
#   current_filename <- names(datasets)[current_dataset]
#   infoLabel <<- tklabel(reviewFrame, text = paste(current_filename, 'trace', current_trace))
#   tkgrid(infoLabel, row = 0, column = 0, columnspan = 2)
#   reviewPlot <<- tkrplot(reviewFrame, fun = function() {
#     ds <- datasets[[current_dataset]]
#     if(current_trace > length(ds$data)){
#       plot.new()
#       text(0.5, 0.5, paste('No more recordings in', current_filename))
#     } else {
#       trace_matrix <- ds$data[[current_trace]]
#       data_column <- as.numeric(tclvalue(dataColVar))
#       if(is.na(data_column) || data_column < 1 || data_column > ncol(trace_matrix))
#         data_column <- 1
#       dt_val <- ds$samplingIntervalInSec * 1000
#       time <- seq(0, by = dt_val, length.out = nrow(trace_matrix))
#       trace <- trace_matrix[, data_column]
      
#       cex <- 0.6
#       par(cex.lab = cex, cex.axis = cex, cex.main = cex)
#       plot(time, trace, col = 'darkgray', xlab = 'time (ms)',
#            xlim = smart_axis_limits(time),
#            ylim = smart_axis_limits(trace),
#            ylab = tclvalue(unitVar), type = 'l', bty = 'l',
#            axes = FALSE, main = paste(current_filename, 'trace', current_trace))
#       axis(1, tcl = -0.2)
#       axis(2, las = 1, tcl = -0.2)

#     }
#   }, hscale = 1, vscale = 1)
#   tkgrid(reviewPlot, row = 1, column = 0, columnspan = 2)
#   acceptButton <- tkbutton(reviewFrame, text = 'Accept', command = function() {
#     traces2average[[current_dataset]] <<- c(traces2average[[current_dataset]], current_trace)
#     tkconfigure(acceptButton, state = 'disabled', relief = 'sunken')
#     tkconfigure(rejectButton, state = 'normal', relief = 'raised')
#   })
#   tkgrid(acceptButton, row = 2, column = 0)
#   rejectButton <- tkbutton(reviewFrame, text = 'Reject', command = function() {
#     tkconfigure(rejectButton, state = 'disabled', relief = 'sunken')
#     tkconfigure(acceptButton, state = 'normal', relief = 'raised')
#   })
#   tkgrid(rejectButton, row = 2, column = 1)
#   nextTraceButton <- tkbutton(reviewFrame, text = 'Next Recording', command = function() {
#     tkconfigure(acceptButton, state = 'normal', relief = 'raised')
#     tkconfigure(rejectButton, state = 'normal', relief = 'raised')
#     ds <- datasets[[current_dataset]]
#     numRecordings <- length(ds$data)
#     if(current_trace < numRecordings){
#       current_trace <<- current_trace + 1
#     } else {
#       tkmessageBox(message = paste('Finished reviewing', current_filename))
#       if(current_dataset < length(datasets)){
#         current_dataset <<- current_dataset + 1
#         current_trace <<- 1
#       } else {
#         tkmessageBox(message = "Review complete. Approved recordings are in 'traces2average'.")
#         return()
#       }
#     }
#     current_filename <- names(datasets)[current_dataset]
#     tkconfigure(infoLabel, text = paste(current_filename, 'trace', current_trace))
#     tkrreplot(reviewPlot)
#   })
#   tkgrid(nextTraceButton, row = 3, column = 0, columnspan = 2)
# }


# # averaging Functions
# # function to average selected groups for concatenated mode.
# average_selected_groups <- function() {
#   if (length(groups_list) == 0) {
#     tkmessageBox(message = "No groups available for averaging. Please select groups first.")
#     return()
#   }

#   dt_val <- master_abf$samplingIntervalInSec * 1000
#   stim_time <- as.numeric(tclvalue(stimTimeVar))
#   base_val <- as.numeric(tclvalue(baselineVar))
#   data_column <- as.numeric(tclvalue(dataColVar))
#   if (is.na(data_column) || data_column < 1) data_column <- 1

#   baseline2zero <- function(y, dt, stimulation_time, baseline) {
#     idx1 <- (stimulation_time - baseline) / dt
#     idx2 <- baseline / dt
#     y1 <- y[idx1:length(y)]
#     y1 <- y1 - mean(y1[1:idx2])
#     y1 - mean(y1[1:idx2])
#   }

#   group_corrected_mean <- lapply(groups_list, function(group_indices) {
#     traces_corrected <- lapply(group_indices, function(i) {
#       trace <- master_abf$data[[i]][, data_column]
#       baseline2zero(trace, dt = dt_val, stimulation_time = stim_time, baseline = base_val)
#     })
#     trace_mat <- do.call(cbind, traces_corrected)
#     rowMeans(trace_mat)
#   })

#   averaged_data <<- group_corrected_mean

#   # Ensure that this update is reflected properly when cycling
#   current_avg_index <<- 1

#   # Code for clearing existing UI elements
#   children <- as.character(tkwinfo('children', plotPanel))
#   for (child in children) {
#     tryCatch({ tkdestroy(.Tk.ID[[child]]) }, error = function(e) {}, silent = TRUE)
#   }

#   # Create a new frame for displaying averages
#   avgFrame <<- tkframe(plotPanel)
#   tkgrid(avgFrame, row = 0, column = 0, sticky = 'nsew')

#   # Function to draw the current average
#   drawAvgPlot <- function() {
#     num_groups <- length(group_corrected_mean)
#     if (num_groups < 1) {
#       plot.new()
#       text(0.5, 0.5, 'No averaged data available')
#       return()
#     }

#     dt_val <- master_abf$samplingIntervalInSec * 1000
#     stim_time <- as.numeric(tclvalue(stimTimeVar))
#     all_y <- unlist(group_corrected_mean)
#     shared_ylim <- range(all_y)
#     max_time <- max(sapply(group_corrected_mean, function(avg) dt_val * (length(avg) - 1)))
#     shared_xlim <- c(0, max_time)

#     par(mfrow = c(1, 1))
#     cex <- 0.6
#     par(cex.lab = cex, cex.axis = cex, cex.main = cex)

#     avg_trace <- group_corrected_mean[[current_avg_index]]
#     time <- seq(from = stim_time - base_val, by = dt_val, length.out = length(avg_trace))
#     stim_y <- avg_trace[which.min(abs(time - stim_time))]

#     egs_plot(x = time, y = avg_trace, color = 'darkgray',
#              show_bar = TRUE, show_text = TRUE,
#              xbar = as.numeric(tclvalue(xbarVar)),
#              ybar = as.numeric(tclvalue(ybarVar)),
#              xlim = shared_xlim, ylim = shared_ylim,
#              cex = cex)
#     points(stim_time, stim_y, pch = 8, col = 'black')
#     text(stim_time, stim_y, labels = 'stim', pos = 3, cex = 0.6)
#   }

#   avgPlot <<- tkrplot(avgFrame, fun = drawAvgPlot, hscale = 1, vscale = 1)
#   tkgrid(avgPlot, row = 0, column = 0)

#   navFrame <- tkframe(avgFrame)
#   tkgrid(navFrame, row = 1, column = 0)

#   tkgrid(tklabel(navFrame, text = 'Average:'), row = 0, column = 0, padx = 5)
#   avgLabel <- tklabel(navFrame, text = paste(current_avg_index, 'of', length(averaged_data)))
#   tkgrid(avgLabel, row = 0, column = 1, padx = 5)

#   nextAvgButton <- tkbutton(navFrame, text = 'Next', command = function(){
#     if(current_avg_index < length(averaged_data)){
#       current_avg_index <<- current_avg_index + 1
#     } else {
#       current_avg_index <<- 1
#     }
#     tkconfigure(avgLabel, text = paste(current_avg_index, 'of', length(averaged_data)))
#     tkrreplot(reviewPlot)
#   })
#   tkgrid(nextAvgButton, row = 0, column = 2, padx = 5)

#   tkdelete(consoleText, '1.0', 'end')
#   tkinsert(consoleText, 'end', 'Averaging complete. Check the updated plot.')
# }

# # for separate mode (non-concatenated):
# averageApprovedTraces_sep <- function() {
#   if(length(traces2average) == 0 || all(sapply(traces2average, length) == 0)){
#     tkmessageBox(message = "No approved traces available. Please review recordings first.")
#     return()
#   }
#   folderPath <- tclvalue(folderPathVar)
#   if(nchar(folderPath) == 0){
#     tkmessageBox(message = "Please select an ABF folder first.")
#     return()
#   }
#   selIndices <- as.integer(tkcurselection(abfListBox))
#   allFiles <- as.character(tkget(abfListBox, 0, 'end'))
#   abf_files <- if(length(selIndices)==0) allFiles else allFiles[selIndices+1]
#   if(length(abf_files)==0){
#     tkmessageBox(message = "No ABF files selected.")
#     return()
#   }
#   baseline <- as.numeric(tclvalue(baselineVar))
#   stimTime <- as.numeric(tclvalue(stimTimeVar))
#   xbar <- as.numeric(tclvalue(xbarVar))
#   ybar <- as.numeric(tclvalue(ybarVar))
  
#   result <- tryCatch({
#     abf_out <<- abf_averages(
#       datasets = abf_analysis_result$datasets,
#       traces2average = traces2average,
#       baseline = baseline,
#       stimulation_time = stimTime,
#       dataCol = as.numeric(tclvalue(dataColVar)),
#       xlim = NULL, ylim = NULL,
#       color = 'darkgray',
#       xbar = xbar, ybar = ybar,
#       width = 5.25, height = 2.75,
#       plotIt = FALSE
#     )
#     abf_out
#   }, error = function(e){
#     tkmessageBox(message = paste('Error during averaging of approved traces:', e$message))
#     NULL
#   })
#   if(!is.null(result)){
#     tkdelete(consoleText, '1.0', 'end')
#     msg <- sprintf("Averaging on approved traces complete.\nProcessed %d file(s).", length(abf_files))
#     tkinsert(consoleText, 'end', paste0(msg, '\n'))
#     abf_analysis_result <<- result
    
#     averaged_data <<- result$baseline_corrected_mean_data
#     datasets <- result$datasets
#     current_avg_index <<- 1
    
#     children <- as.character(tkwinfo('children', plotPanel))
#     for(child in children){
#       tryCatch({ tkdestroy(.Tk.ID[[child]]) }, error = function(e){}, silent = TRUE)
#     }
    
#     avgFrame <- tkframe(plotPanel)
#     tkgrid(avgFrame, row = 0, column = 0, sticky = 'nsew')
    
#     drawSingleAvg <- function(){
#       if(length(averaged_data) == 0){
#         plot.new()
#         text(0.5, 0.5, 'No averaged data')
#         return()
#       }
#       cex <- 0.6
#       par(cex.lab = cex, cex.axis = cex, cex.main = cex)
#       dt_val <- datasets[[current_avg_index]]$samplingIntervalInSec * 1000
#       time <- seq(from = stimTime - baseline, by = dt_val, length.out = length(averaged_data[[current_avg_index]]))
#       all_y <- unlist(averaged_data)
#       ylim <- range(all_y)
#       xlim <- c(0, max(sapply(averaged_data, function(y) length(y) * dt_val)))
#       stim_time <- as.numeric(tclvalue(stimTimeVar))
#       trace_y <- averaged_data[[current_avg_index]]
#       stim_y <- trace_y[which.min(abs(time - stim_time))]

#       egs_plot(x = time, y = trace_y, color = 'darkgray',
#                show_bar = TRUE, show_text = TRUE,
#                xbar = xbar, ybar = ybar,
#                xlim = xlim, ylim = ylim, cex = cex)

#       points(stim_time, stim_y, pch = 8, col = 'black')
#       text(stim_time, stim_y, labels = 'stim', pos = 3, cex = 0.6)
#     }
    
#     reviewPlot <<- tkrplot(avgFrame, fun = drawSingleAvg, hscale = 1, vscale = 1)
#     tkgrid(reviewPlot, row = 0, column = 0, columnspan = 3)
    
#     tkgrid(tklabel(avgFrame, text = 'Average:'), row = 1, column = 0)
#     avgLabel <- tklabel(avgFrame, text = paste(current_avg_index, 'of', length(averaged_data)))
#     tkgrid(avgLabel, row = 1, column = 1)
    
#     nextAvgButton <- tkbutton(avgFrame, text = 'Next', command = function(){
#       if(current_avg_index < length(averaged_data)){
#         current_avg_index <<- current_avg_index + 1
#       } else {
#         current_avg_index <<- 1
#       }
#       tkconfigure(avgLabel, text = paste(current_avg_index, 'of', length(averaged_data)))
#       tkrreplot(reviewPlot)
#     })
#     tkgrid(nextAvgButton, row = 1, column = 2)
#   }
# }

# # UI Setup
# ABF_analysis_tk <- function() {
#   tt <- tktoplevel()
#   tkwm.title(tt, 'ABF Analysis')
#   sidebarFrame <- tkframe(tt)
#   mainFrame <- tkframe(tt)
#   tkgrid(sidebarFrame, row = 0, column = 0, sticky = 'ns')
#   tkgrid(mainFrame, row = 0, column = 1, sticky = 'nsew')
#   tkgrid.rowconfigure(tt, 0, weight = 1)
#   tkgrid.columnconfigure(tt, 1, weight = 1)
  
#   # save mainFrame as the global plot panel
#   plotPanel <<- mainFrame
  
#   folderLabel <- tklabel(sidebarFrame, text = 'Select ABF Folder:')
#   tkgrid(folderLabel, row = 0, column = 0, sticky = 'w')
#   folderPathVar <<- tclVar('')
#   folderEntry <- tkentry(sidebarFrame, textvariable = folderPathVar, width = 30)
#   tkgrid(folderEntry, row = 0, column = 1, sticky = 'w')
  
#   browseFolderButton <- tkbutton(sidebarFrame, text = 'Browse', command = function(){
#     folderPath <- tclvalue(tkchooseDirectory())
#     if(nchar(folderPath) > 0){
#       tclvalue(folderPathVar) <<- folderPath
#       abf_list <- list.files(path = folderPath, pattern = '\\.abf$', ignore.case = TRUE)
#       if(length(abf_list) == 0){
#         tkmessageBox(message = "No ABF files found in the selected folder.")
#       } else {
#         tkdelete(abfListBox, 0, 'end')
#         for(f in abf_list){ tkinsert(abfListBox, 'end', f) }
#         firstFilePath <- file.path(folderPath, abf_list[1])
#         ds <- readABF(firstFilePath)
#         dummy_result <- list(metadata = list(extract_metadata(ds)))
#         updateAdditionalParams(dummy_result)
#       }
#     }
#   })
#   tkgrid(browseFolderButton, row = 0, column = 2, padx = 5)
  
#   abfListLabel <- tklabel(sidebarFrame, text = 'ABF Files:')
#   tkgrid(abfListLabel, row = 1, column = 0, sticky = 'w', pady = 5)
#   abfListBox <<- tklistbox(sidebarFrame, height = 5, width = 40, selectmode = 'multiple')
#   tkgrid(abfListBox, row = 2, column = 0, columnspan = 3, pady = 5, sticky = 'we')
  
#   paramFrame <- tkframe(sidebarFrame)
#   tkgrid(paramFrame, row = 3, column = 0, columnspan = 3, sticky = 'w')
  
#   experimentVar <<- tclVar('voltage clamp')
#   unitVar <<- tclVar('')
#   dataColVar <<- tclVar('')
#   dtVar <<- tclVar('')
#   ntracesVar <<- tclVar('')
  
#   tkgrid(tklabel(paramFrame, text = 'Experiment:'), row = 0, column = 0, sticky = 'w')
#   experimentCombo <- ttkcombobox(paramFrame, textvariable = experimentVar, values = c('voltage clamp', 'current clamp'), width = 15)
#   tkgrid(experimentCombo, row = 0, column = 1, sticky = 'w')
  
#   tkgrid(tklabel(paramFrame, text = 'Units:'), row = 1, column = 0, sticky = 'w')
#   unitEntry <- tkentry(paramFrame, textvariable = unitVar, width = 10)
#   tkgrid(unitEntry, row = 1, column = 1, sticky = 'w')
  
#   tkgrid(tklabel(paramFrame, text = 'Data Column:'), row = 2, column = 0, sticky = 'w')
#   dataColEntry <- tkentry(paramFrame, textvariable = dataColVar, width = 10)
#   tkgrid(dataColEntry, row = 2, column = 1, sticky = 'w')
  
#   tkgrid(tklabel(paramFrame, text = 'dt (ms):'), row = 3, column = 0, sticky = 'w')
#   dtEntry <- tkentry(paramFrame, textvariable = dtVar, width = 10)
#   tkgrid(dtEntry, row = 3, column = 1, sticky = 'w')
  
#   tkgrid(tklabel(paramFrame, text = '# traces:'), row = 4, column = 0, sticky = 'w')
#   ntracesEntry <- tkentry(paramFrame, textvariable = ntracesVar, width = 10)
#   tkgrid(ntracesEntry, row = 4, column = 1, sticky = 'w')
  
#   baselineVar <<- tclVar('100')
#   stimTimeVar <<- tclVar('150')
#   xbarVar <<- tclVar('100')
#   ybarVar <<- tclVar('50')
  
#   tkgrid(tklabel(sidebarFrame, text = 'Baseline:'), row = 4, column = 0, sticky = 'w')
#   tkgrid(tkentry(sidebarFrame, textvariable = baselineVar, width = 10), row = 4, column = 1, sticky = 'w')
#   tkgrid(tklabel(sidebarFrame, text = 'Stimulation Time:'), row = 5, column = 0, sticky = 'w')
#   tkgrid(tkentry(sidebarFrame, textvariable = stimTimeVar, width = 10), row = 5, column = 1, sticky = 'w')
#   tkgrid(tklabel(sidebarFrame, text = 'x-bar length:'), row = 6, column = 0, sticky = 'w')
#   tkgrid(tkentry(sidebarFrame, textvariable = xbarVar, width = 10), row = 6, column = 1, sticky = 'w')
#   tkgrid(tklabel(sidebarFrame, text = 'y-bar length:'), row = 7, column = 0, sticky = 'w')
#   tkgrid(tkentry(sidebarFrame, textvariable = ybarVar, width = 10), row = 7, column = 1, sticky = 'w')
  
#   # Checkbutton for concatenation mode
#   concatMode <<- tclVar('0')   # '0' means separate (default); '1' means concatenate.
#   concatButton <- tkcheckbutton(sidebarFrame, variable = concatMode, text = 'Concatenate Imported ABFs')
#   tkgrid(concatButton, row = 8, column = 0, columnspan = 3)
  
#   consoleText <<- tktext(sidebarFrame, width = 40, height = 4)
#   tkgrid(consoleText, row = 9, column = 0, columnspan = 3, pady = 5)
  
#   updateAdditionalParams <<- function(result) {
#     if (!is.null(result) && length(result$metadata) >= 1) {
#       meta1 <- result$metadata[[1]]
#       tclvalue(dtVar) <<- as.character(meta1$samplingIntervalInSec * 1000)
#       if (as.character(tclvalue(concatMode)) != '1') {
#         if (!is.null(meta1$header$lActualEpisodes)) {
#           tclvalue(ntracesVar) <<- as.character(meta1$header$lActualEpisodes)
#         } else {
#           tclvalue(ntracesVar) <<- 'N/A'
#         }
#       }
#       expType <- tclvalue(experimentVar)
#       col_idx <- choose_data_column(meta1$channelUnits, expType)
#       if (!is.na(col_idx)) {
#         tclvalue(unitVar) <<- meta1$channelUnits[col_idx]
#         tclvalue(dataColVar) <<- as.character(col_idx)
#       } else {
#         tclvalue(unitVar) <<- 'N/A'
#         tclvalue(dataColVar) <<- 'N/A'
#       }
#     }
#   }
  
#   tkbind(experimentCombo, '<<ComboboxSelected>>', function() {
#     if (exists('abf_analysis_result', envir = .GlobalEnv)) {
#       result <- get('abf_analysis_result', envir = .GlobalEnv)
#       updateAdditionalParams(result)
#     }
#   })
  
#   runAnalysis <<- function() {
#     folderPath <- tclvalue(folderPathVar)
#     if (nchar(folderPath) == 0) {
#       tkmessageBox(message = "Please select an ABF folder first.")
#       return()
#     }
#     selIndices <- as.integer(tkcurselection(abfListBox))
#     allFiles <- as.character(tkget(abfListBox, 0, 'end'))
#     abf_files <- if (length(selIndices) == 0) allFiles else allFiles[selIndices + 1]
#     if (length(abf_files) == 0) {
#       tkmessageBox(message = "No ABF files selected.")
#       return()
#     }
#     result <- tryCatch({
#       load_abf_data(abf_files = abf_files, abf_path = folderPath)
#     }, error = function(e) {
#       tkmessageBox(message = paste('Error during data loading:', e$message))
#       return(NULL)
#     })
#     if (!is.null(result)) {
#       tkdelete(consoleText, '1.0', 'end')
#       tkinsert(consoleText, 'end', paste('Data loaded. Processed', length(abf_files), 'file(s).'))
#       assign('abf_analysis_result', result, envir = .GlobalEnv)
#       updateAdditionalParams(result)
      
#       cons_msg <- check_consistency(result$metadata)
#       if (cons_msg == 'Data is consistent') {
#         tkmessageBox(message = cons_msg)
#         # If concatenation is selected, combine; otherwise, use the original separate structure.
#         if (as.character(tclvalue(concatMode)) == '1') {
#           master_abf <<- combine_abf_data(result)
#         } else {
#           master_abf <<- result
#         }

#         if (as.character(tclvalue(concatMode)) == '1') {
#           tclvalue(ntracesVar) <<- as.character(length(master_abf$data))
#         }

#       } else {
#         tkmessageBox(message = paste('ERROR:', cons_msg))
#       }
#       tkconfigure(runAnalysisButton, text = 'Load Data')
#     }
#   }
  
#   runAnalysisButton <<- tkbutton(sidebarFrame, text = 'Load Data', command = runAnalysis)
#   tkgrid(runAnalysisButton, row = 10, column = 0, columnspan = 3, pady = 5)
  
#   reviewButton <<- tkbutton(sidebarFrame, text = 'Review Recordings', command = function() {
#     if (as.character(tclvalue(concatMode)) == '1') {
#       review_master_recordings()
#     } else {
#       review_recordings()
#     }
#   })
#   tkgrid(reviewButton, row = 11, column = 0, columnspan = 3, pady = 5)
  
#   avgApprovedTracesButton <<- tkbutton(sidebarFrame, text = 'Average Approved Traces', 
#                                         command = function() {
#     if (as.character(tclvalue(concatMode)) == '1') {
#       average_selected_groups()
#     } else {
#       averageApprovedTraces_sep()
#     }
#   })
#   tkgrid(avgApprovedTracesButton, row = 12, column = 0, columnspan = 3, pady = 5)
  
#   downloadButton <<- tkbutton(sidebarFrame, text = 'Download Data', command = download_data)
#   tkgrid(downloadButton, row = 13, column = 0, columnspan = 3, pady = 5)
  
#   tkfocus(tt)
# }

# # launch UI
# ABF_analysis_tk()

# # this version concatenates the imported abfs if files are consistent
# # Remove all objects from the environment
# rm(list = ls(all = TRUE))

# # -----------------------------------------------------------------------------
# # Package loading function and required packages:
# # -----------------------------------------------------------------------------
# load_required_packages <- function(packages) {
#   new.packages <- packages[!(packages %in% installed.packages()[, 'Package'])]
#   if (length(new.packages)) install.packages(new.packages)
#   invisible(lapply(packages, library, character.only = TRUE))
# }
# required.packages <- c('robustbase', 'minpack.lm', 'Rcpp', 'signal', 'readABF', 'tcltk', 'tkrplot')
# load_required_packages(required.packages)

# # -----------------------------------------------------------------------------
# # Helper Functions
# # -----------------------------------------------------------------------------

# extract_metadata <- function(abf_dataset) {
#   list(
#     path                  = abf_dataset$path,
#     formatVersion         = abf_dataset$formatVersion,
#     channelNames          = abf_dataset$channelNames,
#     channelUnits          = abf_dataset$channelUnits,
#     samplingIntervalInSec = abf_dataset$samplingIntervalInSec,
#     header                = abf_dataset$header,
#     tags                  = abf_dataset$tags,
#     sections              = abf_dataset$sections
#   )
# }

# choose_data_column <- function(channelUnits, experiment) {
#   if (experiment == 'Voltage Clamp') {
#     idx <- grep('A', channelUnits, ignore.case = TRUE)
#   } else if (experiment == 'Current Clamp') {
#     idx <- grep('V', channelUnits, ignore.case = TRUE)
#   } else {
#     idx <- integer(0)
#   }
#   if (length(idx) > 0) return(idx[1])
#   else return(NA)
# }

# check_consistency <- function(metadata) {
#   dt_values <- sapply(metadata, function(meta) meta$samplingIntervalInSec * 1000)
#   traces_values <- sapply(metadata, function(meta) meta$header$lActualEpisodes)
#   expType <- tclvalue(experimentVar)
#   unit_values <- sapply(metadata, function(meta) {
#     col_idx <- choose_data_column(meta$channelUnits, expType)
#     if (!is.na(col_idx)) meta$channelUnits[col_idx] else NA_character_
#   })
#   dt_good <- (length(unique(dt_values)) == 1)
#   traces_good <- (length(unique(traces_values)) == 1)
#   unit_good <- (length(unique(unit_values)) == 1)
#   if (dt_good && traces_good && unit_good) {
#     return('Data is consistent')
#   } else {
#     error_msgs <- c()
#     if (!dt_good) error_msgs <- c(error_msgs, paste('Inconsistent dt values:', paste(dt_values, collapse = ', ')))
#     if (!unit_good) error_msgs <- c(error_msgs, paste('Inconsistent Units:', paste(unit_values, collapse = ', ')))
#     if (!traces_good) error_msgs <- c(error_msgs, paste('Inconsistent Traces:', paste(traces_values, collapse = ', ')))
#     return(paste(error_msgs, collapse = '; '))
#   }
# }

# # Original plotting function (same as in your original code)
# egs_plot <- function(x, y, sign = -1, xlim = NULL, ylim = NULL, lwd = 1, show_text = FALSE, 
#                      xbar = 100, ybar = 50, color = '#4C77BB', show_bar = FALSE, cex = 0.6) {
#   if (is.null(ylim))
#     ylim <- if (sign == 1) c(0, max(y)) else c(-max(-y), 0)
#   if (is.null(xlim))
#     xlim <- c(min(x), max(x))
#   idx1 <- which.min(abs(x - xlim[1]))
#   idx2 <- which.min(abs(x - xlim[2]))
#   plot(x[idx1:idx2], y[idx1:idx2], type = 'l', col = color,
#        xlim = xlim, ylim = ylim, bty = 'n', lwd = lwd, lty = 1,
#        axes = FALSE, frame = FALSE, xlab = '', ylab = '')
#   if (show_bar) {
#     ybar_start <- min(ylim) + (max(ylim) - min(ylim)) / 20
#     x_start <- max(xlim) - xbar - 50
#     y_start <- ybar_start
#     x_end <- x_start + xbar
#     y_end <- y_start + ybar
#     segments(x_start, y_start, x_end, y_start, lwd = lwd, col = 'black')
#     segments(x_start, y_start, x_start, y_end, lwd = lwd, col = 'black')
#     if (show_text) {
#       text(x = (x_start + x_end) / 2, y = y_start - ybar / 20, 
#            labels = paste(xbar, 'ms'), adj = c(0.5, 1), cex = cex)
#       text(x = x_start - xbar / 4, y = (y_start + y_end) / 2,  
#            labels = paste(ybar, 'pA'), adj = c(0.5, 0.5), srt = 90, cex = cex)
#     }
#   }
# }

# load_abf_data <- function(abf_files = NULL, abf_path = NULL) {
#   abf_path <- if (is.null(abf_path)) getwd() else abf_path
#   setwd(abf_path)
#   N <- length(abf_files)
#   datasets <- lapply(seq_len(N), function(ii) readABF(abf_files[ii]))
#   names(datasets) <- abf_files
#   metadata <- lapply(datasets, extract_metadata)
#   return(list(datasets = datasets, metadata = metadata))
# }

# # -----------------------------------------------------------------------------
# # Download Data Function: Write out the baseline_corrected_mean_data stored
# # in the global averaged_data variable.
# # -----------------------------------------------------------------------------
# download_data <- function() {
#   if (is.null(averaged_data) || length(averaged_data) == 0) {
#     tkmessageBox(message = 'No averaged data available.')
#     return()
#   }
#   dt_val <- master_abf$samplingIntervalInSec * 1000
#   avg_length <- length(averaged_data[[1]])
#   time_vec <- seq(0, by = dt_val, length.out = avg_length)
#   avg_mat <- sapply(averaged_data, function(x) x)
#   df <- data.frame(Time = time_vec, avg_mat)
#   download_folder <- tclvalue(folderPathVar)
#   if(nchar(download_folder) == 0) {
#     tkmessageBox(message = 'No folder selected for download.')
#     return()
#   }
#   file_path <- file.path(download_folder, 'averaged_data.csv')
#   write.csv(df, file = file_path, row.names = FALSE)
#   tkmessageBox(message = paste('Averaged data saved to', file_path))
# }

# # -----------------------------------------------------------------------------
# # Combine ABF Files into a Master Object
# # -----------------------------------------------------------------------------
# combine_abf_data <- function(result) {
#   master_abf <- list()
#   master_abf$data <- list()
#   master_abf$source_files <- c()
#   master_abf$samplingIntervalInSec <- result$datasets[[1]]$samplingIntervalInSec
  
#   for (i in seq_along(result$datasets)) {
#     ds <- result$datasets[[i]]
#     n_traces <- length(ds$data)
#     master_abf$data <- c(master_abf$data, ds$data)
#     master_abf$source_files <- c(master_abf$source_files, rep(names(result$datasets)[i], n_traces))
#   }
#   return(master_abf)
# }

# # -----------------------------------------------------------------------------
# # Global Variables for Master ABF, Review Process, and Averaged Data:
# # -----------------------------------------------------------------------------
# master_abf <<- NULL      # Combined ABF data.
# averaged_data <<- NULL   # Will hold baseline_corrected_mean_data.
# current_trace <<- 1      # Currently displayed trace.
# total_traces <<- 0       # Total number of traces.
# current_group_selected <<- integer(0)  # Indices approved for the current group.
# groups_list <<- list()   # List of groups (each group is a vector of trace indices).

# # -----------------------------------------------------------------------------
# # Revised Review Function for Master ABF
# # -----------------------------------------------------------------------------
# review_master_recordings <- function() {
#   if (is.null(master_abf)) {
#     tkmessageBox(message = 'No master ABF data available for review. Please load ABF files first.')
#     return()
#   }
  
#   total_traces <<- length(master_abf$data)
#   current_trace <<- 1
#   current_group_selected <<- integer(0)
#   groups_list <<- list()
  
#   # Clear previous UI plot panel children.
#   children <- as.character(tkwinfo('children', plotPanel))
#   for (child in children) {
#     tryCatch({ tkdestroy(.Tk.ID[[child]]) }, error = function(e) {}, silent = TRUE)
#   }
  
#   reviewFrame <<- tkframe(plotPanel)
#   tkgrid(reviewFrame, row = 0, column = 0, sticky = 'nsew')
  
#   infoLabel <<- tklabel(reviewFrame, text = paste('Trace', current_trace, 'of', total_traces))
#   tkgrid(infoLabel, row = 0, column = 0, columnspan = 3)
  
#   reviewPlot <<- tkrplot(reviewFrame, fun = function() {
#     if (current_trace > total_traces) {
#       plot.new()
#       text(0.5, 0.5, 'No more traces to review.')
#     } else {
#       trace_matrix <- master_abf$data[[current_trace]]
#       dt_val <- master_abf$samplingIntervalInSec * 1000
#       time <- seq(0, by = dt_val, length.out = nrow(trace_matrix))
#       data_column <- as.numeric(tclvalue(dataColVar))
#       if (is.na(data_column) || data_column < 1 || data_column > ncol(trace_matrix))
#         data_column <- 1
#       trace <- trace_matrix[, data_column]
#       plot(time, trace, col = 'darkgray', xlab = 'Time (ms)',
#            ylab = tclvalue(unitVar), type = 'l', bty = 'l',
#            axes = FALSE, main = paste('Trace', current_trace))
#       axis(1)
#       axis(2, las = 1)
#     }
#   }, hscale = 1, vscale = 1)
#   tkgrid(reviewPlot, row = 1, column = 0, columnspan = 3)
  
#   acceptButton <<- tkbutton(reviewFrame, text = 'Accept', command = function() {
#     current_group_selected <<- c(current_group_selected, current_trace)
#     tkconfigure(acceptButton, state = 'disabled')
#     tkconfigure(rejectButton, state = 'normal')
#   })
#   tkgrid(acceptButton, row = 2, column = 0)
  
#   rejectButton <<- tkbutton(reviewFrame, text = 'Reject', command = function() {
#     tkconfigure(rejectButton, state = 'disabled')
#     tkconfigure(acceptButton, state = 'normal')
#   })
#   tkgrid(rejectButton, row = 2, column = 1)
  
#   nextTraceButton <<- tkbutton(reviewFrame, text = 'Next Trace', command = function() {
#     tkconfigure(acceptButton, state = 'normal')
#     tkconfigure(rejectButton, state = 'normal')
#     if (current_trace < total_traces) {
#       current_trace <<- current_trace + 1
#       tkconfigure(infoLabel, text = paste('Trace', current_trace, 'of', total_traces))
#       tkrreplot(reviewPlot)
#     } else {
#       tkmessageBox(message = 'Review complete for all traces.')
#     }
#   })
#   tkgrid(nextTraceButton, row = 2, column = 2)
  
#   # Button to add the current approved trace group.
#   averageGroupButton <<- tkbutton(reviewFrame, text = 'Average Selected Traces', 
#                                    command = function() {
#     if (length(current_group_selected) == 0) {
#       tkmessageBox(message = 'No traces selected in current group.')
#     } else {
#       groups_list[[length(groups_list) + 1]] <<- current_group_selected
#       tkmessageBox(message = paste('Group', length(groups_list), 
#                                    'selected with traces:', paste(current_group_selected, collapse = ', ')))
#       current_group_selected <<- integer(0)
#     }
#   })
#   tkgrid(averageGroupButton, row = 3, column = 0, columnspan = 3)
  
#   # Replace the averaging step with a 'Selection Complete' button.
#   selectionCompleteButton <<- tkbutton(reviewFrame, text = 'Selection Complete', 
#                                         command = function() {
#     tkmessageBox(message = 'Review complete. Approved traces have been stored.')
#     # Optionally, you could close the review window here.
#   })
#   tkgrid(selectionCompleteButton, row = 4, column = 0, columnspan = 3)
# }

# # -----------------------------------------------------------------------------
# # New Averaging Function: Compute baseline-corrected means (mimicking your original abf_averages)
# # -----------------------------------------------------------------------------
# average_selected_groups <- function() {
#   if (length(groups_list) == 0) {
#     tkmessageBox(message = 'No groups available for averaging. Please select groups first.')
#     return()
#   }
  
#   dt_val <- master_abf$samplingIntervalInSec * 1000
#   stim_time <- as.numeric(tclvalue(stimTimeVar))
#   base_val <- as.numeric(tclvalue(baselineVar))
#   data_column <- as.numeric(tclvalue(dataColVar))
#   if (is.na(data_column) || data_column < 1) data_column <- 1
  
#   # baseline2zero function same as your original.
#   baseline2zero <- function(y, dt, stimulation_time, baseline) {
#     idx1 <- (stimulation_time - baseline) / dt
#     idx2 <- baseline / dt
#     y1 <- y[idx1:length(y)]
#     y1 <- y1 - mean(y1[1:idx2])
#     y1 - mean(y1[1:idx2])
#   }
  
#   # For each group, extract the trace data, apply baseline correction and compute the row–wise mean.
#   group_corrected_mean <- lapply(groups_list, function(group_indices) {
#     traces_corrected <- lapply(group_indices, function(i) {
#       trace <- master_abf$data[[i]][, data_column]
#       baseline2zero(trace, dt = dt_val, stimulation_time = stim_time, baseline = base_val)
#     })
#     trace_mat <- do.call(cbind, traces_corrected)
#     rowMeans(trace_mat)
#   })
  
#   # Save baseline_corrected_mean_data globally.
#   averaged_data <<- group_corrected_mean
  
#   # Plot the averaged traces in the UI.
#   children <- as.character(tkwinfo('children', plotPanel))
#   for (child in children) {
#     tryCatch({ tkdestroy(.Tk.ID[[child]]) }, error = function(e) {}, silent = TRUE)
#   }
#   avgFrame <<- tkframe(plotPanel)
#   tkgrid(avgFrame, row = 0, column = 0, sticky = 'nsew')
  
#   drawAvgPlot <- function() {
#     num_groups <- length(group_corrected_mean)
#     if (num_groups < 1) {
#       plot.new()
#       text(0.5, 0.5, 'No averaged data available')
#       return()
#     }
#     all_y <- unlist(group_corrected_mean)
#     shared_ylim <- range(all_y)
#     max_time <- max(sapply(group_corrected_mean, function(avg) { dt_val * (length(avg) - 1) }))
#     shared_xlim <- c(0, max_time)
    
#     par(mfrow = c(1, num_groups))
#     for (i in seq_along(group_corrected_mean)) {
#       avg_trace <- group_corrected_mean[[i]]
#       time <- seq(0, by = dt_val, length.out = length(avg_trace))
#       show_bar <- (i == num_groups)
#       egs_plot(x = time, y = avg_trace, color = 'darkgray',
#                show_bar = show_bar, show_text = show_bar,
#                xbar = as.numeric(tclvalue(xbarVar)), ybar = as.numeric(tclvalue(ybarVar)),
#                xlim = shared_xlim, ylim = shared_ylim)
#     }
#   }
  
#   avgPlot <<- tkrplot(avgFrame, fun = drawAvgPlot, hscale = 1, vscale = 1)
#   tkgrid(avgPlot, row = 0, column = 0)
#   tkmessageBox(message = 'Averaging complete. Check the updated plot.')
# }

# # -----------------------------------------------------------------------------
# # Main UI Setup
# # -----------------------------------------------------------------------------
# ABF_analysis_tk_A <- function() {
#   tt <- tktoplevel()
#   tkwm.title(tt, 'ABF Analysis')
#   sidebarFrame <- tkframe(tt)
#   mainFrame <- tkframe(tt)
#   tkgrid(sidebarFrame, row = 0, column = 0, sticky = 'ns')
#   tkgrid(mainFrame, row = 0, column = 1, sticky = 'nsew')
#   tkgrid.rowconfigure(tt, 0, weight = 1)
#   tkgrid.columnconfigure(tt, 1, weight = 1)
  
#   # Save mainFrame as the global UI plot panel.
#   plotPanel <<- mainFrame
  
#   folderLabel <- tklabel(sidebarFrame, text = 'Select ABF Folder:')
#   tkgrid(folderLabel, row = 0, column = 0, sticky = 'w')
#   folderPathVar <<- tclVar('')
#   folderEntry <- tkentry(sidebarFrame, textvariable = folderPathVar, width = 30)
#   tkgrid(folderEntry, row = 0, column = 1, sticky = 'w')
  
#   browseFolderButton <- tkbutton(sidebarFrame, text = 'Browse', command = function() {
#     folderPath <- tclvalue(tkchooseDirectory())
#     if (nchar(folderPath) > 0) {
#       tclvalue(folderPathVar) <<- folderPath
#       abf_list <- list.files(path = folderPath, pattern = '\\.abf$', ignore.case = TRUE)
#       if (length(abf_list) == 0) {
#         tkmessageBox(message = 'No ABF files found in the selected folder.')
#       } else {
#         tkdelete(abfListBox, 0, 'end')
#         for (f in abf_list) { tkinsert(abfListBox, 'end', f) }
#         firstFilePath <- file.path(folderPath, abf_list[1])
#         ds <- readABF(firstFilePath)
#         dummy_result <- list(metadata = list(extract_metadata(ds)))
#         updateAdditionalParams(dummy_result)
#       }
#     }
#   })
#   tkgrid(browseFolderButton, row = 0, column = 2, padx = 5)
  
#   abfListLabel <- tklabel(sidebarFrame, text = 'ABF Files:')
#   tkgrid(abfListLabel, row = 1, column = 0, sticky = 'w', pady = 5)
#   abfListBox <<- tklistbox(sidebarFrame, height = 5, width = 40, selectmode = 'multiple')
#   tkgrid(abfListBox, row = 2, column = 0, columnspan = 3, pady = 5, sticky = 'we')
  
#   paramFrame <- tkframe(sidebarFrame)
#   tkgrid(paramFrame, row = 3, column = 0, columnspan = 3, sticky = 'w')
  
#   experimentVar <<- tclVar('Voltage Clamp')
#   unitVar <<- tclVar('')
#   dataColVar <<- tclVar('')
#   dtVar <<- tclVar('')
#   ntracesVar <<- tclVar('')
  
#   tkgrid(tklabel(paramFrame, text = 'Experiment:'), row = 0, column = 0, sticky = 'w')
#   experimentCombo <- ttkcombobox(paramFrame, textvariable = experimentVar, values = c('Voltage Clamp', 'Current Clamp'), width = 15)
#   tkgrid(experimentCombo, row = 0, column = 1, sticky = 'w')
  
#   tkgrid(tklabel(paramFrame, text = 'Units:'), row = 1, column = 0, sticky = 'w')
#   unitEntry <- tkentry(paramFrame, textvariable = unitVar, width = 10, state = 'readonly')
#   tkgrid(unitEntry, row = 1, column = 1, sticky = 'w')
  
#   tkgrid(tklabel(paramFrame, text = 'Data Column:'), row = 2, column = 0, sticky = 'w')
#   dataColEntry <- tkentry(paramFrame, textvariable = dataColVar, width = 10, state = 'readonly')
#   tkgrid(dataColEntry, row = 2, column = 1, sticky = 'w')
  
#   tkgrid(tklabel(paramFrame, text = 'dt (ms):'), row = 3, column = 0, sticky = 'w')
#   dtEntry <- tkentry(paramFrame, textvariable = dtVar, width = 10, state = 'readonly')
#   tkgrid(dtEntry, row = 3, column = 1, sticky = 'w')
  
#   tkgrid(tklabel(paramFrame, text = 'Traces:'), row = 4, column = 0, sticky = 'w')
#   ntracesEntry <- tkentry(paramFrame, textvariable = ntracesVar, width = 10, state = 'readonly')
#   tkgrid(ntracesEntry, row = 4, column = 1, sticky = 'w')
  
#   baselineVar <<- tclVar('100')
#   stimTimeVar <<- tclVar('350')
#   xbarVar <<- tclVar('100')
#   ybarVar <<- tclVar('50')
  
#   tkgrid(tklabel(sidebarFrame, text = 'Baseline:'), row = 4, column = 0, sticky = 'w')
#   tkgrid(tkentry(sidebarFrame, textvariable = baselineVar, width = 10), row = 4, column = 1, sticky = 'w')
#   tkgrid(tklabel(sidebarFrame, text = 'Stimulation Time:'), row = 5, column = 0, sticky = 'w')
#   tkgrid(tkentry(sidebarFrame, textvariable = stimTimeVar, width = 10), row = 5, column = 1, sticky = 'w')
#   tkgrid(tklabel(sidebarFrame, text = 'x-bar length:'), row = 6, column = 0, sticky = 'w')
#   tkgrid(tkentry(sidebarFrame, textvariable = xbarVar, width = 10), row = 6, column = 1, sticky = 'w')
#   tkgrid(tklabel(sidebarFrame, text = 'y-bar length:'), row = 7, column = 0, sticky = 'w')
#   tkgrid(tkentry(sidebarFrame, textvariable = ybarVar, width = 10), row = 7, column = 1, sticky = 'w')
  
#   consoleText <<- tktext(sidebarFrame, width = 40, height = 4)
#   tkgrid(consoleText, row = 8, column = 0, columnspan = 3, pady = 5)
  
#   updateAdditionalParams <<- function(result) {
#     if (!is.null(result) && length(result$metadata) >= 1) {
#       meta1 <- result$metadata[[1]]
#       tclvalue(dtVar) <<- as.character(meta1$samplingIntervalInSec * 1000)
#       if (!is.null(meta1$header$lActualEpisodes)) {
#         tclvalue(ntracesVar) <<- as.character(meta1$header$lActualEpisodes)
#       } else {
#         tclvalue(ntracesVar) <<- 'N/A'
#       }
#       expType <- tclvalue(experimentVar)
#       col_idx <- choose_data_column(meta1$channelUnits, expType)
#       if (!is.na(col_idx)) {
#         tclvalue(unitVar) <<- meta1$channelUnits[col_idx]
#         tclvalue(dataColVar) <<- as.character(col_idx)
#       } else {
#         tclvalue(unitVar) <<- 'N/A'
#         tclvalue(dataColVar) <<- 'N/A'
#       }
#     }
#   }
  
#   tkbind(experimentCombo, '<<ComboboxSelected>>', function() {
#     if (exists('abf_analysis_result', envir = .GlobalEnv)) {
#       result <- get('abf_analysis_result', envir = .GlobalEnv)
#       updateAdditionalParams(result)
#     }
#   })
  
#   runAnalysis <<- function() {
#     folderPath <- tclvalue(folderPathVar)
#     if (nchar(folderPath) == 0) {
#       tkmessageBox(message = 'Please select an ABF folder first.')
#       return()
#     }
#     selIndices <- as.integer(tkcurselection(abfListBox))
#     allFiles <- as.character(tkget(abfListBox, 0, 'end'))
#     abf_files <- if (length(selIndices) == 0) allFiles else allFiles[selIndices + 1]
#     if (length(abf_files) == 0) {
#       tkmessageBox(message = 'No ABF files selected.')
#       return()
#     }
#     result <- tryCatch({
#       load_abf_data(abf_files = abf_files, abf_path = folderPath)
#     }, error = function(e) {
#       tkmessageBox(message = paste('Error during data loading:', e$message))
#       return(NULL)
#     })
#     if (!is.null(result)) {
#       tkdelete(consoleText, '1.0', 'end')
#       tkinsert(consoleText, 'end', paste('Data loaded. Processed', length(abf_files), 'file(s).'))
#       assign('abf_analysis_result', result, envir = .GlobalEnv)
#       updateAdditionalParams(result)
#       cons_msg <- check_consistency(result$metadata)
#       if (cons_msg == 'Data is consistent') {
#         tkmessageBox(message = cons_msg)
#         master_abf <<- combine_abf_data(result)
#       } else {
#         tkmessageBox(message = paste('ERROR:', cons_msg))
#       }
#       tkconfigure(runAnalysisButton, text = 'Load Data')
#     }
#   }
  
#   runAnalysisButton <<- tkbutton(sidebarFrame, text = 'Load Data', command = runAnalysis)
#   tkgrid(runAnalysisButton, row = 9, column = 0, columnspan = 3, pady = 5)
  
#   reviewButton <<- tkbutton(sidebarFrame, text = 'Review Recordings', command = review_master_recordings)
#   tkgrid(reviewButton, row = 10, column = 0, columnspan = 3, pady = 5)
  
#   # New button below 'Review Recordings' for averaging approved traces.
#   avgApprovedTracesButton <<- tkbutton(sidebarFrame, text = 'Average Approved Traces', command = average_selected_groups)
#   tkgrid(avgApprovedTracesButton, row = 11, column = 0, columnspan = 3, pady = 5)
  
#   downloadButton <<- tkbutton(sidebarFrame, text = 'Download Data', command = download_data)
#   tkgrid(downloadButton, row = 12, column = 0, columnspan = 3, pady = 5)
  
#   tkfocus(tt)
# }

# # Launch the UI.
# ABF_analysis_tk_A()



# # This version loads as separated abfs and does not allow averaging between separated loaded abfs

# # Remove all objects from the environment
# rm(list = ls(all = TRUE))

# # -----------------------------------------------------------------------------
# # Package loading function and required packages:
# # -----------------------------------------------------------------------------
# load_required_packages <- function(packages) {
#   new.packages <- packages[!(packages %in% installed.packages()[, 'Package'])]
#   if (length(new.packages)) install.packages(new.packages)
#   invisible(lapply(packages, library, character.only = TRUE))
# }
# required.packages <- c('robustbase', 'minpack.lm', 'Rcpp', 'signal', 'readABF', 'tcltk', 'tkrplot')
# load_required_packages(required.packages)

# # -----------------------------------------------------------------------------
# # Helper Functions
# # -----------------------------------------------------------------------------

# extract_metadata <- function(abf_dataset) {
#   list(
#     path                  = abf_dataset$path,
#     formatVersion         = abf_dataset$formatVersion,
#     channelNames          = abf_dataset$channelNames,
#     channelUnits          = abf_dataset$channelUnits,
#     samplingIntervalInSec = abf_dataset$samplingIntervalInSec,
#     header                = abf_dataset$header,
#     tags                  = abf_dataset$tags,
#     sections              = abf_dataset$sections
#   )
# }

# choose_data_column <- function(channelUnits, experiment) {
#   if (experiment == 'Voltage Clamp') {
#     idx <- grep('A', channelUnits, ignore.case = TRUE)
#   } else if (experiment == 'Current Clamp') {
#     idx <- grep('V', channelUnits, ignore.case = TRUE)
#   } else {
#     idx <- integer(0)
#   }
#   if (length(idx) > 0) return(idx[1])
#   else return(NA)
# }

# check_consistency <- function(metadata) {
#   dt_values <- sapply(metadata, function(meta) meta$samplingIntervalInSec * 1000)
#   traces_values <- sapply(metadata, function(meta) meta$header$lActualEpisodes)
#   expType <- tclvalue(experimentVar)
#   unit_values <- sapply(metadata, function(meta) {
#     col_idx <- choose_data_column(meta$channelUnits, expType)
#     if (!is.na(col_idx)) meta$channelUnits[col_idx] else NA_character_
#   })
#   dt_good <- (length(unique(dt_values)) == 1)
#   traces_good <- (length(unique(traces_values)) == 1)
#   unit_good <- (length(unique(unit_values)) == 1)
#   if (dt_good && traces_good && unit_good) {
#     return('Data is consistent')
#   } else {
#     error_msgs <- c()
#     if (!dt_good) error_msgs <- c(error_msgs, paste('Inconsistent dt values:', paste(dt_values, collapse = ', ')))
#     if (!unit_good) error_msgs <- c(error_msgs, paste('Inconsistent Units:', paste(unit_values, collapse = ', ')))
#     if (!traces_good) error_msgs <- c(error_msgs, paste('Inconsistent Traces:', paste(traces_values, collapse = ', ')))
#     return(paste(error_msgs, collapse = '; '))
#   }
# }

# egs_plot <- function(x, y, sign = -1, xlim = NULL, ylim = NULL, lwd = 1, show_text = FALSE, 
#                      xbar = 100, ybar = 50, color = '#4C77BB', show_bar = FALSE, cex = 0.6) {
#   if (is.null(ylim))
#     ylim <- if (sign == 1) c(0, max(y)) else c(-max(-y), 0)
#   if (is.null(xlim))
#     xlim <- c(min(x), max(x))
#   idx1 <- which.min(abs(x - xlim[1]))
#   idx2 <- which.min(abs(x - xlim[2]))
#   plot(x[idx1:idx2], y[idx1:idx2], type = 'l', col = color,
#        xlim = xlim, ylim = ylim, bty = 'n', lwd = lwd, lty = 1,
#        axes = FALSE, frame = FALSE, xlab = '', ylab = '')
#   if (show_bar) {
#     ybar_start <- min(ylim) + (max(ylim) - min(ylim)) / 20
#     x_start <- max(xlim) - xbar - 50
#     y_start <- ybar_start
#     x_end <- x_start + xbar
#     y_end <- y_start + ybar
#     segments(x_start, y_start, x_end, y_start, lwd = lwd, col = 'black')
#     segments(x_start, y_start, x_start, y_end, lwd = lwd, col = 'black')
#     if (show_text) {
#       text(x = (x_start + x_end) / 2, y = y_start - ybar / 20, 
#            labels = paste(xbar, 'ms'), adj = c(0.5, 1), cex = cex)
#       text(x = x_start - xbar / 4, y = (y_start + y_end) / 2,  
#            labels = paste(ybar, 'pA'), adj = c(0.5, 0.5), srt = 90, cex = cex)
#     }
#   }
# }

# load_abf_data <- function(abf_files = NULL, abf_path = NULL) {
#   abf_path <- if (is.null(abf_path)) getwd() else abf_path
#   setwd(abf_path)
#   N <- length(abf_files)
#   datasets <- lapply(seq_len(N), function(ii) readABF(abf_files[ii]))
#   names(datasets) <- abf_files
#   metadata <- lapply(datasets, extract_metadata)
#   return(list(datasets = datasets, metadata = metadata))
# }

# abf_averages <- function(datasets, 
#                          baseline = 100, 
#                          stimulation_time = 350, 
#                          traces2average = NULL,
#                          dataCol = 1, 
#                          ylim = NULL, 
#                          xlim = NULL, 
#                          color = '#4C77BB', 
#                          xbar = 100, 
#                          ybar = 50, 
#                          width = 5.25, 
#                          height = 2.75, 
#                          save = FALSE) {
#   N <- length(datasets)
#   sampling_intervals <- sapply(datasets, function(ds) ds$samplingIntervalInSec * 1000)
#   responses <- lapply(seq_len(N), function(iii) {
#     sapply(seq_along(datasets[[iii]]$data), function(ii) {
#       datasets[[iii]]$data[[ii]][, dataCol]
#     })
#   })
#   names(responses) <- names(datasets)
#   baseline2zero <- function(y, dt, stimulation_time, baseline) {
#     idx1 <- (stimulation_time - baseline) / dt
#     idx2 <- baseline / dt
#     y1 <- y[idx1:length(y)]
#     y1 <- y1 - mean(y1[1:idx2])
#     y1 - mean(y1[1:idx2])
#   }
#   responses0 <- lapply(seq_len(N), function(iii) {
#     sapply(seq_len(ncol(responses[[iii]])), function(jj) {
#       baseline2zero(responses[[iii]][, jj],
#                     dt = sampling_intervals[iii],
#                     stimulation_time = stimulation_time,
#                     baseline = baseline)
#     })
#   })
#   names(responses0) <- names(responses)
#   responses0_mean <- if (is.null(traces2average)) {
#     lapply(seq_len(N), function(iii) apply(responses0[[iii]], 1, mean))
#   } else {
#     lapply(seq_len(N), function(iii)
#       apply(responses0[[iii]][, traces2average[[iii]], drop = FALSE], 1, mean))
#   }
#   time <- lapply(seq_len(N), function(iii) {
#     dt_val <- sampling_intervals[iii]
#     seq(0, by = dt_val, length.out = length(responses0_mean[[iii]]))
#   })
#   par(mfrow = c(1, N))
#   show_bar <- rep(FALSE, N)
#   if (N > 0) show_bar[N] <- TRUE
#   for (ii in seq_len(N)) {
#     for (ii in seq_len(N)) {
#       egs_plot(x = time[[ii]], y = responses0_mean[[ii]], color = 'darkgray',
#                show_bar = FALSE, show_text = FALSE)
#     }
#   }
#   if (save) {
#     warning('save_graph not implemented in this UI example')
#   }
#   return(list(raw_data = responses,
#               baseline_corrected_data = responses0,
#               baseline_corrected_mean_data = responses0_mean,
#               datasets = datasets))
# }

# traces2average <<- list()

# review_recordings <- function() {
#   if (!exists('abf_analysis_result', envir = .GlobalEnv)) {
#     tkmessageBox(message = 'No analysis result available for review.')
#     return()
#   }
#   result <- get('abf_analysis_result', envir = .GlobalEnv)
#   datasets <- result$datasets
#   traces2average <<- vector('list', length = length(datasets))
#   for (i in seq_along(datasets)) { 
#     traces2average[[i]] <<- integer(0) 
#   }
#   current_dataset <<- 1
#   current_trace <<- 1
#   children <- as.character(tkwinfo('children', plotPanel))
#   for (child in children) {
#     tryCatch({
#       tkdestroy(.Tk.ID[[child]])
#     }, error = function(e) {}, silent = TRUE)
#   }

#   reviewFrame <<- tkframe(plotPanel)
#   tkgrid(reviewFrame, row = 0, column = 0, sticky = 'nsew')
#   current_filename <- names(datasets)[current_dataset]
#   infoLabel <<- tklabel(reviewFrame, text = paste(current_filename, 'trace', current_trace))
#   tkgrid(infoLabel, row = 0, column = 0, columnspan = 2)
#   reviewPlot <<- tkrplot(reviewFrame, fun = function() {
#     ds <- datasets[[current_dataset]]
#     if (current_trace > length(ds$data)) {
#       plot.new()
#       text(0.5, 0.5, paste('No more recordings in', current_filename))
#     } else {
#       trace_matrix <- ds$data[[current_trace]]
#       data_column <- as.numeric(tclvalue(dataColVar))
#       if (is.na(data_column) || data_column < 1 || data_column > ncol(trace_matrix)) {
#         data_column <- 1
#       }
#       dt_val <- ds$samplingIntervalInSec * 1000
#       time <- seq(0, by = dt_val, length.out = nrow(trace_matrix))
#       trace <- trace_matrix[, data_column]
#       plot(time, trace, col = 'darkgray', xlab = 'Time (ms)', ylab = tclvalue(unitVar),
#            type = 'l', bty = 'l', axes = FALSE,
#            main = paste(current_filename, 'trace', current_trace))
#       axis(1)
#       axis(2, las = 1)
#     }
#   }, hscale = 1, vscale = 1)
#   tkgrid(reviewPlot, row = 1, column = 0, columnspan = 2)
#   acceptButton <- tkbutton(reviewFrame, text = 'Accept', command = function() {
#     traces2average[[current_dataset]] <<- c(traces2average[[current_dataset]], current_trace)
#     tkconfigure(acceptButton, state = 'disabled', relief = 'sunken')
#     tkconfigure(rejectButton, state = 'normal', relief = 'raised')
#   })
#   tkgrid(acceptButton, row = 2, column = 0)
#   rejectButton <- tkbutton(reviewFrame, text = 'Reject', command = function() {
#     tkconfigure(rejectButton, state = 'disabled', relief = 'sunken')
#     tkconfigure(acceptButton, state = 'normal', relief = 'raised')
#   })
#   tkgrid(rejectButton, row = 2, column = 1)
#   nextTraceButton <- tkbutton(reviewFrame, text = 'Next Recording', command = function() {
#     tkconfigure(acceptButton, state = 'normal', relief = 'raised')
#     tkconfigure(rejectButton, state = 'normal', relief = 'raised')
#     ds <- datasets[[current_dataset]]
#     numRecordings <- length(ds$data)
#     if (current_trace < numRecordings) {
#       current_trace <<- current_trace + 1
#     } else {
#       tkmessageBox(message = paste('Finished reviewing', current_filename))
#       if (current_dataset < length(datasets)) {
#         current_dataset <<- current_dataset + 1
#         current_trace <<- 1
#       } else {
#         tkmessageBox(message = 'Review complete. Approved recordings are in 'traces2average'.')
#         return()
#       }
#     }
#     current_filename <- names(datasets)[current_dataset]
#     tkconfigure(infoLabel, text = paste(current_filename, 'trace', current_trace))
#     tkrreplot(reviewPlot)
#   })
#   tkgrid(nextTraceButton, row = 3, column = 0, columnspan = 2)
# }

# download_data <- function() {
#   if (!exists('abf_analysis_result', envir = .GlobalEnv)) {
#     tkmessageBox(message = 'No analysis result available for download.')
#     return()
#   }
#   result <- get('abf_analysis_result', envir = .GlobalEnv)
#   abf_out <- result  # structure containing baseline_corrected_mean_data
#   dataList <- abf_out$baseline_corrected_mean_data
#   if (length(dataList) == 0) {
#     tkmessageBox(message = 'No averaged data available.')
#     return()
#   }
#   dt_val <- result$datasets[[1]]$samplingIntervalInSec * 1000
#   time_vec <- seq(0, by = dt_val, length.out = length(dataList[[1]]))
#   # Construct a data frame: first column is Time; next columns (one per ABF file) are the averaged traces,
#   # gathered by sapply as requested.
#   avg_mat <- sapply(1:length(dataList), function(ii) dataList[[ii]])
#   df <- data.frame(Time = time_vec, avg_mat)
#   download_folder <- tclvalue(folderPathVar)
#   if(nchar(download_folder) == 0) {
#     tkmessageBox(message = 'No folder selected for download.')
#     return()
#   }
#   file_path <- file.path(download_folder, 'averaged_data.csv')
#   write.csv(df, file = file_path, row.names = FALSE)
#   tkmessageBox(message = paste('Averaged data saved to', file_path))
# }

# ABF_analysis_tk_B <- function() {
#   tt <- tktoplevel()
#   tkwm.title(tt, 'ABF Analysis')
#   sidebarFrame <- tkframe(tt)
#   mainFrame <- tkframe(tt)
#   tkgrid(sidebarFrame, row = 0, column = 0, sticky = 'ns')
#   tkgrid(mainFrame, row = 0, column = 1, sticky = 'nsew')
#   tkgrid.rowconfigure(tt, 0, weight = 1)
#   tkgrid.columnconfigure(tt, 1, weight = 1)
  
#   # Save mainFrame as the global plot panel.
#   plotPanel <<- mainFrame
  
#   folderLabel <- tklabel(sidebarFrame, text = 'Select ABF Folder:')
#   tkgrid(folderLabel, row = 0, column = 0, sticky = 'w')
#   folderPathVar <<- tclVar('')
#   folderEntry <- tkentry(sidebarFrame, textvariable = folderPathVar, width = 30)
#   tkgrid(folderEntry, row = 0, column = 1, sticky = 'w')
  
#   browseFolderButton <- tkbutton(sidebarFrame, text = 'Browse', command = function() {
#     folderPath <- tclvalue(tkchooseDirectory())
#     if (nchar(folderPath) > 0) {
#       tclvalue(folderPathVar) <<- folderPath
#       abf_list <- list.files(path = folderPath, pattern = '\\.abf$', ignore.case = TRUE)
#       if (length(abf_list) == 0) {
#         tkmessageBox(message = 'No ABF files found in the selected folder.')
#       } else {
#         tkdelete(abfListBox, 0, 'end')
#         for (f in abf_list) { tkinsert(abfListBox, 'end', f) }
#         firstFilePath <- file.path(folderPath, abf_list[1])
#         ds <- readABF(firstFilePath)
#         dummy_result <- list(metadata = list(extract_metadata(ds)))
#         updateAdditionalParams(dummy_result)
#       }
#     }
#   })
#   tkgrid(browseFolderButton, row = 0, column = 2, padx = 5)
  
#   abfListLabel <- tklabel(sidebarFrame, text = 'ABF Files:')
#   tkgrid(abfListLabel, row = 1, column = 0, sticky = 'w', pady = 5)
#   abfListBox <<- tklistbox(sidebarFrame, height = 5, width = 40, selectmode = 'multiple')
#   tkgrid(abfListBox, row = 2, column = 0, columnspan = 3, pady = 5, sticky = 'we')
  
#   paramFrame <- tkframe(sidebarFrame)
#   tkgrid(paramFrame, row = 3, column = 0, columnspan = 3, sticky = 'w')
  
#   experimentVar <<- tclVar('Voltage Clamp')
#   unitVar <<- tclVar('')
#   dataColVar <<- tclVar('')
#   dtVar <<- tclVar('')
#   ntracesVar <<- tclVar('')
  
#   tkgrid(tklabel(paramFrame, text = 'Experiment:'), row = 0, column = 0, sticky = 'w')
#   experimentCombo <- ttkcombobox(paramFrame, textvariable = experimentVar, values = c('Voltage Clamp', 'Current Clamp'), width = 15)
#   tkgrid(experimentCombo, row = 0, column = 1, sticky = 'w')
  
#   tkgrid(tklabel(paramFrame, text = 'Units:'), row = 1, column = 0, sticky = 'w')
#   unitEntry <- tkentry(paramFrame, textvariable = unitVar, width = 10, state = 'readonly')
#   tkgrid(unitEntry, row = 1, column = 1, sticky = 'w')
  
#   tkgrid(tklabel(paramFrame, text = 'Data Column:'), row = 2, column = 0, sticky = 'w')
#   dataColEntry <- tkentry(paramFrame, textvariable = dataColVar, width = 10, state = 'readonly')
#   tkgrid(dataColEntry, row = 2, column = 1, sticky = 'w')
  
#   tkgrid(tklabel(paramFrame, text = 'dt (ms):'), row = 3, column = 0, sticky = 'w')
#   dtEntry <- tkentry(paramFrame, textvariable = dtVar, width = 10, state = 'readonly')
#   tkgrid(dtEntry, row = 3, column = 1, sticky = 'w')
  
#   tkgrid(tklabel(paramFrame, text = 'Traces:'), row = 4, column = 0, sticky = 'w')
#   ntracesEntry <- tkentry(paramFrame, textvariable = ntracesVar, width = 10, state = 'readonly')
#   tkgrid(ntracesEntry, row = 4, column = 1, sticky = 'w')
  
#   baselineVar <<- tclVar('100')
#   stimTimeVar <<- tclVar('350')
#   xbarVar <<- tclVar('100')
#   ybarVar <<- tclVar('50')
  
#   tkgrid(tklabel(sidebarFrame, text = 'Baseline:'), row = 4, column = 0, sticky = 'w')
#   tkgrid(tkentry(sidebarFrame, textvariable = baselineVar, width = 10), row = 4, column = 1, sticky = 'w')
#   tkgrid(tklabel(sidebarFrame, text = 'Stimulation Time:'), row = 5, column = 0, sticky = 'w')
#   tkgrid(tkentry(sidebarFrame, textvariable = stimTimeVar, width = 10), row = 5, column = 1, sticky = 'w')
#   tkgrid(tklabel(sidebarFrame, text = 'x-bar length:'), row = 6, column = 0, sticky = 'w')
#   tkgrid(tkentry(sidebarFrame, textvariable = xbarVar, width = 10), row = 6, column = 1, sticky = 'w')
#   tkgrid(tklabel(sidebarFrame, text = 'y-bar length:'), row = 7, column = 0, sticky = 'w')
#   tkgrid(tkentry(sidebarFrame, textvariable = ybarVar, width = 10), row = 7, column = 1, sticky = 'w')
  
#   consoleText <<- tktext(sidebarFrame, width = 40, height = 4)
#   tkgrid(consoleText, row = 8, column = 0, columnspan = 3, pady = 5)
  
#   updateAdditionalParams <<- function(result) {
#     if (!is.null(result) && length(result$metadata) >= 1) {
#       meta1 <- result$metadata[[1]]
#       tclvalue(dtVar) <<- as.character(meta1$samplingIntervalInSec * 1000)
#       if (!is.null(meta1$header$lActualEpisodes)) {
#         tclvalue(ntracesVar) <<- as.character(meta1$header$lActualEpisodes)
#       } else {
#         tclvalue(ntracesVar) <<- 'N/A'
#       }
#       expType <- tclvalue(experimentVar)
#       col_idx <- choose_data_column(meta1$channelUnits, expType)
#       if (!is.na(col_idx)) {
#         tclvalue(unitVar) <<- meta1$channelUnits[col_idx]
#         tclvalue(dataColVar) <<- as.character(col_idx)
#       } else {
#         tclvalue(unitVar) <<- 'N/A'
#         tclvalue(dataColVar) <<- 'N/A'
#       }
#     }
#   }
  
#   tkbind(experimentCombo, '<<ComboboxSelected>>', function() {
#     if (exists('abf_analysis_result', envir = .GlobalEnv)) {
#       result <- get('abf_analysis_result', envir = .GlobalEnv)
#       updateAdditionalParams(result)
#     }
#   })
  
#   runAnalysis <<- function() {
#     folderPath <- tclvalue(folderPathVar)
#     if (nchar(folderPath) == 0) {
#       tkmessageBox(message = 'Please select an ABF folder first.')
#       return()
#     }
#     selIndices <- as.integer(tkcurselection(abfListBox))
#     allFiles <- as.character(tkget(abfListBox, 0, 'end'))
#     abf_files <- if (length(selIndices) == 0) allFiles else allFiles[selIndices + 1]
#     if (length(abf_files) == 0) {
#       tkmessageBox(message = 'No ABF files selected.')
#       return()
#     }
#     result <- tryCatch({
#       load_abf_data(abf_files = abf_files, abf_path = folderPath)
#     }, error = function(e) {
#       tkmessageBox(message = paste('Error during data loading:', e$message))
#       return(NULL)
#     })
#     result <<- result
#     if (!is.null(result)) {
#       tkdelete(consoleText, '1.0', 'end')
#       tkinsert(consoleText, 'end', paste('Data loaded. Processed', length(abf_files), 'file(s).'))
#       assign('abf_analysis_result', result, envir = .GlobalEnv)
#       updateAdditionalParams(result)
#       cons_msg <- check_consistency(result$metadata)
#       if (cons_msg == 'Data is consistent') {
#         tkmessageBox(message = cons_msg)
#       } else {
#         tkmessageBox(message = paste('ERROR:', cons_msg))
#       }
#       tkconfigure(runAnalysisButton, text = 'Load Data')
#     }
#   }
  
#   runAnalysisButton <<- tkbutton(sidebarFrame, text = 'Load Data', command = runAnalysis)
#   tkgrid(runAnalysisButton, row = 9, column = 0, columnspan = 3, pady = 5)
  
#   reviewButton <<- tkbutton(sidebarFrame, text = 'Review Recordings', command = review_recordings)
#   tkgrid(reviewButton, row = 10, column = 0, columnspan = 3, pady = 5)
  
#   averageApprovedTraces <<- function() {
#     if (length(traces2average) == 0 || all(sapply(traces2average, length) == 0)) {
#       tkmessageBox(message = 'No approved traces available. Please review recordings first.')
#       return()
#     }
#     folderPath <- tclvalue(folderPathVar)
#     if (nchar(folderPath) == 0) {
#       tkmessageBox(message = 'Please select an ABF folder first.')
#       return()
#     }
#     selIndices <- as.integer(tkcurselection(abfListBox))
#     allFiles <- as.character(tkget(abfListBox, 0, 'end'))
#     abf_files <- if (length(selIndices) == 0) allFiles else allFiles[selIndices + 1]
#     if (length(abf_files) == 0) {
#       tkmessageBox(message = 'No ABF files selected.')
#       return()
#     }
#     baseline <- as.numeric(tclvalue(baselineVar))
#     stimTime <- as.numeric(tclvalue(stimTimeVar))
#     xbar <- as.numeric(tclvalue(xbarVar))
#     ybar <- as.numeric(tclvalue(ybarVar))
    
#     result <- tryCatch({
#       abf_out <<- abf_averages(
#         datasets = abf_analysis_result$datasets,
#         traces2average = traces2average,
#         baseline = baseline,
#         stimulation_time = stimTime,
#         dataCol = as.numeric(tclvalue(dataColVar)),
#         xlim = NULL, ylim = NULL,
#         color = 'darkgray', 
#         xbar = xbar, ybar = ybar,
#         width = 5.25, height = 2.75
#       )
#       abf_out
#     }, error = function(e) {
#       tkmessageBox(message = paste('Error during averaging of approved traces:', e$message))
#       NULL
#     })

#     if (!is.null(result)) {
#       tkdelete(consoleText, '1.0', 'end')
#       tkinsert(consoleText, 'end', paste('Averaging on approved traces complete. Processed', length(abf_files), 'file(s).'))
      
#       # Update the global analysis result
#       abf_analysis_result <<- result
      
#       # Define drawPlot to display the averaged trace using scale bars only (no axes or title)
#       drawPlot <- function() {
#         if (exists('abf_analysis_result', envir = .GlobalEnv)) {
#           result <- get('abf_analysis_result', envir = .GlobalEnv)
#           datasets <- result$datasets
#           traces <- result$baseline_corrected_mean_data
#           if (length(datasets) > 0 && length(traces) > 0) {
#             par(mfrow = c(1, length(traces)))
#             # Determine shared limits
#             all_y <- unlist(traces)
#             shared_ylim <- range(all_y)
#             shared_xlim <- range(unlist(lapply(seq_along(traces), function(i) {
#               dt <- datasets[[i]]$samplingIntervalInSec * 1000
#               seq(0, by = dt, length.out = length(traces[[i]]))
#             })))
#             for (i in seq_along(traces)) {
#               dt_val <- datasets[[i]]$samplingIntervalInSec * 1000
#               time <- seq(0, by = dt_val, length.out = length(traces[[i]]))
#               show_bar <- (i == length(traces))
#               egs_plot(x = time, y = traces[[i]], color = 'darkgray',
#                        show_bar = show_bar, show_text = show_bar, 
#                        xbar = xbar, ybar = ybar,
#                        xlim = shared_xlim, ylim = shared_ylim)
#             }
#           } else {
#             plot.new()
#             text(0.5, 0.5, 'No data available')
#           }
#         } else {
#           plot.new()
#           text(0.5, 0.5, 'No analysis result to display')
#         }
#       }
      
#       # If a plot widget already exists, destroy it.
#       children <- as.character(tkwinfo('children', plotPanel))
#       for (child in children) {
#         tryCatch({
#           tkdestroy(.Tk.ID[[child]])
#         }, error = function(e) {}, silent = TRUE)
#       }

#       # Create a new frame to hold the averaged plot
#       avgFrame <- tkframe(plotPanel)
#       tkgrid(avgFrame, row = 0, column = 0, sticky = 'nsew')

#       # Use the same drawPlot function defined above
#       reviewPlot <<- tkrplot(avgFrame, fun = drawPlot, hscale = 1, vscale = 1)
#       tkgrid(reviewPlot, row = 0, column = 0)

#     }
#   }
  
#   averageApprovedButton <<- tkbutton(sidebarFrame, text = 'Average Approved Traces', command = averageApprovedTraces)
#   tkgrid(averageApprovedButton, row = 11, column = 0, columnspan = 3, pady = 5)
  
#   downloadButton <<- tkbutton(sidebarFrame, text = 'Download Data', command = download_data)
#   tkgrid(downloadButton, row = 12, column = 0, columnspan = 3, pady = 5)
  
#   tkfocus(tt)
# }

# ABF_analysis_tk_B()


