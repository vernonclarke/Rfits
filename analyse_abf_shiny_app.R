# this version combines versions A and B below:
# Remove all objects from the environment
rm(list = ls(all = TRUE))

# -----------------------------------------------------------------------------
# Package loading function and required packages:
# -----------------------------------------------------------------------------
load_required_packages <- function(packages) {
  new.packages <- packages[!(packages %in% installed.packages()[, "Package"])]
  if (length(new.packages)) install.packages(new.packages)
  invisible(lapply(packages, library, character.only = TRUE))
}
required.packages <- c("robustbase", "minpack.lm", "Rcpp", "signal", "readABF", "tcltk", "tkrplot")
load_required_packages(required.packages)

# -----------------------------------------------------------------------------
# Helper Functions
# -----------------------------------------------------------------------------

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
  if (experiment == "Voltage Clamp") {
    idx <- grep("A", channelUnits, ignore.case = TRUE)
  } else if (experiment == "Current Clamp") {
    idx <- grep("V", channelUnits, ignore.case = TRUE)
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
    return("Data is consistent")
  } else {
    error_msgs <- c()
    if (!dt_good) error_msgs <- c(error_msgs, paste("Inconsistent dt values:", paste(dt_values, collapse = ", ")))
    if (!unit_good) error_msgs <- c(error_msgs, paste("Inconsistent Units:", paste(unit_values, collapse = ", ")))
    if (!traces_good) error_msgs <- c(error_msgs, paste("Inconsistent Traces:", paste(traces_values, collapse = ", ")))
    return(paste(error_msgs, collapse = "; "))
  }
}

# Original plotting function – unchanged.
egs_plot <- function(x, y, sign = -1, xlim = NULL, ylim = NULL, lwd = 1, show_text = FALSE, 
                     xbar = 100, ybar = 50, color = "#4C77BB", show_bar = FALSE, cex = 0.6) {
  if (is.null(ylim))
    ylim <- if (sign == 1) c(0, max(y)) else c(-max(-y), 0)
  if (is.null(xlim))
    xlim <- c(min(x), max(x))
  idx1 <- which.min(abs(x - xlim[1]))
  idx2 <- which.min(abs(x - xlim[2]))
  plot(x[idx1:idx2], y[idx1:idx2], type = "l", col = color,
       xlim = xlim, ylim = ylim, bty = "n", lwd = lwd, lty = 1,
       axes = FALSE, frame = FALSE, xlab = "", ylab = "")
  if (show_bar) {
    ybar_start <- min(ylim) + (max(ylim) - min(ylim)) / 20
    x_start <- max(xlim) - xbar - 50
    y_start <- ybar_start
    x_end <- x_start + xbar
    y_end <- y_start + ybar
    segments(x_start, y_start, x_end, y_start, lwd = lwd, col = "black")
    segments(x_start, y_start, x_start, y_end, lwd = lwd, col = "black")
    if (show_text) {
      text(x = (x_start + x_end) / 2, y = y_start - ybar / 20, labels = paste(xbar, "ms"), 
           adj = c(0.5, 1), cex = cex)
      text(x = x_start - xbar / 4, y = (y_start + y_end) / 2, labels = paste(ybar, "pA"), 
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

# -----------------------------------------------------------------------------
# Download Data Function
# -----------------------------------------------------------------------------
download_data <- function() {
  if (is.null(averaged_data) || length(averaged_data) == 0) {
    tkmessageBox(message = "No averaged data available.")
    return()
  }
  dt_val <- master_abf$samplingIntervalInSec * 1000
  avg_length <- length(averaged_data[[1]])
  time_vec <- seq(0, by = dt_val, length.out = avg_length)
  avg_mat <- sapply(averaged_data, function(x) x)
  df <- data.frame(Time = time_vec, avg_mat)
  download_folder <- tclvalue(folderPathVar)
  if (nchar(download_folder) == 0) {
    tkmessageBox(message = "No folder selected for download.")
    return()
  }
  file_path <- file.path(download_folder, "averaged_data.csv")
  write.csv(df, file = file_path, row.names = FALSE)
  tkmessageBox(message = paste("Averaged data saved to", file_path))
}

# -----------------------------------------------------------------------------
# Original averaging routine for separate (non-concatenated) mode.
# We add an extra argument plotIt (default TRUE).
abf_averages <- function(datasets, 
                         baseline = 100, 
                         stimulation_time = 350, 
                         traces2average = NULL,
                         dataCol = 1, 
                         ylim = NULL, 
                         xlim = NULL, 
                         color = "#4C77BB", 
                         xbar = 100, 
                         ybar = 50, 
                         width = 5.25, 
                         height = 2.75, 
                         save = FALSE,
                         plotIt = TRUE) {
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
    seq(0, by = dt_val, length.out = length(responses0_mean[[iii]]))
  })
  if(plotIt){
    par(mfrow = c(1, N))
    show_bar <- rep(FALSE, N)
    if (N > 0) show_bar[N] <- TRUE
    for(ii in seq_len(N)) {
      for(ii in seq_len(N)) {
        egs_plot(x = time[[ii]], y = responses0_mean[[ii]], color = "darkgrey",
                 show_bar = FALSE, show_text = FALSE)
      }
    }
  }
  return(list(raw_data = responses,
              baseline_corrected_data = responses0,
              baseline_corrected_mean_data = responses0_mean,
              datasets = datasets))
}

# -----------------------------------------------------------------------------
# Function to concatenate ABF files into one master object.
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

# -----------------------------------------------------------------------------
# Global Variables
# -----------------------------------------------------------------------------
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

# -----------------------------------------------------------------------------
# Review Functions
# -----------------------------------------------------------------------------
# --- For Concatenated Mode ---
review_master_recordings <- function() {
  if(is.null(master_abf)) {
    tkmessageBox(message = "No master ABF data available. Please load data first.")
    return()
  }
  total_traces <<- length(master_abf$data)
  current_trace <<- 1
  current_group_selected <<- integer(0)
  groups_list <<- list()
  
  children <- as.character(tkwinfo("children", plotPanel))
  for(child in children) {
    tryCatch({ tkdestroy(.Tk.ID[[child]]) }, error = function(e) {}, silent = TRUE)
  }
  
  reviewFrame <<- tkframe(plotPanel)
  tkgrid(reviewFrame, row = 0, column = 0, sticky = "nsew")
  
  infoLabel <<- tklabel(reviewFrame, text = paste("Trace", current_trace, "of", total_traces))
  tkgrid(infoLabel, row = 0, column = 0, columnspan = 3)
  
  reviewPlot <<- tkrplot(reviewFrame, fun = function() {
    if(current_trace > total_traces){
      plot.new()
      text(0.5, 0.5, "No more traces to review.")
    } else {
      trace_matrix <- master_abf$data[[current_trace]]
      dt_val <- master_abf$samplingIntervalInSec * 1000
      time <- seq(0, by = dt_val, length.out = nrow(trace_matrix))
      data_column <- as.numeric(tclvalue(dataColVar))
      if(is.na(data_column) || data_column < 1 || data_column > ncol(trace_matrix))
        data_column <- 1
      trace <- trace_matrix[, data_column]
      plot(time, trace, col = "darkgrey", xlab = "Time (ms)",
           ylab = tclvalue(unitVar), type = "l", bty = "l",
           axes = FALSE, main = paste("Trace", current_trace))
      axis(1)
      axis(2, las = 1)
    }
  }, hscale = 1, vscale = 1)
  tkgrid(reviewPlot, row = 1, column = 0, columnspan = 3)
  
  acceptButton <<- tkbutton(reviewFrame, text = "Accept", command = function() {
    current_group_selected <<- c(current_group_selected, current_trace)
    tkconfigure(acceptButton, state = "disabled")
    tkconfigure(rejectButton, state = "normal")
  })
  tkgrid(acceptButton, row = 2, column = 0)
  
  rejectButton <<- tkbutton(reviewFrame, text = "Reject", command = function() {
    tkconfigure(rejectButton, state = "disabled")
    tkconfigure(acceptButton, state = "normal")
  })
  tkgrid(rejectButton, row = 2, column = 1)
  
  nextTraceButton <<- tkbutton(reviewFrame, text = "Next Trace", command = function() {
    tkconfigure(acceptButton, state = "normal")
    tkconfigure(rejectButton, state = "normal")
    if(current_trace < total_traces){
      current_trace <<- current_trace + 1
      tkconfigure(infoLabel, text = paste("Trace", current_trace, "of", total_traces))
      tkrreplot(reviewPlot)
    } else {
      tkmessageBox(message = "Review complete for all traces.")
    }
  })
  tkgrid(nextTraceButton, row = 2, column = 2)
  
  averageGroupButton <<- tkbutton(reviewFrame, text = "Add Selected Group", 
                                   command = function() {
    if(length(current_group_selected) == 0) {
      tkmessageBox(message = "No traces selected in current group.")
    } else {
      groups_list[[length(groups_list) + 1]] <<- current_group_selected
      tkmessageBox(message = paste("Group", length(groups_list), "selected with traces:",
                                   paste(current_group_selected, collapse = ", ")))
      current_group_selected <<- integer(0)
    }
  })
  tkgrid(averageGroupButton, row = 3, column = 0, columnspan = 3)
  
  selectionCompleteButton <<- tkbutton(reviewFrame, text = "Selection Complete", 
                                        command = function() {
    tkmessageBox(message = "Review complete. Approved traces have been stored.")
  })
  tkgrid(selectionCompleteButton, row = 4, column = 0, columnspan = 3)
}

# --- For Separate (Non-concatenated) Mode ---
review_recordings <- function() {
  if(!exists("abf_analysis_result", envir = .GlobalEnv)){
    tkmessageBox(message = "No analysis result available for review.")
    return()
  }
  result <- get("abf_analysis_result", envir = .GlobalEnv)
  datasets <- result$datasets
  traces2average <<- vector("list", length = length(datasets))
  for(i in seq_along(datasets)){
    traces2average[[i]] <<- integer(0)
  }
  current_dataset <<- 1
  current_trace <<- 1
  children <- as.character(tkwinfo("children", plotPanel))
  for(child in children){
    tryCatch({ tkdestroy(.Tk.ID[[child]]) }, error = function(e){}, silent = TRUE)
  }
  reviewFrame <<- tkframe(plotPanel)
  tkgrid(reviewFrame, row = 0, column = 0, sticky = "nsew")
  current_filename <- names(datasets)[current_dataset]
  infoLabel <<- tklabel(reviewFrame, text = paste(current_filename, "trace", current_trace))
  tkgrid(infoLabel, row = 0, column = 0, columnspan = 2)
  reviewPlot <<- tkrplot(reviewFrame, fun = function() {
    ds <- datasets[[current_dataset]]
    if(current_trace > length(ds$data)){
      plot.new()
      text(0.5, 0.5, paste("No more recordings in", current_filename))
    } else {
      trace_matrix <- ds$data[[current_trace]]
      data_column <- as.numeric(tclvalue(dataColVar))
      if(is.na(data_column) || data_column < 1 || data_column > ncol(trace_matrix))
        data_column <- 1
      dt_val <- ds$samplingIntervalInSec * 1000
      time <- seq(0, by = dt_val, length.out = nrow(trace_matrix))
      trace <- trace_matrix[, data_column]
      plot(time, trace, col = "darkgrey", xlab = "Time (ms)", ylab = tclvalue(unitVar),
           type = "l", bty = "l", axes = FALSE,
           main = paste(current_filename, "trace", current_trace))
      axis(1); axis(2, las = 1)
    }
  }, hscale = 1, vscale = 1)
  tkgrid(reviewPlot, row = 1, column = 0, columnspan = 2)
  acceptButton <- tkbutton(reviewFrame, text = "Accept", command = function() {
    traces2average[[current_dataset]] <<- c(traces2average[[current_dataset]], current_trace)
    tkconfigure(acceptButton, state = "disabled", relief = "sunken")
    tkconfigure(rejectButton, state = "normal", relief = "raised")
  })
  tkgrid(acceptButton, row = 2, column = 0)
  rejectButton <- tkbutton(reviewFrame, text = "Reject", command = function() {
    tkconfigure(rejectButton, state = "disabled", relief = "sunken")
    tkconfigure(acceptButton, state = "normal", relief = "raised")
  })
  tkgrid(rejectButton, row = 2, column = 1)
  nextTraceButton <- tkbutton(reviewFrame, text = "Next Recording", command = function() {
    tkconfigure(acceptButton, state = "normal", relief = "raised")
    tkconfigure(rejectButton, state = "normal", relief = "raised")
    ds <- datasets[[current_dataset]]
    numRecordings <- length(ds$data)
    if(current_trace < numRecordings){
      current_trace <<- current_trace + 1
    } else {
      tkmessageBox(message = paste("Finished reviewing", current_filename))
      if(current_dataset < length(datasets)){
        current_dataset <<- current_dataset + 1
        current_trace <<- 1
      } else {
        tkmessageBox(message = "Review complete. Approved recordings are in 'traces2average'.")
        return()
      }
    }
    current_filename <- names(datasets)[current_dataset]
    tkconfigure(infoLabel, text = paste(current_filename, "trace", current_trace))
    tkrreplot(reviewPlot)
  })
  tkgrid(nextTraceButton, row = 3, column = 0, columnspan = 2)
}

# -----------------------------------------------------------------------------
# Averaging Functions
# -----------------------------------------------------------------------------
# For concatenated mode:
average_selected_groups <- function() {
  if(length(groups_list) == 0){
    tkmessageBox(message = "No groups available for averaging. Please select groups first.")
    return()
  }
  dt_val <- master_abf$samplingIntervalInSec * 1000
  stim_time <- as.numeric(tclvalue(stimTimeVar))
  base_val <- as.numeric(tclvalue(baselineVar))
  data_column <- as.numeric(tclvalue(dataColVar))
  if(is.na(data_column) || data_column < 1) data_column <- 1
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
  children <- as.character(tkwinfo("children", plotPanel))
  for(child in children){
    tryCatch({ tkdestroy(.Tk.ID[[child]]) }, error = function(e){}, silent = TRUE)
  }
  avgFrame <<- tkframe(plotPanel)
  tkgrid(avgFrame, row = 0, column = 0, sticky = "nsew")
  drawAvgPlot <- function(){
    num_groups <- length(group_corrected_mean)
    if(num_groups < 1){
      plot.new()
      text(0.5, 0.5, "No averaged data available")
      return()
    }
    all_y <- unlist(group_corrected_mean)
    shared_ylim <- range(all_y)
    max_time <- max(sapply(group_corrected_mean, function(avg){ dt_val * (length(avg) - 1) }))
    shared_xlim <- c(0, max_time)
    par(mfrow = c(1, num_groups))
    for(i in seq_along(group_corrected_mean)){
      avg_trace <- group_corrected_mean[[i]]
      time <- seq(0, by = dt_val, length.out = length(avg_trace))
      show_bar <- (i == num_groups)
      egs_plot(x = time, y = avg_trace, color = "darkgrey",
               show_bar = show_bar, show_text = show_bar,
               xbar = as.numeric(tclvalue(xbarVar)), ybar = as.numeric(tclvalue(ybarVar)),
               xlim = shared_xlim, ylim = shared_ylim)
    }
  }
  avgPlot <<- tkrplot(avgFrame, fun = drawAvgPlot, hscale = 1, vscale = 1)
  tkgrid(avgPlot, row = 0, column = 0)
  tkmessageBox(message = "Averaging complete. Check the updated plot.")
}

# For separate mode (non-concatenated):
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
  allFiles <- as.character(tkget(abfListBox, 0, "end"))
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
    # Note: set plotIt = FALSE so abf_averages does not plot on the default device.
    abf_out <<- abf_averages(
      datasets = abf_analysis_result$datasets,
      traces2average = traces2average,
      baseline = baseline,
      stimulation_time = stimTime,
      dataCol = as.numeric(tclvalue(dataColVar)),
      xlim = NULL, ylim = NULL,
      color = "darkgrey",
      xbar = xbar, ybar = ybar,
      width = 5.25, height = 2.75,
      plotIt = FALSE
    )
    abf_out
  }, error = function(e){
    tkmessageBox(message = paste("Error during averaging of approved traces:", e$message))
    NULL
  })
  if(!is.null(result)){
    tkdelete(consoleText, "1.0", "end")
    tkinsert(consoleText, "end", paste("Averaging on approved traces complete. Processed", length(abf_files), "file(s)."))
    abf_analysis_result <<- result
    # Now, create a new tkrplot widget to display the averaged plot in the UI.
    drawPlot <- function(){
      if(exists("abf_analysis_result", envir = .GlobalEnv)){
        result <- get("abf_analysis_result", envir = .GlobalEnv)
        datasets <- result$datasets
        traces <- result$baseline_corrected_mean_data
        if(length(datasets) > 0 && length(traces) > 0){
          par(mfrow = c(1, length(traces)))
          all_y <- unlist(traces)
          shared_ylim <- range(all_y)
          shared_xlim <- range(unlist(lapply(seq_along(traces), function(i){
            dt <- datasets[[i]]$samplingIntervalInSec * 1000
            seq(0, by = dt, length.out = length(traces[[i]]))
          })))
          for(i in seq_along(traces)){
            dt_val <- datasets[[i]]$samplingIntervalInSec * 1000
            time <- seq(0, by = dt_val, length.out = length(traces[[i]]))
            show_bar <- (i == length(traces))
            egs_plot(x = time, y = traces[[i]], color = "darkgrey",
                     show_bar = show_bar, show_text = show_bar,
                     xbar = xbar, ybar = ybar,
                     xlim = shared_xlim, ylim = shared_ylim)
          }
        } else {
          plot.new()
          text(0.5, 0.5, "No data available")
        }
      } else {
        plot.new()
        text(0.5, 0.5, "No analysis result to display")
      }
    }
    
    children <- as.character(tkwinfo("children", plotPanel))
    for(child in children){
      tryCatch({ tkdestroy(.Tk.ID[[child]]) }, error = function(e){}, silent = TRUE)
    }
    
    avgFrame <- tkframe(plotPanel)
    tkgrid(avgFrame, row = 0, column = 0, sticky = "nsew")
    reviewPlot <<- tkrplot(avgFrame, fun = drawPlot, hscale = 1, vscale = 1)
    tkgrid(reviewPlot, row = 0, column = 0)
  }
}

# -----------------------------------------------------------------------------
# Main UI Setup
# -----------------------------------------------------------------------------
ABF_analysis_tk <- function() {
  tt <- tktoplevel()
  tkwm.title(tt, "ABF Analysis")
  sidebarFrame <- tkframe(tt)
  mainFrame <- tkframe(tt)
  tkgrid(sidebarFrame, row = 0, column = 0, sticky = "ns")
  tkgrid(mainFrame, row = 0, column = 1, sticky = "nsew")
  tkgrid.rowconfigure(tt, 0, weight = 1)
  tkgrid.columnconfigure(tt, 1, weight = 1)
  
  # Save mainFrame as the global plot panel.
  plotPanel <<- mainFrame
  
  folderLabel <- tklabel(sidebarFrame, text = "Select ABF Folder:")
  tkgrid(folderLabel, row = 0, column = 0, sticky = "w")
  folderPathVar <<- tclVar("")
  folderEntry <- tkentry(sidebarFrame, textvariable = folderPathVar, width = 30)
  tkgrid(folderEntry, row = 0, column = 1, sticky = "w")
  
  browseFolderButton <- tkbutton(sidebarFrame, text = "Browse", command = function(){
    folderPath <- tclvalue(tkchooseDirectory())
    if(nchar(folderPath) > 0){
      tclvalue(folderPathVar) <<- folderPath
      abf_list <- list.files(path = folderPath, pattern = "\\.abf$", ignore.case = TRUE)
      if(length(abf_list) == 0){
        tkmessageBox(message = "No ABF files found in the selected folder.")
      } else {
        tkdelete(abfListBox, 0, "end")
        for(f in abf_list){ tkinsert(abfListBox, "end", f) }
        firstFilePath <- file.path(folderPath, abf_list[1])
        ds <- readABF(firstFilePath)
        dummy_result <- list(metadata = list(extract_metadata(ds)))
        updateAdditionalParams(dummy_result)
      }
    }
  })
  tkgrid(browseFolderButton, row = 0, column = 2, padx = 5)
  
  abfListLabel <- tklabel(sidebarFrame, text = "ABF Files:")
  tkgrid(abfListLabel, row = 1, column = 0, sticky = "w", pady = 5)
  abfListBox <<- tklistbox(sidebarFrame, height = 5, width = 40, selectmode = "multiple")
  tkgrid(abfListBox, row = 2, column = 0, columnspan = 3, pady = 5, sticky = "we")
  
  paramFrame <- tkframe(sidebarFrame)
  tkgrid(paramFrame, row = 3, column = 0, columnspan = 3, sticky = "w")
  
  experimentVar <<- tclVar("Voltage Clamp")
  unitVar <<- tclVar("")
  dataColVar <<- tclVar("")
  dtVar <<- tclVar("")
  ntracesVar <<- tclVar("")
  
  tkgrid(tklabel(paramFrame, text = "Experiment:"), row = 0, column = 0, sticky = "w")
  experimentCombo <- ttkcombobox(paramFrame, textvariable = experimentVar, values = c("Voltage Clamp", "Current Clamp"), width = 15)
  tkgrid(experimentCombo, row = 0, column = 1, sticky = "w")
  
  tkgrid(tklabel(paramFrame, text = "Units:"), row = 1, column = 0, sticky = "w")
  unitEntry <- tkentry(paramFrame, textvariable = unitVar, width = 10, state = "readonly")
  tkgrid(unitEntry, row = 1, column = 1, sticky = "w")
  
  tkgrid(tklabel(paramFrame, text = "Data Column:"), row = 2, column = 0, sticky = "w")
  dataColEntry <- tkentry(paramFrame, textvariable = dataColVar, width = 10, state = "readonly")
  tkgrid(dataColEntry, row = 2, column = 1, sticky = "w")
  
  tkgrid(tklabel(paramFrame, text = "dt (ms):"), row = 3, column = 0, sticky = "w")
  dtEntry <- tkentry(paramFrame, textvariable = dtVar, width = 10, state = "readonly")
  tkgrid(dtEntry, row = 3, column = 1, sticky = "w")
  
  tkgrid(tklabel(paramFrame, text = "Traces:"), row = 4, column = 0, sticky = "w")
  ntracesEntry <- tkentry(paramFrame, textvariable = ntracesVar, width = 10, state = "readonly")
  tkgrid(ntracesEntry, row = 4, column = 1, sticky = "w")
  
  baselineVar <<- tclVar("100")
  stimTimeVar <<- tclVar("350")
  xbarVar <<- tclVar("100")
  ybarVar <<- tclVar("50")
  
  tkgrid(tklabel(sidebarFrame, text = "Baseline:"), row = 4, column = 0, sticky = "w")
  tkgrid(tkentry(sidebarFrame, textvariable = baselineVar, width = 10), row = 4, column = 1, sticky = "w")
  tkgrid(tklabel(sidebarFrame, text = "Stimulation Time:"), row = 5, column = 0, sticky = "w")
  tkgrid(tkentry(sidebarFrame, textvariable = stimTimeVar, width = 10), row = 5, column = 1, sticky = "w")
  tkgrid(tklabel(sidebarFrame, text = "x-bar length:"), row = 6, column = 0, sticky = "w")
  tkgrid(tkentry(sidebarFrame, textvariable = xbarVar, width = 10), row = 6, column = 1, sticky = "w")
  tkgrid(tklabel(sidebarFrame, text = "y-bar length:"), row = 7, column = 0, sticky = "w")
  tkgrid(tkentry(sidebarFrame, textvariable = ybarVar, width = 10), row = 7, column = 1, sticky = "w")
  
  # NEW: Checkbutton for concatenation mode.
  concatMode <<- tclVar("0")   # "0" means separate (default); "1" means concatenate.
  concatButton <- tkcheckbutton(sidebarFrame, variable = concatMode, text = "Concatenate Imported ABFs")
  tkgrid(concatButton, row = 8, column = 0, columnspan = 3)
  
  consoleText <<- tktext(sidebarFrame, width = 40, height = 4)
  tkgrid(consoleText, row = 9, column = 0, columnspan = 3, pady = 5)
  
  updateAdditionalParams <<- function(result) {
    if (!is.null(result) && length(result$metadata) >= 1) {
      meta1 <- result$metadata[[1]]
      tclvalue(dtVar) <<- as.character(meta1$samplingIntervalInSec * 1000)
      if (!is.null(meta1$header$lActualEpisodes)) {
        tclvalue(ntracesVar) <<- as.character(meta1$header$lActualEpisodes)
      } else {
        tclvalue(ntracesVar) <<- "N/A"
      }
      expType <- tclvalue(experimentVar)
      col_idx <- choose_data_column(meta1$channelUnits, expType)
      if (!is.na(col_idx)) {
        tclvalue(unitVar) <<- meta1$channelUnits[col_idx]
        tclvalue(dataColVar) <<- as.character(col_idx)
      } else {
        tclvalue(unitVar) <<- "N/A"
        tclvalue(dataColVar) <<- "N/A"
      }
    }
  }
  
  tkbind(experimentCombo, "<<ComboboxSelected>>", function() {
    if (exists("abf_analysis_result", envir = .GlobalEnv)) {
      result <- get("abf_analysis_result", envir = .GlobalEnv)
      updateAdditionalParams(result)
    }
  })
  
  runAnalysis <<- function() {
    folderPath <- tclvalue(folderPathVar)
    if (nchar(folderPath) == 0) {
      tkmessageBox(message = "Please select an ABF folder first.")
      return()
    }
    selIndices <- as.integer(tkcurselection(abfListBox))
    allFiles <- as.character(tkget(abfListBox, 0, "end"))
    abf_files <- if (length(selIndices) == 0) allFiles else allFiles[selIndices + 1]
    if (length(abf_files) == 0) {
      tkmessageBox(message = "No ABF files selected.")
      return()
    }
    result <- tryCatch({
      load_abf_data(abf_files = abf_files, abf_path = folderPath)
    }, error = function(e) {
      tkmessageBox(message = paste("Error during data loading:", e$message))
      return(NULL)
    })
    if (!is.null(result)) {
      tkdelete(consoleText, "1.0", "end")
      tkinsert(consoleText, "end", paste("Data loaded. Processed", length(abf_files), "file(s)."))
      assign("abf_analysis_result", result, envir = .GlobalEnv)
      updateAdditionalParams(result)
      
      cons_msg <- check_consistency(result$metadata)
      if (cons_msg == "Data is consistent") {
        tkmessageBox(message = cons_msg)
        # If concatenation is selected, combine; otherwise, use the original separate structure.
        if (as.character(tclvalue(concatMode)) == "1") {
          master_abf <<- combine_abf_data(result)
        } else {
          master_abf <<- result
        }
      } else {
        tkmessageBox(message = paste("ERROR:", cons_msg))
      }
      tkconfigure(runAnalysisButton, text = "Load Data")
    }
  }
  
  runAnalysisButton <<- tkbutton(sidebarFrame, text = "Load Data", command = runAnalysis)
  tkgrid(runAnalysisButton, row = 10, column = 0, columnspan = 3, pady = 5)
  
  reviewButton <<- tkbutton(sidebarFrame, text = "Review Recordings", command = function() {
    if (as.character(tclvalue(concatMode)) == "1") {
      review_master_recordings()
    } else {
      review_recordings()
    }
  })
  tkgrid(reviewButton, row = 11, column = 0, columnspan = 3, pady = 5)
  
  avgApprovedTracesButton <<- tkbutton(sidebarFrame, text = "Average Approved Traces", 
                                        command = function() {
    if (as.character(tclvalue(concatMode)) == "1") {
      average_selected_groups()
    } else {
      averageApprovedTraces_sep()
    }
  })
  tkgrid(avgApprovedTracesButton, row = 12, column = 0, columnspan = 3, pady = 5)
  
  downloadButton <<- tkbutton(sidebarFrame, text = "Download Data", command = download_data)
  tkgrid(downloadButton, row = 13, column = 0, columnspan = 3, pady = 5)
  
  tkfocus(tt)
}

# Launch the UI.
ABF_analysis_tk()

# this version concatenates the imported abfs if files are consistent
# Remove all objects from the environment
rm(list = ls(all = TRUE))

# -----------------------------------------------------------------------------
# Package loading function and required packages:
# -----------------------------------------------------------------------------
load_required_packages <- function(packages) {
  new.packages <- packages[!(packages %in% installed.packages()[, "Package"])]
  if (length(new.packages)) install.packages(new.packages)
  invisible(lapply(packages, library, character.only = TRUE))
}
required.packages <- c("robustbase", "minpack.lm", "Rcpp", "signal", "readABF", "tcltk", "tkrplot")
load_required_packages(required.packages)

# -----------------------------------------------------------------------------
# Helper Functions
# -----------------------------------------------------------------------------

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
  if (experiment == "Voltage Clamp") {
    idx <- grep("A", channelUnits, ignore.case = TRUE)
  } else if (experiment == "Current Clamp") {
    idx <- grep("V", channelUnits, ignore.case = TRUE)
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
    return("Data is consistent")
  } else {
    error_msgs <- c()
    if (!dt_good) error_msgs <- c(error_msgs, paste("Inconsistent dt values:", paste(dt_values, collapse = ", ")))
    if (!unit_good) error_msgs <- c(error_msgs, paste("Inconsistent Units:", paste(unit_values, collapse = ", ")))
    if (!traces_good) error_msgs <- c(error_msgs, paste("Inconsistent Traces:", paste(traces_values, collapse = ", ")))
    return(paste(error_msgs, collapse = "; "))
  }
}

# Original plotting function (same as in your original code)
egs_plot <- function(x, y, sign = -1, xlim = NULL, ylim = NULL, lwd = 1, show_text = FALSE, 
                     xbar = 100, ybar = 50, color = "#4C77BB", show_bar = FALSE, cex = 0.6) {
  if (is.null(ylim))
    ylim <- if (sign == 1) c(0, max(y)) else c(-max(-y), 0)
  if (is.null(xlim))
    xlim <- c(min(x), max(x))
  idx1 <- which.min(abs(x - xlim[1]))
  idx2 <- which.min(abs(x - xlim[2]))
  plot(x[idx1:idx2], y[idx1:idx2], type = "l", col = color,
       xlim = xlim, ylim = ylim, bty = "n", lwd = lwd, lty = 1,
       axes = FALSE, frame = FALSE, xlab = "", ylab = "")
  if (show_bar) {
    ybar_start <- min(ylim) + (max(ylim) - min(ylim)) / 20
    x_start <- max(xlim) - xbar - 50
    y_start <- ybar_start
    x_end <- x_start + xbar
    y_end <- y_start + ybar
    segments(x_start, y_start, x_end, y_start, lwd = lwd, col = "black")
    segments(x_start, y_start, x_start, y_end, lwd = lwd, col = "black")
    if (show_text) {
      text(x = (x_start + x_end) / 2, y = y_start - ybar / 20, 
           labels = paste(xbar, "ms"), adj = c(0.5, 1), cex = cex)
      text(x = x_start - xbar / 4, y = (y_start + y_end) / 2,  
           labels = paste(ybar, "pA"), adj = c(0.5, 0.5), srt = 90, cex = cex)
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

# -----------------------------------------------------------------------------
# Download Data Function: Write out the baseline_corrected_mean_data stored
# in the global averaged_data variable.
# -----------------------------------------------------------------------------
download_data <- function() {
  if (is.null(averaged_data) || length(averaged_data) == 0) {
    tkmessageBox(message = "No averaged data available.")
    return()
  }
  dt_val <- master_abf$samplingIntervalInSec * 1000
  avg_length <- length(averaged_data[[1]])
  time_vec <- seq(0, by = dt_val, length.out = avg_length)
  avg_mat <- sapply(averaged_data, function(x) x)
  df <- data.frame(Time = time_vec, avg_mat)
  download_folder <- tclvalue(folderPathVar)
  if(nchar(download_folder) == 0) {
    tkmessageBox(message = "No folder selected for download.")
    return()
  }
  file_path <- file.path(download_folder, "averaged_data.csv")
  write.csv(df, file = file_path, row.names = FALSE)
  tkmessageBox(message = paste("Averaged data saved to", file_path))
}

# -----------------------------------------------------------------------------
# Combine ABF Files into a Master Object
# -----------------------------------------------------------------------------
combine_abf_data <- function(result) {
  master_abf <- list()
  master_abf$data <- list()
  master_abf$source_files <- c()
  master_abf$samplingIntervalInSec <- result$datasets[[1]]$samplingIntervalInSec
  
  for (i in seq_along(result$datasets)) {
    ds <- result$datasets[[i]]
    n_traces <- length(ds$data)
    master_abf$data <- c(master_abf$data, ds$data)
    master_abf$source_files <- c(master_abf$source_files, rep(names(result$datasets)[i], n_traces))
  }
  return(master_abf)
}

# -----------------------------------------------------------------------------
# Global Variables for Master ABF, Review Process, and Averaged Data:
# -----------------------------------------------------------------------------
master_abf <<- NULL      # Combined ABF data.
averaged_data <<- NULL   # Will hold baseline_corrected_mean_data.
current_trace <<- 1      # Currently displayed trace.
total_traces <<- 0       # Total number of traces.
current_group_selected <<- integer(0)  # Indices approved for the current group.
groups_list <<- list()   # List of groups (each group is a vector of trace indices).

# -----------------------------------------------------------------------------
# Revised Review Function for Master ABF
# -----------------------------------------------------------------------------
review_master_recordings <- function() {
  if (is.null(master_abf)) {
    tkmessageBox(message = "No master ABF data available for review. Please load ABF files first.")
    return()
  }
  
  total_traces <<- length(master_abf$data)
  current_trace <<- 1
  current_group_selected <<- integer(0)
  groups_list <<- list()
  
  # Clear previous UI plot panel children.
  children <- as.character(tkwinfo("children", plotPanel))
  for (child in children) {
    tryCatch({ tkdestroy(.Tk.ID[[child]]) }, error = function(e) {}, silent = TRUE)
  }
  
  reviewFrame <<- tkframe(plotPanel)
  tkgrid(reviewFrame, row = 0, column = 0, sticky = "nsew")
  
  infoLabel <<- tklabel(reviewFrame, text = paste("Trace", current_trace, "of", total_traces))
  tkgrid(infoLabel, row = 0, column = 0, columnspan = 3)
  
  reviewPlot <<- tkrplot(reviewFrame, fun = function() {
    if (current_trace > total_traces) {
      plot.new()
      text(0.5, 0.5, "No more traces to review.")
    } else {
      trace_matrix <- master_abf$data[[current_trace]]
      dt_val <- master_abf$samplingIntervalInSec * 1000
      time <- seq(0, by = dt_val, length.out = nrow(trace_matrix))
      data_column <- as.numeric(tclvalue(dataColVar))
      if (is.na(data_column) || data_column < 1 || data_column > ncol(trace_matrix))
        data_column <- 1
      trace <- trace_matrix[, data_column]
      plot(time, trace, col = "darkgrey", xlab = "Time (ms)",
           ylab = tclvalue(unitVar), type = "l", bty = "l",
           axes = FALSE, main = paste("Trace", current_trace))
      axis(1)
      axis(2, las = 1)
    }
  }, hscale = 1, vscale = 1)
  tkgrid(reviewPlot, row = 1, column = 0, columnspan = 3)
  
  acceptButton <<- tkbutton(reviewFrame, text = "Accept", command = function() {
    current_group_selected <<- c(current_group_selected, current_trace)
    tkconfigure(acceptButton, state = "disabled")
    tkconfigure(rejectButton, state = "normal")
  })
  tkgrid(acceptButton, row = 2, column = 0)
  
  rejectButton <<- tkbutton(reviewFrame, text = "Reject", command = function() {
    tkconfigure(rejectButton, state = "disabled")
    tkconfigure(acceptButton, state = "normal")
  })
  tkgrid(rejectButton, row = 2, column = 1)
  
  nextTraceButton <<- tkbutton(reviewFrame, text = "Next Trace", command = function() {
    tkconfigure(acceptButton, state = "normal")
    tkconfigure(rejectButton, state = "normal")
    if (current_trace < total_traces) {
      current_trace <<- current_trace + 1
      tkconfigure(infoLabel, text = paste("Trace", current_trace, "of", total_traces))
      tkrreplot(reviewPlot)
    } else {
      tkmessageBox(message = "Review complete for all traces.")
    }
  })
  tkgrid(nextTraceButton, row = 2, column = 2)
  
  # Button to add the current approved trace group.
  averageGroupButton <<- tkbutton(reviewFrame, text = "Average Selected Traces", 
                                   command = function() {
    if (length(current_group_selected) == 0) {
      tkmessageBox(message = "No traces selected in current group.")
    } else {
      groups_list[[length(groups_list) + 1]] <<- current_group_selected
      tkmessageBox(message = paste("Group", length(groups_list), 
                                   "selected with traces:", paste(current_group_selected, collapse = ", ")))
      current_group_selected <<- integer(0)
    }
  })
  tkgrid(averageGroupButton, row = 3, column = 0, columnspan = 3)
  
  # Replace the averaging step with a "Selection Complete" button.
  selectionCompleteButton <<- tkbutton(reviewFrame, text = "Selection Complete", 
                                        command = function() {
    tkmessageBox(message = "Review complete. Approved traces have been stored.")
    # Optionally, you could close the review window here.
  })
  tkgrid(selectionCompleteButton, row = 4, column = 0, columnspan = 3)
}

# -----------------------------------------------------------------------------
# New Averaging Function: Compute baseline-corrected means (mimicking your original abf_averages)
# -----------------------------------------------------------------------------
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
  
  # baseline2zero function same as your original.
  baseline2zero <- function(y, dt, stimulation_time, baseline) {
    idx1 <- (stimulation_time - baseline) / dt
    idx2 <- baseline / dt
    y1 <- y[idx1:length(y)]
    y1 <- y1 - mean(y1[1:idx2])
    y1 - mean(y1[1:idx2])
  }
  
  # For each group, extract the trace data, apply baseline correction and compute the row–wise mean.
  group_corrected_mean <- lapply(groups_list, function(group_indices) {
    traces_corrected <- lapply(group_indices, function(i) {
      trace <- master_abf$data[[i]][, data_column]
      baseline2zero(trace, dt = dt_val, stimulation_time = stim_time, baseline = base_val)
    })
    trace_mat <- do.call(cbind, traces_corrected)
    rowMeans(trace_mat)
  })
  
  # Save baseline_corrected_mean_data globally.
  averaged_data <<- group_corrected_mean
  
  # Plot the averaged traces in the UI.
  children <- as.character(tkwinfo("children", plotPanel))
  for (child in children) {
    tryCatch({ tkdestroy(.Tk.ID[[child]]) }, error = function(e) {}, silent = TRUE)
  }
  avgFrame <<- tkframe(plotPanel)
  tkgrid(avgFrame, row = 0, column = 0, sticky = "nsew")
  
  drawAvgPlot <- function() {
    num_groups <- length(group_corrected_mean)
    if (num_groups < 1) {
      plot.new()
      text(0.5, 0.5, "No averaged data available")
      return()
    }
    all_y <- unlist(group_corrected_mean)
    shared_ylim <- range(all_y)
    max_time <- max(sapply(group_corrected_mean, function(avg) { dt_val * (length(avg) - 1) }))
    shared_xlim <- c(0, max_time)
    
    par(mfrow = c(1, num_groups))
    for (i in seq_along(group_corrected_mean)) {
      avg_trace <- group_corrected_mean[[i]]
      time <- seq(0, by = dt_val, length.out = length(avg_trace))
      show_bar <- (i == num_groups)
      egs_plot(x = time, y = avg_trace, color = "darkgrey",
               show_bar = show_bar, show_text = show_bar,
               xbar = as.numeric(tclvalue(xbarVar)), ybar = as.numeric(tclvalue(ybarVar)),
               xlim = shared_xlim, ylim = shared_ylim)
    }
  }
  
  avgPlot <<- tkrplot(avgFrame, fun = drawAvgPlot, hscale = 1, vscale = 1)
  tkgrid(avgPlot, row = 0, column = 0)
  tkmessageBox(message = "Averaging complete. Check the updated plot.")
}

# -----------------------------------------------------------------------------
# Main UI Setup
# -----------------------------------------------------------------------------
ABF_analysis_tk_A <- function() {
  tt <- tktoplevel()
  tkwm.title(tt, "ABF Analysis")
  sidebarFrame <- tkframe(tt)
  mainFrame <- tkframe(tt)
  tkgrid(sidebarFrame, row = 0, column = 0, sticky = "ns")
  tkgrid(mainFrame, row = 0, column = 1, sticky = "nsew")
  tkgrid.rowconfigure(tt, 0, weight = 1)
  tkgrid.columnconfigure(tt, 1, weight = 1)
  
  # Save mainFrame as the global UI plot panel.
  plotPanel <<- mainFrame
  
  folderLabel <- tklabel(sidebarFrame, text = "Select ABF Folder:")
  tkgrid(folderLabel, row = 0, column = 0, sticky = "w")
  folderPathVar <<- tclVar("")
  folderEntry <- tkentry(sidebarFrame, textvariable = folderPathVar, width = 30)
  tkgrid(folderEntry, row = 0, column = 1, sticky = "w")
  
  browseFolderButton <- tkbutton(sidebarFrame, text = "Browse", command = function() {
    folderPath <- tclvalue(tkchooseDirectory())
    if (nchar(folderPath) > 0) {
      tclvalue(folderPathVar) <<- folderPath
      abf_list <- list.files(path = folderPath, pattern = "\\.abf$", ignore.case = TRUE)
      if (length(abf_list) == 0) {
        tkmessageBox(message = "No ABF files found in the selected folder.")
      } else {
        tkdelete(abfListBox, 0, "end")
        for (f in abf_list) { tkinsert(abfListBox, "end", f) }
        firstFilePath <- file.path(folderPath, abf_list[1])
        ds <- readABF(firstFilePath)
        dummy_result <- list(metadata = list(extract_metadata(ds)))
        updateAdditionalParams(dummy_result)
      }
    }
  })
  tkgrid(browseFolderButton, row = 0, column = 2, padx = 5)
  
  abfListLabel <- tklabel(sidebarFrame, text = "ABF Files:")
  tkgrid(abfListLabel, row = 1, column = 0, sticky = "w", pady = 5)
  abfListBox <<- tklistbox(sidebarFrame, height = 5, width = 40, selectmode = "multiple")
  tkgrid(abfListBox, row = 2, column = 0, columnspan = 3, pady = 5, sticky = "we")
  
  paramFrame <- tkframe(sidebarFrame)
  tkgrid(paramFrame, row = 3, column = 0, columnspan = 3, sticky = "w")
  
  experimentVar <<- tclVar("Voltage Clamp")
  unitVar <<- tclVar("")
  dataColVar <<- tclVar("")
  dtVar <<- tclVar("")
  ntracesVar <<- tclVar("")
  
  tkgrid(tklabel(paramFrame, text = "Experiment:"), row = 0, column = 0, sticky = "w")
  experimentCombo <- ttkcombobox(paramFrame, textvariable = experimentVar, values = c("Voltage Clamp", "Current Clamp"), width = 15)
  tkgrid(experimentCombo, row = 0, column = 1, sticky = "w")
  
  tkgrid(tklabel(paramFrame, text = "Units:"), row = 1, column = 0, sticky = "w")
  unitEntry <- tkentry(paramFrame, textvariable = unitVar, width = 10, state = "readonly")
  tkgrid(unitEntry, row = 1, column = 1, sticky = "w")
  
  tkgrid(tklabel(paramFrame, text = "Data Column:"), row = 2, column = 0, sticky = "w")
  dataColEntry <- tkentry(paramFrame, textvariable = dataColVar, width = 10, state = "readonly")
  tkgrid(dataColEntry, row = 2, column = 1, sticky = "w")
  
  tkgrid(tklabel(paramFrame, text = "dt (ms):"), row = 3, column = 0, sticky = "w")
  dtEntry <- tkentry(paramFrame, textvariable = dtVar, width = 10, state = "readonly")
  tkgrid(dtEntry, row = 3, column = 1, sticky = "w")
  
  tkgrid(tklabel(paramFrame, text = "Traces:"), row = 4, column = 0, sticky = "w")
  ntracesEntry <- tkentry(paramFrame, textvariable = ntracesVar, width = 10, state = "readonly")
  tkgrid(ntracesEntry, row = 4, column = 1, sticky = "w")
  
  baselineVar <<- tclVar("100")
  stimTimeVar <<- tclVar("350")
  xbarVar <<- tclVar("100")
  ybarVar <<- tclVar("50")
  
  tkgrid(tklabel(sidebarFrame, text = "Baseline:"), row = 4, column = 0, sticky = "w")
  tkgrid(tkentry(sidebarFrame, textvariable = baselineVar, width = 10), row = 4, column = 1, sticky = "w")
  tkgrid(tklabel(sidebarFrame, text = "Stimulation Time:"), row = 5, column = 0, sticky = "w")
  tkgrid(tkentry(sidebarFrame, textvariable = stimTimeVar, width = 10), row = 5, column = 1, sticky = "w")
  tkgrid(tklabel(sidebarFrame, text = "x-bar length:"), row = 6, column = 0, sticky = "w")
  tkgrid(tkentry(sidebarFrame, textvariable = xbarVar, width = 10), row = 6, column = 1, sticky = "w")
  tkgrid(tklabel(sidebarFrame, text = "y-bar length:"), row = 7, column = 0, sticky = "w")
  tkgrid(tkentry(sidebarFrame, textvariable = ybarVar, width = 10), row = 7, column = 1, sticky = "w")
  
  consoleText <<- tktext(sidebarFrame, width = 40, height = 4)
  tkgrid(consoleText, row = 8, column = 0, columnspan = 3, pady = 5)
  
  updateAdditionalParams <<- function(result) {
    if (!is.null(result) && length(result$metadata) >= 1) {
      meta1 <- result$metadata[[1]]
      tclvalue(dtVar) <<- as.character(meta1$samplingIntervalInSec * 1000)
      if (!is.null(meta1$header$lActualEpisodes)) {
        tclvalue(ntracesVar) <<- as.character(meta1$header$lActualEpisodes)
      } else {
        tclvalue(ntracesVar) <<- "N/A"
      }
      expType <- tclvalue(experimentVar)
      col_idx <- choose_data_column(meta1$channelUnits, expType)
      if (!is.na(col_idx)) {
        tclvalue(unitVar) <<- meta1$channelUnits[col_idx]
        tclvalue(dataColVar) <<- as.character(col_idx)
      } else {
        tclvalue(unitVar) <<- "N/A"
        tclvalue(dataColVar) <<- "N/A"
      }
    }
  }
  
  tkbind(experimentCombo, "<<ComboboxSelected>>", function() {
    if (exists("abf_analysis_result", envir = .GlobalEnv)) {
      result <- get("abf_analysis_result", envir = .GlobalEnv)
      updateAdditionalParams(result)
    }
  })
  
  runAnalysis <<- function() {
    folderPath <- tclvalue(folderPathVar)
    if (nchar(folderPath) == 0) {
      tkmessageBox(message = "Please select an ABF folder first.")
      return()
    }
    selIndices <- as.integer(tkcurselection(abfListBox))
    allFiles <- as.character(tkget(abfListBox, 0, "end"))
    abf_files <- if (length(selIndices) == 0) allFiles else allFiles[selIndices + 1]
    if (length(abf_files) == 0) {
      tkmessageBox(message = "No ABF files selected.")
      return()
    }
    result <- tryCatch({
      load_abf_data(abf_files = abf_files, abf_path = folderPath)
    }, error = function(e) {
      tkmessageBox(message = paste("Error during data loading:", e$message))
      return(NULL)
    })
    if (!is.null(result)) {
      tkdelete(consoleText, "1.0", "end")
      tkinsert(consoleText, "end", paste("Data loaded. Processed", length(abf_files), "file(s)."))
      assign("abf_analysis_result", result, envir = .GlobalEnv)
      updateAdditionalParams(result)
      cons_msg <- check_consistency(result$metadata)
      if (cons_msg == "Data is consistent") {
        tkmessageBox(message = cons_msg)
        master_abf <<- combine_abf_data(result)
      } else {
        tkmessageBox(message = paste("ERROR:", cons_msg))
      }
      tkconfigure(runAnalysisButton, text = "Load Data")
    }
  }
  
  runAnalysisButton <<- tkbutton(sidebarFrame, text = "Load Data", command = runAnalysis)
  tkgrid(runAnalysisButton, row = 9, column = 0, columnspan = 3, pady = 5)
  
  reviewButton <<- tkbutton(sidebarFrame, text = "Review Recordings", command = review_master_recordings)
  tkgrid(reviewButton, row = 10, column = 0, columnspan = 3, pady = 5)
  
  # New button below "Review Recordings" for averaging approved traces.
  avgApprovedTracesButton <<- tkbutton(sidebarFrame, text = "Average Approved Traces", command = average_selected_groups)
  tkgrid(avgApprovedTracesButton, row = 11, column = 0, columnspan = 3, pady = 5)
  
  downloadButton <<- tkbutton(sidebarFrame, text = "Download Data", command = download_data)
  tkgrid(downloadButton, row = 12, column = 0, columnspan = 3, pady = 5)
  
  tkfocus(tt)
}

# Launch the UI.
ABF_analysis_tk_A()



# This version loads as separated abfs and does not allow averaging between separated loaded abfs

# Remove all objects from the environment
rm(list = ls(all = TRUE))

# -----------------------------------------------------------------------------
# Package loading function and required packages:
# -----------------------------------------------------------------------------
load_required_packages <- function(packages) {
  new.packages <- packages[!(packages %in% installed.packages()[, "Package"])]
  if (length(new.packages)) install.packages(new.packages)
  invisible(lapply(packages, library, character.only = TRUE))
}
required.packages <- c("robustbase", "minpack.lm", "Rcpp", "signal", "readABF", "tcltk", "tkrplot")
load_required_packages(required.packages)

# -----------------------------------------------------------------------------
# Helper Functions
# -----------------------------------------------------------------------------

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
  if (experiment == "Voltage Clamp") {
    idx <- grep("A", channelUnits, ignore.case = TRUE)
  } else if (experiment == "Current Clamp") {
    idx <- grep("V", channelUnits, ignore.case = TRUE)
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
    return("Data is consistent")
  } else {
    error_msgs <- c()
    if (!dt_good) error_msgs <- c(error_msgs, paste("Inconsistent dt values:", paste(dt_values, collapse = ", ")))
    if (!unit_good) error_msgs <- c(error_msgs, paste("Inconsistent Units:", paste(unit_values, collapse = ", ")))
    if (!traces_good) error_msgs <- c(error_msgs, paste("Inconsistent Traces:", paste(traces_values, collapse = ", ")))
    return(paste(error_msgs, collapse = "; "))
  }
}

egs_plot <- function(x, y, sign = -1, xlim = NULL, ylim = NULL, lwd = 1, show_text = FALSE, 
                     xbar = 100, ybar = 50, color = "#4C77BB", show_bar = FALSE, cex = 0.6) {
  if (is.null(ylim))
    ylim <- if (sign == 1) c(0, max(y)) else c(-max(-y), 0)
  if (is.null(xlim))
    xlim <- c(min(x), max(x))
  idx1 <- which.min(abs(x - xlim[1]))
  idx2 <- which.min(abs(x - xlim[2]))
  plot(x[idx1:idx2], y[idx1:idx2], type = "l", col = color,
       xlim = xlim, ylim = ylim, bty = "n", lwd = lwd, lty = 1,
       axes = FALSE, frame = FALSE, xlab = "", ylab = "")
  if (show_bar) {
    ybar_start <- min(ylim) + (max(ylim) - min(ylim)) / 20
    x_start <- max(xlim) - xbar - 50
    y_start <- ybar_start
    x_end <- x_start + xbar
    y_end <- y_start + ybar
    segments(x_start, y_start, x_end, y_start, lwd = lwd, col = "black")
    segments(x_start, y_start, x_start, y_end, lwd = lwd, col = "black")
    if (show_text) {
      text(x = (x_start + x_end) / 2, y = y_start - ybar / 20, 
           labels = paste(xbar, "ms"), adj = c(0.5, 1), cex = cex)
      text(x = x_start - xbar / 4, y = (y_start + y_end) / 2,  
           labels = paste(ybar, "pA"), adj = c(0.5, 0.5), srt = 90, cex = cex)
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

abf_averages <- function(datasets, 
                         baseline = 100, 
                         stimulation_time = 350, 
                         traces2average = NULL,
                         dataCol = 1, 
                         ylim = NULL, 
                         xlim = NULL, 
                         color = "#4C77BB", 
                         xbar = 100, 
                         ybar = 50, 
                         width = 5.25, 
                         height = 2.75, 
                         save = FALSE) {
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
  responses0_mean <- if (is.null(traces2average)) {
    lapply(seq_len(N), function(iii) apply(responses0[[iii]], 1, mean))
  } else {
    lapply(seq_len(N), function(iii)
      apply(responses0[[iii]][, traces2average[[iii]], drop = FALSE], 1, mean))
  }
  time <- lapply(seq_len(N), function(iii) {
    dt_val <- sampling_intervals[iii]
    seq(0, by = dt_val, length.out = length(responses0_mean[[iii]]))
  })
  par(mfrow = c(1, N))
  show_bar <- rep(FALSE, N)
  if (N > 0) show_bar[N] <- TRUE
  for (ii in seq_len(N)) {
    for (ii in seq_len(N)) {
      egs_plot(x = time[[ii]], y = responses0_mean[[ii]], color = "darkgrey",
               show_bar = FALSE, show_text = FALSE)
    }
  }
  if (save) {
    warning("save_graph not implemented in this UI example")
  }
  return(list(raw_data = responses,
              baseline_corrected_data = responses0,
              baseline_corrected_mean_data = responses0_mean,
              datasets = datasets))
}

traces2average <<- list()

review_recordings <- function() {
  if (!exists("abf_analysis_result", envir = .GlobalEnv)) {
    tkmessageBox(message = "No analysis result available for review.")
    return()
  }
  result <- get("abf_analysis_result", envir = .GlobalEnv)
  datasets <- result$datasets
  traces2average <<- vector("list", length = length(datasets))
  for (i in seq_along(datasets)) { 
    traces2average[[i]] <<- integer(0) 
  }
  current_dataset <<- 1
  current_trace <<- 1
  children <- as.character(tkwinfo("children", plotPanel))
  for (child in children) {
    tryCatch({
      tkdestroy(.Tk.ID[[child]])
    }, error = function(e) {}, silent = TRUE)
  }

  reviewFrame <<- tkframe(plotPanel)
  tkgrid(reviewFrame, row = 0, column = 0, sticky = "nsew")
  current_filename <- names(datasets)[current_dataset]
  infoLabel <<- tklabel(reviewFrame, text = paste(current_filename, "trace", current_trace))
  tkgrid(infoLabel, row = 0, column = 0, columnspan = 2)
  reviewPlot <<- tkrplot(reviewFrame, fun = function() {
    ds <- datasets[[current_dataset]]
    if (current_trace > length(ds$data)) {
      plot.new()
      text(0.5, 0.5, paste("No more recordings in", current_filename))
    } else {
      trace_matrix <- ds$data[[current_trace]]
      data_column <- as.numeric(tclvalue(dataColVar))
      if (is.na(data_column) || data_column < 1 || data_column > ncol(trace_matrix)) {
        data_column <- 1
      }
      dt_val <- ds$samplingIntervalInSec * 1000
      time <- seq(0, by = dt_val, length.out = nrow(trace_matrix))
      trace <- trace_matrix[, data_column]
      plot(time, trace, col = "darkgrey", xlab = "Time (ms)", ylab = tclvalue(unitVar),
           type = "l", bty = "l", axes = FALSE,
           main = paste(current_filename, "trace", current_trace))
      axis(1)
      axis(2, las = 1)
    }
  }, hscale = 1, vscale = 1)
  tkgrid(reviewPlot, row = 1, column = 0, columnspan = 2)
  acceptButton <- tkbutton(reviewFrame, text = "Accept", command = function() {
    traces2average[[current_dataset]] <<- c(traces2average[[current_dataset]], current_trace)
    tkconfigure(acceptButton, state = "disabled", relief = "sunken")
    tkconfigure(rejectButton, state = "normal", relief = "raised")
  })
  tkgrid(acceptButton, row = 2, column = 0)
  rejectButton <- tkbutton(reviewFrame, text = "Reject", command = function() {
    tkconfigure(rejectButton, state = "disabled", relief = "sunken")
    tkconfigure(acceptButton, state = "normal", relief = "raised")
  })
  tkgrid(rejectButton, row = 2, column = 1)
  nextTraceButton <- tkbutton(reviewFrame, text = "Next Recording", command = function() {
    tkconfigure(acceptButton, state = "normal", relief = "raised")
    tkconfigure(rejectButton, state = "normal", relief = "raised")
    ds <- datasets[[current_dataset]]
    numRecordings <- length(ds$data)
    if (current_trace < numRecordings) {
      current_trace <<- current_trace + 1
    } else {
      tkmessageBox(message = paste("Finished reviewing", current_filename))
      if (current_dataset < length(datasets)) {
        current_dataset <<- current_dataset + 1
        current_trace <<- 1
      } else {
        tkmessageBox(message = "Review complete. Approved recordings are in 'traces2average'.")
        return()
      }
    }
    current_filename <- names(datasets)[current_dataset]
    tkconfigure(infoLabel, text = paste(current_filename, "trace", current_trace))
    tkrreplot(reviewPlot)
  })
  tkgrid(nextTraceButton, row = 3, column = 0, columnspan = 2)
}

download_data <- function() {
  if (!exists("abf_analysis_result", envir = .GlobalEnv)) {
    tkmessageBox(message = "No analysis result available for download.")
    return()
  }
  result <- get("abf_analysis_result", envir = .GlobalEnv)
  abf_out <- result  # structure containing baseline_corrected_mean_data
  dataList <- abf_out$baseline_corrected_mean_data
  if (length(dataList) == 0) {
    tkmessageBox(message = "No averaged data available.")
    return()
  }
  dt_val <- result$datasets[[1]]$samplingIntervalInSec * 1000
  time_vec <- seq(0, by = dt_val, length.out = length(dataList[[1]]))
  # Construct a data frame: first column is Time; next columns (one per ABF file) are the averaged traces,
  # gathered by sapply as requested.
  avg_mat <- sapply(1:length(dataList), function(ii) dataList[[ii]])
  df <- data.frame(Time = time_vec, avg_mat)
  download_folder <- tclvalue(folderPathVar)
  if(nchar(download_folder) == 0) {
    tkmessageBox(message = "No folder selected for download.")
    return()
  }
  file_path <- file.path(download_folder, "averaged_data.csv")
  write.csv(df, file = file_path, row.names = FALSE)
  tkmessageBox(message = paste("Averaged data saved to", file_path))
}

ABF_analysis_tk_B <- function() {
  tt <- tktoplevel()
  tkwm.title(tt, "ABF Analysis")
  sidebarFrame <- tkframe(tt)
  mainFrame <- tkframe(tt)
  tkgrid(sidebarFrame, row = 0, column = 0, sticky = "ns")
  tkgrid(mainFrame, row = 0, column = 1, sticky = "nsew")
  tkgrid.rowconfigure(tt, 0, weight = 1)
  tkgrid.columnconfigure(tt, 1, weight = 1)
  
  # Save mainFrame as the global plot panel.
  plotPanel <<- mainFrame
  
  folderLabel <- tklabel(sidebarFrame, text = "Select ABF Folder:")
  tkgrid(folderLabel, row = 0, column = 0, sticky = "w")
  folderPathVar <<- tclVar("")
  folderEntry <- tkentry(sidebarFrame, textvariable = folderPathVar, width = 30)
  tkgrid(folderEntry, row = 0, column = 1, sticky = "w")
  
  browseFolderButton <- tkbutton(sidebarFrame, text = "Browse", command = function() {
    folderPath <- tclvalue(tkchooseDirectory())
    if (nchar(folderPath) > 0) {
      tclvalue(folderPathVar) <<- folderPath
      abf_list <- list.files(path = folderPath, pattern = "\\.abf$", ignore.case = TRUE)
      if (length(abf_list) == 0) {
        tkmessageBox(message = "No ABF files found in the selected folder.")
      } else {
        tkdelete(abfListBox, 0, "end")
        for (f in abf_list) { tkinsert(abfListBox, "end", f) }
        firstFilePath <- file.path(folderPath, abf_list[1])
        ds <- readABF(firstFilePath)
        dummy_result <- list(metadata = list(extract_metadata(ds)))
        updateAdditionalParams(dummy_result)
      }
    }
  })
  tkgrid(browseFolderButton, row = 0, column = 2, padx = 5)
  
  abfListLabel <- tklabel(sidebarFrame, text = "ABF Files:")
  tkgrid(abfListLabel, row = 1, column = 0, sticky = "w", pady = 5)
  abfListBox <<- tklistbox(sidebarFrame, height = 5, width = 40, selectmode = "multiple")
  tkgrid(abfListBox, row = 2, column = 0, columnspan = 3, pady = 5, sticky = "we")
  
  paramFrame <- tkframe(sidebarFrame)
  tkgrid(paramFrame, row = 3, column = 0, columnspan = 3, sticky = "w")
  
  experimentVar <<- tclVar("Voltage Clamp")
  unitVar <<- tclVar("")
  dataColVar <<- tclVar("")
  dtVar <<- tclVar("")
  ntracesVar <<- tclVar("")
  
  tkgrid(tklabel(paramFrame, text = "Experiment:"), row = 0, column = 0, sticky = "w")
  experimentCombo <- ttkcombobox(paramFrame, textvariable = experimentVar, values = c("Voltage Clamp", "Current Clamp"), width = 15)
  tkgrid(experimentCombo, row = 0, column = 1, sticky = "w")
  
  tkgrid(tklabel(paramFrame, text = "Units:"), row = 1, column = 0, sticky = "w")
  unitEntry <- tkentry(paramFrame, textvariable = unitVar, width = 10, state = "readonly")
  tkgrid(unitEntry, row = 1, column = 1, sticky = "w")
  
  tkgrid(tklabel(paramFrame, text = "Data Column:"), row = 2, column = 0, sticky = "w")
  dataColEntry <- tkentry(paramFrame, textvariable = dataColVar, width = 10, state = "readonly")
  tkgrid(dataColEntry, row = 2, column = 1, sticky = "w")
  
  tkgrid(tklabel(paramFrame, text = "dt (ms):"), row = 3, column = 0, sticky = "w")
  dtEntry <- tkentry(paramFrame, textvariable = dtVar, width = 10, state = "readonly")
  tkgrid(dtEntry, row = 3, column = 1, sticky = "w")
  
  tkgrid(tklabel(paramFrame, text = "Traces:"), row = 4, column = 0, sticky = "w")
  ntracesEntry <- tkentry(paramFrame, textvariable = ntracesVar, width = 10, state = "readonly")
  tkgrid(ntracesEntry, row = 4, column = 1, sticky = "w")
  
  baselineVar <<- tclVar("100")
  stimTimeVar <<- tclVar("350")
  xbarVar <<- tclVar("100")
  ybarVar <<- tclVar("50")
  
  tkgrid(tklabel(sidebarFrame, text = "Baseline:"), row = 4, column = 0, sticky = "w")
  tkgrid(tkentry(sidebarFrame, textvariable = baselineVar, width = 10), row = 4, column = 1, sticky = "w")
  tkgrid(tklabel(sidebarFrame, text = "Stimulation Time:"), row = 5, column = 0, sticky = "w")
  tkgrid(tkentry(sidebarFrame, textvariable = stimTimeVar, width = 10), row = 5, column = 1, sticky = "w")
  tkgrid(tklabel(sidebarFrame, text = "x-bar length:"), row = 6, column = 0, sticky = "w")
  tkgrid(tkentry(sidebarFrame, textvariable = xbarVar, width = 10), row = 6, column = 1, sticky = "w")
  tkgrid(tklabel(sidebarFrame, text = "y-bar length:"), row = 7, column = 0, sticky = "w")
  tkgrid(tkentry(sidebarFrame, textvariable = ybarVar, width = 10), row = 7, column = 1, sticky = "w")
  
  consoleText <<- tktext(sidebarFrame, width = 40, height = 4)
  tkgrid(consoleText, row = 8, column = 0, columnspan = 3, pady = 5)
  
  updateAdditionalParams <<- function(result) {
    if (!is.null(result) && length(result$metadata) >= 1) {
      meta1 <- result$metadata[[1]]
      tclvalue(dtVar) <<- as.character(meta1$samplingIntervalInSec * 1000)
      if (!is.null(meta1$header$lActualEpisodes)) {
        tclvalue(ntracesVar) <<- as.character(meta1$header$lActualEpisodes)
      } else {
        tclvalue(ntracesVar) <<- "N/A"
      }
      expType <- tclvalue(experimentVar)
      col_idx <- choose_data_column(meta1$channelUnits, expType)
      if (!is.na(col_idx)) {
        tclvalue(unitVar) <<- meta1$channelUnits[col_idx]
        tclvalue(dataColVar) <<- as.character(col_idx)
      } else {
        tclvalue(unitVar) <<- "N/A"
        tclvalue(dataColVar) <<- "N/A"
      }
    }
  }
  
  tkbind(experimentCombo, "<<ComboboxSelected>>", function() {
    if (exists("abf_analysis_result", envir = .GlobalEnv)) {
      result <- get("abf_analysis_result", envir = .GlobalEnv)
      updateAdditionalParams(result)
    }
  })
  
  runAnalysis <<- function() {
    folderPath <- tclvalue(folderPathVar)
    if (nchar(folderPath) == 0) {
      tkmessageBox(message = "Please select an ABF folder first.")
      return()
    }
    selIndices <- as.integer(tkcurselection(abfListBox))
    allFiles <- as.character(tkget(abfListBox, 0, "end"))
    abf_files <- if (length(selIndices) == 0) allFiles else allFiles[selIndices + 1]
    if (length(abf_files) == 0) {
      tkmessageBox(message = "No ABF files selected.")
      return()
    }
    result <- tryCatch({
      load_abf_data(abf_files = abf_files, abf_path = folderPath)
    }, error = function(e) {
      tkmessageBox(message = paste("Error during data loading:", e$message))
      return(NULL)
    })
    result <<- result
    if (!is.null(result)) {
      tkdelete(consoleText, "1.0", "end")
      tkinsert(consoleText, "end", paste("Data loaded. Processed", length(abf_files), "file(s)."))
      assign("abf_analysis_result", result, envir = .GlobalEnv)
      updateAdditionalParams(result)
      cons_msg <- check_consistency(result$metadata)
      if (cons_msg == "Data is consistent") {
        tkmessageBox(message = cons_msg)
      } else {
        tkmessageBox(message = paste("ERROR:", cons_msg))
      }
      tkconfigure(runAnalysisButton, text = "Load Data")
    }
  }
  
  runAnalysisButton <<- tkbutton(sidebarFrame, text = "Load Data", command = runAnalysis)
  tkgrid(runAnalysisButton, row = 9, column = 0, columnspan = 3, pady = 5)
  
  reviewButton <<- tkbutton(sidebarFrame, text = "Review Recordings", command = review_recordings)
  tkgrid(reviewButton, row = 10, column = 0, columnspan = 3, pady = 5)
  
  averageApprovedTraces <<- function() {
    if (length(traces2average) == 0 || all(sapply(traces2average, length) == 0)) {
      tkmessageBox(message = "No approved traces available. Please review recordings first.")
      return()
    }
    folderPath <- tclvalue(folderPathVar)
    if (nchar(folderPath) == 0) {
      tkmessageBox(message = "Please select an ABF folder first.")
      return()
    }
    selIndices <- as.integer(tkcurselection(abfListBox))
    allFiles <- as.character(tkget(abfListBox, 0, "end"))
    abf_files <- if (length(selIndices) == 0) allFiles else allFiles[selIndices + 1]
    if (length(abf_files) == 0) {
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
        color = "darkgrey", 
        xbar = xbar, ybar = ybar,
        width = 5.25, height = 2.75
      )
      abf_out
    }, error = function(e) {
      tkmessageBox(message = paste("Error during averaging of approved traces:", e$message))
      NULL
    })

    if (!is.null(result)) {
      tkdelete(consoleText, "1.0", "end")
      tkinsert(consoleText, "end", paste("Averaging on approved traces complete. Processed", length(abf_files), "file(s)."))
      
      # Update the global analysis result
      abf_analysis_result <<- result
      
      # Define drawPlot to display the averaged trace using scale bars only (no axes or title)
      drawPlot <- function() {
        if (exists("abf_analysis_result", envir = .GlobalEnv)) {
          result <- get("abf_analysis_result", envir = .GlobalEnv)
          datasets <- result$datasets
          traces <- result$baseline_corrected_mean_data
          if (length(datasets) > 0 && length(traces) > 0) {
            par(mfrow = c(1, length(traces)))
            # Determine shared limits
            all_y <- unlist(traces)
            shared_ylim <- range(all_y)
            shared_xlim <- range(unlist(lapply(seq_along(traces), function(i) {
              dt <- datasets[[i]]$samplingIntervalInSec * 1000
              seq(0, by = dt, length.out = length(traces[[i]]))
            })))
            for (i in seq_along(traces)) {
              dt_val <- datasets[[i]]$samplingIntervalInSec * 1000
              time <- seq(0, by = dt_val, length.out = length(traces[[i]]))
              show_bar <- (i == length(traces))
              egs_plot(x = time, y = traces[[i]], color = "darkgrey",
                       show_bar = show_bar, show_text = show_bar, 
                       xbar = xbar, ybar = ybar,
                       xlim = shared_xlim, ylim = shared_ylim)
            }
          } else {
            plot.new()
            text(0.5, 0.5, "No data available")
          }
        } else {
          plot.new()
          text(0.5, 0.5, "No analysis result to display")
        }
      }
      
      # If a plot widget already exists, destroy it.
      children <- as.character(tkwinfo("children", plotPanel))
      for (child in children) {
        tryCatch({
          tkdestroy(.Tk.ID[[child]])
        }, error = function(e) {}, silent = TRUE)
      }

      # Create a new frame to hold the averaged plot
      avgFrame <- tkframe(plotPanel)
      tkgrid(avgFrame, row = 0, column = 0, sticky = "nsew")

      # Use the same drawPlot function defined above
      reviewPlot <<- tkrplot(avgFrame, fun = drawPlot, hscale = 1, vscale = 1)
      tkgrid(reviewPlot, row = 0, column = 0)

    }
  }
  
  averageApprovedButton <<- tkbutton(sidebarFrame, text = "Average Approved Traces", command = averageApprovedTraces)
  tkgrid(averageApprovedButton, row = 11, column = 0, columnspan = 3, pady = 5)
  
  downloadButton <<- tkbutton(sidebarFrame, text = "Download Data", command = download_data)
  tkgrid(downloadButton, row = 12, column = 0, columnspan = 3, pady = 5)
  
  tkfocus(tt)
}

ABF_analysis_tk_B()




# # Remove all objects from the environment
# rm(list = ls(all = TRUE))

# # -----------------------------------------------------------------------------
# # Package loading function and required packages:
# # -----------------------------------------------------------------------------
# load_required_packages <- function(packages) {
#   new.packages <- packages[!(packages %in% installed.packages()[, "Package"])]
#   if (length(new.packages)) install.packages(new.packages)
#   invisible(lapply(packages, library, character.only = TRUE))
# }
# required.packages <- c("robustbase", "minpack.lm", "Rcpp", "signal", "readABF", "tcltk", "tkrplot")
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
#   if (experiment == "Voltage Clamp") {
#     idx <- grep("A", channelUnits, ignore.case = TRUE)
#   } else if (experiment == "Current Clamp") {
#     idx <- grep("V", channelUnits, ignore.case = TRUE)
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
#     return("Data is consistent")
#   } else {
#     error_msgs <- c()
#     if (!dt_good) error_msgs <- c(error_msgs, paste("Inconsistent dt values:", paste(dt_values, collapse = ", ")))
#     if (!unit_good) error_msgs <- c(error_msgs, paste("Inconsistent Units:", paste(unit_values, collapse = ", ")))
#     if (!traces_good) error_msgs <- c(error_msgs, paste("Inconsistent Traces:", paste(traces_values, collapse = ", ")))
#     return(paste(error_msgs, collapse = "; "))
#   }
# }

# egs_plot <- function(x, y, sign = -1, xlim = NULL, ylim = NULL, lwd = 1, show_text = FALSE, 
#                       xbar = 100, ybar = 50, color = "#4C77BB", show_bar = FALSE) {
#   if (is.null(ylim))
#     ylim <- if (sign == 1) c(0, max(y)) else c(-max(-y), 0)
#   if (is.null(xlim))
#     xlim <- c(min(x), max(x))
#   idx1 <- which.min(abs(x - xlim[1]))
#   idx2 <- which.min(abs(x - xlim[2]))
#   plot(x[idx1:idx2], y[idx1:idx2], type = "l", col = color,
#        xlim = xlim, ylim = ylim, bty = "n", lwd = lwd, lty = 1,
#        axes = FALSE, frame = FALSE, xlab = "", ylab = "")
#   if (show_bar) {
#     ybar_start <- min(ylim) + (max(ylim) - min(ylim)) / 20
#     x_start <- max(xlim) - xbar - 50
#     y_start <- ybar_start
#     x_end <- x_start + xbar
#     y_end <- y_start + ybar
#     segments(x_start, y_start, x_end, y_start, lwd = lwd, col = "black")
#     segments(x_start, y_start, x_start, y_end, lwd = lwd, col = "black")
#     if (show_text) {
#       text(x = (x_start + x_end) / 2, y = y_start - ybar / 20, 
#            labels = paste(xbar, "ms"), adj = c(0.5, 1))
#       text(x = x_start - xbar / 4, y = (y_start + y_end) / 2,  
#            labels = paste(ybar, "pA"), adj = c(0.5, 0.5), srt = 90)
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
#                          color = "#4C77BB", 
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
#     egs_plot(x = time[[ii]], y = responses0_mean[[ii]], ylim = ylim,
#              xlim = xlim, color = color, show_text = FALSE, show_bar = show_bar[ii])
#     axis(1)
#     axis(2)
#   }
#   if (save) {
#     warning("save_graph not implemented in this UI example")
#   }
#   return(list(raw_data = responses,
#               baseline_corrected_data = responses0,
#               baseline_corrected_mean_data = responses0_mean,
#               datasets = datasets))
# }

# traces2average <<- list()

# review_recordings <- function() {
#   if (!exists("abf_analysis_result", envir = .GlobalEnv)) {
#     tkmessageBox(message = "No analysis result available for review.")
#     return()
#   }
#   result <- get("abf_analysis_result", envir = .GlobalEnv)
#   datasets <- result$datasets
#   traces2average <<- vector("list", length = length(datasets))
#   for (i in seq_along(datasets)) { 
#     traces2average[[i]] <<- integer(0) 
#   }
#   current_dataset <<- 1
#   current_trace <<- 1
#   children <- as.character(tkwinfo("children", plotPanel))
#   if (length(children) > 0) {
#     for (child in children) tkdestroy(child)
#   }
#   reviewFrame <<- tkframe(plotPanel)
#   tkgrid(reviewFrame, row = 0, column = 0, sticky = "nsew")
#   current_filename <- names(datasets)[current_dataset]
#   infoLabel <<- tklabel(reviewFrame, text = paste(current_filename, "trace", current_trace))
#   tkgrid(infoLabel, row = 0, column = 0, columnspan = 2)
#   reviewPlot <<- tkrplot(reviewFrame, fun = function() {
#     ds <- datasets[[current_dataset]]
#     if (current_trace > length(ds$data)) {
#       plot.new()
#       text(0.5, 0.5, paste("No more recordings in", current_filename))
#     } else {
#       trace_matrix <- ds$data[[current_trace]]
#       data_column <- as.numeric(tclvalue(dataColVar))
#       if (is.na(data_column) || data_column < 1 || data_column > ncol(trace_matrix)) {
#         data_column <- 1
#       }
#       dt_val <- ds$samplingIntervalInSec * 1000
#       time <- seq(0, by = dt_val, length.out = nrow(trace_matrix))
#       trace <- trace_matrix[, data_column]
#       plot(time, trace, col = "darkgrey", xlab = "Time (ms)", ylab = tclvalue(unitVar),
#            type = "l", bty = "l", axes = FALSE,
#            main = paste(current_filename, "trace", current_trace))
#       axis(1)
#       axis(2, las = 1)
#     }
#   }, hscale = 1, vscale = 1)
#   tkgrid(reviewPlot, row = 1, column = 0, columnspan = 2)
#   acceptButton <- tkbutton(reviewFrame, text = "Accept", command = function() {
#     traces2average[[current_dataset]] <<- c(traces2average[[current_dataset]], current_trace)
#     tkconfigure(acceptButton, state = "disabled", relief = "sunken")
#     tkconfigure(rejectButton, state = "normal", relief = "raised")
#   })
#   tkgrid(acceptButton, row = 2, column = 0)
#   rejectButton <- tkbutton(reviewFrame, text = "Reject", command = function() {
#     tkconfigure(rejectButton, state = "disabled", relief = "sunken")
#     tkconfigure(acceptButton, state = "normal", relief = "raised")
#   })
#   tkgrid(rejectButton, row = 2, column = 1)
#   nextTraceButton <- tkbutton(reviewFrame, text = "Next Recording", command = function() {
#     tkconfigure(acceptButton, state = "normal", relief = "raised")
#     tkconfigure(rejectButton, state = "normal", relief = "raised")
#     ds <- datasets[[current_dataset]]
#     numRecordings <- length(ds$data)
#     if (current_trace < numRecordings) {
#       current_trace <<- current_trace + 1
#     } else {
#       tkmessageBox(message = paste("Finished reviewing", current_filename))
#       if (current_dataset < length(datasets)) {
#         current_dataset <<- current_dataset + 1
#         current_trace <<- 1
#       } else {
#         tkmessageBox(message = "Review complete. Approved recordings are in 'traces2average'.")
#         return()
#       }
#     }
#     current_filename <- names(datasets)[current_dataset]
#     tkconfigure(infoLabel, text = paste(current_filename, "trace", current_trace))
#     tkrreplot(reviewPlot)
#   })
#   tkgrid(nextTraceButton, row = 3, column = 0, columnspan = 2)
# }

# download_data <- function() {
#   if (!exists("abf_analysis_result", envir = .GlobalEnv)) {
#     tkmessageBox(message = "No analysis result available for download.")
#     return()
#   }
#   result <- get("abf_analysis_result", envir = .GlobalEnv)
#   abf_out <- result  # structure containing baseline_corrected_mean_data
#   dataList <- abf_out$baseline_corrected_mean_data
#   if (length(dataList) == 0) {
#     tkmessageBox(message = "No averaged data available.")
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
#     tkmessageBox(message = "No folder selected for download.")
#     return()
#   }
#   file_path <- file.path(download_folder, "averaged_data.csv")
#   write.csv(df, file = file_path, row.names = FALSE)
#   tkmessageBox(message = paste("Averaged data saved to", file_path))
# }

# ABF_analysis_tk <- function() {
#   tt <- tktoplevel()
#   tkwm.title(tt, "ABF Analysis")
#   sidebarFrame <- tkframe(tt)
#   mainFrame <- tkframe(tt)
#   tkgrid(sidebarFrame, row = 0, column = 0, sticky = "ns")
#   tkgrid(mainFrame, row = 0, column = 1, sticky = "nsew")
#   tkgrid.rowconfigure(tt, 0, weight = 1)
#   tkgrid.columnconfigure(tt, 1, weight = 1)
  
#   # Save mainFrame as the global plot panel.
#   plotPanel <<- mainFrame
  
#   folderLabel <- tklabel(sidebarFrame, text = "Select ABF Folder:")
#   tkgrid(folderLabel, row = 0, column = 0, sticky = "w")
#   folderPathVar <<- tclVar("")
#   folderEntry <- tkentry(sidebarFrame, textvariable = folderPathVar, width = 30)
#   tkgrid(folderEntry, row = 0, column = 1, sticky = "w")
  
#   browseFolderButton <- tkbutton(sidebarFrame, text = "Browse", command = function() {
#     folderPath <- tclvalue(tkchooseDirectory())
#     if (nchar(folderPath) > 0) {
#       tclvalue(folderPathVar) <<- folderPath
#       abf_list <- list.files(path = folderPath, pattern = "\\.abf$", ignore.case = TRUE)
#       if (length(abf_list) == 0) {
#         tkmessageBox(message = "No ABF files found in the selected folder.")
#       } else {
#         tkdelete(abfListBox, 0, "end")
#         for (f in abf_list) { tkinsert(abfListBox, "end", f) }
#         firstFilePath <- file.path(folderPath, abf_list[1])
#         ds <- readABF(firstFilePath)
#         dummy_result <- list(metadata = list(extract_metadata(ds)))
#         updateAdditionalParams(dummy_result)
#       }
#     }
#   })
#   tkgrid(browseFolderButton, row = 0, column = 2, padx = 5)
  
#   abfListLabel <- tklabel(sidebarFrame, text = "ABF Files:")
#   tkgrid(abfListLabel, row = 1, column = 0, sticky = "w", pady = 5)
#   abfListBox <<- tklistbox(sidebarFrame, height = 5, width = 40, selectmode = "multiple")
#   tkgrid(abfListBox, row = 2, column = 0, columnspan = 3, pady = 5, sticky = "we")
  
#   paramFrame <- tkframe(sidebarFrame)
#   tkgrid(paramFrame, row = 3, column = 0, columnspan = 3, sticky = "w")
  
#   experimentVar <<- tclVar("Voltage Clamp")
#   unitVar <<- tclVar("")
#   dataColVar <<- tclVar("")
#   dtVar <<- tclVar("")
#   ntracesVar <<- tclVar("")
  
#   tkgrid(tklabel(paramFrame, text = "Experiment:"), row = 0, column = 0, sticky = "w")
#   experimentCombo <- ttkcombobox(paramFrame, textvariable = experimentVar, values = c("Voltage Clamp", "Current Clamp"), width = 15)
#   tkgrid(experimentCombo, row = 0, column = 1, sticky = "w")
  
#   tkgrid(tklabel(paramFrame, text = "Units:"), row = 1, column = 0, sticky = "w")
#   unitEntry <- tkentry(paramFrame, textvariable = unitVar, width = 10, state = "readonly")
#   tkgrid(unitEntry, row = 1, column = 1, sticky = "w")
  
#   tkgrid(tklabel(paramFrame, text = "Data Column:"), row = 2, column = 0, sticky = "w")
#   dataColEntry <- tkentry(paramFrame, textvariable = dataColVar, width = 10, state = "readonly")
#   tkgrid(dataColEntry, row = 2, column = 1, sticky = "w")
  
#   tkgrid(tklabel(paramFrame, text = "dt (ms):"), row = 3, column = 0, sticky = "w")
#   dtEntry <- tkentry(paramFrame, textvariable = dtVar, width = 10, state = "readonly")
#   tkgrid(dtEntry, row = 3, column = 1, sticky = "w")
  
#   tkgrid(tklabel(paramFrame, text = "Traces:"), row = 4, column = 0, sticky = "w")
#   ntracesEntry <- tkentry(paramFrame, textvariable = ntracesVar, width = 10, state = "readonly")
#   tkgrid(ntracesEntry, row = 4, column = 1, sticky = "w")
  
#   baselineVar <<- tclVar("100")
#   stimTimeVar <<- tclVar("350")
#   xbarVar <<- tclVar("100")
#   ybarVar <<- tclVar("50")
  
#   tkgrid(tklabel(sidebarFrame, text = "Baseline:"), row = 4, column = 0, sticky = "w")
#   tkgrid(tkentry(sidebarFrame, textvariable = baselineVar, width = 10), row = 4, column = 1, sticky = "w")
#   tkgrid(tklabel(sidebarFrame, text = "Stimulation Time:"), row = 5, column = 0, sticky = "w")
#   tkgrid(tkentry(sidebarFrame, textvariable = stimTimeVar, width = 10), row = 5, column = 1, sticky = "w")
#   tkgrid(tklabel(sidebarFrame, text = "x-bar length:"), row = 6, column = 0, sticky = "w")
#   tkgrid(tkentry(sidebarFrame, textvariable = xbarVar, width = 10), row = 6, column = 1, sticky = "w")
#   tkgrid(tklabel(sidebarFrame, text = "y-bar length:"), row = 7, column = 0, sticky = "w")
#   tkgrid(tkentry(sidebarFrame, textvariable = ybarVar, width = 10), row = 7, column = 1, sticky = "w")
  
#   consoleText <<- tktext(sidebarFrame, width = 40, height = 4)
#   tkgrid(consoleText, row = 8, column = 0, columnspan = 3, pady = 5)
  
#   updateAdditionalParams <<- function(result) {
#     if (!is.null(result) && length(result$metadata) >= 1) {
#       meta1 <- result$metadata[[1]]
#       tclvalue(dtVar) <<- as.character(meta1$samplingIntervalInSec * 1000)
#       if (!is.null(meta1$header$lActualEpisodes)) {
#         tclvalue(ntracesVar) <<- as.character(meta1$header$lActualEpisodes)
#       } else {
#         tclvalue(ntracesVar) <<- "N/A"
#       }
#       expType <- tclvalue(experimentVar)
#       col_idx <- choose_data_column(meta1$channelUnits, expType)
#       if (!is.na(col_idx)) {
#         tclvalue(unitVar) <<- meta1$channelUnits[col_idx]
#         tclvalue(dataColVar) <<- as.character(col_idx)
#       } else {
#         tclvalue(unitVar) <<- "N/A"
#         tclvalue(dataColVar) <<- "N/A"
#       }
#     }
#   }
  
#   tkbind(experimentCombo, "<<ComboboxSelected>>", function() {
#     if (exists("abf_analysis_result", envir = .GlobalEnv)) {
#       result <- get("abf_analysis_result", envir = .GlobalEnv)
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
#     allFiles <- as.character(tkget(abfListBox, 0, "end"))
#     abf_files <- if (length(selIndices) == 0) allFiles else allFiles[selIndices + 1]
#     if (length(abf_files) == 0) {
#       tkmessageBox(message = "No ABF files selected.")
#       return()
#     }
#     result <- tryCatch({
#       load_abf_data(abf_files = abf_files, abf_path = folderPath)
#     }, error = function(e) {
#       tkmessageBox(message = paste("Error during data loading:", e$message))
#       return(NULL)
#     })
#     if (!is.null(result)) {
#       tkdelete(consoleText, "1.0", "end")
#       tkinsert(consoleText, "end", paste("Data loaded. Processed", length(abf_files), "file(s)."))
#       assign("abf_analysis_result", result, envir = .GlobalEnv)
#       updateAdditionalParams(result)
#       cons_msg <- check_consistency(result$metadata)
#       if (cons_msg == "Data is consistent") {
#         tkmessageBox(message = cons_msg)
#       } else {
#         tkmessageBox(message = paste("ERROR:", cons_msg))
#       }
#       tkconfigure(runAnalysisButton, text = "Load Data")
#     }
#   }
  
#   runAnalysisButton <<- tkbutton(sidebarFrame, text = "Load Data", command = runAnalysis)
#   tkgrid(runAnalysisButton, row = 9, column = 0, columnspan = 3, pady = 5)
  
#   reviewButton <<- tkbutton(sidebarFrame, text = "Review Recordings", command = review_recordings)
#   tkgrid(reviewButton, row = 10, column = 0, columnspan = 3, pady = 5)
  
#   averageApprovedTraces <<- function() {
#     if (length(traces2average) == 0 || all(sapply(traces2average, length) == 0)) {
#       tkmessageBox(message = "No approved traces available. Please review recordings first.")
#       return()
#     }
#     folderPath <- tclvalue(folderPathVar)
#     if (nchar(folderPath) == 0) {
#       tkmessageBox(message = "Please select an ABF folder first.")
#       return()
#     }
#     selIndices <- as.integer(tkcurselection(abfListBox))
#     allFiles <- as.character(tkget(abfListBox, 0, "end"))
#     abf_files <- if (length(selIndices) == 0) allFiles else allFiles[selIndices + 1]
#     if (length(abf_files) == 0) {
#       tkmessageBox(message = "No ABF files selected.")
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
#         color = "darkgrey", 
#         xbar = xbar, ybar = ybar,
#         width = 5.25, height = 2.75
#       )
#       abf_out
#     }, error = function(e) {
#       tkmessageBox(message = paste("Error during averaging of approved traces:", e$message))
#       NULL
#     })

#     if (!is.null(result)) {
#       tkdelete(consoleText, "1.0", "end")
#       tkinsert(consoleText, "end", paste("Averaging on approved traces complete. Processed", length(abf_files), "file(s)."))
      
#       # Update the global analysis result
#       abf_analysis_result <<- result
      
#       # Define drawPlot to display the averaged trace using scale bars only (no axes or title)
#       drawPlot <- function() {
#         if (exists("abf_analysis_result", envir = .GlobalEnv)) {
#           result <- get("abf_analysis_result", envir = .GlobalEnv)
#           datasets <- result$datasets
#           if (length(datasets) >= 1) {
#             iii <- 1
#             dt_val <- datasets[[iii]]$samplingIntervalInSec * 1000
#             time <- seq(0, by = dt_val, length.out = nrow(datasets[[iii]]$data[[1]]))
#             channel_index <- as.numeric(tclvalue(dataColVar))
#             if (is.na(channel_index) || channel_index < 1 || channel_index > length(datasets[[iii]]$data)) {
#               channel_index <- 1
#             }
#             trace <- datasets[[iii]]$data[[channel_index]][, 1]
#             # Draw the averaged trace using egs_plot with scale bars and labels (axes are not added)
#             egs_plot(x = time, y = trace, show_bar = TRUE, show_text = TRUE, color = "darkgrey")
#           } else {
#             plot.new()
#             text(0.5, 0.5, "No data available")
#           }
#         } else {
#           plot.new()
#           text(0.5, 0.5, "No analysis result to display")
#         }
#       }
      
#       # If a plot widget already exists, destroy it.
#       if (exists("plotWidget", envir = .GlobalEnv)) tkdestroy(plotWidget)
      
#       # Embed the averaged plot into the right-hand panel
#       plotWidget <<- tkrplot(plotPanel, fun = drawPlot, hscale = 1, vscale = 1)
#       tkgrid(plotWidget, row = 0, column = 0, sticky = "nsew")
#     }
#   }
  
#   averageApprovedButton <<- tkbutton(sidebarFrame, text = "Average Approved Traces", command = averageApprovedTraces)
#   tkgrid(averageApprovedButton, row = 11, column = 0, columnspan = 3, pady = 5)
  
#   downloadButton <<- tkbutton(sidebarFrame, text = "Download Data", command = download_data)
#   tkgrid(downloadButton, row = 12, column = 0, columnspan = 3, pady = 5)
  
#   tkfocus(tt)
# }

# ABF_analysis_tk()

# # Remove all objects from the environment
# rm(list = ls(all = TRUE))

# # -----------------------------------------------------------------------------
# # Package loading function and required packages:
# # -----------------------------------------------------------------------------
# load_required_packages <- function(packages) {
#   new.packages <- packages[!(packages %in% installed.packages()[, "Package"])]
#   if (length(new.packages)) install.packages(new.packages)
#   invisible(lapply(packages, library, character.only = TRUE))
# }
# required.packages <- c("robustbase", "minpack.lm", "Rcpp", "signal", "readABF", "tcltk", "tkrplot")
# load_required_packages(required.packages)

# # -----------------------------------------------------------------------------
# # Helper Functions
# # -----------------------------------------------------------------------------

# # Extract key metadata from an ABF dataset.
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

# # Choose data column index based on experiment type.
# choose_data_column <- function(channelUnits, experiment) {
#   if (experiment == "Voltage Clamp") {
#     idx <- grep("A", channelUnits, ignore.case = TRUE)
#   } else if (experiment == "Current Clamp") {
#     idx <- grep("V", channelUnits, ignore.case = TRUE)
#   } else {
#     idx <- integer(0)
#   }
#   if (length(idx) > 0) return(idx[1])
#   else return(NA)
# }

# # Consistency check: Verify that dt, chosen unit, and number of traces are identical.
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
#     return("Data is consistent")
#   } else {
#     error_msgs <- c()
#     if (!dt_good) error_msgs <- c(error_msgs, paste("Inconsistent dt values:", paste(dt_values, collapse = ", ")))
#     if (!unit_good) error_msgs <- c(error_msgs, paste("Inconsistent Units:", paste(unit_values, collapse = ", ")))
#     if (!traces_good) error_msgs <- c(error_msgs, paste("Inconsistent Traces:", paste(traces_values, collapse = ", ")))
#     return(paste(error_msgs, collapse = "; "))
#   }
# }

# # Helper function for drawing scale bars (legacy; used by egs_plot)
# egs_plot <- function(x, y, sign = -1, xlim = NULL, ylim = NULL, lwd = 1, 
#                      show_text = FALSE, height = 4, width = 2.5, 
#                      xbar = 100, ybar = 50, color = "#4C77BB", show_bar = FALSE) {
#   if (is.null(ylim))
#     ylim <- if (sign == 1) c(0, max(y)) else c(-max(-y), 0)
#   if (is.null(xlim))
#     xlim <- c(min(x), max(x))
  
#   idx1 <- which.min(abs(x - xlim[1]))
#   idx2 <- which.min(abs(x - xlim[2]))
  
#   plot(x[idx1:idx2], y[idx1:idx2], type = "l", col = color,
#        xlim = xlim, ylim = ylim, bty = "n", lwd = lwd, lty = 1,
#        axes = FALSE, frame = FALSE, xlab = "", ylab = "")
  
#   if (show_bar) {
#     ybar_start <- min(ylim) + (max(ylim) - min(ylim)) / 20
#     x_start <- max(xlim) - xbar - 50
#     y_start <- ybar_start
#     x_end <- x_start + xbar
#     y_end <- y_start + ybar
#     segments(x_start, y_start, x_end, y_start, lwd = lwd, col = "black")
#     segments(x_start, y_start, x_start, y_end, lwd = lwd, col = "black")
#     if (show_text) {
#       text(x = (x_start + x_end) / 2, y = y_start - ybar / 20, 
#            labels = paste(xbar, "ms"), adj = c(0.5, 1))
#       text(x = x_start - xbar / 4, y = (y_start + y_end) / 2,  
#            labels = paste(ybar, "pA"), adj = c(0.5, 0.5), srt = 90)
#     }
#   }
# }

# # New function to load ABF data.
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
# # Revised abf_averages Function
# # -----------------------------------------------------------------------------
# abf_averages <- function(datasets, 
#                          baseline = 100, 
#                          stimulation_time = 350, 
#                          traces2average = NULL,
#                          dataCol = 1, 
#                          ylim = NULL, 
#                          xlim = NULL, 
#                          color = "#4C77BB", 
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
#     egs_plot(x = time[[ii]], y = responses0_mean[[ii]], ylim = ylim,
#              xlim = xlim, color = color, width = width, height = height,
#              show_text = FALSE, show_bar = show_bar[ii])
#     axis(1)
#     axis(2)
#   }
  
#   if (save) {
#     warning("no save_graph in this UI example yet")
#   }
  
#   return(list(raw_data = responses,
#               baseline_corrected_data = responses0,
#               baseline_corrected_mean_data = responses0_mean,
#               datasets = datasets))
# }

# # -----------------------------------------------------------------------------
# # Global variable for approved trace indices (one per dataset)
# # -----------------------------------------------------------------------------
# traces2average <<- list()

# # -----------------------------------------------------------------------------
# # Revised review_recordings Function (embedded into the right-hand plot panel)
# # -----------------------------------------------------------------------------
# review_recordings <- function() {
#   if (!exists("abf_analysis_result", envir = .GlobalEnv)) {
#     tkmessageBox(message = "No analysis result available for review.")
#     return()
#   }
#   result <- get("abf_analysis_result", envir = .GlobalEnv)
#   datasets <- result$datasets
#   traces2average <<- vector("list", length = length(datasets))
#   for (i in seq_along(datasets)) { 
#     traces2average[[i]] <<- integer(0) 
#   }
  
#   current_dataset <<- 1
#   current_trace <<- 1
  
#   # Clear existing content from the plot panel.
#   children <- as.character(tkwinfo("children", plotPanel))
#   if (length(children) > 0) {
#     for (child in children) tkdestroy(child)
#   }
  
#   # Create a review frame inside the plot panel.
#   reviewFrame <<- tkframe(plotPanel)
#   tkgrid(reviewFrame, row = 0, column = 0, sticky = "nsew")
  
#   current_filename <- names(datasets)[current_dataset]
#   infoLabel <<- tklabel(reviewFrame, text = paste(current_filename, "trace", current_trace))
#   tkgrid(infoLabel, row = 0, column = 0, columnspan = 2)
  
#   # Create a review plot using standard plot() with axes added.
#   reviewPlot <<- tkrplot(reviewFrame, fun = function() {
#     ds <- datasets[[current_dataset]]
#     if (current_trace > length(ds$data)) {
#       plot.new()
#       text(0.5, 0.5, paste("No more recordings in", current_filename))
#     } else {
#       trace_matrix <- ds$data[[current_trace]]
#       data_column <- as.numeric(tclvalue(dataColVar))
#       if (is.na(data_column) || data_column < 1 || data_column > ncol(trace_matrix)) {
#         data_column <- 1
#       }
#       dt_val <- ds$samplingIntervalInSec * 1000
#       time <- seq(0, by = dt_val, length.out = nrow(trace_matrix))
#       trace <- trace_matrix[, data_column]
#       plot(time, trace, col = "darkgrey", xlab = "Time (ms)", ylab = tclvalue(unitVar),
#            type = "l", bty = "l", axes = FALSE, 
#            main = paste(current_filename, "trace", current_trace))
#       axis(1)
#       axis(2, las = 1)
#     }
#   }, hscale = 1, vscale = 1)
#   tkgrid(reviewPlot, row = 1, column = 0, columnspan = 2)
  
#   acceptButton <- tkbutton(reviewFrame, text = "Accept", command = function() {
#     traces2average[[current_dataset]] <<- c(traces2average[[current_dataset]], current_trace)
#     tkconfigure(acceptButton, state = "disabled", relief = "sunken")
#     tkconfigure(rejectButton, state = "normal", relief = "raised")
#   })
#   tkgrid(acceptButton, row = 2, column = 0)
  
#   rejectButton <- tkbutton(reviewFrame, text = "Reject", command = function() {
#     tkconfigure(rejectButton, state = "disabled", relief = "sunken")
#     tkconfigure(acceptButton, state = "normal", relief = "raised")
#   })
#   tkgrid(rejectButton, row = 2, column = 1)
  
#   nextTraceButton <- tkbutton(reviewFrame, text = "Next Recording", command = function() {
#     tkconfigure(acceptButton, state = "normal", relief = "raised")
#     tkconfigure(rejectButton, state = "normal", relief = "raised")
#     ds <- datasets[[current_dataset]]
#     numRecordings <- length(ds$data)
#     if (current_trace < numRecordings) {
#       current_trace <<- current_trace + 1
#     } else {
#       tkmessageBox(message = paste("Finished reviewing", current_filename))
#       if (current_dataset < length(datasets)) {
#         current_dataset <<- current_dataset + 1
#         current_trace <<- 1
#       } else {
#         tkmessageBox(message = "Review complete. Approved recordings are in 'traces2average'.")
#         return()
#       }
#     }
#     current_filename <- names(datasets)[current_dataset]
#     tkconfigure(infoLabel, text = paste(current_filename, "trace", current_trace))
#     tkrreplot(reviewPlot)
#   })
#   tkgrid(nextTraceButton, row = 3, column = 0, columnspan = 2)
# }

# # -----------------------------------------------------------------------------
# # Download function to save averaged responses as CSV (all ABFs in one sheet)
# # -----------------------------------------------------------------------------
# download_data <- function() {
#   if (!exists("abf_analysis_result", envir = .GlobalEnv)) {
#     tkmessageBox(message = "No analysis result available for download.")
#     return()
#   }
#   result <- get("abf_analysis_result", envir = .GlobalEnv)
#   dataList <- result$baseline_corrected_mean_data
#   if (length(dataList) == 0) {
#     tkmessageBox(message = "No averaged data available.")
#     return()
#   }
#   # Assume all vectors are the same length (by consistency check).
#   dt_val <- result$datasets[[1]]$samplingIntervalInSec * 1000
#   time_vec <- seq(0, by = dt_val, length.out = length(dataList[[1]]))
#   df <- data.frame(Time = time_vec)
#   for (fname in names(dataList)) {
#     df[[fname]] <- dataList[[fname]]
#   }
#   download_folder <- tclvalue(folderPathVar)
#   if(nchar(download_folder) == 0) {
#     tkmessageBox(message = "No folder selected for download.")
#     return()
#   }
#   file_path <- file.path(download_folder, "averaged_data.csv")
#   write.csv(df, file = file_path, row.names = FALSE)
#   tkmessageBox(message = paste("Averaged data saved to", file_path))
# }

# # -----------------------------------------------------------------------------
# # Main ABF Analysis UI
# # -----------------------------------------------------------------------------
# ABF_analysis_tk <- function() {
#   tt <- tktoplevel()
#   tkwm.title(tt, "ABF Analysis")
  
#   sidebarFrame <- tkframe(tt)
#   mainFrame <- tkframe(tt)
#   tkgrid(sidebarFrame, row = 0, column = 0, sticky = "ns")
#   tkgrid(mainFrame, row = 0, column = 1, sticky = "nsew")
#   tkgrid.rowconfigure(tt, 0, weight = 1)
#   tkgrid.columnconfigure(tt, 1, weight = 1)
  
#   # Save mainFrame as the global plot panel.
#   plotPanel <<- mainFrame
  
#   # --- Sidebar Controls ---
#   folderLabel <- tklabel(sidebarFrame, text = "Select ABF Folder:")
#   tkgrid(folderLabel, row = 0, column = 0, sticky = "w")
#   folderPathVar <<- tclVar("")
#   folderEntry <- tkentry(sidebarFrame, textvariable = folderPathVar, width = 30)
#   tkgrid(folderEntry, row = 0, column = 1, sticky = "w")
  
#   browseFolderButton <- tkbutton(sidebarFrame, text = "Browse", command = function() {
#     folderPath <- tclvalue(tkchooseDirectory())
#     if (nchar(folderPath) > 0) {
#       tclvalue(folderPathVar) <<- folderPath
#       abf_list <- list.files(path = folderPath, pattern = "\\.abf$", ignore.case = TRUE)
#       if (length(abf_list) == 0) {
#         tkmessageBox(message = "No ABF files found in the selected folder.")
#       } else {
#         tkdelete(abfListBox, 0, "end")
#         for (f in abf_list) { tkinsert(abfListBox, "end", f) }
#         firstFilePath <- file.path(folderPath, abf_list[1])
#         ds <- readABF(firstFilePath)
#         dummy_result <- list(metadata = list(extract_metadata(ds)))
#         updateAdditionalParams(dummy_result)
#       }
#     }
#   })
#   tkgrid(browseFolderButton, row = 0, column = 2, padx = 5)
  
#   abfListLabel <- tklabel(sidebarFrame, text = "ABF Files:")
#   tkgrid(abfListLabel, row = 1, column = 0, sticky = "w", pady = 5)
#   abfListBox <<- tklistbox(sidebarFrame, height = 5, width = 40, selectmode = "multiple")
#   tkgrid(abfListBox, row = 2, column = 0, columnspan = 3, pady = 5, sticky = "we")
  
#   # --- Additional Parameters (above Baseline) ---
#   paramFrame <- tkframe(sidebarFrame)
#   tkgrid(paramFrame, row = 3, column = 0, columnspan = 3, sticky = "w")
  
#   experimentVar <<- tclVar("Voltage Clamp")
#   unitVar <<- tclVar("")
#   dataColVar <<- tclVar("")
#   dtVar <<- tclVar("")
#   ntracesVar <<- tclVar("")
  
#   tkgrid(tklabel(paramFrame, text = "Experiment:"), row = 0, column = 0, sticky = "w")
#   experimentCombo <- ttkcombobox(paramFrame, textvariable = experimentVar, 
#                                  values = c("Voltage Clamp", "Current Clamp"), width = 15)
#   tkgrid(experimentCombo, row = 0, column = 1, sticky = "w")
  
#   tkgrid(tklabel(paramFrame, text = "Units:"), row = 1, column = 0, sticky = "w")
#   unitEntry <- tkentry(paramFrame, textvariable = unitVar, width = 10, state = "readonly")
#   tkgrid(unitEntry, row = 1, column = 1, sticky = "w")
  
#   tkgrid(tklabel(paramFrame, text = "Data Column:"), row = 2, column = 0, sticky = "w")
#   dataColEntry <- tkentry(paramFrame, textvariable = dataColVar, width = 10, state = "readonly")
#   tkgrid(dataColEntry, row = 2, column = 1, sticky = "w")
  
#   tkgrid(tklabel(paramFrame, text = "dt (ms):"), row = 3, column = 0, sticky = "w")
#   dtEntry <- tkentry(paramFrame, textvariable = dtVar, width = 10, state = "readonly")
#   tkgrid(dtEntry, row = 3, column = 1, sticky = "w")
  
#   tkgrid(tklabel(paramFrame, text = "Traces:"), row = 4, column = 0, sticky = "w")
#   ntracesEntry <- tkentry(paramFrame, textvariable = ntracesVar, width = 10, state = "readonly")
#   tkgrid(ntracesEntry, row = 4, column = 1, sticky = "w")
  
#   # --- Baseline and Other Parameter Inputs ---
#   baselineVar <<- tclVar("100")
#   stimTimeVar <<- tclVar("350")
#   xbarVar <<- tclVar("100")
#   ybarVar <<- tclVar("50")
  
#   tkgrid(tklabel(sidebarFrame, text = "Baseline:"), row = 4, column = 0, sticky = "w")
#   tkgrid(tkentry(sidebarFrame, textvariable = baselineVar, width = 10), row = 4, column = 1, sticky = "w")
#   tkgrid(tklabel(sidebarFrame, text = "Stimulation Time:"), row = 5, column = 0, sticky = "w")
#   tkgrid(tkentry(sidebarFrame, textvariable = stimTimeVar, width = 10), row = 5, column = 1, sticky = "w")
#   tkgrid(tklabel(sidebarFrame, text = "x-bar length:"), row = 6, column = 0, sticky = "w")
#   tkgrid(tkentry(sidebarFrame, textvariable = xbarVar, width = 10), row = 6, column = 1, sticky = "w")
#   tkgrid(tklabel(sidebarFrame, text = "y-bar length:"), row = 7, column = 0, sticky = "w")
#   tkgrid(tkentry(sidebarFrame, textvariable = ybarVar, width = 10), row = 7, column = 1, sticky = "w")
  
#   consoleText <<- tktext(sidebarFrame, width = 40, height = 4)
#   tkgrid(consoleText, row = 8, column = 0, columnspan = 3, pady = 5)
  
#   # --- Update Function ---
#   updateAdditionalParams <<- function(result) {
#     if (!is.null(result) && length(result$metadata) >= 1) {
#       meta1 <- result$metadata[[1]]
#       tclvalue(dtVar) <<- as.character(meta1$samplingIntervalInSec * 1000)
#       if (!is.null(meta1$header$lActualEpisodes)) {
#         tclvalue(ntracesVar) <<- as.character(meta1$header$lActualEpisodes)
#       } else {
#         tclvalue(ntracesVar) <<- "N/A"
#       }
#       expType <- tclvalue(experimentVar)
#       col_idx <- choose_data_column(meta1$channelUnits, expType)
#       if (!is.na(col_idx)) {
#         tclvalue(unitVar) <<- meta1$channelUnits[col_idx]
#         tclvalue(dataColVar) <<- as.character(col_idx)
#       } else {
#         tclvalue(unitVar) <<- "N/A"
#         tclvalue(dataColVar) <<- "N/A"
#       }
#     }
#   }
  
#   tkbind(experimentCombo, "<<ComboboxSelected>>", function() {
#     if (exists("abf_analysis_result", envir = .GlobalEnv)) {
#       result <- get("abf_analysis_result", envir = .GlobalEnv)
#       updateAdditionalParams(result)
#     }
#   })
  
#   # --- Analysis Action: Load Data button (loads data only, no plotting) ---
#   runAnalysis <<- function() {
#     folderPath <- tclvalue(folderPathVar)
#     if (nchar(folderPath) == 0) {
#       tkmessageBox(message = "Please select an ABF folder first.")
#       return()
#     }
#     selIndices <- as.integer(tkcurselection(abfListBox))
#     allFiles <- as.character(tkget(abfListBox, 0, "end"))
#     abf_files <- if (length(selIndices) == 0) allFiles else allFiles[selIndices + 1]
#     if (length(abf_files) == 0) {
#       tkmessageBox(message = "No ABF files selected.")
#       return()
#     }
#     result <- tryCatch({
#       load_abf_data(abf_files = abf_files, abf_path = folderPath)
#     }, error = function(e) {
#       tkmessageBox(message = paste("Error during data loading:", e$message))
#       return(NULL)
#     })
#     if (!is.null(result)) {
#       tkdelete(consoleText, "1.0", "end")
#       tkinsert(consoleText, "end", paste("Data loaded. Processed", length(abf_files), "file(s)."))
#       assign("abf_analysis_result", result, envir = .GlobalEnv)
#       updateAdditionalParams(result)
#       cons_msg <- check_consistency(result$metadata)
#       if (cons_msg == "Data is consistent") {
#         tkmessageBox(message = cons_msg)
#       } else {
#         tkmessageBox(message = paste("ERROR:", cons_msg))
#       }
#       tkconfigure(runAnalysisButton, text = "Load Data")
#     }
#   }
  
#   runAnalysisButton <<- tkbutton(sidebarFrame, text = "Load Data", command = runAnalysis)
#   tkgrid(runAnalysisButton, row = 9, column = 0, columnspan = 3, pady = 5)
  
#   # --- New Button: Review Recordings (placed above Average Approved Traces) ---
#   reviewButton <<- tkbutton(sidebarFrame, text = "Review Recordings", command = review_recordings)
#   tkgrid(reviewButton, row = 10, column = 0, columnspan = 3, pady = 5)
  
#   # --- New Button: Average Approved Traces ---
#   averageApprovedTraces <<- function() {
#     if (length(traces2average) == 0 || all(sapply(traces2average, length) == 0)) {
#       tkmessageBox(message = "No approved traces available. Please review recordings first.")
#       return()
#     }
#     folderPath <- tclvalue(folderPathVar)
#     if (nchar(folderPath) == 0) {
#       tkmessageBox(message = "Please select an ABF folder first.")
#       return()
#     }
#     selIndices <- as.integer(tkcurselection(abfListBox))
#     allFiles <- as.character(tkget(abfListBox, 0, "end"))
#     abf_files <- if (length(selIndices) == 0) allFiles else allFiles[selIndices + 1]
#     if (length(abf_files) == 0) {
#       tkmessageBox(message = "No ABF files selected.")
#       return()
#     }
#     baseline <- as.numeric(tclvalue(baselineVar))
#     stimTime <- as.numeric(tclvalue(stimTimeVar))
#     xbar <- as.numeric(tclvalue(xbarVar))
#     ybar <- as.numeric(tclvalue(ybarVar))
    
#     result <- tryCatch({
#       abf_out <<- abf_averages(datasets = abf_analysis_result$datasets,
#                               traces2average = traces2average,
#                               baseline = baseline, stimulation_time = stimTime,
#                               dataCol = as.numeric(tclvalue(dataColVar)),
#                               xlim = NULL, ylim = NULL,
#                               color = "darkgrey", xbar = xbar, ybar = ybar,
#                               width = 5.25, height = 2.75)
#       return(abf_out)
#     }, error = function(e) {
#       tkmessageBox(message = paste("Error during averaging of approved traces:", e$message))
#       return(NULL)
#     })
#     if (!is.null(result)) {
#       tkdelete(consoleText, "1.0", "end")
#       tkinsert(consoleText, "end", paste("Averaging on approved traces complete. Processed", length(abf_files), "file(s)."))
#       assign("abf_analysis_result", result, envir = .GlobalEnv)
      
#       # Define the drawing function for the averaged plot.
#       drawPlot <- function() {
#         if (exists("abf_analysis_result", envir = .GlobalEnv)) {
#           result <- get("abf_analysis_result", envir = .GlobalEnv)
#           datasets <- result$datasets
#           if (length(datasets) >= 1) {
#             iii <- 1
#             dt_val <- datasets[[iii]]$samplingIntervalInSec * 1000
#             time <- seq(0, by = dt_val, length.out = nrow(datasets[[iii]]$data[[1]]))
#             channel_index <- as.numeric(tclvalue(dataColVar))
#             if (is.na(channel_index) || channel_index < 1 || channel_index > length(datasets[[iii]]$data)) {
#               channel_index <- 1
#             }
#             trace <- datasets[[iii]]$data[[channel_index]][, 1]
#             # Draw the trace using egs_plot with scale bars only (no axes)
#             egs_plot(x = time, y = trace, show_bar = TRUE, show_text = TRUE, color = "darkgrey")
#           } else {
#             plot.new()
#             text(0.5, 0.5, "No data available")
#           }
#         } else {
#           plot.new()
#           text(0.5, 0.5, "No analysis result to display")
#         }
#       }
#       # Destroy any previous plotWidget and embed the averaged plot into the right-hand panel.
#       if (exists("plotWidget", envir = .GlobalEnv)) tkdestroy(plotWidget)
#       plotWidget <<- tkrplot(plotPanel, fun = drawPlot, hscale = 1, vscale = 1, width = 600, height = 400)
#       tkgrid(plotWidget, row = 0, column = 0, sticky = "nsew")
#     }
#   }
  
#   averageApprovedButton <<- tkbutton(sidebarFrame, text = "Average Approved Traces", command = averageApprovedTraces)
#   tkgrid(averageApprovedButton, row = 11, column = 0, columnspan = 3, pady = 5)
  
#   # --- New Button: Download Data (CSV) ---
#   sapply(1:length(abf_out$baseline_corrected_mean_data), function(ii) abf_out$baseline_corrected_mean_data[[ii]])

#   downloadButton <<- tkbutton(sidebarFrame, text = "Download Data", command = download_data)
#   tkgrid(downloadButton, row = 12, column = 0, columnspan = 3, pady = 5)
  
#   tkfocus(tt)
# }

# # -----------------------------------------------------------------------------
# # Launch the ABF Analysis interface:
# # -----------------------------------------------------------------------------
# ABF_analysis_tk()


# # Remove all objects from the environment
# rm(list = ls(all = TRUE))

# # -----------------------------------------------------------------------------
# # Package loading function and required packages:
# # -----------------------------------------------------------------------------
# load_required_packages <- function(packages) {
#   new.packages <- packages[!(packages %in% installed.packages()[, "Package"])]
#   if (length(new.packages)) install.packages(new.packages)
#   invisible(lapply(packages, library, character.only = TRUE))
# }
# required.packages <- c("robustbase", "minpack.lm", "Rcpp", "signal", "readABF", "tcltk", "tkrplot")
# load_required_packages(required.packages)

# # -----------------------------------------------------------------------------
# # Helper Functions
# # -----------------------------------------------------------------------------

# # Extract key metadata from an ABF dataset.
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

# # Choose data column index based on experiment type.
# choose_data_column <- function(channelUnits, experiment) {
#   if (experiment == "Voltage Clamp") {
#     idx <- grep("A", channelUnits, ignore.case = TRUE)
#   } else if (experiment == "Current Clamp") {
#     idx <- grep("V", channelUnits, ignore.case = TRUE)
#   } else {
#     idx <- integer(0)
#   }
#   if (length(idx) > 0) return(idx[1])
#   else return(NA)
# }

# # Consistency check: Verify that dt, chosen unit, and number of traces are identical.
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
#     return("Data is consistent")
#   } else {
#     error_msgs <- c()
#     if (!dt_good) error_msgs <- c(error_msgs, paste("Inconsistent dt values:", paste(dt_values, collapse = ", ")))
#     if (!unit_good) error_msgs <- c(error_msgs, paste("Inconsistent Units:", paste(unit_values, collapse = ", ")))
#     if (!traces_good) error_msgs <- c(error_msgs, paste("Inconsistent Traces:", paste(traces_values, collapse = ", ")))
#     return(paste(error_msgs, collapse = "; "))
#   }
# }

# # (The egs_plot function is kept here for legacy use.)
# egs_plot <- function(x, y, sign = -1, xlim = NULL, ylim = NULL, lwd = 1, 
#                      show_text = FALSE, height = 4, width = 2.5, 
#                      xbar = 100, ybar = 50, color = "#4C77BB", show_bar = FALSE) {
#   if (is.null(ylim))
#     ylim <- if (sign == 1) c(0, max(y)) else c(-max(-y), 0)
#   if (is.null(xlim))
#     xlim <- c(min(x), max(x))
  
#   idx1 <- which.min(abs(x - xlim[1]))
#   idx2 <- which.min(abs(x - xlim[2]))
  
#   plot(x[idx1:idx2], y[idx1:idx2], type = "l", col = color,
#        xlim = xlim, ylim = ylim, bty = "n", lwd = lwd, lty = 1,
#        axes = FALSE, frame = FALSE, xlab = "", ylab = "")
  
#   if (show_bar) {
#     ybar_start <- min(ylim) + (max(ylim) - min(ylim)) / 20
#     x_start <- max(xlim) - xbar - 50
#     y_start <- ybar_start
#     x_end <- x_start + xbar
#     y_end <- y_start + ybar
#     segments(x_start, y_start, x_end, y_start, lwd = lwd, col = "black")
#     segments(x_start, y_start, x_start, y_end, lwd = lwd, col = "black")
#     if (show_text) {
#       text(x = (x_start + x_end) / 2, y = y_start - ybar / 20, 
#            labels = paste(xbar, "ms"), adj = c(0.5, 1))
#       text(x = x_start - xbar / 4, y = (y_start + y_end) / 2,  
#            labels = paste(ybar, "pA"), adj = c(0.5, 0.5), srt = 90)
#     }
#   }
# }

# # New function to load ABF data.
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
# # Revised abf_averages Function
# # -----------------------------------------------------------------------------
# abf_averages <- function(datasets, 
#                          baseline = 100, 
#                          stimulation_time = 350, 
#                          traces2average = NULL,
#                          dataCol = 1, 
#                          ylim = NULL, 
#                          xlim = NULL, 
#                          color = "#4C77BB", 
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
#     egs_plot(x = time[[ii]], y = responses0_mean[[ii]], ylim = ylim,
#              xlim = xlim, color = color, width = width, height = height,
#              show_text = FALSE, show_bar = show_bar[ii])
#     axis(1)
#     axis(2)
#   }
  
#   if (save) {
#     warning("save_graph not implemented in this UI example")
#   }
  
#   return(list(raw_data = responses,
#               baseline_corrected_data = responses0,
#               baseline_corrected_mean_data = responses0_mean,
#               datasets = datasets))
# }

# # -----------------------------------------------------------------------------
# # Global variable for approved trace indices (one per dataset)
# # -----------------------------------------------------------------------------
# traces2average <<- list()

# # -----------------------------------------------------------------------------
# # Revised review_recordings Function (embedded into the right-hand plot panel)
# # -----------------------------------------------------------------------------
# review_recordings <- function() {
#   if (!exists("abf_analysis_result", envir = .GlobalEnv)) {
#     tkmessageBox(message = "No analysis result available for review.")
#     return()
#   }
#   result <- get("abf_analysis_result", envir = .GlobalEnv)
#   datasets <- result$datasets
#   traces2average <<- vector("list", length = length(datasets))
#   for (i in seq_along(datasets)) { 
#     traces2average[[i]] <<- integer(0) 
#   }
  
#   current_dataset <<- 1
#   current_trace <<- 1
  
#   # Clear existing content from the plot panel.
#   children <- as.character(tkwinfo("children", plotPanel))
#   if (length(children) > 0) {
#     for (child in children) tkdestroy(child)
#   }
  
#   # Create a review frame inside the plot panel.
#   reviewFrame <<- tkframe(plotPanel)
#   tkgrid(reviewFrame, row = 0, column = 0, sticky = "nsew")
  
#   # Use the current file name as the header.
#   current_filename <- names(datasets)[current_dataset]
#   infoLabel <<- tklabel(reviewFrame, text = paste(current_filename, "trace", current_trace))
#   tkgrid(infoLabel, row = 0, column = 0, columnspan = 2)
  
#   # Create a review plot using standard plot() with axes suppressed, then add axes manually.
#   reviewPlot <<- tkrplot(reviewFrame, fun = function() {
#     ds <- datasets[[current_dataset]]
#     if (current_trace > length(ds$data)) {
#       plot.new()
#       text(0.5, 0.5, paste("No more recordings in", current_filename))
#     } else {
#       trace_matrix <- ds$data[[current_trace]]
#       data_column <- as.numeric(tclvalue(dataColVar))
#       if (is.na(data_column) || data_column < 1 || data_column > ncol(trace_matrix)) {
#         data_column <- 1
#       }
#       dt_val <- ds$samplingIntervalInSec * 1000
#       time <- seq(0, by = dt_val, length.out = nrow(trace_matrix))
#       trace <- trace_matrix[, data_column]
#       # Suppress default axes and then add only horizontal tick labels.
#       plot(time, trace, col = "darkgrey", xlab = "Time (ms)", ylab = tclvalue(unitVar),
#            type = "l", bty = "l", axes = FALSE, 
#            main = paste(current_filename, "trace", current_trace))
#       axis(1)
#       axis(2, las = 1)
#     }
#   }, hscale = 1, vscale = 1)
#   tkgrid(reviewPlot, row = 1, column = 0, columnspan = 2)
  
#   acceptButton <- tkbutton(reviewFrame, text = "Accept", command = function() {
#     traces2average[[current_dataset]] <<- c(traces2average[[current_dataset]], current_trace)
#     tkconfigure(acceptButton, state = "disabled", relief = "sunken")
#     tkconfigure(rejectButton, state = "normal", relief = "raised")
#   })
#   tkgrid(acceptButton, row = 2, column = 0)
  
#   rejectButton <- tkbutton(reviewFrame, text = "Reject", command = function() {
#     tkconfigure(rejectButton, state = "disabled", relief = "sunken")
#     tkconfigure(acceptButton, state = "normal", relief = "raised")
#   })
#   tkgrid(rejectButton, row = 2, column = 1)
  
#   nextTraceButton <- tkbutton(reviewFrame, text = "Next Recording", command = function() {
#     tkconfigure(acceptButton, state = "normal", relief = "raised")
#     tkconfigure(rejectButton, state = "normal", relief = "raised")
#     ds <- datasets[[current_dataset]]
#     numRecordings <- length(ds$data)
#     if (current_trace < numRecordings) {
#       current_trace <<- current_trace + 1
#     } else {
#       tkmessageBox(message = paste("Finished reviewing", current_filename))
#       if (current_dataset < length(datasets)) {
#         current_dataset <<- current_dataset + 1
#         current_trace <<- 1
#       } else {
#         tkmessageBox(message = "Review complete. Approved recordings are in 'traces2average'.")
#         return()
#       }
#     }
#     current_filename <- names(datasets)[current_dataset]
#     tkconfigure(infoLabel, text = paste(current_filename, "trace", current_trace))
#     tkrreplot(reviewPlot)
#   })
#   tkgrid(nextTraceButton, row = 3, column = 0, columnspan = 2)
# }

# # -----------------------------------------------------------------------------
# # Main ABF Analysis UI
# # -----------------------------------------------------------------------------
# ABF_analysis_tk <- function() {
#   tt <- tktoplevel()
#   tkwm.title(tt, "ABF Analysis")
  
#   sidebarFrame <- tkframe(tt)
#   mainFrame <- tkframe(tt)
#   tkgrid(sidebarFrame, row = 0, column = 0, sticky = "ns")
#   tkgrid(mainFrame, row = 0, column = 1, sticky = "nsew")
#   tkgrid.rowconfigure(tt, 0, weight = 1)
#   tkgrid.columnconfigure(tt, 1, weight = 1)
  
#   # Save mainFrame as the global plot panel.
#   plotPanel <<- mainFrame
  
#   # --- Sidebar Controls ---
#   folderLabel <- tklabel(sidebarFrame, text = "Select ABF Folder:")
#   tkgrid(folderLabel, row = 0, column = 0, sticky = "w")
#   folderPathVar <<- tclVar("")
#   folderEntry <- tkentry(sidebarFrame, textvariable = folderPathVar, width = 30)
#   tkgrid(folderEntry, row = 0, column = 1, sticky = "w")
  
#   browseFolderButton <- tkbutton(sidebarFrame, text = "Browse", command = function() {
#     folderPath <- tclvalue(tkchooseDirectory())
#     if (nchar(folderPath) > 0) {
#       tclvalue(folderPathVar) <<- folderPath
#       abf_list <- list.files(path = folderPath, pattern = "\\.abf$", ignore.case = TRUE)
#       if (length(abf_list) == 0) {
#         tkmessageBox(message = "No ABF files found in the selected folder.")
#       } else {
#         tkdelete(abfListBox, 0, "end")
#         for (f in abf_list) { tkinsert(abfListBox, "end", f) }
#         firstFilePath <- file.path(folderPath, abf_list[1])
#         ds <- readABF(firstFilePath)
#         dummy_result <- list(metadata = list(extract_metadata(ds)))
#         updateAdditionalParams(dummy_result)
#       }
#     }
#   })
#   tkgrid(browseFolderButton, row = 0, column = 2, padx = 5)
  
#   abfListLabel <- tklabel(sidebarFrame, text = "ABF Files:")
#   tkgrid(abfListLabel, row = 1, column = 0, sticky = "w", pady = 5)
#   abfListBox <<- tklistbox(sidebarFrame, height = 5, width = 40, selectmode = "multiple")
#   tkgrid(abfListBox, row = 2, column = 0, columnspan = 3, pady = 5, sticky = "we")
  
#   # --- Additional Parameters (above Baseline) ---
#   paramFrame <- tkframe(sidebarFrame)
#   tkgrid(paramFrame, row = 3, column = 0, columnspan = 3, sticky = "w")
  
#   experimentVar <<- tclVar("Voltage Clamp")
#   unitVar <<- tclVar("")
#   dataColVar <<- tclVar("")
#   dtVar <<- tclVar("")
#   ntracesVar <<- tclVar("")
  
#   tkgrid(tklabel(paramFrame, text = "Experiment:"), row = 0, column = 0, sticky = "w")
#   experimentCombo <- ttkcombobox(paramFrame, textvariable = experimentVar, 
#                                  values = c("Voltage Clamp", "Current Clamp"), width = 15)
#   tkgrid(experimentCombo, row = 0, column = 1, sticky = "w")
  
#   tkgrid(tklabel(paramFrame, text = "Units:"), row = 1, column = 0, sticky = "w")
#   unitEntry <- tkentry(paramFrame, textvariable = unitVar, width = 10, state = "readonly")
#   tkgrid(unitEntry, row = 1, column = 1, sticky = "w")
  
#   tkgrid(tklabel(paramFrame, text = "Data Column:"), row = 2, column = 0, sticky = "w")
#   dataColEntry <- tkentry(paramFrame, textvariable = dataColVar, width = 10, state = "readonly")
#   tkgrid(dataColEntry, row = 2, column = 1, sticky = "w")
  
#   tkgrid(tklabel(paramFrame, text = "dt (ms):"), row = 3, column = 0, sticky = "w")
#   dtEntry <- tkentry(paramFrame, textvariable = dtVar, width = 10, state = "readonly")
#   tkgrid(dtEntry, row = 3, column = 1, sticky = "w")
  
#   tkgrid(tklabel(paramFrame, text = "Traces:"), row = 4, column = 0, sticky = "w")
#   ntracesEntry <- tkentry(paramFrame, textvariable = ntracesVar, width = 10, state = "readonly")
#   tkgrid(ntracesEntry, row = 4, column = 1, sticky = "w")
  
#   # --- Baseline and Other Parameter Inputs ---
#   baselineVar <<- tclVar("100")
#   stimTimeVar <<- tclVar("350")
#   xbarVar <<- tclVar("100")
#   ybarVar <<- tclVar("50")
  
#   tkgrid(tklabel(sidebarFrame, text = "Baseline:"), row = 4, column = 0, sticky = "w")
#   tkgrid(tkentry(sidebarFrame, textvariable = baselineVar, width = 10), row = 4, column = 1, sticky = "w")
#   tkgrid(tklabel(sidebarFrame, text = "Stimulation Time:"), row = 5, column = 0, sticky = "w")
#   tkgrid(tkentry(sidebarFrame, textvariable = stimTimeVar, width = 10), row = 5, column = 1, sticky = "w")
#   tkgrid(tklabel(sidebarFrame, text = "x-bar length:"), row = 6, column = 0, sticky = "w")
#   tkgrid(tkentry(sidebarFrame, textvariable = xbarVar, width = 10), row = 6, column = 1, sticky = "w")
#   tkgrid(tklabel(sidebarFrame, text = "y-bar length:"), row = 7, column = 0, sticky = "w")
#   tkgrid(tkentry(sidebarFrame, textvariable = ybarVar, width = 10), row = 7, column = 1, sticky = "w")
  
#   consoleText <<- tktext(sidebarFrame, width = 40, height = 4)
#   tkgrid(consoleText, row = 8, column = 0, columnspan = 3, pady = 5)
  
#   # --- Update Function ---
#   updateAdditionalParams <<- function(result) {
#     if (!is.null(result) && length(result$metadata) >= 1) {
#       meta1 <- result$metadata[[1]]
#       tclvalue(dtVar) <<- as.character(meta1$samplingIntervalInSec * 1000)
#       if (!is.null(meta1$header$lActualEpisodes)) {
#         tclvalue(ntracesVar) <<- as.character(meta1$header$lActualEpisodes)
#       } else {
#         tclvalue(ntracesVar) <<- "N/A"
#       }
#       expType <- tclvalue(experimentVar)
#       col_idx <- choose_data_column(meta1$channelUnits, expType)
#       if (!is.na(col_idx)) {
#         tclvalue(unitVar) <<- meta1$channelUnits[col_idx]
#         tclvalue(dataColVar) <<- as.character(col_idx)
#       } else {
#         tclvalue(unitVar) <<- "N/A"
#         tclvalue(dataColVar) <<- "N/A"
#       }
#     }
#   }
  
#   tkbind(experimentCombo, "<<ComboboxSelected>>", function() {
#     if (exists("abf_analysis_result", envir = .GlobalEnv)) {
#       result <- get("abf_analysis_result", envir = .GlobalEnv)
#       updateAdditionalParams(result)
#     }
#   })
  
#   # --- Analysis Action: Load Data button (loads data only, no plotting) ---
#   runAnalysis <<- function() {
#     folderPath <- tclvalue(folderPathVar)
#     if (nchar(folderPath) == 0) {
#       tkmessageBox(message = "Please select an ABF folder first.")
#       return()
#     }
#     selIndices <- as.integer(tkcurselection(abfListBox))
#     allFiles <- as.character(tkget(abfListBox, 0, "end"))
#     abf_files <- if (length(selIndices) == 0) allFiles else allFiles[selIndices + 1]
#     if (length(abf_files) == 0) {
#       tkmessageBox(message = "No ABF files selected.")
#       return()
#     }
#     result <- tryCatch({
#       load_abf_data(abf_files = abf_files, abf_path = folderPath)
#     }, error = function(e) {
#       tkmessageBox(message = paste("Error during data loading:", e$message))
#       return(NULL)
#     })
#     if (!is.null(result)) {
#       tkdelete(consoleText, "1.0", "end")
#       tkinsert(consoleText, "end", paste("Data loaded. Processed", length(abf_files), "file(s)."))
#       assign("abf_analysis_result", result, envir = .GlobalEnv)
#       updateAdditionalParams(result)
#       cons_msg <- check_consistency(result$metadata)
#       if (cons_msg == "Data is consistent") {
#         tkmessageBox(message = cons_msg)
#       } else {
#         tkmessageBox(message = paste("ERROR:", cons_msg))
#       }
#       tkconfigure(runAnalysisButton, text = "Load Data")
#     }
#   }
  
#   runAnalysisButton <<- tkbutton(sidebarFrame, text = "Load Data", command = runAnalysis)
#   tkgrid(runAnalysisButton, row = 9, column = 0, columnspan = 3, pady = 5)
  
#   # --- New Button: Review Recordings (placed above Average Approved Traces) ---
#   reviewButton <<- tkbutton(sidebarFrame, text = "Review Recordings", command = review_recordings)
#   tkgrid(reviewButton, row = 10, column = 0, columnspan = 3, pady = 5)
  
#   # --- New Button: Average Approved Traces ---
#   averageApprovedTraces <<- function() {
#     if (length(traces2average) == 0 || all(sapply(traces2average, length) == 0)) {
#       tkmessageBox(message = "No approved traces available. Please review recordings first.")
#       return()
#     }
#     folderPath <- tclvalue(folderPathVar)
#     if (nchar(folderPath) == 0) {
#       tkmessageBox(message = "Please select an ABF folder first.")
#       return()
#     }
#     selIndices <- as.integer(tkcurselection(abfListBox))
#     allFiles <- as.character(tkget(abfListBox, 0, "end"))
#     abf_files <- if (length(selIndices) == 0) allFiles else allFiles[selIndices + 1]
#     if (length(abf_files) == 0) {
#       tkmessageBox(message = "No ABF files selected.")
#       return()
#     }
#     baseline <- as.numeric(tclvalue(baselineVar))
#     stimTime <- as.numeric(tclvalue(stimTimeVar))
#     xbar <- as.numeric(tclvalue(xbarVar))
#     ybar <- as.numeric(tclvalue(ybarVar))
    
#     result <- tryCatch({
#       abf_out <- abf_averages(datasets = abf_analysis_result$datasets,
#                               traces2average = traces2average,
#                               baseline = baseline, stimulation_time = stimTime,
#                               dataCol = as.numeric(tclvalue(dataColVar)),
#                               xlim = NULL, ylim = NULL,
#                               color = "darkgrey", xbar = xbar, ybar = ybar,
#                               width = 5.25, height = 2.75)
#       return(abf_out)
#     }, error = function(e) {
#       tkmessageBox(message = paste("Error during averaging of approved traces:", e$message))
#       return(NULL)
#     })
#     if (!is.null(result)) {
#       tkdelete(consoleText, "1.0", "end")
#       tkinsert(consoleText, "end", paste("Averaging on approved traces complete. Processed", length(abf_files), "file(s)."))
#       assign("abf_analysis_result", result, envir = .GlobalEnv)
#       if (exists("plotWidget", envir = .GlobalEnv)) tkdestroy(plotWidget)
  
#       drawPlot <- function() {
#         if (exists("abf_analysis_result", envir = .GlobalEnv)) {
#           result <- get("abf_analysis_result", envir = .GlobalEnv)
#           datasets <- result$datasets
#           if (length(datasets) >= 1) {
#             iii <- 1
#             dt_val <- datasets[[iii]]$samplingIntervalInSec * 1000
#             time <- seq(0, by = dt_val, length.out = nrow(datasets[[iii]]$data[[1]]))
#             channel_index <- as.numeric(tclvalue(dataColVar))
#             if (is.na(channel_index) || channel_index < 1 || channel_index > length(datasets[[iii]]$data)) {
#               channel_index <- 1
#             }
#             trace <- datasets[[iii]]$data[[channel_index]][, 1]
            
#             # Draw the trace using egs_plot, which draws no axes,
#             # but with scale bars and labels (via show_bar & show_text)
#             egs_plot(x = time, y = trace, show_bar = TRUE, show_text = TRUE, color = "darkgrey")
#             # Do not add axes or title.
#           } else {
#             plot.new()
#             text(0.5, 0.5, "No data available")
#           }
#         } else {
#           plot.new()
#           text(0.5, 0.5, "No analysis result to display")
#         }
#       }

#       # Then embed the averaged plot into the right-hand panel by specifying the parent,
#       # and providing explicit width and height so that it appears within the UI rather than in a separate Quartz window:
#       plotWidget <<- tkrplot(plotPanel, fun = drawPlot, hscale = 1, vscale = 1, width = 600, height = 400)
#       tkgrid(plotWidget, row = 0, column = 0, sticky = "nsew")

#     }
#   }
  
#   averageApprovedButton <<- tkbutton(sidebarFrame, text = "Average Approved Traces", command = averageApprovedTraces)
#   tkgrid(averageApprovedButton, row = 11, column = 0, columnspan = 3, pady = 5)
  
#   tkfocus(tt)
# }

# # -----------------------------------------------------------------------------
# # Launch the ABF Analysis interface:
# # -----------------------------------------------------------------------------
# ABF_analysis_tk()





# # Remove all objects from the environment
# rm(list = ls(all = TRUE))

# # -----------------------------------------------------------------------------
# # Package loading function and required packages:
# # -----------------------------------------------------------------------------
# load_required_packages <- function(packages) {
#   new.packages <- packages[!(packages %in% installed.packages()[, "Package"])]
#   if (length(new.packages)) install.packages(new.packages)
#   invisible(lapply(packages, library, character.only = TRUE))
# }
# required.packages <- c("robustbase", "minpack.lm", "Rcpp", "signal", "readABF", "tcltk", "tkrplot")
# load_required_packages(required.packages)

# # -----------------------------------------------------------------------------
# # Helper Functions
# # -----------------------------------------------------------------------------

# # Extract key metadata from an ABF dataset.
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

# # Choose data column index based on experiment type.
# choose_data_column <- function(channelUnits, experiment) {
#   if (experiment == "Voltage Clamp") {
#     idx <- grep("A", channelUnits, ignore.case = TRUE)
#   } else if (experiment == "Current Clamp") {
#     idx <- grep("V", channelUnits, ignore.case = TRUE)
#   } else {
#     idx <- integer(0)
#   }
#   if (length(idx) > 0) return(idx[1])
#   else return(NA)
# }

# # Consistency check: Verify dt, chosen unit, and number of traces are identical.
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
#     return("Data is consistent")
#   } else {
#     error_msgs <- c()
#     if (!dt_good) error_msgs <- c(error_msgs, paste("Inconsistent dt values:", paste(dt_values, collapse = ", ")))
#     if (!unit_good) error_msgs <- c(error_msgs, paste("Inconsistent Units:", paste(unit_values, collapse = ", ")))
#     if (!traces_good) error_msgs <- c(error_msgs, paste("Inconsistent Traces:", paste(traces_values, collapse = ", ")))
#     return(paste(error_msgs, collapse = "; "))
#   }
# }

# # Custom plotting function (egs_plot)
# egs_plot <- function(x, y, sign = -1, xlim = NULL, ylim = NULL, lwd = 1, 
#                      show_text = FALSE, height = 4, width = 2.5, 
#                      xbar = 100, ybar = 50, color = "#4C77BB", show_bar = FALSE) {
#   if (is.null(ylim))
#     ylim <- if (sign == 1) c(0, max(y)) else c(-max(-y), 0)
#   if (is.null(xlim))
#     xlim <- c(min(x), max(x))
  
#   idx1 <- which.min(abs(x - xlim[1]))
#   idx2 <- which.min(abs(x - xlim[2]))
  
#   plot(x[idx1:idx2], y[idx1:idx2], type = "l", col = color,
#        xlim = xlim, ylim = ylim, bty = "n", lwd = lwd, lty = 1,
#        axes = FALSE, frame = FALSE, xlab = "", ylab = "")
  
#   if (show_bar) {
#     ybar_start <- min(ylim) + (max(ylim) - min(ylim)) / 20
#     x_start <- max(xlim) - xbar - 50
#     y_start <- ybar_start
#     x_end <- x_start + xbar
#     y_end <- y_start + ybar
#     segments(x_start, y_start, x_end, y_start, lwd = lwd, col = "black")
#     segments(x_start, y_start, x_start, y_end, lwd = lwd, col = "black")
#     if (show_text) {
#       text(x = (x_start + x_end) / 2, y = y_start - ybar / 20, 
#            labels = paste(xbar, "ms"), adj = c(0.5, 1))
#       text(x = x_start - xbar / 4, y = (y_start + y_end) / 2,  
#            labels = paste(ybar, "pA"), adj = c(0.5, 0.5), srt = 90)
#     }
#   }
# }

# # New function to load ABF data.
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
# # Revised abf_averages Function
# # -----------------------------------------------------------------------------
# abf_averages <- function(datasets, 
#                          baseline = 100, 
#                          stimulation_time = 350, 
#                          traces2average = NULL,
#                          dataCol = 1, 
#                          ylim = NULL, 
#                          xlim = NULL, 
#                          color = "#4C77BB", 
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
#     egs_plot(x = time[[ii]], y = responses0_mean[[ii]], ylim = ylim,
#              xlim = xlim, color = color, width = width, height = height,
#              show_text = FALSE, show_bar = show_bar[ii])
#     axis(1)
#     axis(2)
#   }
  
#   if (save) {
#     warning("save_graph not implemented in this UI example")
#   }
  
#   return(list(raw_data = responses,
#               baseline_corrected_data = responses0,
#               baseline_corrected_mean_data = responses0_mean,
#               datasets = datasets))
# }

# # -----------------------------------------------------------------------------
# # Global variable for approved trace indices (one per dataset)
# # -----------------------------------------------------------------------------
# traces2average <<- list()

# # -----------------------------------------------------------------------------
# # Revised review_recordings Function (embedded into the right-hand plot panel)
# # -----------------------------------------------------------------------------
# review_recordings <- function() {
#   if (!exists("abf_analysis_result", envir = .GlobalEnv)) {
#     tkmessageBox(message = "No analysis result available for review.")
#     return()
#   }
#   result <- get("abf_analysis_result", envir = .GlobalEnv)
#   datasets <- result$datasets
#   traces2average <<- vector("list", length = length(datasets))
#   for (i in seq_along(datasets)) { 
#     traces2average[[i]] <<- integer(0) 
#   }
  
#   current_dataset <<- 1
#   current_trace <<- 1
  
#   # Clear existing content from the plot panel.
#   children <- as.character(tkwinfo("children", plotPanel))
#   if (length(children) > 0) {
#     for (child in children) tkdestroy(child)
#   }
  
#   # Create a review frame inside the plot panel.
#   reviewFrame <<- tkframe(plotPanel)
#   tkgrid(reviewFrame, row = 0, column = 0, sticky = "nsew")
  
#   # Use the file name from the datasets list as the header.
#   current_filename <- names(datasets)[current_dataset]
#   infoLabel <<- tklabel(reviewFrame, text = paste(current_filename, "trace", current_trace))
#   tkgrid(infoLabel, row = 0, column = 0, columnspan = 2)
  
#   reviewPlot <<- tkrplot(reviewFrame, fun = function() {
#     ds <- datasets[[current_dataset]]
#     if (current_trace > length(ds$data)) {
#       plot.new()
#       text(0.5, 0.5, paste("No more recordings in", current_filename))
#     } else {
#       trace_matrix <- ds$data[[current_trace]]
#       data_column <- as.numeric(tclvalue(dataColVar))
#       if (is.na(data_column) || data_column < 1 || data_column > ncol(trace_matrix)) {
#         data_column <- 1
#       }
#       dt_val <- ds$samplingIntervalInSec * 1000
#       time <- seq(0, by = dt_val, length.out = nrow(trace_matrix))
#       trace <- trace_matrix[, data_column]
#       # Draw the trace with egs_plot (which draws scale bars and no axes)
#       egs_plot(x = time, y = trace, show_bar = TRUE, show_text = TRUE)
#       # Add a title indicating filename and trace number.
#       title(main = paste(current_filename, "trace", current_trace), line = 0.5)
#     }
#   }, hscale = 1, vscale = 1)
#   tkgrid(reviewPlot, row = 1, column = 0, columnspan = 2)
  
#   acceptButton <- tkbutton(reviewFrame, text = "Accept", command = function() {
#     traces2average[[current_dataset]] <<- c(traces2average[[current_dataset]], current_trace)
#     tkconfigure(acceptButton, state = "disabled", relief = "sunken")
#     tkconfigure(rejectButton, state = "normal", relief = "raised")
#   })
#   tkgrid(acceptButton, row = 2, column = 0)
  
#   rejectButton <- tkbutton(reviewFrame, text = "Reject", command = function() {
#     tkconfigure(rejectButton, state = "disabled", relief = "sunken")
#     tkconfigure(acceptButton, state = "normal", relief = "raised")
#   })
#   tkgrid(rejectButton, row = 2, column = 1)
  
#   nextTraceButton <- tkbutton(reviewFrame, text = "Next Recording", command = function() {
#     tkconfigure(acceptButton, state = "normal", relief = "raised")
#     tkconfigure(rejectButton, state = "normal", relief = "raised")
#     ds <- datasets[[current_dataset]]
#     numRecordings <- length(ds$data)
#     if (current_trace < numRecordings) {
#       current_trace <<- current_trace + 1
#     } else {
#       tkmessageBox(message = paste("Finished reviewing", current_filename))
#       if (current_dataset < length(datasets)) {
#         current_dataset <<- current_dataset + 1
#         current_trace <<- 1
#       } else {
#         tkmessageBox(message = "Review complete. Approved recordings are in 'traces2average'.")
#         return()
#       }
#     }
#     current_filename <- names(datasets)[current_dataset]
#     tkconfigure(infoLabel, text = paste(current_filename, "trace", current_trace))
#     tkrreplot(reviewPlot)
#   })
#   tkgrid(nextTraceButton, row = 3, column = 0, columnspan = 2)
# }

# # -----------------------------------------------------------------------------
# # Main ABF Analysis UI
# # -----------------------------------------------------------------------------
# ABF_analysis_tk <- function() {
#   tt <- tktoplevel()
#   tkwm.title(tt, "ABF Analysis")
  
#   sidebarFrame <- tkframe(tt)
#   mainFrame <- tkframe(tt)
#   tkgrid(sidebarFrame, row = 0, column = 0, sticky = "ns")
#   tkgrid(mainFrame, row = 0, column = 1, sticky = "nsew")
#   tkgrid.rowconfigure(tt, 0, weight = 1)
#   tkgrid.columnconfigure(tt, 1, weight = 1)
  
#   # Save mainFrame as the global plot panel.
#   plotPanel <<- mainFrame
  
#   # --- Sidebar Controls ---
#   folderLabel <- tklabel(sidebarFrame, text = "Select ABF Folder:")
#   tkgrid(folderLabel, row = 0, column = 0, sticky = "w")
#   folderPathVar <<- tclVar("")
#   folderEntry <- tkentry(sidebarFrame, textvariable = folderPathVar, width = 30)
#   tkgrid(folderEntry, row = 0, column = 1, sticky = "w")
  
#   browseFolderButton <- tkbutton(sidebarFrame, text = "Browse", command = function() {
#     folderPath <- tclvalue(tkchooseDirectory())
#     if (nchar(folderPath) > 0) {
#       tclvalue(folderPathVar) <<- folderPath
#       abf_list <- list.files(path = folderPath, pattern = "\\.abf$", ignore.case = TRUE)
#       if (length(abf_list) == 0) {
#         tkmessageBox(message = "No ABF files found in the selected folder.")
#       } else {
#         tkdelete(abfListBox, 0, "end")
#         for (f in abf_list) { tkinsert(abfListBox, "end", f) }
#         firstFilePath <- file.path(folderPath, abf_list[1])
#         ds <- readABF(firstFilePath)
#         dummy_result <- list(metadata = list(extract_metadata(ds)))
#         updateAdditionalParams(dummy_result)
#       }
#     }
#   })
#   tkgrid(browseFolderButton, row = 0, column = 2, padx = 5)
  
#   abfListLabel <- tklabel(sidebarFrame, text = "ABF Files:")
#   tkgrid(abfListLabel, row = 1, column = 0, sticky = "w", pady = 5)
#   abfListBox <<- tklistbox(sidebarFrame, height = 5, width = 40, selectmode = "multiple")
#   tkgrid(abfListBox, row = 2, column = 0, columnspan = 3, pady = 5, sticky = "we")
  
#   # --- Additional Parameters (above Baseline) ---
#   paramFrame <- tkframe(sidebarFrame)
#   tkgrid(paramFrame, row = 3, column = 0, columnspan = 3, sticky = "w")
  
#   experimentVar <<- tclVar("Voltage Clamp")
#   unitVar <<- tclVar("")
#   dataColVar <<- tclVar("")
#   dtVar <<- tclVar("")
#   ntracesVar <<- tclVar("")
  
#   tkgrid(tklabel(paramFrame, text = "Experiment:"), row = 0, column = 0, sticky = "w")
#   experimentCombo <- ttkcombobox(paramFrame, textvariable = experimentVar, 
#                                  values = c("Voltage Clamp", "Current Clamp"), width = 15)
#   tkgrid(experimentCombo, row = 0, column = 1, sticky = "w")
  
#   tkgrid(tklabel(paramFrame, text = "Units:"), row = 1, column = 0, sticky = "w")
#   unitEntry <- tkentry(paramFrame, textvariable = unitVar, width = 10, state = "readonly")
#   tkgrid(unitEntry, row = 1, column = 1, sticky = "w")
  
#   tkgrid(tklabel(paramFrame, text = "Data Column:"), row = 2, column = 0, sticky = "w")
#   dataColEntry <- tkentry(paramFrame, textvariable = dataColVar, width = 10, state = "readonly")
#   tkgrid(dataColEntry, row = 2, column = 1, sticky = "w")
  
#   tkgrid(tklabel(paramFrame, text = "dt (ms):"), row = 3, column = 0, sticky = "w")
#   dtEntry <- tkentry(paramFrame, textvariable = dtVar, width = 10, state = "readonly")
#   tkgrid(dtEntry, row = 3, column = 1, sticky = "w")
  
#   tkgrid(tklabel(paramFrame, text = "Traces:"), row = 4, column = 0, sticky = "w")
#   ntracesEntry <- tkentry(paramFrame, textvariable = ntracesVar, width = 10, state = "readonly")
#   tkgrid(ntracesEntry, row = 4, column = 1, sticky = "w")
  
#   # --- Baseline and Other Parameter Inputs ---
#   baselineVar <<- tclVar("100")
#   stimTimeVar <<- tclVar("350")
#   xbarVar <<- tclVar("100")
#   ybarVar <<- tclVar("50")
  
#   tkgrid(tklabel(sidebarFrame, text = "Baseline:"), row = 4, column = 0, sticky = "w")
#   tkgrid(tkentry(sidebarFrame, textvariable = baselineVar, width = 10), row = 4, column = 1, sticky = "w")
#   tkgrid(tklabel(sidebarFrame, text = "Stimulation Time:"), row = 5, column = 0, sticky = "w")
#   tkgrid(tkentry(sidebarFrame, textvariable = stimTimeVar, width = 10), row = 5, column = 1, sticky = "w")
#   tkgrid(tklabel(sidebarFrame, text = "x-bar length:"), row = 6, column = 0, sticky = "w")
#   tkgrid(tkentry(sidebarFrame, textvariable = xbarVar, width = 10), row = 6, column = 1, sticky = "w")
#   tkgrid(tklabel(sidebarFrame, text = "y-bar length:"), row = 7, column = 0, sticky = "w")
#   tkgrid(tkentry(sidebarFrame, textvariable = ybarVar, width = 10), row = 7, column = 1, sticky = "w")
  
#   consoleText <<- tktext(sidebarFrame, width = 40, height = 4)
#   tkgrid(consoleText, row = 8, column = 0, columnspan = 3, pady = 5)
  
#   # --- Update Function ---
#   updateAdditionalParams <<- function(result) {
#     if (!is.null(result) && length(result$metadata) >= 1) {
#       meta1 <- result$metadata[[1]]
#       tclvalue(dtVar) <<- as.character(meta1$samplingIntervalInSec * 1000)
#       if (!is.null(meta1$header$lActualEpisodes)) {
#         tclvalue(ntracesVar) <<- as.character(meta1$header$lActualEpisodes)
#       } else {
#         tclvalue(ntracesVar) <<- "N/A"
#       }
#       expType <- tclvalue(experimentVar)
#       col_idx <- choose_data_column(meta1$channelUnits, expType)
#       if (!is.na(col_idx)) {
#         tclvalue(unitVar) <<- meta1$channelUnits[col_idx]
#         tclvalue(dataColVar) <<- as.character(col_idx)
#       } else {
#         tclvalue(unitVar) <<- "N/A"
#         tclvalue(dataColVar) <<- "N/A"
#       }
#     }
#   }
  
#   tkbind(experimentCombo, "<<ComboboxSelected>>", function() {
#     if (exists("abf_analysis_result", envir = .GlobalEnv)) {
#       result <- get("abf_analysis_result", envir = .GlobalEnv)
#       updateAdditionalParams(result)
#     }
#   })
  
#   # --- Analysis Action: Load Data button (loads data only, no plotting) ---
#   runAnalysis <<- function() {
#     folderPath <- tclvalue(folderPathVar)
#     if (nchar(folderPath) == 0) {
#       tkmessageBox(message = "Please select an ABF folder first.")
#       return()
#     }
#     selIndices <- as.integer(tkcurselection(abfListBox))
#     allFiles <- as.character(tkget(abfListBox, 0, "end"))
#     abf_files <- if (length(selIndices) == 0) allFiles else allFiles[selIndices + 1]
#     if (length(abf_files) == 0) {
#       tkmessageBox(message = "No ABF files selected.")
#       return()
#     }
#     result <- tryCatch({
#       load_abf_data(abf_files = abf_files, abf_path = folderPath)
#     }, error = function(e) {
#       tkmessageBox(message = paste("Error during data loading:", e$message))
#       return(NULL)
#     })
#     if (!is.null(result)) {
#       tkdelete(consoleText, "1.0", "end")
#       tkinsert(consoleText, "end", paste("Data loaded. Processed", length(abf_files), "file(s)."))
#       assign("abf_analysis_result", result, envir = .GlobalEnv)
#       updateAdditionalParams(result)
#       cons_msg <- check_consistency(result$metadata)
#       if (cons_msg == "Data is consistent") {
#         tkmessageBox(message = cons_msg)
#       } else {
#         tkmessageBox(message = paste("ERROR:", cons_msg))
#       }
#       tkconfigure(runAnalysisButton, text = "Load Data")
#     }
#   }
  
#   runAnalysisButton <<- tkbutton(sidebarFrame, text = "Load Data", command = runAnalysis)
#   tkgrid(runAnalysisButton, row = 9, column = 0, columnspan = 3, pady = 5)
  
#   # --- New Button: Review Recordings (placed above Average Approved Traces) ---
#   reviewButton <<- tkbutton(sidebarFrame, text = "Review Recordings", command = review_recordings)
#   tkgrid(reviewButton, row = 10, column = 0, columnspan = 3, pady = 5)
  
#   # --- New Button: Average Approved Traces ---
#   averageApprovedTraces <<- function() {
#     if (length(traces2average) == 0 || all(sapply(traces2average, length) == 0)) {
#       tkmessageBox(message = "No approved traces available. Please review recordings first.")
#       return()
#     }
#     folderPath <- tclvalue(folderPathVar)
#     if (nchar(folderPath) == 0) {
#       tkmessageBox(message = "Please select an ABF folder first.")
#       return()
#     }
#     selIndices <- as.integer(tkcurselection(abfListBox))
#     allFiles <- as.character(tkget(abfListBox, 0, "end"))
#     abf_files <- if (length(selIndices) == 0) allFiles else allFiles[selIndices + 1]
#     if (length(abf_files) == 0) {
#       tkmessageBox(message = "No ABF files selected.")
#       return()
#     }
#     baseline <- as.numeric(tclvalue(baselineVar))
#     stimTime <- as.numeric(tclvalue(stimTimeVar))
#     xbar <- as.numeric(tclvalue(xbarVar))
#     ybar <- as.numeric(tclvalue(ybarVar))
    
#     result <- tryCatch({
#       abf_out <- abf_averages(datasets = abf_analysis_result$datasets,
#                               traces2average = traces2average,
#                               baseline = baseline, stimulation_time = stimTime,
#                               dataCol = as.numeric(tclvalue(dataColVar)),
#                               xlim = NULL, ylim = NULL,
#                               color = "darkgrey", xbar = xbar, ybar = ybar,
#                               width = 5.25, height = 2.75)
#       return(abf_out)
#     }, error = function(e) {
#       tkmessageBox(message = paste("Error during averaging of approved traces:", e$message))
#       return(NULL)
#     })
#     if (!is.null(result)) {
#       tkdelete(consoleText, "1.0", "end")
#       tkinsert(consoleText, "end", paste("Averaging on approved traces complete. Processed", length(abf_files), "file(s)."))
#       assign("abf_analysis_result", result, envir = .GlobalEnv)
#       if (exists("plotWidget", envir = .GlobalEnv)) tkdestroy(plotWidget)
#       drawPlot <- function() {
#         if (exists("abf_analysis_result", envir = .GlobalEnv)) {
#           result <- get("abf_analysis_result", envir = .GlobalEnv)
#           datasets <- result$datasets
#           if (length(datasets) >= 1) {
#             iii <- 1
#             dt_val <- datasets[[iii]]$samplingIntervalInSec * 1000
#             time <- seq(0, by = dt_val, length.out = nrow(datasets[[iii]]$data[[1]]))
#             channel_index <- as.numeric(tclvalue(dataColVar))
#             if (is.na(channel_index) || channel_index < 1 || channel_index > length(datasets[[iii]]$data)) {
#               channel_index <- 1
#             }
#             trace <- datasets[[iii]]$data[[channel_index]][, 1]
#             egs_plot(x=time, y=trace, show_bar=TRUE, show_text=TRUE)
#             title(main = paste(names(datasets)[iii], "trace 1"), line = 0.5)
#           } else {
#             plot.new()
#             text(0.5, 0.5, "No data available")
#           }
#         } else {
#           plot.new()
#           text(0.5, 0.5, "No analysis result to display")
#         }
#       }
#       plotWidget <<- tkrplot(plotPanel, fun = drawPlot, hscale = 1, vscale = 1)
#       tkgrid(plotWidget, row = 0, column = 0, sticky = "nsew")
#     }
#   }
  
#   averageApprovedButton <<- tkbutton(sidebarFrame, text = "Average Approved Traces", command = averageApprovedTraces)
#   tkgrid(averageApprovedButton, row = 11, column = 0, columnspan = 3, pady = 5)
  
#   tkfocus(tt)
# }

# # -----------------------------------------------------------------------------
# # Launch the ABF Analysis interface:
# # -----------------------------------------------------------------------------
# ABF_analysis_tk()








# changes to code:

# 1. button changes name from Run analysis to load data ... I want name to always be load data
# 2. when loaded chick consistency says  ERROR: data is consistent (should say data is consistent (not ERROR))
# 3. put review recordings button above average selected traces 
# 4. accept reject traces has too many buttons to click I do not need to ok after I press accept otr reject (perhaps just make relevant button appear 'pressed')
# 5. all plotted graphs should appera in  a panel onright side of UI currently they plot to QUARTZ

# # Remove all objects from the environment
# rm(list = ls(all = TRUE))

# # -----------------------------------------------------------------------------
# # Package loading function and required packages:
# # -----------------------------------------------------------------------------
# load_required_packages <- function(packages) {
#   new.packages <- packages[!(packages %in% installed.packages()[, "Package"])]
#   if (length(new.packages)) install.packages(new.packages)
#   invisible(lapply(packages, library, character.only = TRUE))
# }
# required.packages <- c("robustbase", "minpack.lm", "Rcpp", "signal", "readABF", "tcltk", "tkrplot")
# load_required_packages(required.packages)

# # -----------------------------------------------------------------------------
# # Helper Functions
# # -----------------------------------------------------------------------------

# # Extract key metadata from an ABF dataset.
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

# # Choose data column index based on experiment type.
# # For Voltage Clamp, choose the first channel whose unit contains "A" (e.g., pA, nA);
# # For Current Clamp, choose the first channel whose unit contains "V" (e.g., mV, nV).
# choose_data_column <- function(channelUnits, experiment) {
#   if (experiment == "Voltage Clamp") {
#     idx <- grep("A", channelUnits, ignore.case = TRUE)
#   } else if (experiment == "Current Clamp") {
#     idx <- grep("V", channelUnits, ignore.case = TRUE)
#   } else {
#     idx <- integer(0)
#   }
#   if (length(idx) > 0) return(idx[1])
#   else return(NA)
# }

# # Consistency check: Verify that dt, chosen unit, and number of traces are identical across ABF files.
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
#     return("Data is consistent")
#   } else {
#     error_msgs <- c()
#     if (!dt_good) error_msgs <- c(error_msgs, paste("Inconsistent dt values:", paste(dt_values, collapse = ", ")))
#     if (!unit_good) error_msgs <- c(error_msgs, paste("Inconsistent Units:", paste(unit_values, collapse = ", ")))
#     if (!traces_good) error_msgs <- c(error_msgs, paste("Inconsistent Traces:", paste(traces_values, collapse = ", ")))
#     return(paste(error_msgs, collapse = "; "))
#   }
# }

# # Custom plotting function (egs_plot)
# egs_plot <- function(x, y, sign = -1, xlim = NULL, ylim = NULL, lwd = 1, 
#                      show_text = FALSE, height = 4, width = 2.5, 
#                      xbar = 100, ybar = 50, color = "#4C77BB", show_bar = FALSE) {
#   if (is.null(ylim))
#     ylim <- if (sign == 1) c(0, max(y)) else c(-max(-y), 0)
#   if (is.null(xlim))
#     xlim <- c(min(x), max(x))
  
#   idx1 <- which.min(abs(x - xlim[1]))
#   idx2 <- which.min(abs(x - xlim[2]))
  
#   plot(x[idx1:idx2], y[idx1:idx2], type = "l", col = color,
#        xlim = xlim, ylim = ylim, bty = "n", lwd = lwd, lty = 1,
#        axes = FALSE, frame = FALSE, xlab = "", ylab = "")
  
#   if (show_bar) {
#     ybar_start <- min(ylim) + (max(ylim) - min(ylim)) / 20
#     x_start <- max(xlim) - xbar - 50
#     y_start <- ybar_start
#     x_end <- x_start + xbar
#     y_end <- y_start + ybar
#     segments(x_start, y_start, x_end, y_start, lwd = lwd, col = "black")
#     segments(x_start, y_start, x_start, y_end, lwd = lwd, col = "black")
#     if (show_text) {
#       text(x = (x_start + x_end) / 2, y = y_start - ybar / 20, 
#            labels = paste(xbar, "ms"), adj = c(0.5, 1))
#       text(x = x_start - xbar / 4, y = (y_start + y_end) / 2,  
#            labels = paste(ybar, "pA"), adj = c(0.5, 0.5), srt = 90)
#     }
#   }
# }

# # New function to load ABF data (loads files and metadata, no processing)
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
# # Revised abf_averages Function
# # This version expects the preloaded 'datasets', and extracts data from each dataset
# # as follows: for each dataset (iii), for each recording (ii) in datasets[[iii]]$data,
# # extract the column specified by 'dataCol' (which comes from dataColVar).
# # If traces2average is provided, only average those recordings; otherwise, average all.
# # -----------------------------------------------------------------------------
# abf_averages <- function(datasets, 
#                          baseline = 100, 
#                          stimulation_time = 350, 
#                          traces2average = NULL,
#                          dataCol = 1, 
#                          ylim = NULL, 
#                          xlim = NULL, 
#                          color = "#4C77BB", 
#                          xbar = 100, 
#                          ybar = 50, 
#                          width = 5.25, 
#                          height = 2.75, 
#                          save = FALSE) {
#   N <- length(datasets)
#   # Calculate sampling intervals (in ms) for each dataset.
#   sampling_intervals <- sapply(datasets, function(ds) ds$samplingIntervalInSec * 1000)
  
#   # Extract raw responses:
#   # For each dataset (iii), for each recording (ii) in datasets[[iii]]$data,
#   # extract the column specified by 'dataCol' from that matrix.
#   responses <- lapply(seq_len(N), function(iii) {
#     sapply(seq_along(datasets[[iii]]$data), function(ii) {
#       datasets[[iii]]$data[[ii]][, dataCol]
#     })
#   })
#   names(responses) <- names(datasets)
  
#   # Baseline correction helper.
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
  
#   # Average the baseline-corrected responses.
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
#     egs_plot(x = time[[ii]], y = responses0_mean[[ii]], ylim = ylim,
#              xlim = xlim, color = color, width = width, height = height,
#              show_text = FALSE, show_bar = show_bar[ii])
#     axis(1)
#     axis(2)
#   }
  
#   if (save) {
#     warning("save_graph not implemented in this UI example")
#   }
  
#   return(list(raw_data = responses,
#               baseline_corrected_data = responses0,
#               baseline_corrected_mean_data = responses0_mean,
#               datasets = datasets))
# }

# # -----------------------------------------------------------------------------
# # Global variable for approved trace indices (one per dataset)
# # -----------------------------------------------------------------------------
# traces2average <<- list()

# # -----------------------------------------------------------------------------
# # Revised review_recordings Function
# # For each dataset, we now use each element of ds$data (i.e. each recording)
# # and extract the column specified by dataColVar.
# # -----------------------------------------------------------------------------
# review_recordings <- function() {
#   if (!exists("abf_analysis_result", envir = .GlobalEnv)) {
#     tkmessageBox(message = "No analysis result available for review.")
#     return()
#   }
#   result <- get("abf_analysis_result", envir = .GlobalEnv)
#   datasets <- result$datasets
#   traces2average <<- vector("list", length = length(datasets))
#   for(i in seq_along(datasets)) { traces2average[[i]] <<- integer(0) }
  
#   current_dataset <<- 1
#   current_trace <<- 1
  
#   reviewWin <<- tktoplevel()
#   tkwm.title(reviewWin, "Review Recordings")
  
#   infoLabel <<- tklabel(reviewWin, text = paste("Dataset", current_dataset, "Trace", current_trace))
#   tkgrid(infoLabel, row = 0, column = 0, columnspan = 2)
  
#   # In this revised version, for the current dataset,
#   # we treat each element of ds$data as a recording.
#   # From each recording, we extract the trace using the dataCol index.
#   reviewPlot <<- tkrplot(reviewWin, fun = function() {
#     ds <- datasets[[current_dataset]]
#     # We assume that for review, we loop over the recordings (i.e. the list elements in ds$data)
#     if (current_trace > length(ds$data)) {
#       plot.new()
#       text(0.5, 0.5, paste("No more recordings in dataset", current_dataset))
#     } else {
#       # For the current recording, extract the matrix.
#       trace_matrix <- ds$data[[current_trace]]
#       data_column <- as.numeric(tclvalue(dataColVar))
#       if (is.na(data_column) || data_column < 1 || data_column > ncol(trace_matrix)) {
#         data_column <- 1
#       }
#       dt_val <- ds$samplingIntervalInSec * 1000
#       time <- seq(0, by = dt_val, length.out = nrow(trace_matrix))
#       trace <- trace_matrix[, data_column]
#       plot(time, trace, col = "indianred", xlab = "Time (ms)", type = "l",
#            bty = "l", las = 1, main = paste("Dataset", current_dataset, "Recording", current_trace))
#       axis(1)
#       axis(2)
#     }
#   }, hscale = 1.5, vscale = 1.5)
#   tkgrid(reviewPlot, row = 1, column = 0, columnspan = 2)
  
#   acceptButton <- tkbutton(reviewWin, text = "Accept", command = function() {
#     traces2average[[current_dataset]] <<- c(traces2average[[current_dataset]], current_trace)
#     tkmessageBox(message = paste("Accepted dataset", current_dataset, "recording", current_trace))
#   })
#   tkgrid(acceptButton, row = 2, column = 0)
  
#   rejectButton <- tkbutton(reviewWin, text = "Reject", command = function() {
#     tkmessageBox(message = paste("Rejected dataset", current_dataset, "recording", current_trace))
#   })
#   tkgrid(rejectButton, row = 2, column = 1)
  
#   nextTraceButton <- tkbutton(reviewWin, text = "Next Recording", command = function() {
#     ds <- datasets[[current_dataset]]
#     numRecordings <- length(ds$data)
#     if (current_trace < numRecordings) {
#       current_trace <<- current_trace + 1
#     } else {
#       tkmessageBox(message = paste("Finished reviewing dataset", current_dataset))
#       if (current_dataset < length(datasets)) {
#         current_dataset <<- current_dataset + 1
#         current_trace <<- 1
#       } else {
#         tkmessageBox(message = "Review complete. Approved recordings are in 'traces2average'.")
#         tkdestroy(reviewWin)
#         return()
#       }
#     }
#     tkconfigure(infoLabel, text = paste("Dataset", current_dataset, "Recording", current_trace))
#     tkrreplot(reviewPlot)
#   })
#   tkgrid(nextTraceButton, row = 3, column = 0, columnspan = 2)
# }

# # -----------------------------------------------------------------------------
# # Main ABF Analysis UI
# # -----------------------------------------------------------------------------
# ABF_analysis_tk <- function() {
#   tt <- tktoplevel()
#   tkwm.title(tt, "ABF Analysis")
  
#   sidebarFrame <- tkframe(tt)
#   mainFrame <- tkframe(tt)
#   tkgrid(sidebarFrame, row = 0, column = 0, sticky = "ns")
#   tkgrid(mainFrame, row = 0, column = 1, sticky = "nsew")
#   tkgrid.rowconfigure(tt, 0, weight = 1)
#   tkgrid.columnconfigure(tt, 1, weight = 1)
  
#   # --- Sidebar Controls ---
#   folderLabel <- tklabel(sidebarFrame, text = "Select ABF Folder:")
#   tkgrid(folderLabel, row = 0, column = 0, sticky = "w")
#   folderPathVar <<- tclVar("")
#   folderEntry <- tkentry(sidebarFrame, textvariable = folderPathVar, width = 30)
#   tkgrid(folderEntry, row = 0, column = 1, sticky = "w")
  
#   browseFolderButton <- tkbutton(sidebarFrame, text = "Browse", command = function() {
#     folderPath <- tclvalue(tkchooseDirectory())
#     if (nchar(folderPath) > 0) {
#       tclvalue(folderPathVar) <<- folderPath
#       abf_list <- list.files(path = folderPath, pattern = "\\.abf$", ignore.case = TRUE)
#       if (length(abf_list) == 0) {
#         tkmessageBox(message = "No ABF files found in the selected folder.")
#       } else {
#         tkdelete(abfListBox, 0, "end")
#         for (f in abf_list) { tkinsert(abfListBox, "end", f) }
#         firstFilePath <- file.path(folderPath, abf_list[1])
#         ds <- readABF(firstFilePath)
#         dummy_result <- list(metadata = list(extract_metadata(ds)))
#         updateAdditionalParams(dummy_result)
#       }
#     }
#   })
#   tkgrid(browseFolderButton, row = 0, column = 2, padx = 5)
  
#   abfListLabel <- tklabel(sidebarFrame, text = "ABF Files:")
#   tkgrid(abfListLabel, row = 1, column = 0, sticky = "w", pady = 5)
#   abfListBox <<- tklistbox(sidebarFrame, height = 5, width = 40, selectmode = "multiple")
#   tkgrid(abfListBox, row = 2, column = 0, columnspan = 3, pady = 5, sticky = "we")
  
#   # --- Additional Parameters (above Baseline) ---
#   paramFrame <- tkframe(sidebarFrame)
#   tkgrid(paramFrame, row = 3, column = 0, columnspan = 3, sticky = "w")
  
#   experimentVar <<- tclVar("Voltage Clamp")
#   unitVar <<- tclVar("")
#   dataColVar <<- tclVar("")
#   dtVar <<- tclVar("")
#   ntracesVar <<- tclVar("")
  
#   tkgrid(tklabel(paramFrame, text = "Experiment:"), row = 0, column = 0, sticky = "w")
#   experimentCombo <- ttkcombobox(paramFrame, textvariable = experimentVar, 
#                                  values = c("Voltage Clamp", "Current Clamp"), width = 15)
#   tkgrid(experimentCombo, row = 0, column = 1, sticky = "w")
  
#   tkgrid(tklabel(paramFrame, text = "Units:"), row = 1, column = 0, sticky = "w")
#   unitEntry <- tkentry(paramFrame, textvariable = unitVar, width = 10, state = "readonly")
#   tkgrid(unitEntry, row = 1, column = 1, sticky = "w")
  
#   tkgrid(tklabel(paramFrame, text = "Data Column:"), row = 2, column = 0, sticky = "w")
#   dataColEntry <- tkentry(paramFrame, textvariable = dataColVar, width = 10, state = "readonly")
#   tkgrid(dataColEntry, row = 2, column = 1, sticky = "w")
  
#   tkgrid(tklabel(paramFrame, text = "dt (ms):"), row = 3, column = 0, sticky = "w")
#   dtEntry <- tkentry(paramFrame, textvariable = dtVar, width = 10, state = "readonly")
#   tkgrid(dtEntry, row = 3, column = 1, sticky = "w")
  
#   tkgrid(tklabel(paramFrame, text = "Traces:"), row = 4, column = 0, sticky = "w")
#   ntracesEntry <- tkentry(paramFrame, textvariable = ntracesVar, width = 10, state = "readonly")
#   tkgrid(ntracesEntry, row = 4, column = 1, sticky = "w")
  
#   # --- Baseline and Other Parameter Inputs ---
#   baselineVar <<- tclVar("100")
#   stimTimeVar <<- tclVar("350")
#   xbarVar <<- tclVar("100")
#   ybarVar <<- tclVar("50")
  
#   tkgrid(tklabel(sidebarFrame, text = "Baseline:"), row = 4, column = 0, sticky = "w")
#   tkgrid(tkentry(sidebarFrame, textvariable = baselineVar, width = 10), row = 4, column = 1, sticky = "w")
#   tkgrid(tklabel(sidebarFrame, text = "Stimulation Time:"), row = 5, column = 0, sticky = "w")
#   tkgrid(tkentry(sidebarFrame, textvariable = stimTimeVar, width = 10), row = 5, column = 1, sticky = "w")
#   tkgrid(tklabel(sidebarFrame, text = "x-bar length:"), row = 6, column = 0, sticky = "w")
#   tkgrid(tkentry(sidebarFrame, textvariable = xbarVar, width = 10), row = 6, column = 1, sticky = "w")
#   tkgrid(tklabel(sidebarFrame, text = "y-bar length:"), row = 7, column = 0, sticky = "w")
#   tkgrid(tkentry(sidebarFrame, textvariable = ybarVar, width = 10), row = 7, column = 1, sticky = "w")
  
#   consoleText <<- tktext(sidebarFrame, width = 40, height = 4)
#   tkgrid(consoleText, row = 8, column = 0, columnspan = 3, pady = 5)
  
#   # --- Update Function ---
#   updateAdditionalParams <<- function(result) {
#     if (!is.null(result) && length(result$metadata) >= 1) {
#       meta1 <- result$metadata[[1]]
#       tclvalue(dtVar) <<- as.character(meta1$samplingIntervalInSec * 1000)
#       if (!is.null(meta1$header$lActualEpisodes)) {
#         tclvalue(ntracesVar) <<- as.character(meta1$header$lActualEpisodes)
#       } else {
#         tclvalue(ntracesVar) <<- "N/A"
#       }
#       expType <- tclvalue(experimentVar)
#       col_idx <- choose_data_column(meta1$channelUnits, expType)
#       if (!is.na(col_idx)) {
#         tclvalue(unitVar) <<- meta1$channelUnits[col_idx]
#         tclvalue(dataColVar) <<- as.character(col_idx)
#       } else {
#         tclvalue(unitVar) <<- "N/A"
#         tclvalue(dataColVar) <<- "N/A"
#       }
#     }
#   }
  
#   tkbind(experimentCombo, "<<ComboboxSelected>>", function() {
#     if (exists("abf_analysis_result", envir = .GlobalEnv)) {
#       result <- get("abf_analysis_result", envir = .GlobalEnv)
#       updateAdditionalParams(result)
#     }
#   })
  
#   # --- Analysis Action: Run Analysis button (loads data only, no plotting) ---
#   runAnalysis <<- function() {
#     folderPath <- tclvalue(folderPathVar)
#     if (nchar(folderPath) == 0) {
#       tkmessageBox(message = "Please select an ABF folder first.")
#       return()
#     }
#     selIndices <- as.integer(tkcurselection(abfListBox))
#     allFiles <- as.character(tkget(abfListBox, 0, "end"))
#     abf_files <- if (length(selIndices) == 0) allFiles else allFiles[selIndices + 1]
#     if (length(abf_files) == 0) {
#       tkmessageBox(message = "No ABF files selected.")
#       return()
#     }
#     result <- tryCatch({
#       load_abf_data(abf_files = abf_files, abf_path = folderPath)
#     }, error = function(e) {
#       tkmessageBox(message = paste("Error during data loading:", e$message))
#       return(NULL)
#     })
#     if (!is.null(result)) {
#       tkdelete(consoleText, "1.0", "end")
#       tkinsert(consoleText, "end", paste("Data loaded. Processed", length(abf_files), "file(s)."))
#       assign("abf_analysis_result", result, envir = .GlobalEnv)
#       updateAdditionalParams(result)
#       cons_msg <- check_consistency(result$metadata)
#       if (cons_msg != "GOOD") {
#         tkmessageBox(message = paste("ERROR:", cons_msg))
#       } else {
#         tkmessageBox(message = "GOOD: All parameters are consistent.")
#       }
#       tkconfigure(runAnalysisButton, text = "Load Data")
#       tkgrid.forget(plotWidget)  # Remove the main plot (if any)
#     }
#   }
  
#   runAnalysisButton <<- tkbutton(sidebarFrame, text = "Run Analysis", command = runAnalysis)
#   tkgrid(runAnalysisButton, row = 9, column = 0, columnspan = 3, pady = 5)
  
#   # --- New Button: Average Approved Traces ---
#   averageApprovedTraces <<- function() {
#     if (length(traces2average) == 0 || all(sapply(traces2average, length) == 0)) {
#       tkmessageBox(message = "No approved traces available. Please review recordings first.")
#       return()
#     }
#     folderPath <- tclvalue(folderPathVar)
#     if (nchar(folderPath) == 0) {
#       tkmessageBox(message = "Please select an ABF folder first.")
#       return()
#     }
#     selIndices <- as.integer(tkcurselection(abfListBox))
#     allFiles <- as.character(tkget(abfListBox, 0, "end"))
#     abf_files <- if (length(selIndices) == 0) allFiles else allFiles[selIndices + 1]
#     if (length(abf_files) == 0) {
#       tkmessageBox(message = "No ABF files selected.")
#       return()
#     }
#     baseline <- as.numeric(tclvalue(baselineVar))
#     stimTime <- as.numeric(tclvalue(stimTimeVar))
#     xbar <- as.numeric(tclvalue(xbarVar))
#     ybar <- as.numeric(tclvalue(ybarVar))
    
#     result <- tryCatch({
#       abf_out <- abf_averages(datasets = abf_analysis_result$datasets,
#                               traces2average = traces2average,
#                               baseline = baseline, stimulation_time = stimTime,
#                               dataCol = as.numeric(tclvalue(dataColVar)),
#                               xlim = NULL, ylim = NULL,
#                               color = "darkgrey", xbar = xbar, ybar = ybar,
#                               width = 5.25, height = 2.75)
#       return(abf_out)
#     }, error = function(e) {
#       tkmessageBox(message = paste("Error during averaging of approved traces:", e$message))
#       return(NULL)
#     })
#     if (!is.null(result)) {
#       tkdelete(consoleText, "1.0", "end")
#       tkinsert(consoleText, "end", paste("Averaging on approved traces complete. Processed", length(abf_files), "file(s)."))
#       assign("abf_analysis_result", result, envir = .GlobalEnv)
#       tkrreplot(plotWidget)
#     }
#   }
  
#   averageApprovedButton <<- tkbutton(sidebarFrame, text = "Average Approved Traces", command = averageApprovedTraces)
#   tkgrid(averageApprovedButton, row = 10, column = 0, columnspan = 3, pady = 5)
  
#   # --- New Button: Review Recordings ---
#   reviewButton <<- tkbutton(sidebarFrame, text = "Review Recordings", command = review_recordings)
#   tkgrid(reviewButton, row = 11, column = 0, columnspan = 3, pady = 5)
  
#   # --- Main Panel: Plot display via tkrplot ---
#   drawPlot <- function() {
#     if (exists("abf_analysis_result", envir = .GlobalEnv)) {
#       result <- get("abf_analysis_result", envir = .GlobalEnv)
#       datasets <- result$datasets
#       if (length(datasets) >= 1) {
#         iii <- 1
#         dt_val <- datasets[[iii]]$samplingIntervalInSec * 1000
#         time <- seq(0, by = dt_val, length.out = nrow(datasets[[iii]]$data[[1]]))
#         channel_index <- as.numeric(tclvalue(dataColVar))
#         if (is.na(channel_index) || channel_index < 1 || channel_index > length(datasets[[iii]]$data)) {
#           channel_index <- 1
#         }
#         # In the main panel, we simply display the first recording of the chosen channel.
#         trace <- datasets[[iii]]$data[[channel_index]][, 1]
#         plot(time, trace, col = "indianred", xlab = "Time (ms)", type = "l",
#              bty = "l", las = 1, main = paste("Trace", iii))
#         axis(1)
#         axis(2)
#       } else {
#         plot.new()
#         text(0.5, 0.5, "No data available")
#       }
#     } else {
#       plot.new()
#       text(0.5, 0.5, "No analysis result to display")
#     }
#   }
  
#   plotWidget <<- tkrplot(mainFrame, fun = drawPlot, hscale = 1.5, vscale = 1.5)
#   tkgrid(plotWidget, row = 0, column = 0, sticky = "nsew")
  
#   tkfocus(tt)
# }

# # -----------------------------------------------------------------------------
# # Launch the ABF Analysis interface:
# # -----------------------------------------------------------------------------
# ABF_analysis_tk()






# # Remove all objects from the environment
# rm(list = ls(all = TRUE))

# # -----------------------------------------------------------------------------
# # Package loading function and required packages:
# # -----------------------------------------------------------------------------
# load_required_packages <- function(packages) {
#   new.packages <- packages[!(packages %in% installed.packages()[, "Package"])]
#   if (length(new.packages)) install.packages(new.packages)
#   invisible(lapply(packages, library, character.only = TRUE))
# }
# required.packages <- c("robustbase", "minpack.lm", "Rcpp", "signal", "readABF", "tcltk", "tkrplot")
# load_required_packages(required.packages)

# # -----------------------------------------------------------------------------
# # Helper functions
# # -----------------------------------------------------------------------------

# # Function to extract key metadata from an ABF dataset
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

# # Function to choose data column index based on experiment type.
# # For Voltage Clamp, choose the first channel whose unit contains "A"
# # (e.g., pA, nA). For Current Clamp, choose the first channel whose unit
# # contains "V" (e.g., mV, nV).
# choose_data_column <- function(channelUnits, experiment) {
#   if (experiment == "Voltage Clamp") {
#     idx <- grep("A", channelUnits, ignore.case = TRUE)
#   } else if (experiment == "Current Clamp") {
#     idx <- grep("V", channelUnits, ignore.case = TRUE)
#   } else {
#     idx <- integer(0)
#   }
#   if (length(idx) > 0) return(idx[1])
#   else return(NA)
# }

# # A consistency check: verify that dt, the chosen unit, and the number of traces are identical across all ABF files.
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
#     return("GOOD")
#   } else {
#     error_msgs <- c()
#     if (!dt_good) error_msgs <- c(error_msgs, paste("Inconsistent dt values:", paste(dt_values, collapse = ", ")))
#     if (!unit_good) error_msgs <- c(error_msgs, paste("Inconsistent Units:", paste(unit_values, collapse = ", ")))
#     if (!traces_good) error_msgs <- c(error_msgs, paste("Inconsistent Traces:", paste(traces_values, collapse = ", ")))
#     return(paste(error_msgs, collapse = "; "))
#   }
# }

# # -----------------------------------------------------------------------------
# # Custom plotting function (egs_plot)
# # -----------------------------------------------------------------------------
# egs_plot <- function(x, y, sign = -1, xlim = NULL, ylim = NULL, lwd = 1, 
#                      show_text = FALSE, height = 4, width = 2.5, 
#                      xbar = 100, ybar = 50, color = "#4C77BB", show_bar = FALSE) {
#   if (is.null(ylim))
#     ylim <- if (sign == 1) c(0, max(y)) else c(-max(-y), 0)
#   if (is.null(xlim))
#     xlim <- c(min(x), max(x))
  
#   idx1 <- which.min(abs(x - xlim[1]))
#   idx2 <- which.min(abs(x - xlim[2]))
  
#   plot(x[idx1:idx2], y[idx1:idx2], type = "l", col = color,
#        xlim = xlim, ylim = ylim, bty = "n", lwd = lwd, lty = 1,
#        axes = FALSE, frame = FALSE, xlab = "", ylab = "")
  
#   if (show_bar) {
#     ybar_start <- min(ylim) + (max(ylim) - min(ylim)) / 20
#     x_start <- max(xlim) - xbar - 50
#     y_start <- ybar_start
#     x_end <- x_start + xbar
#     y_end <- y_start + ybar
#     segments(x_start, y_start, x_end, y_start, lwd = lwd, col = "black")
#     segments(x_start, y_start, x_start, y_end, lwd = lwd, col = "black")
#     if (show_text) {
#       text(x = (x_start + x_end) / 2, y = y_start - ybar / 20, 
#            labels = paste(xbar, "ms"), adj = c(0.5, 1))
#       text(x = x_start - xbar / 4, y = (y_start + y_end) / 2,  
#            labels = paste(ybar, "pA"), adj = c(0.5, 0.5), srt = 90)
#     }
#   }
# }

# # -----------------------------------------------------------------------------
# # New function to load ABF data without plotting (for Run Analysis)
# # -----------------------------------------------------------------------------
# load_abf_data <- function(abf_files = NULL, abf_path = NULL) {
#   abf_path <- if (is.null(abf_path)) getwd() else abf_path
#   setwd(abf_path)
#   N <- length(abf_files)
#   datasets <- lapply(1:N, function(ii) readABF(abf_files[ii]))
#   names(datasets) <- abf_files
#   metadata <- lapply(datasets, extract_metadata)
#   return(list(datasets = datasets, metadata = metadata))
# }

# # -----------------------------------------------------------------------------
# # Global variable to hold approved trace indices (one element per dataset)
# # -----------------------------------------------------------------------------
# traces2average <<- list()

# # -----------------------------------------------------------------------------
# # Function to review recordings (trace-by-trace) for each dataset.
# # -----------------------------------------------------------------------------
# review_recordings <- function() {
#   if (!exists("abf_analysis_result", envir = .GlobalEnv)) {
#     tkmessageBox(message = "No analysis result available for review.")
#     return()
#   }
#   result <- get("abf_analysis_result", envir = .GlobalEnv)
#   datasets <- result$datasets
#   traces2average <<- vector("list", length = length(datasets))
#   for(i in seq_along(datasets)) { traces2average[[i]] <<- integer(0) }
  
#   current_dataset <<- 1
#   current_trace <<- 1
  
#   reviewWin <<- tktoplevel()
#   tkwm.title(reviewWin, "Review Recordings")
  
#   infoLabel <<- tklabel(reviewWin, text = paste("Dataset", current_dataset, "Trace", current_trace))
#   tkgrid(infoLabel, row = 0, column = 0, columnspan = 2)
  
#   reviewPlot <<- tkrplot(reviewWin, fun = function() {
#     ds <- datasets[[current_dataset]]
#     dt_val <- ds$samplingIntervalInSec * 1000
#     time <- seq(0, by = dt_val, length.out = nrow(ds$data[[1]]))
#     col_num <- as.numeric(tclvalue(dataColVar))
#     if (is.na(col_num) || col_num < 1 || col_num > ncol(ds$data[[1]])) { col_num <- 1 }
#     if (current_trace > ncol(ds$data[[1]])) {
#       plot.new()
#       text(0.5, 0.5, paste("No more traces in dataset", current_dataset))
#     } else {
#       trace <- ds$data[[1]][, current_trace]
#       plot(time, trace, col = "indianred", xlab = "Time (ms)", type = "l",
#            bty = "l", las = 1, main = paste("Dataset", current_dataset, "Trace", current_trace))
#       axis(1)
#       axis(2)
#     }
#   }, hscale = 1.5, vscale = 1.5)
#   tkgrid(reviewPlot, row = 1, column = 0, columnspan = 2)
  
#   acceptButton <- tkbutton(reviewWin, text = "Accept", command = function() {
#     traces2average[[current_dataset]] <<- c(traces2average[[current_dataset]], current_trace)
#     tkmessageBox(message = paste("Accepted dataset", current_dataset, "trace", current_trace))
#   })
#   tkgrid(acceptButton, row = 2, column = 0)
  
#   rejectButton <- tkbutton(reviewWin, text = "Reject", command = function() {
#     tkmessageBox(message = paste("Rejected dataset", current_dataset, "trace", current_trace))
#   })
#   tkgrid(rejectButton, row = 2, column = 1)
  
#   nextTraceButton <- tkbutton(reviewWin, text = "Next Trace", command = function() {
#     ds <- datasets[[current_dataset]]
#     numTraces <- ncol(ds$data[[1]])
#     if (current_trace < numTraces) {
#       current_trace <<- current_trace + 1
#     } else {
#       tkmessageBox(message = paste("Finished reviewing dataset", current_dataset))
#       if (current_dataset < length(datasets)) {
#         current_dataset <<- current_dataset + 1
#         current_trace <<- 1
#       } else {
#         tkmessageBox(message = "Review complete. Final accepted traces are stored in 'traces2average'.")
#         tkdestroy(reviewWin)
#         return()
#       }
#     }
#     tkconfigure(infoLabel, text = paste("Dataset", current_dataset, "Trace", current_trace))
#     tkrreplot(reviewPlot)
#   })
#   tkgrid(nextTraceButton, row = 3, column = 0, columnspan = 2)
# }

# # -----------------------------------------------------------------------------
# # Now build the main ABF Analysis UI
# # -----------------------------------------------------------------------------
# ABF_analysis_tk <- function() {
#   tt <- tktoplevel()
#   tkwm.title(tt, "ABF Analysis")
  
#   sidebarFrame <- tkframe(tt)
#   mainFrame <- tkframe(tt)
#   tkgrid(sidebarFrame, row = 0, column = 0, sticky = "ns")
#   tkgrid(mainFrame, row = 0, column = 1, sticky = "nsew")
#   tkgrid.rowconfigure(tt, 0, weight = 1)
#   tkgrid.columnconfigure(tt, 1, weight = 1)
  
#   # --- Sidebar Controls ---
  
#   folderLabel <- tklabel(sidebarFrame, text = "Select ABF Folder:")
#   tkgrid(folderLabel, row = 0, column = 0, sticky = "w")
#   folderPathVar <<- tclVar("")
#   folderEntry <- tkentry(sidebarFrame, textvariable = folderPathVar, width = 30)
#   tkgrid(folderEntry, row = 0, column = 1, sticky = "w")
  
#   browseFolderButton <- tkbutton(sidebarFrame, text = "Browse", command = function() {
#     folderPath <- tclvalue(tkchooseDirectory())
#     if (nchar(folderPath) > 0) {
#       tclvalue(folderPathVar) <<- folderPath
#       abf_list <- list.files(path = folderPath, pattern = "\\.abf$", ignore.case = TRUE)
#       if (length(abf_list) == 0) {
#         tkmessageBox(message = "No ABF files found in the selected folder.")
#       } else {
#         tkdelete(abfListBox, 0, "end")
#         for (f in abf_list) { tkinsert(abfListBox, "end", f) }
#         firstFilePath <- file.path(folderPath, abf_list[1])
#         ds <- readABF(firstFilePath)
#         dummy_result <- list(metadata = list(extract_metadata(ds)))
#         updateAdditionalParams(dummy_result)
#       }
#     }
#   })
#   tkgrid(browseFolderButton, row = 0, column = 2, padx = 5)
  
#   abfListLabel <- tklabel(sidebarFrame, text = "ABF Files:")
#   tkgrid(abfListLabel, row = 1, column = 0, sticky = "w", pady = 5)
#   abfListBox <<- tklistbox(sidebarFrame, height = 5, width = 40, selectmode = "multiple")
#   tkgrid(abfListBox, row = 2, column = 0, columnspan = 3, pady = 5, sticky = "we")
  
#   # --- Additional Parameters (above Baseline) ---
#   paramFrame <- tkframe(sidebarFrame)
#   tkgrid(paramFrame, row = 3, column = 0, columnspan = 3, sticky = "w")
  
#   experimentVar <<- tclVar("Voltage Clamp")
#   unitVar <<- tclVar("")
#   dataColVar <<- tclVar("")
#   dtVar <<- tclVar("")
#   ntracesVar <<- tclVar("")
  
#   tkgrid(tklabel(paramFrame, text = "Experiment:"), row = 0, column = 0, sticky = "w")
#   experimentCombo <- ttkcombobox(paramFrame, textvariable = experimentVar, 
#                                  values = c("Voltage Clamp", "Current Clamp"), width = 15)
#   tkgrid(experimentCombo, row = 0, column = 1, sticky = "w")
  
#   tkgrid(tklabel(paramFrame, text = "Units:"), row = 1, column = 0, sticky = "w")
#   unitEntry <- tkentry(paramFrame, textvariable = unitVar, width = 10, state = "readonly")
#   tkgrid(unitEntry, row = 1, column = 1, sticky = "w")
  
#   tkgrid(tklabel(paramFrame, text = "Data Column:"), row = 2, column = 0, sticky = "w")
#   dataColEntry <- tkentry(paramFrame, textvariable = dataColVar, width = 10, state = "readonly")
#   tkgrid(dataColEntry, row = 2, column = 1, sticky = "w")
  
#   tkgrid(tklabel(paramFrame, text = "dt (ms):"), row = 3, column = 0, sticky = "w")
#   dtEntry <- tkentry(paramFrame, textvariable = dtVar, width = 10, state = "readonly")
#   tkgrid(dtEntry, row = 3, column = 1, sticky = "w")
  
#   tkgrid(tklabel(paramFrame, text = "Traces:"), row = 4, column = 0, sticky = "w")
#   ntracesEntry <- tkentry(paramFrame, textvariable = ntracesVar, width = 10, state = "readonly")
#   tkgrid(ntracesEntry, row = 4, column = 1, sticky = "w")
  
#   # --- Baseline and Other Parameter Inputs ---
#   baselineVar <<- tclVar("100")
#   stimTimeVar <<- tclVar("350")
#   xbarVar <<- tclVar("100")
#   ybarVar <<- tclVar("50")
  
#   tkgrid(tklabel(sidebarFrame, text = "Baseline:"), row = 4, column = 0, sticky = "w")
#   tkgrid(tkentry(sidebarFrame, textvariable = baselineVar, width = 10), row = 4, column = 1, sticky = "w")
#   tkgrid(tklabel(sidebarFrame, text = "Stimulation Time:"), row = 5, column = 0, sticky = "w")
#   tkgrid(tkentry(sidebarFrame, textvariable = stimTimeVar, width = 10), row = 5, column = 1, sticky = "w")
#   tkgrid(tklabel(sidebarFrame, text = "x-bar length:"), row = 6, column = 0, sticky = "w")
#   tkgrid(tkentry(sidebarFrame, textvariable = xbarVar, width = 10), row = 6, column = 1, sticky = "w")
#   tkgrid(tklabel(sidebarFrame, text = "y-bar length:"), row = 7, column = 0, sticky = "w")
#   tkgrid(tkentry(sidebarFrame, textvariable = ybarVar, width = 10), row = 7, column = 1, sticky = "w")
  
#   consoleText <<- tktext(sidebarFrame, width = 40, height = 4)
#   tkgrid(consoleText, row = 8, column = 0, columnspan = 3, pady = 5)
  
#   # --- Update Function ---
#   updateAdditionalParams <<- function(result) {
#     if (!is.null(result) && length(result$metadata) >= 1) {
#       meta1 <- result$metadata[[1]]
#       tclvalue(dtVar) <<- as.character(meta1$samplingIntervalInSec * 1000)
#       if (!is.null(meta1$header$lActualEpisodes)) {
#         tclvalue(ntracesVar) <<- as.character(meta1$header$lActualEpisodes)
#       } else {
#         tclvalue(ntracesVar) <<- "N/A"
#       }
#       expType <- tclvalue(experimentVar)
#       col_idx <- choose_data_column(meta1$channelUnits, expType)
#       if (!is.na(col_idx)) {
#         tclvalue(unitVar) <<- meta1$channelUnits[col_idx]
#         tclvalue(dataColVar) <<- as.character(col_idx)
#       } else {
#         tclvalue(unitVar) <<- "N/A"
#         tclvalue(dataColVar) <<- "N/A"
#       }
#     }
#   }
  
#   tkbind(experimentCombo, "<<ComboboxSelected>>", function() {
#     if (exists("abf_analysis_result", envir = .GlobalEnv)) {
#       result <- get("abf_analysis_result", envir = .GlobalEnv)
#       updateAdditionalParams(result)
#     }
#   })
  
#   # --- Analysis Action: Run Analysis button (loads data only, no plotting) ---
#   runAnalysis <<- function() {
#     folderPath <- tclvalue(folderPathVar)
#     if (nchar(folderPath) == 0) {
#       tkmessageBox(message = "Please select an ABF folder first.")
#       return()
#     }
#     selIndices <- as.integer(tkcurselection(abfListBox))
#     allFiles <- as.character(tkget(abfListBox, 0, "end"))
#     abf_files <- if (length(selIndices) == 0) allFiles else allFiles[selIndices + 1]
#     if (length(abf_files) == 0) {
#       tkmessageBox(message = "No ABF files selected.")
#       return()
#     }
#     # In this step we load the data only (without averaging/plotting)
#     result <- tryCatch({
#       load_abf_data(abf_files = abf_files, abf_path = folderPath)
#     }, error = function(e) {
#       tkmessageBox(message = paste("Error during data loading:", e$message))
#       return(NULL)
#     })
#     if (!is.null(result)) {
#       tkdelete(consoleText, "1.0", "end")
#       tkinsert(consoleText, "end", paste("Data loaded. Processed", length(abf_files), "file(s)."))
#       assign("abf_analysis_result", result, envir = .GlobalEnv)
#       updateAdditionalParams(result)
#       cons_msg <- check_consistency(result$metadata)
#       if (cons_msg != "GOOD") {
#         tkmessageBox(message = paste("ERROR:", cons_msg))
#       } else {
#         tkmessageBox(message = "GOOD: All parameters are consistent.")
#       }
#       tkconfigure(runAnalysisButton, text = "Load Data")
#       tkgrid.forget(plotWidget)  # Remove the main plot (if any)
#     }
#   }
  
#   runAnalysisButton <<- tkbutton(sidebarFrame, text = "Run Analysis", command = runAnalysis)
#   tkgrid(runAnalysisButton, row = 9, column = 0, columnspan = 3, pady = 5)
  
#   # --- New Button: Average Approved Traces ---
#   averageApprovedTraces <<- function() {
#     if (length(traces2average) == 0 || all(sapply(traces2average, length) == 0)) {
#       tkmessageBox(message = "No approved traces available. Please review recordings first.")
#       return()
#     }
#     folderPath <- tclvalue(folderPathVar)
#     if (nchar(folderPath) == 0) {
#       tkmessageBox(message = "Please select an ABF folder first.")
#       return()
#     }
#     selIndices <- as.integer(tkcurselection(abfListBox))
#     allFiles <- as.character(tkget(abfListBox, 0, "end"))
#     abf_files <- if (length(selIndices) == 0) allFiles else allFiles[selIndices + 1]
#     if (length(abf_files) == 0) {
#       tkmessageBox(message = "No ABF files selected.")
#       return()
#     }
#     baseline <- as.numeric(tclvalue(baselineVar))
#     stimTime <- as.numeric(tclvalue(stimTimeVar))
#     xbar <- as.numeric(tclvalue(xbarVar))
#     ybar <- as.numeric(tclvalue(ybarVar))
    
#     result <- tryCatch({
#       abf_out <- abf_averages(abf_files = abf_files, abf_path = folderPath,
#                               traces2average = traces2average,
#                               baseline = baseline, stimulation_time = stimTime,
#                               xbar = xbar, ybar = ybar)
#       return(abf_out)
#     }, error = function(e) {
#       tkmessageBox(message = paste("Error during averaging of approved traces:", e$message))
#       return(NULL)
#     })
    
#     if (!is.null(result)) {
#       tkdelete(consoleText, "1.0", "end")
#       tkinsert(consoleText, "end", paste("Averaging on approved traces complete. Processed", length(abf_files), "file(s)."))
#       assign("abf_analysis_result", result, envir = .GlobalEnv)
#       tkrreplot(plotWidget)
#     }
#   }
  
#   averageApprovedButton <<- tkbutton(sidebarFrame, text = "Average Approved Traces", command = averageApprovedTraces)
#   tkgrid(averageApprovedButton, row = 10, column = 0, columnspan = 3, pady = 5)
  
#   # --- New Button: Review Recordings ---
#   reviewButton <<- tkbutton(sidebarFrame, text = "Review Recordings", command = review_recordings)
#   tkgrid(reviewButton, row = 11, column = 0, columnspan = 3, pady = 5)
  
#   # --- Main Panel: Plot display via tkrplot ---
#   drawPlot <- function() {
#     if (exists("abf_analysis_result", envir = .GlobalEnv)) {
#       result <- get("abf_analysis_result", envir = .GlobalEnv)
#       datasets <- result$datasets
#       if (length(datasets) >= 1) {
#         iii <- 1
#         dt_val <- datasets[[iii]]$samplingIntervalInSec * 1000
#         time <- seq(0, by = dt_val, length.out = nrow(datasets[[iii]]$data[[1]]))
#         col_num <- as.numeric(tclvalue(dataColVar))
#         if (is.na(col_num) || col_num < 1 || col_num > ncol(datasets[[iii]]$data[[1]])) {
#           col_num <- 1
#         }
#         trace <- datasets[[iii]]$data[[1]][, col_num]
#         plot(time, trace, col = "indianred", xlab = "Time (ms)", type = "l",
#              bty = "l", las = 1, main = paste("Trace", iii))
#         axis(1)
#         axis(2)
#       } else {
#         plot.new()
#         text(0.5, 0.5, "No data available")
#       }
#     } else {
#       plot.new()
#       text(0.5, 0.5, "No analysis result to display")
#     }
#   }
  
#   plotWidget <<- tkrplot(mainFrame, fun = drawPlot, hscale = 1.5, vscale = 1.5)
#   tkgrid(plotWidget, row = 0, column = 0, sticky = "nsew")
#   tkfocus(tt)
# }

# # -----------------------------------------------------------------------------
# # Launch the ABF Analysis interface:
# # -----------------------------------------------------------------------------
# ABF_analysis_tk()



# # Remove all objects from the environment
# rm(list = ls(all = TRUE))

# # -----------------------------------------------------------------------------
# # Package loading function and required packages:
# # -----------------------------------------------------------------------------
# load_required_packages <- function(packages) {
#   new.packages <- packages[!(packages %in% installed.packages()[, "Package"])]
#   if (length(new.packages)) install.packages(new.packages)
#   invisible(lapply(packages, library, character.only = TRUE))
# }
# required.packages <- c("robustbase", "minpack.lm", "Rcpp", "signal", "readABF", "tcltk", "tkrplot")
# load_required_packages(required.packages)

# # -----------------------------------------------------------------------------
# # Helper functions
# # -----------------------------------------------------------------------------

# # Extract key metadata from an ABF dataset
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

# # Choose the data column index based on experiment type.
# # For Voltage Clamp choose the first channel whose unit contains "A";
# # For Current Clamp choose the first channel whose unit contains "V".
# choose_data_column <- function(channelUnits, experiment) {
#   if (experiment == "Voltage Clamp") {
#     idx <- grep("A", channelUnits, ignore.case = TRUE)
#   } else if (experiment == "Current Clamp") {
#     idx <- grep("V", channelUnits, ignore.case = TRUE)
#   } else {
#     idx <- integer(0)
#   }
#   if (length(idx) > 0) return(idx[1])
#   else return(NA)
# }

# # -----------------------------------------------------------------------------
# # Custom plotting function (egs_plot)
# # -----------------------------------------------------------------------------
# egs_plot <- function(x, y, sign = -1, xlim = NULL, ylim = NULL, lwd = 1, 
#                      show_text = FALSE, height = 4, width = 2.5, 
#                      xbar = 100, ybar = 50, color = "#4C77BB", show_bar = FALSE) {
#   if (is.null(ylim))
#     ylim <- if (sign == 1) c(0, max(y)) else c(-max(-y), 0)
#   if (is.null(xlim))
#     xlim <- c(min(x), max(x))
  
#   idx1 <- which.min(abs(x - xlim[1]))
#   idx2 <- which.min(abs(x - xlim[2]))
#   plot(x[idx1:idx2], y[idx1:idx2], type = "l", col = color,
#        xlim = xlim, ylim = ylim, bty = "n", lwd = lwd, lty = 1,
#        axes = FALSE, frame = FALSE, xlab = "", ylab = "")
  
#   if (show_bar) {
#     ybar_start <- min(ylim) + (max(ylim) - min(ylim)) / 20
#     x_start <- max(xlim) - xbar - 50
#     y_start <- ybar_start
#     x_end <- x_start + xbar
#     y_end <- y_start + ybar
#     segments(x_start, y_start, x_end, y_start, lwd = lwd, col = "black")
#     segments(x_start, y_start, x_start, y_end, lwd = lwd, col = "black")
#     if (show_text) {
#       text(x = (x_start + x_end) / 2, y = y_start - ybar / 20, 
#            labels = paste(xbar, "ms"), adj = c(0.5, 1))
#       text(x = x_start - xbar / 4, y = (y_start + y_end) / 2,  
#            labels = paste(ybar, "pA"), adj = c(0.5, 0.5), srt = 90)
#     }
#   }
# }

# # -----------------------------------------------------------------------------
# # Modified abf_averages function (adapted for UI use)
# # -----------------------------------------------------------------------------
# abf_averages <- function(abf_files = NULL, abf_path = NULL, traces2average = NULL,
#                          baseline = 100, stimulation_time = 350, ylim = NULL, 
#                          xlim = NULL, color = "#4C77BB", xbar = 100, ybar = 50, 
#                          width = 5.25, height = 2.75, save = FALSE) {
  
#   if (!is.null(traces2average) && length(traces2average) != length(abf_files)) {
#     stop("Error: if specified, the traces2average list must have the same length as the input abf_files list.")
#   }
  
#   baseline2zero <- function(y, dt, stimulation_time, baseline) {
#     idx1 <- (stimulation_time - baseline) / dt
#     idx2 <- baseline / dt
#     y1 <- y[idx1:length(y)]
#     y1 <- y1 - mean(y1[1:idx2])
#     return(y1 - mean(y1[1:idx2]))
#   }
  
#   abf_path <- if (is.null(abf_path)) getwd() else abf_path
#   setwd(abf_path)
#   N <- length(abf_files)
  
#   datasets <- lapply(1:N, function(ii) readABF(abf_files[ii]))
#   names(datasets) <- abf_files
  
#   metadata <- lapply(datasets, extract_metadata)
#   headers <- lapply(1:N, function(iii) datasets[[iii]][names(datasets[[iii]]) != "data"])
#   names(headers) <- abf_files
#   sampling_intervals <- sapply(1:N, function(ii) datasets[[ii]]$samplingIntervalInSec * 1000)
  
#   responses <- lapply(1:N, function(iii) {
#     sapply(1:length(datasets[[iii]]$data), function(ii) 
#       datasets[[iii]]$data[[ii]][, 1])
#   })
#   names(responses) <- abf_files
  
#   responses0 <- lapply(1:N, function(iii) {
#     sapply(1:dim(responses[[iii]])[2], function(ii)
#       baseline2zero(responses[[iii]][, ii],
#                     dt = sampling_intervals[iii],
#                     stimulation_time = stimulation_time,
#                     baseline = baseline))
#   })
#   names(responses0) <- abf_files
  
#   if (is.null(traces2average)) {
#     responses0_mean <- lapply(1:N, function(iii) apply(responses0[[iii]], 1, mean))
#   } else {
#     responses0_mean <- lapply(1:N, function(iii)
#       apply(responses0[[iii]][, traces2average[[iii]]], 1, mean))
#   }
  
#   time <- lapply(1:N, function(iii) {
#     dt_val <- datasets[[iii]]$samplingIntervalInSec * 1000
#     seq(0, by = dt_val, length.out = length(responses0_mean[[iii]]))
#   })
  
#   par(mfrow = c(1, N))
#   show_bar <- rep(FALSE, N)
#   show_bar[N] <- TRUE
#   for (ii in 1:N) {
#     egs_plot(x = time[[ii]], y = responses0_mean[[ii]], ylim = ylim,
#              xlim = xlim, color = "darkgrey", width = width, height = height,
#              show_text = FALSE, show_bar = show_bar[ii])
#     axis(1)
#     axis(2)
#   }
  
#   if (save) {
#     warning("save_graph not implemented in this UI example")
#   }
  
#   return(list("headers" = headers,
#               "raw_data" = responses,
#               "baseline_corrected_data" = responses0,
#               "baseline_corrected_mean_data" = responses0_mean,
#               "datasets" = datasets,
#               "metadata" = metadata))
# }

# -----------------------------------------------------------------------------
# Global variable to hold accepted trace indices (one element per dataset)
# -----------------------------------------------------------------------------
traces2average <<- list()

# -----------------------------------------------------------------------------
# Function to review recordings (trace-by-trace) for each dataset.
# -----------------------------------------------------------------------------
review_recordings <- function() {
  if (!exists("abf_analysis_result", envir = .GlobalEnv)) {
    tkmessageBox(message = "No analysis result available for review.")
    return()
  }
  result <- get("abf_analysis_result", envir = .GlobalEnv)
  datasets <- result$datasets
  # Initialize global list for accepted traces
  traces2average <<- vector("list", length = length(datasets))
  for(i in seq_along(datasets)) { traces2average[[i]] <<- integer(0) }
  
  current_dataset <<- 1
  current_trace <<- 1
  
  reviewWin <<- tktoplevel()
  tkwm.title(reviewWin, "Review Recordings")
  
  infoLabel <<- tklabel(reviewWin, text = paste("Dataset", current_dataset, "Trace", current_trace))
  tkgrid(infoLabel, row = 0, column = 0, columnspan = 2)
  
  reviewPlot <<- tkrplot(reviewWin, fun = function() {
    ds <- datasets[[current_dataset]]
    dt_val <- ds$samplingIntervalInSec * 1000
    time <- seq(0, by = dt_val, length.out = nrow(ds$data[[1]]))
    col_num <- as.numeric(tclvalue(dataColVar))
    if (is.na(col_num) || col_num < 1 || col_num > ncol(ds$data[[1]])) { col_num <- 1 }
    if (current_trace > ncol(ds$data[[1]])) {
      plot.new()
      text(0.5, 0.5, paste("No more traces in dataset", current_dataset))
    } else {
      trace <- ds$data[[1]][, current_trace]
      plot(time, trace, col = "indianred", xlab = "Time (ms)", type = "l",
           bty = "l", las = 1, main = paste("Dataset", current_dataset, "Trace", current_trace))
      axis(1)
      axis(2)
    }
  }, hscale = 1.5, vscale = 1.5)
  tkgrid(reviewPlot, row = 1, column = 0, columnspan = 2)
  
  acceptButton <- tkbutton(reviewWin, text = "Accept", command = function() {
    traces2average[[current_dataset]] <<- c(traces2average[[current_dataset]], current_trace)
    tkmessageBox(message = paste("Accepted dataset", current_dataset, "trace", current_trace))
  })
  tkgrid(acceptButton, row = 2, column = 0)
  
  rejectButton <- tkbutton(reviewWin, text = "Reject", command = function() {
    tkmessageBox(message = paste("Rejected dataset", current_dataset, "trace", current_trace))
  })
  tkgrid(rejectButton, row = 2, column = 1)
  
  nextTraceButton <- tkbutton(reviewWin, text = "Next Trace", command = function() {
    ds <- datasets[[current_dataset]]
    numTraces <- ncol(ds$data[[1]])
    if (current_trace < numTraces) {
      current_trace <<- current_trace + 1
    } else {
      tkmessageBox(message = paste("Finished reviewing dataset", current_dataset))
      if (current_dataset < length(datasets)) {
        current_dataset <<- current_dataset + 1
        current_trace <<- 1
      } else {
        tkmessageBox(message = "Review complete. Final accepted traces are stored in 'traces2average'.")
        tkdestroy(reviewWin)
        return()
      }
    }
    tkconfigure(infoLabel, text = paste("Dataset", current_dataset, "Trace", current_trace))
    tkrreplot(reviewPlot)
  })
  tkgrid(nextTraceButton, row = 3, column = 0, columnspan = 2)
}

# -----------------------------------------------------------------------------
# Function to average only the approved traces.
# This re-runs abf_averages with traces2average set to the approved indices.
# -----------------------------------------------------------------------------
averageApprovedTraces <- function() {
  folderPath <- tclvalue(folderPathVar)
  if (nchar(folderPath) == 0) {
    tkmessageBox(message = "Please select an ABF folder first.")
    return()
  }
  selIndices <- as.integer(tkcurselection(abfListBox))
  allFiles <- as.character(tkget(abfListBox, 0, "end"))
  abf_files <- if (length(selIndices) == 0) allFiles else allFiles[selIndices + 1]
  if (length(abf_files) == 0) {
    tkmessageBox(message = "No ABF files selected.")
    return()
  }
  baseline <- as.numeric(tclvalue(baselineVar))
  stimTime <- as.numeric(tclvalue(stimTimeVar))
  xbar     <- as.numeric(tclvalue(xbarVar))
  ybar     <- as.numeric(tclvalue(ybarVar))
  
  if (length(traces2average) != length(abf_files)) {
    tkmessageBox(message = "Number of datasets does not match the reviewed results.")
    return()
  }
  
  result <- tryCatch({
    abf_out <- abf_averages(abf_files = abf_files, abf_path = folderPath, 
                            traces2average = traces2average,
                            baseline = baseline, stimulation_time = stimTime, 
                            xbar = xbar, ybar = ybar)
    return(abf_out)
  }, error = function(e) {
    tkmessageBox(message = paste("Error during averaging of approved traces:", e$message))
    return(NULL)
  })
  
  if (!is.null(result)) {
    tkdelete(consoleText, "1.0", "end")
    tkinsert(consoleText, "end", paste("Averaging on approved traces complete. Processed", length(abf_files), "file(s)."))
    assign("abf_analysis_result", result, envir = .GlobalEnv)
    tkrreplot(plotWidget)
  }
}

# -----------------------------------------------------------------------------
# Build the ABF Analysis user interface using Tcl/Tk
# -----------------------------------------------------------------------------
ABF_analysis_tk <- function() {
  tt <- tktoplevel()
  tkwm.title(tt, "ABF Analysis")
  
  sidebarFrame <- tkframe(tt)
  mainFrame <- tkframe(tt)
  tkgrid(sidebarFrame, row = 0, column = 0, sticky = "ns")
  tkgrid(mainFrame, row = 0, column = 1, sticky = "nsew")
  tkgrid.rowconfigure(tt, 0, weight = 1)
  tkgrid.columnconfigure(tt, 1, weight = 1)
  
  # --- Sidebar Controls ---
  folderLabel <- tklabel(sidebarFrame, text = "Select ABF Folder:")
  tkgrid(folderLabel, row = 0, column = 0, sticky = "w")
  folderPathVar <<- tclVar("")
  folderEntry <- tkentry(sidebarFrame, textvariable = folderPathVar, width = 30)
  tkgrid(folderEntry, row = 0, column = 1, sticky = "w")
  
  browseFolderButton <- tkbutton(sidebarFrame, text = "Browse", command = function() {
    folderPath <- tclvalue(tkchooseDirectory())
    if (nchar(folderPath) > 0) {
      tclvalue(folderPathVar) <<- folderPath
      abf_list <- list.files(path = folderPath, pattern = "\\.abf$", ignore.case = TRUE)
      if (length(abf_list) == 0) {
        tkmessageBox(message = "No ABF files found in the selected folder.")
      } else {
        tkdelete(abfListBox, 0, "end")
        for (f in abf_list) { tkinsert(abfListBox, "end", f) }
        firstFilePath <- file.path(folderPath, abf_list[1])
        ds <- readABF(firstFilePath)
        dummy_result <- list(metadata = list(extract_metadata(ds)))
        updateAdditionalParams(dummy_result)
      }
    }
  })
  tkgrid(browseFolderButton, row = 0, column = 2, padx = 5)
  
  abfListLabel <- tklabel(sidebarFrame, text = "ABF Files:")
  tkgrid(abfListLabel, row = 1, column = 0, sticky = "w", pady = 5)
  abfListBox <<- tklistbox(sidebarFrame, height = 5, width = 40, selectmode = "multiple")
  tkgrid(abfListBox, row = 2, column = 0, columnspan = 3, pady = 5, sticky = "we")
  
  # --- Additional Parameters (above Baseline) ---
  paramFrame <- tkframe(sidebarFrame)
  tkgrid(paramFrame, row = 3, column = 0, columnspan = 3, sticky = "w")
  
  experimentVar <<- tclVar("Voltage Clamp")
  unitVar <<- tclVar("")
  dataColVar <<- tclVar("")
  dtVar <<- tclVar("")
  ntracesVar <<- tclVar("")
  
  tkgrid(tklabel(paramFrame, text = "Experiment:"), row = 0, column = 0, sticky = "w")
  experimentCombo <- ttkcombobox(paramFrame, textvariable = experimentVar, 
                                 values = c("Voltage Clamp", "Current Clamp"), width = 15)
  tkgrid(experimentCombo, row = 0, column = 1, sticky = "w")
  
  tkgrid(tklabel(paramFrame, text = "Units:"), row = 1, column = 0, sticky = "w")
  unitEntry <- tkentry(paramFrame, textvariable = unitVar, width = 10, state = "readonly")
  tkgrid(unitEntry, row = 1, column = 1, sticky = "w")
  
  tkgrid(tklabel(paramFrame, text = "Data Column:"), row = 2, column = 0, sticky = "w")
  dataColEntry <- tkentry(paramFrame, textvariable = dataColVar, width = 10, state = "readonly")
  tkgrid(dataColEntry, row = 2, column = 1, sticky = "w")
  
  tkgrid(tklabel(paramFrame, text = "dt (ms):"), row = 3, column = 0, sticky = "w")
  dtEntry <- tkentry(paramFrame, textvariable = dtVar, width = 10, state = "readonly")
  tkgrid(dtEntry, row = 3, column = 1, sticky = "w")
  
  tkgrid(tklabel(paramFrame, text = "Traces:"), row = 4, column = 0, sticky = "w")
  ntracesEntry <- tkentry(paramFrame, textvariable = ntracesVar, width = 10, state = "readonly")
  tkgrid(ntracesEntry, row = 4, column = 1, sticky = "w")
  
  # --- Baseline and Other Parameter Inputs ---
  baselineVar <<- tclVar("100")
  stimTimeVar <<- tclVar("350")
  xbarVar <<- tclVar("100")
  ybarVar <<- tclVar("50")
  
  tkgrid(tklabel(sidebarFrame, text = "Baseline:"), row = 4, column = 0, sticky = "w")
  tkgrid(tkentry(sidebarFrame, textvariable = baselineVar, width = 10), row = 4, column = 1, sticky = "w")
  tkgrid(tklabel(sidebarFrame, text = "Stimulation Time:"), row = 5, column = 0, sticky = "w")
  tkgrid(tkentry(sidebarFrame, textvariable = stimTimeVar, width = 10), row = 5, column = 1, sticky = "w")
  tkgrid(tklabel(sidebarFrame, text = "x-bar length:"), row = 6, column = 0, sticky = "w")
  tkgrid(tkentry(sidebarFrame, textvariable = xbarVar, width = 10), row = 6, column = 1, sticky = "w")
  tkgrid(tklabel(sidebarFrame, text = "y-bar length:"), row = 7, column = 0, sticky = "w")
  tkgrid(tkentry(sidebarFrame, textvariable = ybarVar, width = 10), row = 7, column = 1, sticky = "w")
  
  consoleText <<- tktext(sidebarFrame, width = 40, height = 4)
  tkgrid(consoleText, row = 8, column = 0, columnspan = 3, pady = 5)
  
  # --- Update Function ---
  updateAdditionalParams <<- function(result) {
    if (!is.null(result) && length(result$metadata) >= 1) {
      meta1 <- result$metadata[[1]]
      tclvalue(dtVar) <<- as.character(meta1$samplingIntervalInSec * 1000)
      if (!is.null(meta1$header$lActualEpisodes)) {
        tclvalue(ntracesVar) <<- as.character(meta1$header$lActualEpisodes)
      } else {
        tclvalue(ntracesVar) <<- "N/A"
      }
      expType <- tclvalue(experimentVar)
      col_idx <- choose_data_column(meta1$channelUnits, expType)
      if (!is.na(col_idx)) {
        tclvalue(unitVar) <<- meta1$channelUnits[col_idx]
        tclvalue(dataColVar) <<- as.character(col_idx)
      } else {
        tclvalue(unitVar) <<- "N/A"
        tclvalue(dataColVar) <<- "N/A"
      }
    }
  }
  
  tkbind(experimentCombo, "<<ComboboxSelected>>", function() {
    if (exists("abf_analysis_result", envir = .GlobalEnv)) {
      result <- get("abf_analysis_result", envir = .GlobalEnv)
      updateAdditionalParams(result)
    }
  })
  
  # --- Analysis Action: Run Analysis button (which becomes "Load Data")
  runAnalysis <<- function() {
    folderPath <- tclvalue(folderPathVar)
    if (nchar(folderPath) == 0) {
      tkmessageBox(message = "Please select an ABF folder first.")
      return()
    }
    selIndices <- as.integer(tkcurselection(abfListBox))
    allFiles <- as.character(tkget(abfListBox, 0, "end"))
    abf_files <- if (length(selIndices) == 0) allFiles else allFiles[selIndices + 1]
    if (length(abf_files) == 0) {
      tkmessageBox(message = "No ABF files selected.")
      return()
    }
    baseline <- as.numeric(tclvalue(baselineVar))
    stimTime <- as.numeric(tclvalue(stimTimeVar))
    xbar <- as.numeric(tclvalue(xbarVar))
    ybar <- as.numeric(tclvalue(ybarVar))
    
    result <- tryCatch({
      abf_out <- abf_averages(abf_files = abf_files, abf_path = folderPath, 
                              baseline = baseline, stimulation_time = stimTime, 
                              xbar = xbar, ybar = ybar)
      return(abf_out)
    }, error = function(e) {
      tkmessageBox(message = paste("Error during analysis:", e$message))
      return(NULL)
    })
    
    if (!is.null(result)) {
      tkdelete(consoleText, "1.0", "end")
      tkinsert(consoleText, "end", paste("Analysis complete. Processed", length(abf_files), "file(s)."))
      assign("abf_analysis_result", result, envir = .GlobalEnv)
      updateAdditionalParams(result)
      cons_msg <- check_consistency(result$metadata)
      if (cons_msg != "GOOD") {
        tkmessageBox(message = paste("ERROR:", cons_msg))
      } else {
        tkmessageBox(message = "GOOD: All parameters are consistent.")
      }
      tkconfigure(runAnalysisButton, text = "Load Data")
      tkgrid.forget(plotWidget)
    }
  }
  
  runAnalysisButton <<- tkbutton(sidebarFrame, text = "Run Analysis", command = runAnalysis)
  tkgrid(runAnalysisButton, row = 9, column = 0, columnspan = 3, pady = 5)
  
  # --- New Button: Average Approved Traces ---
  averageApprovedTraces <<- function() {
    if (length(traces2average) == 0 || all(sapply(traces2average, length) == 0)) {
      tkmessageBox(message = "No approved traces available. Please review recordings first.")
      return()
    }
    folderPath <- tclvalue(folderPathVar)
    if (nchar(folderPath) == 0) {
      tkmessageBox(message = "Please select an ABF folder first.")
      return()
    }
    selIndices <- as.integer(tkcurselection(abfListBox))
    allFiles <- as.character(tkget(abfListBox, 0, "end"))
    abf_files <- if (length(selIndices) == 0) allFiles else allFiles[selIndices + 1]
    if (length(abf_files) == 0) {
      tkmessageBox(message = "No ABF files selected.")
      return()
    }
    baseline <- as.numeric(tclvalue(baselineVar))
    stimTime <- as.numeric(tclvalue(stimTimeVar))
    xbar <- as.numeric(tclvalue(xbarVar))
    ybar <- as.numeric(tclvalue(ybarVar))
    
    result <- tryCatch({
      abf_out <- abf_averages(abf_files = abf_files, abf_path = folderPath,
                              traces2average = traces2average,
                              baseline = baseline, stimulation_time = stimTime,
                              xbar = xbar, ybar = ybar)
      return(abf_out)
    }, error = function(e) {
      tkmessageBox(message = paste("Error during averaging of approved traces:", e$message))
      return(NULL)
    })
    
    if (!is.null(result)) {
      tkdelete(consoleText, "1.0", "end")
      tkinsert(consoleText, "end", paste("Averaging on approved traces complete. Processed", length(abf_files), "file(s)."))
      assign("abf_analysis_result", result, envir = .GlobalEnv)
      tkrreplot(plotWidget)
    }
  }
  
  averageApprovedButton <<- tkbutton(sidebarFrame, text = "Average Approved Traces", command = averageApprovedTraces)
  tkgrid(averageApprovedButton, row = 10, column = 0, columnspan = 3, pady = 5)
  
  # --- New Button: Review Recordings ---
  reviewButton <<- tkbutton(sidebarFrame, text = "Review Recordings", command = review_recordings)
  tkgrid(reviewButton, row = 11, column = 0, columnspan = 3, pady = 5)
  
  # --- Main Panel: Plot display via tkrplot ---
  drawPlot <- function() {
    if (exists("abf_analysis_result", envir = .GlobalEnv)) {
      result <- get("abf_analysis_result", envir = .GlobalEnv)
      datasets <- result$datasets
      if (length(datasets) >= 1) {
        iii <- 1
        dt_val <- datasets[[iii]]$samplingIntervalInSec * 1000
        time <- seq(0, by = dt_val, length.out = nrow(datasets[[iii]]$data[[1]]))
        col_num <- as.numeric(tclvalue(dataColVar))
        if (is.na(col_num) || col_num < 1 || col_num > ncol(datasets[[iii]]$data[[1]])) {
          col_num <- 1
        }
        trace <- datasets[[iii]]$data[[1]][, col_num]
        plot(time, trace, col = "indianred", xlab = "Time (ms)", type = "l",
             bty = "l", las = 1, main = paste("Trace", iii))
        axis(1)
        axis(2)
      } else {
        plot.new()
        text(0.5, 0.5, "No data available")
      }
    } else {
      plot.new()
      text(0.5, 0.5, "No analysis result to display")
    }
  }
  
  plotWidget <<- tkrplot(mainFrame, fun = drawPlot, hscale = 1.5, vscale = 1.5)
  tkgrid(plotWidget, row = 0, column = 0, sticky = "nsew")
  tkfocus(tt)
}

# -----------------------------------------------------------------------------
# Launch the ABF Analysis interface:
# -----------------------------------------------------------------------------
ABF_analysis_tk()



# # Remove all objects from the environment
# rm(list = ls(all = TRUE))

# # -----------------------------------------------------------------------------
# # Package loading function and required packages:
# # -----------------------------------------------------------------------------
# load_required_packages <- function(packages) {
#   new.packages <- packages[!(packages %in% installed.packages()[, "Package"])]
#   if (length(new.packages)) install.packages(new.packages)
#   invisible(lapply(packages, library, character.only = TRUE))
# }

# required.packages <- c("robustbase", "minpack.lm", "Rcpp", "signal", "readABF", "tcltk", "tkrplot")
# load_required_packages(required.packages)

# # -----------------------------------------------------------------------------
# # Helper functions
# # -----------------------------------------------------------------------------

# # Function to extract key metadata from an ABF dataset
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

# # Function to choose data column index based on experiment type.
# # For Voltage Clamp, choose the first channel whose unit contains "A" (e.g., pA, nA).
# # For Current Clamp, choose the first channel whose unit contains "V" (e.g., mV, nV).
# choose_data_column <- function(channelUnits, experiment) {
#   if (experiment == "Voltage Clamp") {
#     idx <- grep("A", channelUnits, ignore.case = TRUE)
#   } else if (experiment == "Current Clamp") {
#     idx <- grep("V", channelUnits, ignore.case = TRUE)
#   } else {
#     idx <- integer(0)
#   }
#   if (length(idx) > 0) return(idx[1])
#   else return(NA)
# }

# # -----------------------------------------------------------------------------
# # Custom plotting function (egs_plot)
# # -----------------------------------------------------------------------------
# egs_plot <- function(x, y, sign = -1, xlim = NULL, ylim = NULL, lwd = 1, 
#                      show_text = FALSE, height = 4, width = 2.5, 
#                      xbar = 100, ybar = 50,  color = "#4C77BB", show_bar = FALSE) {
  
#   if (is.null(ylim))
#     ylim <- if (sign == 1) c(0, max(y)) else c(-max(-y), 0)
#   if (is.null(xlim))
#     xlim <- c(min(x), max(x))
  
#   idx1 <- which.min(abs(x - xlim[1]))
#   idx2 <- which.min(abs(x - xlim[2]))
  
#   plot(x[idx1:idx2], y[idx1:idx2], type = "l", col = color,
#        xlim = xlim, ylim = ylim, bty = "n", lwd = lwd, lty = 1,
#        axes = FALSE, frame = FALSE, xlab = "", ylab = "")
  
#   if (show_bar) {
#     # Calculate positions for scale bar and labels
#     ybar_start <- min(ylim) + (max(ylim) - min(ylim)) / 20
#     x_start <- max(xlim) - xbar - 50
#     y_start <- ybar_start
#     x_end <- x_start + xbar
#     y_end <- y_start + ybar
    
#     segments(x_start, y_start, x_end, y_start, lwd = lwd, col = "black")
#     segments(x_start, y_start, x_start, y_end, lwd = lwd, col = "black")
    
#     if (show_text) {
#       text(x = (x_start + x_end) / 2, y = y_start - ybar / 20, 
#            labels = paste(xbar, "ms"), adj = c(0.5, 1))
#       text(x = x_start - xbar / 4, y = (y_start + y_end) / 2,  
#            labels = paste(ybar, "pA"), adj = c(0.5, 0.5), srt = 90)
#     }
#   }
# }

# # -----------------------------------------------------------------------------
# # Modified abf_averages function (adapted for UI use)
# # -----------------------------------------------------------------------------
# abf_averages <- function(abf_files = NULL, abf_path = NULL, traces2average = NULL,
#                          baseline = 100, stimulation_time = 350, ylim = NULL, 
#                          xlim = NULL, color = "#4C77BB", xbar = 100, ybar = 50, 
#                          width = 5.25, height = 2.75, save = FALSE) {
  
#   if (!is.null(traces2average) && length(traces2average) != length(abf_files)) {
#     stop("Error: if specified, the traces2average list must have the same length as the input abf_files list.")
#   }
  
#   baseline2zero <- function(y, dt, stimulation_time, baseline) {
#     idx1 <- (stimulation_time - baseline) / dt
#     idx2 <- baseline / dt
#     y1 <- y[idx1:length(y)]
#     y1 <- y1 - mean(y1[1:idx2])
#     return(y1 - mean(y1[1:idx2]))
#   }
  
#   # Set working folder
#   abf_path <- if (is.null(abf_path)) getwd() else abf_path
#   setwd(abf_path)
  
#   N <- length(abf_files)
  
#   # Read each ABF file into a list of datasets
#   datasets <- lapply(1:N, function(ii) readABF(abf_files[ii]))
#   names(datasets) <- abf_files
  
#   # Extract metadata from each dataset
#   metadata <- lapply(datasets, extract_metadata)
  
#   headers <- lapply(1:N, function(iii) datasets[[iii]][names(datasets[[iii]]) != "data"])
#   names(headers) <- abf_files
  
#   sampling_intervals <- sapply(1:N, function(ii) datasets[[ii]]$samplingIntervalInSec * 1000)
  
#   responses <- lapply(1:N, function(iii) {
#     sapply(1:length(datasets[[iii]]$data), function(ii) 
#       datasets[[iii]]$data[[ii]][, 1])
#   })
#   names(responses) <- abf_files
  
#   responses0 <- lapply(1:N, function(iii) {
#     sapply(1:dim(responses[[iii]])[2], function(ii) 
#       baseline2zero(responses[[iii]][, ii], dt = sampling_intervals[iii], 
#                     stimulation_time = stimulation_time, baseline = baseline))
#   })
#   names(responses0) <- abf_files
  
#   if (is.null(traces2average)) {
#     responses0_mean <- lapply(1:N, function(iii) apply(responses0[[iii]], 1, mean))
#   } else {
#     responses0_mean <- lapply(1:N, function(iii) 
#       apply(responses0[[iii]][, traces2average[[iii]]], 1, mean))
#   }
  
#   time <- lapply(1:N, function(iii) {
#     dt_val <- datasets[[iii]]$samplingIntervalInSec * 1000
#     seq(0, by = dt_val, length.out = length(responses0_mean[[iii]]))
#   })
  
#   # Plot onto the current device (e.g., the tkrplot window)
#   par(mfrow = c(1, N))
#   show_bar <- rep(FALSE, N)
#   show_bar[N] <- TRUE
  
#   for (ii in 1:N) {
#     egs_plot(x = time[[ii]], y = responses0_mean[[ii]], ylim = ylim, 
#              xlim = xlim, color = "darkgrey", width = width, height = height, 
#              show_text = FALSE, show_bar = show_bar[ii])
#     axis(1)
#     axis(2)
#   }
  
#   if (save) {
#     warning("save_graph not implemented in this UI example")
#   }
  
#   return(list("headers" = headers, 
#               "raw_data" = responses, 
#               "baseline_corrected_data" = responses0, 
#               "baseline_corrected_mean_data" = responses0_mean,
#               "datasets" = datasets,
#               "metadata" = metadata))
# }

# # -----------------------------------------------------------------------------
# # Build the ABF Analysis user interface using Tcl/Tk
# # -----------------------------------------------------------------------------
# ABF_analysis_tk <- function() {
#   tt <- tktoplevel()
#   tkwm.title(tt, "ABF Analysis")
  
#   # Divide the window into sidebar (controls) and main (plot) panels
#   sidebarFrame <- tkframe(tt)
#   mainFrame <- tkframe(tt)
#   tkgrid(sidebarFrame, row = 0, column = 0, sticky = "ns")
#   tkgrid(mainFrame, row = 0, column = 1, sticky = "nsew")
#   tkgrid.rowconfigure(tt, 0, weight = 1)
#   tkgrid.columnconfigure(tt, 1, weight = 1)
  
#   # --- Sidebar Controls ---
  
#   # Folder selection for ABF files
#   folderLabel <- tklabel(sidebarFrame, text = "Select ABF Folder:")
#   tkgrid(folderLabel, row = 0, column = 0, sticky = "w")
#   folderPathVar <- tclVar("")
#   folderEntry <- tkentry(sidebarFrame, textvariable = folderPathVar, width = 30)
#   tkgrid(folderEntry, row = 0, column = 1, sticky = "w")
  
#   browseFolderButton <- tkbutton(sidebarFrame, text = "Browse", command = function() {
#     folderPath <- tclvalue(tkchooseDirectory())
#     if (nchar(folderPath) > 0) {
#       tclvalue(folderPathVar) <- folderPath
#       # List available ABF files in the folder:
#       abf_list <- list.files(path = folderPath, pattern = "\\.abf$", ignore.case = TRUE)
#       if (length(abf_list) == 0) {
#         tkmessageBox(message = "No ABF files found in the selected folder.")
#       } else {
#         tkdelete(abfListBox, 0, "end")
#         for (f in abf_list) {
#           tkinsert(abfListBox, "end", f)
#         }
#         # Auto-load metadata from the first ABF file to update additional parameters:
#         firstFilePath <- file.path(folderPath, abf_list[1])
#         ds <- readABF(firstFilePath)
#         dummy_result <- list(metadata = list(extract_metadata(ds)))
#         updateAdditionalParams(dummy_result)
#       }
#     }
#   })
#   tkgrid(browseFolderButton, row = 0, column = 2, padx = 5)
  
#   # Listbox to show available ABF files
#   abfListLabel <- tklabel(sidebarFrame, text = "ABF Files:")
#   tkgrid(abfListLabel, row = 1, column = 0, sticky = "w", pady = 5)
#   abfListBox <- tklistbox(sidebarFrame, height = 5, width = 40, selectmode = "multiple")
#   tkgrid(abfListBox, row = 2, column = 0, columnspan = 3, pady = 5, sticky = "we")
  
#   # --- Additional Parameters (above Baseline) ---
#   # Create a new frame in sidebarFrame for experiment parameters
#   paramFrame <- tkframe(sidebarFrame)
#   tkgrid(paramFrame, row = 3, column = 0, columnspan = 3, sticky = "w")
  
#   # Create tcl variables for additional parameters:
#   experimentVar <- tclVar("Voltage Clamp")
#   unitVar       <- tclVar("")
#   dataColVar    <- tclVar("")
#   dtVar         <- tclVar("")
#   ntracesVar    <- tclVar("")
  
#   # Experiment type selector (Voltage Clamp vs. Current Clamp)
#   tkgrid(tklabel(paramFrame, text = "Experiment:"), row = 0, column = 0, sticky = "w")
#   experimentCombo <- ttkcombobox(paramFrame, textvariable = experimentVar, 
#                                  values = c("Voltage Clamp", "Current Clamp"), width = 15)
#   tkgrid(experimentCombo, row = 0, column = 1, sticky = "w")
  
#   # Display Units (auto-populated)
#   tkgrid(tklabel(paramFrame, text = "Units:"), row = 1, column = 0, sticky = "w")
#   unitEntry <- tkentry(paramFrame, textvariable = unitVar, width = 10, state = "readonly")
#   tkgrid(unitEntry, row = 1, column = 1, sticky = "w")
  
#   # Data column (auto-populated)
#   tkgrid(tklabel(paramFrame, text = "Data Column:"), row = 2, column = 0, sticky = "w")
#   dataColEntry <- tkentry(paramFrame, textvariable = dataColVar, width = 10, state = "readonly")
#   tkgrid(dataColEntry, row = 2, column = 1, sticky = "w")
  
#   # dt (sampling interval in ms)
#   tkgrid(tklabel(paramFrame, text = "dt (ms):"), row = 3, column = 0, sticky = "w")
#   dtEntry <- tkentry(paramFrame, textvariable = dtVar, width = 10, state = "readonly")
#   tkgrid(dtEntry, row = 3, column = 1, sticky = "w")
  
#   # Number of Traces (from header$lActualEpisodes)
#   tkgrid(tklabel(paramFrame, text = "Traces:"), row = 4, column = 0, sticky = "w")
#   ntracesEntry <- tkentry(paramFrame, textvariable = ntracesVar, width = 10, state = "readonly")
#   tkgrid(ntracesEntry, row = 4, column = 1, sticky = "w")
  
#   # --- Baseline and Other Parameter Inputs ---
#   baselineVar <- tclVar("100")
#   stimTimeVar <- tclVar("350")
#   xbarVar <- tclVar("100")
#   ybarVar <- tclVar("50")
  
#   tkgrid(tklabel(sidebarFrame, text = "Baseline:"), row = 4, column = 0, sticky = "w")
#   tkgrid(tkentry(sidebarFrame, textvariable = baselineVar, width = 10), row = 4, column = 1, sticky = "w")
#   tkgrid(tklabel(sidebarFrame, text = "Stimulation Time:"), row = 5, column = 0, sticky = "w")
#   tkgrid(tkentry(sidebarFrame, textvariable = stimTimeVar, width = 10), row = 5, column = 1, sticky = "w")
#   tkgrid(tklabel(sidebarFrame, text = "x-bar length:"), row = 6, column = 0, sticky = "w")
#   tkgrid(tkentry(sidebarFrame, textvariable = xbarVar, width = 10), row = 6, column = 1, sticky = "w")
#   tkgrid(tklabel(sidebarFrame, text = "y-bar length:"), row = 7, column = 0, sticky = "w")
#   tkgrid(tkentry(sidebarFrame, textvariable = ybarVar, width = 10), row = 7, column = 1, sticky = "w")
  
#   # Text widget for status messages (console)
#   consoleText <- tktext(sidebarFrame, width = 40, height = 4)
#   tkgrid(consoleText, row = 8, column = 0, columnspan = 3, pady = 5)
  
#   # --- Update Function ---
#   # Update Units, Data Column, dt, and Ntraces based on analysis result and experiment type.
#   updateAdditionalParams <- function(result) {
#     if (!is.null(result) && length(result$metadata) >= 1) {
#       meta1 <- result$metadata[[1]]
#       # Update dt (sampling interval in ms)
#       tclvalue(dtVar) <- as.character(meta1$samplingIntervalInSec * 1000)
#       # Update number of traces from header information
#       if (!is.null(meta1$header$lActualEpisodes)) {
#         tclvalue(ntracesVar) <- as.character(meta1$header$lActualEpisodes)
#       } else {
#         tclvalue(ntracesVar) <- "N/A"
#       }
#       # Choose column index based on experiment type:
#       expType <- tclvalue(experimentVar)
#       col_idx <- choose_data_column(meta1$channelUnits, expType)
#       if (!is.na(col_idx)) {
#         tclvalue(unitVar) <- meta1$channelUnits[col_idx]
#         tclvalue(dataColVar) <- as.character(col_idx)
#       } else {
#         tclvalue(unitVar) <- "N/A"
#         tclvalue(dataColVar) <- "N/A"
#       }
#     }
#   }
  
#   # Function to check consistency across all ABF files.
#   # It verifies that dt, the unit for the chosen column, and the number of traces are identical.
#   check_consistency <- function(metadata) {
#     dt_values <- sapply(metadata, function(meta) meta$samplingIntervalInSec * 1000)
#     traces_values <- sapply(metadata, function(meta) meta$header$lActualEpisodes)
#     expType <- tclvalue(experimentVar)
#     unit_values <- sapply(metadata, function(meta) {
#       col_idx <- choose_data_column(meta$channelUnits, expType)
#       if (!is.na(col_idx)) meta$channelUnits[col_idx] else NA_character_
#     })
#     dt_good <- (length(unique(dt_values)) == 1)
#     traces_good <- (length(unique(traces_values)) == 1)
#     unit_good <- (length(unique(unit_values)) == 1)
#     if (dt_good && traces_good && unit_good) {
#       return("GOOD")
#     } else {
#       error_msgs <- c()
#       if (!dt_good) error_msgs <- c(error_msgs, paste("Inconsistent dt values:", paste(dt_values, collapse=", ")))
#       if (!unit_good) error_msgs <- c(error_msgs, paste("Inconsistent Units:", paste(unit_values, collapse=", ")))
#       if (!traces_good) error_msgs <- c(error_msgs, paste("Inconsistent Traces:", paste(traces_values, collapse=", ")))
#       return(paste(error_msgs, collapse="; "))
#     }
#   }
  
#   # Bind event when an item is selected from the experiment combobox
#   tkbind(experimentCombo, "<<ComboboxSelected>>", function() {
#     if (exists("abf_analysis_result", envir = .GlobalEnv)) {
#       result <- get("abf_analysis_result", envir = .GlobalEnv)
#       updateAdditionalParams(result)
#     }
#   })
  
#   # --- Analysis Action: Run Analysis button ---
#   runAnalysis <- function() {
#     folderPath <- tclvalue(folderPathVar)
#     if (nchar(folderPath) == 0) {
#       tkmessageBox(message = "Please select an ABF folder first.")
#       return()
#     }
    
#     selIndices <- as.integer(tkcurselection(abfListBox))
#     allFiles <- as.character(tkget(abfListBox, 0, "end"))
#     abf_files <- if (length(selIndices) == 0) allFiles else allFiles[selIndices + 1]
    
#     if (length(abf_files) == 0) {
#       tkmessageBox(message = "No ABF files selected.")
#       return()
#     }
    
#     baseline <- as.numeric(tclvalue(baselineVar))
#     stimTime <- as.numeric(tclvalue(stimTimeVar))
#     xbar     <- as.numeric(tclvalue(xbarVar))
#     ybar     <- as.numeric(tclvalue(ybarVar))
    
#     result <- tryCatch({
#       abf_out <- abf_averages(abf_files = abf_files, abf_path = folderPath, 
#                               baseline = baseline, stimulation_time = stimTime, 
#                               xbar = xbar, ybar = ybar)
#       return(abf_out)
#     }, error = function(e) {
#       tkmessageBox(message = paste("Error during analysis:", e$message))
#       return(NULL)
#     })
    
#     if (!is.null(result)) {
#       tkdelete(consoleText, "1.0", "end")
#       tkinsert(consoleText, "end", paste("Analysis complete. Processed", length(abf_files), "file(s)."))
#       assign("abf_analysis_result", result, envir = .GlobalEnv)
#       updateAdditionalParams(result)
#       # Check consistency across loaded ABF files
#       consistency_msg <- check_consistency(result$metadata)
#       if (consistency_msg != "GOOD") {
#         tkmessageBox(message = paste("ERROR:", consistency_msg))
#       } else {
#         tkmessageBox(message = "GOOD: All parameters are consistent.")
#       }
#       tkrreplot(plotWidget)
#     }
#   }
  
#   runAnalysisButton <- tkbutton(sidebarFrame, text = "Run Analysis", command = runAnalysis)
#   tkgrid(runAnalysisButton, row = 9, column = 0, columnspan = 3, pady = 5)
  
#   # --- Main Panel: Plot display via tkrplot ---
#   drawPlot <- function() {
#     if (exists("abf_analysis_result", envir = .GlobalEnv)) {
#       result <- get("abf_analysis_result", envir = .GlobalEnv)
#       datasets <- result$datasets
#       if (length(datasets) >= 1) {
#         iii <- 1
#         dt_val <- datasets[[iii]]$samplingIntervalInSec * 1000
#         time <- seq(0, by = dt_val, length.out = nrow(datasets[[iii]]$data[[1]]))
#         col_num <- as.numeric(tclvalue(dataColVar))
#         if (is.na(col_num) || col_num < 1 || col_num > ncol(datasets[[iii]]$data[[1]])) {
#           col_num <- 1
#         }
#         trace <- datasets[[iii]]$data[[1]][, col_num]
#         plot(time, trace, col = "indianred", xlab = "Time (ms)", type = "l",
#              bty = "l", las = 1, main = paste("Trace", iii))
#         axis(1)
#         axis(2)
#       } else {
#         plot.new()
#         text(0.5, 0.5, "No data available")
#       }
#     } else {
#       plot.new()
#       text(0.5, 0.5, "No analysis result to display")
#     }
#   }
  
#   # Create the plot widget using tkrplot() (not tkrreplot())
#   plotWidget <- tkrplot(mainFrame, fun = drawPlot, hscale = 1.5, vscale = 1.5)
#   tkgrid(plotWidget, row = 0, column = 0, sticky = "nsew")
  
#   tkfocus(tt)
# }

# # -----------------------------------------------------------------------------
# # Launch the ABF Analysis interface:
# # -----------------------------------------------------------------------------
# ABF_analysis_tk()





# # Remove all objects from the environment
# rm(list = ls(all = TRUE))

# # -----------------------------------------------------------------------------
# # Package loading function and required packages:
# # -----------------------------------------------------------------------------
# load_required_packages <- function(packages) {
#   new.packages <- packages[!(packages %in% installed.packages()[, "Package"])]
#   if (length(new.packages)) install.packages(new.packages)
#   invisible(lapply(packages, library, character.only = TRUE))
# }

# required.packages <- c("robustbase", "minpack.lm", "Rcpp", "signal", "readABF", "tcltk", "tkrplot")
# load_required_packages(required.packages)

# # -----------------------------------------------------------------------------
# # Helper functions
# # -----------------------------------------------------------------------------

# # Function to extract key metadata from an ABF dataset
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

# # Function to choose data column index based on experiment type.
# # For Voltage Clamp, choose the first channel whose unit contains "A" (e.g., pA, nA).
# # For Current Clamp, choose the first channel whose unit contains "V" (e.g., mV, nV).
# choose_data_column <- function(channelUnits, experiment) {
#   if (experiment == "Voltage Clamp") {
#     idx <- grep("A", channelUnits, ignore.case = TRUE)
#   } else if (experiment == "Current Clamp") {
#     idx <- grep("V", channelUnits, ignore.case = TRUE)
#   } else {
#     idx <- integer(0)
#   }
#   if (length(idx) > 0) return(idx[1])
#   else return(NA)
# }

# # -----------------------------------------------------------------------------
# # Custom plotting function (egs_plot)
# # -----------------------------------------------------------------------------
# egs_plot <- function(x, y, sign = -1, xlim = NULL, ylim = NULL, lwd = 1, 
#                      show_text = FALSE, height = 4, width = 2.5, 
#                      xbar = 100, ybar = 50,  color = "#4C77BB", show_bar = FALSE) {
  
#   if (is.null(ylim))
#     ylim <- if (sign == 1) c(0, max(y)) else c(-max(-y), 0)
#   if (is.null(xlim))
#     xlim <- c(min(x), max(x))
  
#   idx1 <- which.min(abs(x - xlim[1]))
#   idx2 <- which.min(abs(x - xlim[2]))
  
#   plot(x[idx1:idx2], y[idx1:idx2], type = "l", col = color,
#        xlim = xlim, ylim = ylim, bty = "n", lwd = lwd, lty = 1,
#        axes = FALSE, frame = FALSE, xlab = "", ylab = "")
  
#   if (show_bar) {
#     # Calculate positions for scale bar and labels
#     ybar_start <- min(ylim) + (max(ylim) - min(ylim)) / 20
#     x_start <- max(xlim) - xbar - 50
#     y_start <- ybar_start
#     x_end <- x_start + xbar
#     y_end <- y_start + ybar
    
#     segments(x_start, y_start, x_end, y_start, lwd = lwd, col = "black")
#     segments(x_start, y_start, x_start, y_end, lwd = lwd, col = "black")
    
#     if (show_text) {
#       text(x = (x_start + x_end) / 2, y = y_start - ybar / 20, 
#            labels = paste(xbar, "ms"), adj = c(0.5, 1))
#       text(x = x_start - xbar / 4, y = (y_start + y_end) / 2,  
#            labels = paste(ybar, "pA"), adj = c(0.5, 0.5), srt = 90)
#     }
#   }
# }

# # -----------------------------------------------------------------------------
# # Modified abf_averages function (adapted for UI use)
# # -----------------------------------------------------------------------------
# abf_averages <- function(abf_files = NULL, abf_path = NULL, traces2average = NULL,
#                          baseline = 100, stimulation_time = 350, ylim = NULL, 
#                          xlim = NULL, color = "#4C77BB", xbar = 100, ybar = 50, 
#                          width = 5.25, height = 2.75, save = FALSE) {
  
#   if (!is.null(traces2average) && length(traces2average) != length(abf_files)) {
#     stop("Error: if specified, the traces2average list must have the same length as the input abf_files list.")
#   }
  
#   baseline2zero <- function(y, dt, stimulation_time, baseline) {
#     idx1 <- (stimulation_time - baseline) / dt
#     idx2 <- baseline / dt
#     y1 <- y[idx1:length(y)]
#     y1 <- y1 - mean(y1[1:idx2])
#     return(y1 - mean(y1[1:idx2]))
#   }
  
#   # Set working folder
#   abf_path <- if (is.null(abf_path)) getwd() else abf_path
#   setwd(abf_path)
  
#   N <- length(abf_files)
  
#   # Read each ABF file into a list of datasets
#   datasets <- lapply(1:N, function(ii) readABF(abf_files[ii]))
#   names(datasets) <- abf_files
  
#   # Extract metadata from each dataset
#   metadata <- lapply(datasets, extract_metadata)
  
#   headers <- lapply(1:N, function(iii) datasets[[iii]][names(datasets[[iii]]) != "data"])
#   names(headers) <- abf_files
  
#   sampling_intervals <- sapply(1:N, function(ii) datasets[[ii]]$samplingIntervalInSec * 1000)
  
#   responses <- lapply(1:N, function(iii) {
#     sapply(1:length(datasets[[iii]]$data), function(ii) 
#       datasets[[iii]]$data[[ii]][, 1])
#   })
#   names(responses) <- abf_files
  
#   responses0 <- lapply(1:N, function(iii) {
#     sapply(1:dim(responses[[iii]])[2], function(ii) 
#       baseline2zero(responses[[iii]][, ii], dt = sampling_intervals[iii], 
#                     stimulation_time = stimulation_time, baseline = baseline))
#   })
#   names(responses0) <- abf_files
  
#   if (is.null(traces2average)) {
#     responses0_mean <- lapply(1:N, function(iii) apply(responses0[[iii]], 1, mean))
#   } else {
#     responses0_mean <- lapply(1:N, function(iii) 
#       apply(responses0[[iii]][, traces2average[[iii]]], 1, mean))
#   }
  
#   time <- lapply(1:N, function(iii) {
#     dt_val <- datasets[[iii]]$samplingIntervalInSec * 1000
#     seq(0, by = dt_val, length.out = length(responses0_mean[[iii]]))
#   })
  
#   # Plot on the current device (e.g., our tkrplot window)
#   par(mfrow = c(1, N))
#   show_bar <- rep(FALSE, N)
#   show_bar[N] <- TRUE
  
#   for (ii in 1:N) {
#     egs_plot(x = time[[ii]], y = responses0_mean[[ii]], ylim = ylim, 
#              xlim = xlim, color = "darkgrey", width = width, height = height, 
#              show_text = FALSE, show_bar = show_bar[ii])
#     axis(1)
#     axis(2)
#   }
  
#   if (save) {
#     warning("save_graph not implemented in this UI example")
#   }
  
#   return(list("headers" = headers, 
#               "raw_data" = responses, 
#               "baseline_corrected_data" = responses0, 
#               "baseline_corrected_mean_data" = responses0_mean,
#               "datasets" = datasets,
#               "metadata" = metadata))
# }

# # -----------------------------------------------------------------------------
# # Build the ABF Analysis user interface using Tcl/Tk
# # -----------------------------------------------------------------------------
# ABF_analysis_tk <- function() {
#   tt <- tktoplevel()
#   tkwm.title(tt, "ABF Analysis")
  
#   # Divide the window into sidebar (controls) and main (plot) panels
#   sidebarFrame <- tkframe(tt)
#   mainFrame <- tkframe(tt)
#   tkgrid(sidebarFrame, row = 0, column = 0, sticky = "ns")
#   tkgrid(mainFrame, row = 0, column = 1, sticky = "nsew")
#   tkgrid.rowconfigure(tt, 0, weight = 1)
#   tkgrid.columnconfigure(tt, 1, weight = 1)
  
#   # --- Sidebar Controls ---
  
#   # Folder selection for ABF files
#   folderLabel <- tklabel(sidebarFrame, text = "Select ABF Folder:")
#   tkgrid(folderLabel, row = 0, column = 0, sticky = "w")
#   folderPathVar <- tclVar("")
#   folderEntry <- tkentry(sidebarFrame, textvariable = folderPathVar, width = 30)
#   tkgrid(folderEntry, row = 0, column = 1, sticky = "w")
  
#   browseFolderButton <- tkbutton(sidebarFrame, text = "Browse", command = function() {
#     folderPath <- tclvalue(tkchooseDirectory())
#     if (nchar(folderPath) > 0) {
#       tclvalue(folderPathVar) <- folderPath
#       # List available ABF files in the folder:
#       abf_list <- list.files(path = folderPath, pattern = "\\.abf$", ignore.case = TRUE)
#       if (length(abf_list) == 0) {
#         tkmessageBox(message = "No ABF files found in the selected folder.")
#       } else {
#         tkdelete(abfListBox, 0, "end")
#         for (f in abf_list) {
#           tkinsert(abfListBox, "end", f)
#         }
#         # Auto-load metadata from the first ABF file to update additional parameters
#         firstFilePath <- file.path(folderPath, abf_list[1])
#         ds <- readABF(firstFilePath)
#         dummy_result <- list(metadata = list(extract_metadata(ds)))
#         updateAdditionalParams(dummy_result)
#       }
#     }
#   })
#   tkgrid(browseFolderButton, row = 0, column = 2, padx = 5)
  
#   # Listbox to show available ABF files
#   abfListLabel <- tklabel(sidebarFrame, text = "ABF Files:")
#   tkgrid(abfListLabel, row = 1, column = 0, sticky = "w", pady = 5)
#   abfListBox <- tklistbox(sidebarFrame, height = 5, width = 40, selectmode = "multiple")
#   tkgrid(abfListBox, row = 2, column = 0, columnspan = 3, pady = 5, sticky = "we")
  
#   # --- Additional Parameters (above Baseline) ---
#   # Create a new frame in sidebarFrame for experiment parameters
#   paramFrame <- tkframe(sidebarFrame)
#   tkgrid(paramFrame, row = 3, column = 0, columnspan = 3, sticky = "w")
  
#   # Create tcl variables for additional parameters:
#   experimentVar <- tclVar("Voltage Clamp")
#   unitVar       <- tclVar("")
#   dataColVar    <- tclVar("")
#   dtVar         <- tclVar("")
#   ntracesVar    <- tclVar("")  # New variable for number of traces
  
#   # Experiment type selector (Voltage Clamp vs. Current Clamp)
#   tkgrid(tklabel(paramFrame, text = "Experiment:"), row = 0, column = 0, sticky = "w")
#   experimentCombo <- ttkcombobox(paramFrame, textvariable = experimentVar, 
#                                  values = c("Voltage Clamp", "Current Clamp"), width = 15)
#   tkgrid(experimentCombo, row = 0, column = 1, sticky = "w")
  
#   # Display Units (auto-populated)
#   tkgrid(tklabel(paramFrame, text = "Units:"), row = 1, column = 0, sticky = "w")
#   unitEntry <- tkentry(paramFrame, textvariable = unitVar, width = 10, state = "readonly")
#   tkgrid(unitEntry, row = 1, column = 1, sticky = "w")
  
#   # Data column (auto-populated)
#   tkgrid(tklabel(paramFrame, text = "Data Column:"), row = 2, column = 0, sticky = "w")
#   dataColEntry <- tkentry(paramFrame, textvariable = dataColVar, width = 10, state = "readonly")
#   tkgrid(dataColEntry, row = 2, column = 1, sticky = "w")
  
#   # dt (sampling interval in ms)
#   tkgrid(tklabel(paramFrame, text = "dt (ms):"), row = 3, column = 0, sticky = "w")
#   dtEntry <- tkentry(paramFrame, textvariable = dtVar, width = 10, state = "readonly")
#   tkgrid(dtEntry, row = 3, column = 1, sticky = "w")
  
#   # Ntraces (number of traces)
#   tkgrid(tklabel(paramFrame, text = "Traces:"), row = 4, column = 0, sticky = "w")
#   ntracesEntry <- tkentry(paramFrame, textvariable = ntracesVar, width = 10, state = "readonly")
#   tkgrid(ntracesEntry, row = 4, column = 1, sticky = "w")
  
#   # --- Baseline and Other Parameter Inputs ---
#   baselineVar <- tclVar("100")
#   stimTimeVar <- tclVar("350")
#   xbarVar <- tclVar("100")
#   ybarVar <- tclVar("50")
  
#   tkgrid(tklabel(sidebarFrame, text = "Baseline:"), row = 4, column = 0, sticky = "w")
#   tkgrid(tkentry(sidebarFrame, textvariable = baselineVar, width = 10), row = 4, column = 1, sticky = "w")
#   tkgrid(tklabel(sidebarFrame, text = "Stimulation Time:"), row = 5, column = 0, sticky = "w")
#   tkgrid(tkentry(sidebarFrame, textvariable = stimTimeVar, width = 10), row = 5, column = 1, sticky = "w")
#   tkgrid(tklabel(sidebarFrame, text = "x-bar length:"), row = 6, column = 0, sticky = "w")
#   tkgrid(tkentry(sidebarFrame, textvariable = xbarVar, width = 10), row = 6, column = 1, sticky = "w")
#   tkgrid(tklabel(sidebarFrame, text = "y-bar length:"), row = 7, column = 0, sticky = "w")
#   tkgrid(tkentry(sidebarFrame, textvariable = ybarVar, width = 10), row = 7, column = 1, sticky = "w")
  
#   # Text widget for status messages (console)
#   consoleText <- tktext(sidebarFrame, width = 40, height = 4)
#   tkgrid(consoleText, row = 8, column = 0, columnspan = 3, pady = 5)
  
#   # --- Update Function ---
#   # Update Units, Data Column, dt, and Ntraces based on analysis result and experiment type.
#   updateAdditionalParams <- function(result) {
#     if (!is.null(result) && length(result$metadata) >= 1) {
#       meta1 <- result$metadata[[1]]
#       # Update dt (sampling interval in ms)
#       tclvalue(dtVar) <- as.character(meta1$samplingIntervalInSec * 1000)
#       # Update Ntraces from header information
#       if (!is.null(meta1$header$lActualEpisodes)) {
#         tclvalue(ntracesVar) <- as.character(meta1$header$lActualEpisodes)
#       } else {
#         tclvalue(ntracesVar) <- "N/A"
#       }
#       # Choose column index based on experiment type:
#       expType <- tclvalue(experimentVar)
#       col_idx <- choose_data_column(meta1$channelUnits, expType)
#       if (!is.na(col_idx)) {
#         tclvalue(unitVar) <- meta1$channelUnits[col_idx]
#         tclvalue(dataColVar) <- as.character(col_idx)
#       } else {
#         tclvalue(unitVar) <- "N/A"
#         tclvalue(dataColVar) <- "N/A"
#       }
#     }
#   }
  
#   # Bind event when an item is selected from the experiment combobox
#   tkbind(experimentCombo, "<<ComboboxSelected>>", function() {
#     if (exists("abf_analysis_result", envir = .GlobalEnv)) {
#       result <- get("abf_analysis_result", envir = .GlobalEnv)
#       updateAdditionalParams(result)
#     }
#   })
  
#   # --- Analysis Action: Run Analysis button ---
#   runAnalysis <- function() {
#     folderPath <- tclvalue(folderPathVar)
#     if (nchar(folderPath) == 0) {
#       tkmessageBox(message = "Please select an ABF folder first.")
#       return()
#     }
    
#     selIndices <- as.integer(tkcurselection(abfListBox))
#     allFiles <- as.character(tkget(abfListBox, 0, "end"))
#     abf_files <- if (length(selIndices) == 0) allFiles else allFiles[selIndices + 1]
    
#     if (length(abf_files) == 0) {
#       tkmessageBox(message = "No ABF files selected.")
#       return()
#     }
    
#     baseline <- as.numeric(tclvalue(baselineVar))
#     stimTime <- as.numeric(tclvalue(stimTimeVar))
#     xbar     <- as.numeric(tclvalue(xbarVar))
#     ybar     <- as.numeric(tclvalue(ybarVar))
    
#     result <- tryCatch({
#       abf_out <- abf_averages(abf_files = abf_files, abf_path = folderPath, 
#                               baseline = baseline, stimulation_time = stimTime, 
#                               xbar = xbar, ybar = ybar)
#       return(abf_out)
#     }, error = function(e) {
#       tkmessageBox(message = paste("Error during analysis:", e$message))
#       return(NULL)
#     })
    
#     if (!is.null(result)) {
#       tkdelete(consoleText, "1.0", "end")
#       tkinsert(consoleText, "end", paste("Analysis complete. Processed", length(abf_files), "file(s)."))
#       assign("abf_analysis_result", result, envir = .GlobalEnv)
#       updateAdditionalParams(result)
#       tkrreplot(plotWidget)
#     }
#   }
  
#   runAnalysisButton <- tkbutton(sidebarFrame, text = "Run Analysis", command = runAnalysis)
#   tkgrid(runAnalysisButton, row = 9, column = 0, columnspan = 3, pady = 5)
  
#   # --- Main Panel: Plot display via tkrplot ---
#   drawPlot <- function() {
#     if (exists("abf_analysis_result", envir = .GlobalEnv)) {
#       result <- get("abf_analysis_result", envir = .GlobalEnv)
#       datasets <- result$datasets
#       if (length(datasets) >= 1) {
#         iii <- 1
#         dt_val <- datasets[[iii]]$samplingIntervalInSec * 1000
#         time <- seq(0, by = dt_val, length.out = nrow(datasets[[iii]]$data[[1]]))
#         col_num <- as.numeric(tclvalue(dataColVar))
#         if (is.na(col_num) || col_num < 1 || col_num > ncol(datasets[[iii]]$data[[1]])) {
#           col_num <- 1
#         }
#         trace <- datasets[[iii]]$data[[1]][, col_num]
#         plot(time, trace, col = "indianred", xlab = "Time (ms)", type = "l",
#              bty = "l", las = 1, main = paste("Trace", iii))
#         axis(1)
#         axis(2)
#       } else {
#         plot.new()
#         text(0.5, 0.5, "No data available")
#       }
#     } else {
#       plot.new()
#       text(0.5, 0.5, "No analysis result to display")
#     }
#   }
  
#   # Create the plot widget using tkrplot() (not tkrreplot())
#   plotWidget <- tkrplot(mainFrame, fun = drawPlot, hscale = 1.5, vscale = 1.5)
#   tkgrid(plotWidget, row = 0, column = 0, sticky = "nsew")
  
#   tkfocus(tt)
# }

# # -----------------------------------------------------------------------------
# # Launch the ABF Analysis interface:
# # -----------------------------------------------------------------------------
# ABF_analysis_tk()





# improve this code :


# # Remove all objects from the environment
# rm(list = ls(all = TRUE))

# # -----------------------------------------------------------------------------
# # Package loading function and required packages:
# # -----------------------------------------------------------------------------
# load_required_packages <- function(packages) {
#   new.packages <- packages[!(packages %in% installed.packages()[, "Package"])]
#   if (length(new.packages)) install.packages(new.packages)
#   invisible(lapply(packages, library, character.only = TRUE))
# }

# required.packages <- c("robustbase", "minpack.lm", "Rcpp", "signal", "readABF", "tcltk", "tkrplot")
# load_required_packages(required.packages)

# # -----------------------------------------------------------------------------
# # Function to extract key metadata from an ABF dataset
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

# # -----------------------------------------------------------------------------
# # Custom plotting function (egs_plot)
# # -----------------------------------------------------------------------------
# egs_plot <- function(x, y, sign = -1, xlim = NULL, ylim = NULL, lwd = 1, 
#                      show_text = FALSE, height = 4, width = 2.5, 
#                      xbar = 100, ybar = 50,  color = "#4C77BB", show_bar = FALSE) {
  
#   if (is.null(ylim))
#     ylim <- if (sign == 1) c(0, max(y)) else c(-max(-y), 0)
  
#   if (is.null(xlim))
#     xlim <- c(min(x), max(x))
  
#   idx1 <- which.min(abs(x - xlim[1]))
#   idx2 <- which.min(abs(x - xlim[2]))
  
#   plot(x[idx1:idx2], y[idx1:idx2], type = "l", col = color, xlim = xlim, 
#        ylim = ylim, bty = "n", lwd = lwd, lty = 1, axes = FALSE, 
#        frame = FALSE, xlab = "", ylab = "")
  
#   if (show_bar) {
#     # Calculate positions for scale bar and labels
#     ybar_start <- min(ylim) + (max(ylim) - min(ylim)) / 20
#     x_start <- max(xlim) - xbar - 50
#     y_start <- ybar_start
#     x_end <- x_start + xbar
#     y_end <- y_start + ybar
    
#     segments(x_start, y_start, x_end, y_start, lwd = lwd, col = "black")
#     segments(x_start, y_start, x_start, y_end, lwd = lwd, col = "black")
    
#     if (show_text) {
#       text(x = (x_start + x_end) / 2, y = y_start - ybar / 20, 
#            labels = paste(xbar, "ms"), adj = c(0.5, 1))
#       text(x = x_start - xbar / 4, y = (y_start + y_end) / 2,  
#            labels = paste(ybar, "pA"), adj = c(0.5, 0.5), srt = 90)
#     }
#   }
# }

# # -----------------------------------------------------------------------------
# # Modified abf_averages function (adapted for UI use)
# # -----------------------------------------------------------------------------
# abf_averages <- function(abf_files = NULL, abf_path = NULL, traces2average = NULL,
#                          baseline = 100, stimulation_time = 350, ylim = NULL, 
#                          xlim = NULL, color = "#4C77BB", xbar = 100, ybar = 50, 
#                          width = 5.25, height = 2.75, save = FALSE) {
  
#   if (!is.null(traces2average) && length(traces2average) != length(abf_files)) {
#     stop("error: if specified, traces2average list must have the same length as the input abf_files list")
#   }
  
#   baseline2zero <- function(y, dt, stimulation_time, baseline) {
#     idx1 <- (stimulation_time - baseline) / dt
#     idx2 <- baseline / dt
#     y1 <- y[idx1:length(y)]
#     y1 <- y1 - mean(y1[1:idx2])
#     return(y1 - mean(y1[1:idx2]))
#   }
  
#   # Set working folder
#   abf_path <- if (is.null(abf_path)) getwd() else abf_path
#   setwd(abf_path)
  
#   N <- length(abf_files)
  
#   # Read each ABF file to create a list of datasets
#   datasets <- lapply(1:N, function(ii) readABF(abf_files[ii]))
#   names(datasets) <- abf_files
  
#   # Extract metadata from each dataset
#   metadata <- lapply(datasets, extract_metadata)
  
#   headers <- lapply(1:N, function(iii) datasets[[iii]][names(datasets[[iii]]) != "data"])
#   names(headers) <- abf_files
  
#   sampling_intervals <- sapply(1:N, function(ii) datasets[[ii]]$samplingIntervalInSec * 1000)
  
#   responses <- lapply(1:N, function(iii) {
#     sapply(1:length(datasets[[iii]]$data), function(ii) 
#       datasets[[iii]]$data[[ii]][, 1])
#   })
#   names(responses) <- abf_files
  
#   responses0 <- lapply(1:N, function(iii) {
#     sapply(1:dim(responses[[iii]])[2], function(ii) 
#       baseline2zero(responses[[iii]][, ii], dt = sampling_intervals[iii], 
#                     stimulation_time = stimulation_time, baseline = baseline))
#   })
#   names(responses0) <- abf_files
  
#   if (is.null(traces2average)) {
#     responses0_mean <- lapply(1:N, function(iii) apply(responses0[[iii]], 1, mean))
#   } else {
#     responses0_mean <- lapply(1:N, function(iii) 
#       apply(responses0[[iii]][, traces2average[[iii]]], 1, mean))
#   }
  
#   time <- lapply(1:N, function(iii) {
#     dt <- datasets[[iii]]$samplingIntervalInSec * 1000
#     0:dt:(length(responses0_mean[[iii]]) - 1) * dt
#   })
  
#   # Plot onto the current device (i.e. our tkrplot window)
#   par(mfrow = c(1, N))
#   show_bar <- rep(FALSE, N)
#   show_bar[N] <- TRUE
  
#   for (ii in 1:N) {
#     egs_plot(x = time[[ii]], y = responses0_mean[[ii]], ylim = ylim, 
#              xlim = xlim, color = "darkgrey", width = width, height = height, 
#              show_text = FALSE, show_bar = show_bar[ii])
#     axis(1)
#     axis(2)
#   }
  
#   if (save) {
#     warning("save_graph not implemented in this UI example")
#   }
  
#   return(list("headers" = headers, 
#               "raw_data" = responses, 
#               "baseline_corrected_data" = responses0, 
#               "baseline_corrected_mean_data" = responses0_mean,
#               "datasets" = datasets,
#               "metadata" = metadata))
# }

# # If you have additional functions in your nNLS functions file, source it here:
# # source('/path/to/your/nNLS functions.R')

# # -----------------------------------------------------------------------------
# # Build the ABF Analysis user interface using Tcl/Tk
# # -----------------------------------------------------------------------------
# ABF_analysis_tk <- function() {
#   tt <- tktoplevel()
#   tkwm.title(tt, "ABF Analysis")
  
#   # Divide the window into sidebar (controls) and main (plot) panels
#   sidebarFrame <- tkframe(tt)
#   mainFrame <- tkframe(tt)
#   tkgrid(sidebarFrame, row = 0, column = 0, sticky = "ns")
#   tkgrid(mainFrame, row = 0, column = 1, sticky = "nsew")
#   tkgrid.rowconfigure(tt, 0, weight = 1)
#   tkgrid.columnconfigure(tt, 1, weight = 1)
  
#   # --- Sidebar Controls ---
#   # Folder selection for ABF files
#   folderLabel <- tklabel(sidebarFrame, text = "Select ABF Folder:")
#   tkgrid(folderLabel, row = 0, column = 0, sticky = "w")
#   folderPathVar <- tclVar("")
#   folderEntry <- tkentry(sidebarFrame, textvariable = folderPathVar, width = 30)
#   tkgrid(folderEntry, row = 0, column = 1, sticky = "w")
  
#   browseFolderButton <- tkbutton(sidebarFrame, text = "Browse", command = function() {
#     folderPath <- tclvalue(tkchooseDirectory())
#     if (nchar(folderPath) > 0) {
#       tclvalue(folderPathVar) <- folderPath
#       # List available ABF files in the folder:
#       abf_list <- list.files(path = folderPath, pattern = "\\.abf$", ignore.case = TRUE)
#       if (length(abf_list) == 0) {
#         tkmessageBox(message = "No ABF files found in the selected folder.")
#       } else {
#         tkdelete(abfListBox, 0, "end")
#         for (f in abf_list) {
#           tkinsert(abfListBox, "end", f)
#         }
#       }
#     }
#   })

#   tkgrid(browseFolderButton, row = 0, column = 2, padx = 5)
  
#   # Listbox to show available ABF files; user can select one or more files
#   abfListLabel <- tklabel(sidebarFrame, text = "ABF Files:")
#   tkgrid(abfListLabel, row = 1, column = 0, sticky = "w", pady = 5)
#   abfListBox <- tklistbox(sidebarFrame, height = 5, width = 40, selectmode = "multiple")
#   tkgrid(abfListBox, row = 2, column = 0, columnspan = 3, pady = 5, sticky = "we")
  
#   # Parameter inputs
#   baselineVar <- tclVar("100")
#   stimTimeVar <- tclVar("350")
#   xbarVar <- tclVar("100")
#   ybarVar <- tclVar("50")
  
#   tkgrid(tklabel(sidebarFrame, text = "Baseline:"), row = 3, column = 0, sticky = "w")
#   tkgrid(tkentry(sidebarFrame, textvariable = baselineVar, width = 10), row = 3, column = 1, sticky = "w")
#   tkgrid(tklabel(sidebarFrame, text = "Stimulation Time:"), row = 4, column = 0, sticky = "w")
#   tkgrid(tkentry(sidebarFrame, textvariable = stimTimeVar, width = 10), row = 4, column = 1, sticky = "w")
#   tkgrid(tklabel(sidebarFrame, text = "x-bar length:"), row = 5, column = 0, sticky = "w")
#   tkgrid(tkentry(sidebarFrame, textvariable = xbarVar, width = 10), row = 5, column = 1, sticky = "w")
#   tkgrid(tklabel(sidebarFrame, text = "y-bar length:"), row = 6, column = 0, sticky = "w")
#   tkgrid(tkentry(sidebarFrame, textvariable = ybarVar, width = 10), row = 6, column = 1, sticky = "w")
  
#   # Text widget for status messages (console)
#   consoleText <- tktext(sidebarFrame, width = 40, height = 4)
#   tkgrid(consoleText, row = 7, column = 0, columnspan = 3, pady = 5)
  
#   # --- Analysis action: Run Analysis button
#   runAnalysis <- function() {
#     folderPath <- tclvalue(folderPathVar)
#     if (nchar(folderPath) == 0) {
#       tkmessageBox(message = "Please select an ABF folder first.")
#       return()
#     }
    
#     # Get files selected from the listbox; if none selected, use all available
#     selIndices <- as.integer(tkcurselection(abfListBox))
#     allFiles <- as.character(tkget(abfListBox, 0, "end"))
#     if (length(selIndices) == 0) {
#       abf_files <- allFiles
#     } else {
#       abf_files <- allFiles[selIndices + 1]  # Note: indices are 0-based
#     }
    
#     if (length(abf_files) == 0) {
#       tkmessageBox(message = "No ABF files selected.")
#       return()
#     }
    
#     baseline <- as.numeric(tclvalue(baselineVar))
#     stimTime <- as.numeric(tclvalue(stimTimeVar))
#     xbar <- as.numeric(tclvalue(xbarVar))
#     ybar <- as.numeric(tclvalue(ybarVar))
    
#     # Run the ABF analysis using your abf_averages function
#     result <- tryCatch({
#       abf_out <- abf_averages(abf_files = abf_files, abf_path = folderPath, 
#                               baseline = baseline, stimulation_time = stimTime, 
#                               xbar = xbar, ybar = ybar)
#       return(abf_out)
#     }, error = function(e) {
#       tkmessageBox(message = paste("Error during analysis:", e$message))
#       return(NULL)
#     })
    
#     if (!is.null(result)) {
#       tkdelete(consoleText, "1.0", "end")
#       tkinsert(consoleText, "end", paste("Analysis complete. Processed", length(abf_files), "file(s)."))
#       # Save result to the global environment for further use if needed
#       assign("abf_analysis_result", result, envir = .GlobalEnv)
#       tkrreplot(plotWidget)
#     }
#   }
  
#   runAnalysisButton <- tkbutton(sidebarFrame, text = "Run Analysis", command = runAnalysis)
#   tkgrid(runAnalysisButton, row = 8, column = 0, columnspan = 3, pady = 5)
  
#   # --- Main Panel: Plot display via tkrplot ---
#   drawPlot <- function() {
#     # If the analysis result exists, plot the first trace from the first dataset.
#     if (exists("abf_analysis_result", envir = .GlobalEnv)) {
#       result <- get("abf_analysis_result", envir = .GlobalEnv)
#       datasets <- result$datasets
#       if (length(datasets) >= 1) {
#         iii <- 1
#         dt <- datasets[[iii]]$samplingIntervalInSec * 1000
#         time <- seq(0, by = dt, length.out = nrow(datasets[[iii]]$data[[1]]))
#         trace <- datasets[[iii]]$data[[1]][, 1]
#         plot(time, trace, col = "indianred", xlab = "Time (ms)", type = "l",
#              bty = "l", las = 1, main = paste("Trace", iii))
#         axis(1)
#         axis(2)
#       } else {
#         plot.new()
#         text(0.5, 0.5, "No data available")
#       }
#     } else {
#       plot.new()
#       text(0.5, 0.5, "No analysis result to display")
#     }
#   }
  
#   plotWidget <- tkrplot(mainFrame, fun = drawPlot, hscale = 1.5, vscale = 1.5)
#   tkgrid(plotWidget, row = 0, column = 0, sticky = "nsew")
  
#   tkfocus(tt)
# }

# # -----------------------------------------------------------------------------
# # Launch the ABF Analysis interface:
# # -----------------------------------------------------------------------------
# ABF_analysis_tk()


# ii=1; all_metadata[[ii]]$samplingIntervalInSec * 1000


#         dt <- datasets[[iii]]$samplingIntervalInSec * 1000
#         time <- seq(0, by = dt, length.out = nrow(datasets[[iii]]$data[[1]]))


# Ntraces <- all_metadata[[ii]]$header$lActualEpisodes

# Units <- all_metadata[[ii]]$channelUnits 


# datasets correspond to

# datasets[[ii]]$data[[1]][1:10,]
#            [,1]      [,2]       [,3]
#  [1,] -454.7119 -50.32349 10.9863276
#  [2,] -449.8291 -49.19434  5.4931638
#  [3,] -444.9463 -49.16382  9.7656245
#  [4,] -445.5566 -50.47608 -0.6103515
#  [5,] -449.2187 -50.56763 11.5966791
#  [6,] -454.7119 -49.77417 11.5966791
#  [7,] -454.7119 -49.10278  4.8828123
#  [8,] -452.2705 -49.40796  3.6621092
#  [9,] -449.8291 -50.62866 10.9863276
# [10,] -451.0498 -50.56763  7.3242184

# all_metadata[[ii]]$channelUnits 
# [1] "pA" "mV" "pA"
#  with Ntrace == length(datasets[[ii]]$data)


# to main panel (above baseline) add the following

# 1. 'experiment' (options voltage-clamp or current-clamp)
# 2. 'units' (auto populate from all_metadata[[ii]]$channelUnits  if voltage clamp take value in first col containing A eg pA if cuurent clamp take value in first col containing V eg nV
# 3. 'data column' (auto populate from all_metadata[[ii]]$channelUnits  if voltage clamp take first col containing A (maybe in any unit of A eg nA or pA) if current clamp then first col containing V eg mV)
# use this col number to extract the data datasets[[ii]]$data[[1]][,col_num]
# 4. 'dt' extract all_metadata[[ii]]$samplingIntervalInSec * 1e3 (units are ms)



# # Remove all objects from the environment
# rm(list = ls(all = TRUE))

# # -----------------------------------------------------------------------------
# # Package loading function and required packages:
# # -----------------------------------------------------------------------------
# load_required_packages <- function(packages) {
#   new.packages <- packages[!(packages %in% installed.packages()[, "Package"])]
#   if (length(new.packages)) install.packages(new.packages)
#   invisible(lapply(packages, library, character.only = TRUE))
# }

# required.packages <- c("robustbase", "minpack.lm", "Rcpp", "signal", "readABF", "tcltk", "tkrplot")
# load_required_packages(required.packages)

# # -----------------------------------------------------------------------------
# # Custom plotting function (egs_plot)
# # -----------------------------------------------------------------------------
# egs_plot <- function(x, y, sign = -1, xlim = NULL, ylim = NULL, lwd = 1, 
#                      show_text = FALSE, height = 4, width = 2.5, 
#                      xbar = 100, ybar = 50,  color = "#4C77BB", show_bar = FALSE) {
  
#   if (is.null(ylim))
#     ylim <- if (sign == 1) c(0, max(y)) else c(-max(-y), 0)
  
#   if (is.null(xlim))
#     xlim <- c(min(x), max(x))
  
#   idx1 <- which.min(abs(x - xlim[1]))
#   idx2 <- which.min(abs(x - xlim[2]))
  
#   plot(x[idx1:idx2], y[idx1:idx2], type = "l", col = color, xlim = xlim, 
#        ylim = ylim, bty = "n", lwd = lwd, lty = 1, axes = FALSE, 
#        frame = FALSE, xlab = "", ylab = "")
  
#   if (show_bar) {
#     # Calculate positions for scale bar and labels
#     ybar_start <- min(ylim) + (max(ylim) - min(ylim)) / 20
#     x_start <- max(xlim) - xbar - 50
#     y_start <- ybar_start
#     x_end <- x_start + xbar
#     y_end <- y_start + ybar
    
#     segments(x_start, y_start, x_end, y_start, lwd = lwd, col = "black")
#     segments(x_start, y_start, x_start, y_end, lwd = lwd, col = "black")
    
#     if (show_text) {
#       text(x = (x_start + x_end) / 2, y = y_start - ybar / 20, 
#            labels = paste(xbar, "ms"), adj = c(0.5, 1))
#       text(x = x_start - xbar / 4, y = (y_start + y_end) / 2,  
#            labels = paste(ybar, "pA"), adj = c(0.5, 0.5), srt = 90)
#     }
#   }
# }

# # -----------------------------------------------------------------------------
# # Modified abf_averages function (adapted for UI use)
# # -----------------------------------------------------------------------------
# abf_averages <- function(abf_files = NULL, abf_path = NULL, traces2average = NULL,
#                          baseline = 100, stimulation_time = 350, ylim = NULL, 
#                          xlim = NULL, color = "#4C77BB", xbar = 100, ybar = 50, 
#                          width = 5.25, height = 2.75, save = FALSE) {
  
#   if (!is.null(traces2average) && length(traces2average) != length(abf_files)) {
#     stop("error: if specified, traces2average list must have the same length as the input abf_files list")
#   }
  
#   baseline2zero <- function(y, dt, stimulation_time, baseline) {
#     idx1 <- (stimulation_time - baseline) / dt
#     idx2 <- baseline / dt
#     y1 <- y[idx1:length(y)]
#     y1 <- y1 - mean(y1[1:idx2])
#     return(y1 - mean(y1[1:idx2]))
#   }
  
#   abf_path <- if (is.null(abf_path)) getwd() else abf_path
#   setwd(abf_path)
  
#   N <- length(abf_files)
#   datasets <- lapply(1:N, function(ii) readABF(abf_files[ii]))
#   names(datasets) <- abf_files
  
#   headers <- lapply(1:N, function(iii) datasets[[iii]][names(datasets[[iii]]) != "data"])
#   names(headers) <- abf_files
  
#   sampling_intervals <- sapply(1:N, function(ii) datasets[[ii]]$samplingIntervalInSec * 1000)
  
#   responses <- lapply(1:N, function(iii) {
#     sapply(1:length(datasets[[iii]]$data), function(ii) 
#       datasets[[iii]]$data[[ii]][, 1])
#   })
#   names(responses) <- abf_files
  
#   responses0 <- lapply(1:N, function(iii) {
#     sapply(1:dim(responses[[iii]])[2], function(ii) 
#       baseline2zero(responses[[iii]][, ii], dt = sampling_intervals[iii], 
#                     stimulation_time = stimulation_time, baseline = baseline))
#   })
#   names(responses0) <- abf_files
  
#   if (is.null(traces2average)) {
#     responses0_mean <- lapply(1:N, function(iii) apply(responses0[[iii]], 1, mean))
#   } else {
#     responses0_mean <- lapply(1:N, function(iii) 
#       apply(responses0[[iii]][, traces2average[[iii]]], 1, mean))
#   }
  
#   time <- lapply(1:N, function(iii) {
#     dt <- datasets[[iii]]$samplingIntervalInSec * 1000
#     0:dt:(length(responses0_mean[[iii]]) - 1) * dt
#   })
  
#   # Plot onto the current device (i.e. our tkrplot window)
#   par(mfrow = c(1, N))
#   show_bar <- rep(FALSE, N)
#   show_bar[N] <- TRUE
  
#   for (ii in 1:N) {
#     egs_plot(x = time[[ii]], y = responses0_mean[[ii]], ylim = ylim, 
#              xlim = xlim, color = "darkgrey", width = width, height = height, 
#              show_text = FALSE, show_bar = show_bar[ii])
#     axis(1)
#     axis(2)
#   }
  
#   if (save) {
#     warning("save_graph not implemented in this UI example")
#   }
  
#   return(list("headers" = headers, 
#               "raw_data" = responses, 
#               "baseline_corrected_data" = responses0, 
#               "baseline_corrected_mean_data" = responses0_mean,
#               "datasets" = datasets))
# }

# # If you have additional functions in your nNLS functions file, source it here:
# # source('/path/to/your/nNLS functions.R')


# # -----------------------------------------------------------------------------
# # Build the ABF Analysis user interface using Tcl/Tk
# # -----------------------------------------------------------------------------
# ABF_analysis_tk <- function() {
#   tt <- tktoplevel()
#   tkwm.title(tt, "ABF Analysis")
  
#   # Divide the window into sidebar (controls) and main (plot) panels
#   sidebarFrame <- tkframe(tt)
#   mainFrame <- tkframe(tt)
#   tkgrid(sidebarFrame, row = 0, column = 0, sticky = "ns")
#   tkgrid(mainFrame, row = 0, column = 1, sticky = "nsew")
#   tkgrid.rowconfigure(tt, 0, weight = 1)
#   tkgrid.columnconfigure(tt, 1, weight = 1)
  
#   # --- Sidebar Controls ---
#   # Folder selection for ABF files
#   folderLabel <- tklabel(sidebarFrame, text = "Select ABF Folder:")
#   tkgrid(folderLabel, row = 0, column = 0, sticky = "w")
#   folderPathVar <- tclVar("")
#   folderEntry <- tkentry(sidebarFrame, textvariable = folderPathVar, width = 30)
#   tkgrid(folderEntry, row = 0, column = 1, sticky = "w")
  
#   browseFolderButton <- tkbutton(sidebarFrame, text = "Browse", command = function() {
#     folderPath <- tclvalue(tkchooseDirectory())
#     if (nchar(folderPath) > 0) {
#       tclvalue(folderPathVar) <- folderPath
#       # List available ABF files in the folder:
#       abf_list <- list.files(path = folderPath, pattern = "\\.abf$", ignore.case = TRUE)
#       if (length(abf_list) == 0) {
#         tkmessageBox(message = "No ABF files found in the selected folder.")
#       } else {
#         tkdelete(abfListBox, 0, "end")
#         for (f in abf_list) {
#           tkinsert(abfListBox, "end", f)
#         }
#       }
#     }
#   })

#   tkgrid(browseFolderButton, row = 0, column = 2, padx = 5)
  
#   # Listbox to show available ABF files; user can select one or more files
#   abfListLabel <- tklabel(sidebarFrame, text = "ABF Files:")
#   tkgrid(abfListLabel, row = 1, column = 0, sticky = "w", pady = 5)
#   abfListBox <- tklistbox(sidebarFrame, height = 5, width = 40, selectmode = "multiple")
#   tkgrid(abfListBox, row = 2, column = 0, columnspan = 3, pady = 5, sticky = "we")
  
#   # Parameter inputs
#   baselineVar <- tclVar("100")
#   stimTimeVar <- tclVar("350")
#   xbarVar <- tclVar("100")
#   ybarVar <- tclVar("50")
  
#   tkgrid(tklabel(sidebarFrame, text = "Baseline:"), row = 3, column = 0, sticky = "w")
#   tkgrid(tkentry(sidebarFrame, textvariable = baselineVar, width = 10), row = 3, column = 1, sticky = "w")
#   tkgrid(tklabel(sidebarFrame, text = "Stimulation Time:"), row = 4, column = 0, sticky = "w")
#   tkgrid(tkentry(sidebarFrame, textvariable = stimTimeVar, width = 10), row = 4, column = 1, sticky = "w")
#   tkgrid(tklabel(sidebarFrame, text = "x-bar length:"), row = 5, column = 0, sticky = "w")
#   tkgrid(tkentry(sidebarFrame, textvariable = xbarVar, width = 10), row = 5, column = 1, sticky = "w")
#   tkgrid(tklabel(sidebarFrame, text = "y-bar length:"), row = 6, column = 0, sticky = "w")
#   tkgrid(tkentry(sidebarFrame, textvariable = ybarVar, width = 10), row = 6, column = 1, sticky = "w")
  
#   # Text widget for status messages (console)
#   consoleText <- tktext(sidebarFrame, width = 40, height = 4)
#   tkgrid(consoleText, row = 7, column = 0, columnspan = 3, pady = 5)
  
#   # --- Analysis action: Run Analysis button
#   runAnalysis <- function() {
#     folderPath <- tclvalue(folderPathVar)
#     if (nchar(folderPath) == 0) {
#       tkmessageBox(message = "Please select an ABF folder first.")
#       return()
#     }
    
#     # Get files selected from the listbox; if none selected, use all available
#     selIndices <- as.integer(tkcurselection(abfListBox))
#     allFiles <- as.character(tkget(abfListBox, 0, "end"))
#     if (length(selIndices) == 0) {
#       abf_files <- allFiles
#     } else {
#       abf_files <- allFiles[selIndices + 1]  # Note: indices are 0-based
#     }
    
#     if (length(abf_files) == 0) {
#       tkmessageBox(message = "No ABF files selected.")
#       return()
#     }
    
#     baseline <- as.numeric(tclvalue(baselineVar))
#     stimTime <- as.numeric(tclvalue(stimTimeVar))
#     xbar <- as.numeric(tclvalue(xbarVar))
#     ybar <- as.numeric(tclvalue(ybarVar))
    
#     # Run the ABF analysis using your abf_averages function
#     result <- tryCatch({
#       abf_out <- abf_averages(abf_files = abf_files, abf_path = folderPath, 
#                               baseline = baseline, stimulation_time = stimTime, 
#                               xbar = xbar, ybar = ybar)
#       return(abf_out)
#     }, error = function(e) {
#       tkmessageBox(message = paste("Error during analysis:", e$message))
#       return(NULL)
#     })
    
#     if (!is.null(result)) {
#       tkdelete(consoleText, "1.0", "end")
#       tkinsert(consoleText, "end", paste("Analysis complete. Processed", length(abf_files), "file(s)."))
#       # Save result to the global environment for further use if needed
#       assign("abf_analysis_result", result, envir = .GlobalEnv)
#       tkrreplot(plotWidget)
#     }
#   }
  
#   runAnalysisButton <- tkbutton(sidebarFrame, text = "Run Analysis", command = runAnalysis)
#   tkgrid(runAnalysisButton, row = 8, column = 0, columnspan = 3, pady = 5)
  
#   # --- Main Panel: Plot display via tkrplot ---
#   drawPlot <- function() {
#     # If the analysis result exists, plot the first trace from the first dataset.
#     if (exists("abf_analysis_result", envir = .GlobalEnv)) {
#       result <- get("abf_analysis_result", envir = .GlobalEnv)
#       datasets <- result$datasets
#       if (length(datasets) >= 1) {
#         iii <- 1
#         dt <- datasets[[iii]]$samplingIntervalInSec * 1000
#         time <- seq(0, by = dt, length.out = nrow(datasets[[iii]]$data[[1]]))
#         trace <- datasets[[iii]]$data[[1]][, 1]
#         plot(time, trace, col = "indianred", xlab = "Time (ms)", type = "l",
#              bty = "l", las = 1, main = paste("Trace", iii))
#         axis(1)
#         axis(2)
#       } else {
#         plot.new()
#         text(0.5, 0.5, "No data available")
#       }
#     } else {
#       plot.new()
#       text(0.5, 0.5, "No analysis result to display")
#     }
#   }
  
#   plotWidget <- tkrplot(mainFrame, fun = drawPlot, hscale = 1.5, vscale = 1.5)
#   tkgrid(plotWidget, row = 0, column = 0, sticky = "nsew")
  
#   tkfocus(tt)
# }

# # -----------------------------------------------------------------------------
# # Launch the ABF Analysis interface:
# # -----------------------------------------------------------------------------
# ABF_analysis_tk()