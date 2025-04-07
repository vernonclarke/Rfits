######################################################################################
# for xquartz to work properly in (some) systems open R from terminal:
# open -n -a R
rm(list = ls(all = TRUE))
graphics.off()

# Load and install necessary packages
load_required_packages <- function(packages) {
  new.packages <- packages[!(packages %in% installed.packages()[, 'Package'])]
  if (length(new.packages)) install.packages(new.packages)
  invisible(lapply(packages, library, character.only = TRUE))
}
required.packages <- c('ARTool', 'robustbase', 'minpack.lm', 'Rcpp', 'signal', 'dbscan')
load_required_packages(required.packages)

# insert your username here to define the correct path
username <- 'YourUsername'
username <- 'euo9382'

path_repository <- '/Documents/Repositories/Rfits'
# construct file path
file_path1 <- paste0('/Users/', username, path_repository)

source(paste0(file_path1, '/nNLS functions.R'))

library(shiny)
library(readxl)



# Define the PSC Analysis widget as a Shiny Gadget
PSC_analysis_widget <- function() {
  ui <- fluidPage(
    titlePanel('PSC Analysis'),
    sidebarLayout(
      sidebarPanel(
        fileInput('file', 'upload csv or xlsx', accept = c('.csv', '.xlsx')),
        uiOutput('column_selector'),
        
        tabsetPanel(
          tabPanel('main options',
            numericInput('dt', 'dt (ms):', 0.1),
            numericInput('stimulation_time', 'stimulation time:', 100),
            numericInput('baseline', 'baseline:', 50),
            numericInput('n', 'n:', 30),
            numericInput('y_abline', 'fit cutoff:', 0.1),
            selectInput('func', 'function:', choices = c('product1N', 'product2N', 'product3N'))
          ),
          tabPanel('fit options',
            numericInput('N', 'N:', 1),
            numericInput('IEI', 'IEI:', 50),
            numericInput('smooth', 'smooth:', 5),
            selectInput('method', 'method', c('BF.LM', 'LM', 'GN', 'port', 'robust', 'MLE')),
            selectInput('weight_method', 'weighting', c('none', '~y_sqrt', '~y')),
            checkboxInput('sequential_fit', 'sequential fit', FALSE),
            numericInput('interval_min', 'min interval:', 0.1),
            numericInput('interval_max', 'max interval:', 0.9),
            textInput('lower', 'lower bounds (comma-separated):', ''),
            textInput('upper', 'upper bounds (comma-separated):', ''),
            textInput('latency.limit', 'latency limit:', '')
          ),
          tabPanel('MLE Settings',
            numericInput('iter', 'MLE iterations:', 1000),
            numericInput('metropolis_scale', 'metropolis scale:', 1.5),
            numericInput('fit_attempts', 'fit attempts:', 10),
            checkboxInput('RWm', 'random walk metropolis:', FALSE)
          ),
          tabPanel('Advanced',
            checkboxInput('filter', 'filter', FALSE),
            numericInput('fc', 'filter cutoff (Hz):', 1000),
            numericInput('rel_decay_fit_limit', 'relative decay fit limit:', 0.1),
            numericInput('half_width_fit_limit', 'half-width fit limit:', 500),
            numericInput('seed', 'seed:', 42),
            numericInput('dp', 'decimal points:', 3),
            checkboxInput('fast_constraint', 'fast constraint', FALSE),
            selectInput('fast_constraint_method', 'fast constraint method', c('rise', 'peak')),
            textInput('fast_decay_limit', 'fast decay limit(s) (comma-separated):', ''),
            checkboxInput('first_delay_constraint', 'first delay constraint', FALSE)
          ),
          tabPanel('Graph Settings',
            numericInput('lwd', 'line width:', 1.2),
            textInput('xlab', 'x-axis label:', 'time (ms)'),
            textInput('ylab', 'y-axis label:', 'PSC (pA)'),
            numericInput('plot_height', 'plot height:', 5),
            numericInput('plot_width', 'plot width:', 5)
          )
        ),
        
        actionButton('run_analysis', 'run initial analysis'),
        conditionalPanel(
          condition = 'output.promptTmax',
          numericInput('user_tmax', 'enter maximum time for fitting:', value=NA),
          numericInput('stimulation_time_adj', 'adjust stimulation time:', value=100),
          numericInput('baseline_adj', 'adjust baseline:', value=50),
          actionButton('update_plot', 'update plot'),
          actionButton('run_main_analysis', 'run main analysis')
        ),
        conditionalPanel(
          condition = 'output.promptFastConstraint',
          checkboxInput('repeat_constraint', 'repeat with fast constraint?', FALSE),
          actionButton('run_final_analysis', 'run analysis with fast constraint applied')
        ),
        hr(),
        downloadButton('download_output', 'download output'),
        actionButton("clear_output", "Clear Output")
      ),
      mainPanel(
        plotOutput('plot'),
        verbatimTextOutput('console')
      )
    )
  )
  
  server <- function(input, output, session) {
    analysis_output <- reactiveVal(NULL)
    response_data <- reactiveVal(NULL)
    prompt_tmax <- reactiveVal(FALSE)
    suggested_tmax <- reactiveVal(NA)
    prompt_fast_constraint <- reactiveVal(FALSE)
    
    # A reactive flag to control clearing
    clear_flag <- reactiveVal(FALSE)
  
    data_loaded <- reactive({
      req(input$file)
      ext <- tools::file_ext(input$file$name)
      if(ext == 'csv') read.csv(input$file$datapath) else readxl::read_excel(input$file$datapath)
    })
    
    output$column_selector <- renderUI({
      req(data_loaded())
      selectInput('data_col', 'select column to analyse', choices=names(data_loaded()))
    })
    
    current_params <- reactiveValues(
      stimulation_time = NULL,
      baseline = NULL
    )
    
    observeEvent(input$run_analysis, {
      req(data_loaded(), input$data_col)
      response_data(data_loaded()[[input$data_col]])
      
      current_params$stimulation_time <- input$stimulation_time
      current_params$baseline <- input$baseline
      
      prompt_tmax(TRUE)
      clear_flag(FALSE)  # Ensure outputs are enabled
    })
    
    observeEvent(input$update_plot, {
      current_params$stimulation_time <- input$stimulation_time_adj
      current_params$baseline <- input$baseline_adj
    })
    
    output$plot <- renderPlot({
      if(clear_flag()) return(NULL)  # If cleared, render nothing
      req(prompt_tmax())
      suggested <- determine_tmax(
        y = response_data(), N = input$N, dt = input$dt,
        stimulation_time = current_params$stimulation_time,
        baseline = current_params$baseline, smooth = input$smooth,
        y_abline = input$y_abline, ylab=input$ylab, prompt = FALSE
      )
      suggested_tmax(suggested)
      # Insert your plotting code here (e.g., plot(suggested, ...))
    })
    
    output$console <- renderPrint({
      if(clear_flag()) return(NULL)
      req(analysis_output())
      analysis_output()
    })
    
    observeEvent(input$run_main_analysis, {
      req(input$user_tmax)
      # Reset clear_flag when a new analysis is run
      clear_flag(FALSE)
      
      latency.limit <- if (nzchar(input$latency.limit)) as.numeric(unlist(strsplit(input$latency.limit, ','))) else NULL
      lower <- if (nzchar(input$lower)) as.numeric(unlist(strsplit(input$lower, ','))) else NULL
      upper <- if (nzchar(input$upper)) as.numeric(unlist(strsplit(input$upper, ','))) else NULL
      fast.decay.limit <- if (nzchar(input$fast_decay_limit)) as.numeric(unlist(strsplit(input$fast_decay_limit, ','))) else NULL
  
      result <- analyse_PSC(
        response = response_data(), dt = input$dt, n = input$n, N = input$N, IEI = input$IEI,
        stimulation_time = current_params$stimulation_time,
        baseline = current_params$baseline, smooth = input$smooth,
        fit.limits = input$user_tmax, fast.constraint = FALSE,
        method = input$method, weight_method = input$weight_method, filter = input$filter, fc = input$fc,
        rel.decay.fit.limit = input$rel_decay_fit_limit,
        half_width_fit_limit = input$half_width_fit_limit, seed = input$seed,
        sequential.fit = input$sequential_fit, dp = input$dp,
        func = get(input$func), interval = c(input$interval_min, input$interval_max),
        lower = lower, upper = upper,
        MLEsettings = list(iter = input$iter, metropolis.scale = input$metropolis_scale,
                           fit.attempts = input$fit_attempts, RWm = input$RWm),
        fast.decay.limit = fast.decay.limit, fast.constraint.method = input$fast_constraint_method,
        first.delay.constraint = input$first_delay_constraint, latency.limit = latency.limit,
        lwd = input$lwd, xlab = input$xlab, ylab = input$ylab,
        height = input$plot_height, width = input$plot_width,
        return.output = TRUE
      )
  
      analysis_output(result)
      prompt_tmax(FALSE)
      prompt_fast_constraint(TRUE)
    })
    
    observeEvent(input$run_final_analysis, {
      clear_flag(FALSE)
      
      latency.limit <- if (nzchar(input$latency.limit)) as.numeric(unlist(strsplit(input$latency.limit, ','))) else NULL
      lower <- if (nzchar(input$lower)) as.numeric(unlist(strsplit(input$lower, ','))) else NULL
      upper <- if (nzchar(input$upper)) as.numeric(unlist(strsplit(input$upper, ','))) else NULL
      fast.decay.limit <- if (nzchar(input$fast_decay_limit)) as.numeric(unlist(strsplit(input$fast_decay_limit, ','))) else NULL
  
      result <- analyse_PSC(
        response = response_data(), dt = input$dt, n = input$n, N = input$N, IEI = input$IEI,
        stimulation_time = current_params$stimulation_time,
        baseline = current_params$baseline, smooth = input$smooth,
        fit.limits = input$user_tmax, fast.constraint = input$repeat_constraint,
        method = input$method, weight_method = input$weight_method, filter = input$filter, fc = input$fc,
        rel.decay.fit.limit = input$rel_decay_fit_limit,
        half_width_fit_limit = input$half_width_fit_limit, seed = input$seed,
        sequential.fit = input$sequential_fit, dp = input$dp,
        func = get(input$func), interval = c(input$interval_min, input$interval_max),
        lower = lower, upper = upper,
        MLEsettings = list(iter = input$iter, metropolis.scale = input$metropolis_scale, 
                           fit.attempts = input$fit_attempts, RWm = input$RWm),
        fast.decay.limit = fast.decay.limit, fast.constraint.method = input$fast_constraint_method,
        first.delay.constraint = input$first_delay_constraint, latency.limit = latency.limit,
        lwd = input$lwd, xlab = input$xlab, ylab = input$ylab,
        height = input$plot_height, width = input$plot_width,
        return.output = TRUE
      )
      analysis_output(result)
      prompt_fast_constraint(FALSE)
    })
    
    output$promptTmax <- reactive({ prompt_tmax() })
    output$promptFastConstraint <- reactive({ prompt_fast_constraint() })
    outputOptions(output, 'promptTmax', suspendWhenHidden = FALSE)
    outputOptions(output, 'promptFastConstraint', suspendWhenHidden = FALSE)
    
    output$download_output <- downloadHandler(
      filename = function() {
        paste0(tools::file_path_sans_ext(basename(input$file$name)), '_', input$data_col, '_PSC_analysis.rds')
      },
      content = function(file) {
        saveRDS(analysis_output(), file)
      }
    )
    
    # Clear button observer: only reset reactive data and set clear flag
    observeEvent(input$clear_output, {
      analysis_output(NULL)
      clear_flag(TRUE)
    })
  }
  
  shiny::runGadget(ui, server, viewer = shiny::dialogViewer("PSC Analysis", width = 1200, height = 800))
}

# To launch your PSC analysis widget from the R console, simply call:
PSC_analysis_widget()




# ######################################################################################
# rm(list = ls(all = TRUE))
# graphics.off()

# # Load and install necessary packages
# load_required_packages <- function(packages) {
#   new.packages <- packages[!(packages %in% installed.packages()[, 'Package'])]
#   if (length(new.packages)) install.packages(new.packages)
#   invisible(lapply(packages, library, character.only = TRUE))
# }
# required.packages <- c('ARTool', 'robustbase', 'minpack.lm', 'Rcpp', 'signal', 'dbscan')
# load_required_packages(required.packages)

# # Insert your username here to define the correct path
# username <- 'YourUsername'
# username <- 'euo9382'

# path_repository <- '/Documents/Repositories/Rfits'
# # Construct file path
# file_path1 <- paste0('/Users/', username, path_repository)

# source(paste0(file_path1, '/nNLS functions.R'))

# required.packages <- c('tkrplot', 'tcltk', 'readxl')
# load_required_packages(required.packages)

# determine_tmax2 <- function(y, 
#                             N = 1, 
#                             dt = 0.1, 
#                             stimulation_time = 0, 
#                             baseline = 0, 
#                             smooth = 5, 
#                             tmax = NULL, 
#                             y_abline = 0.1, 
#                             ylab = NULL, 
#                             height = 5, 
#                             width = 5) {
#   if (is.null(tmax)) {
#     # Calculate peak information (assuming peak.fun and abline_fun are defined)
#     peak <- peak.fun(y = y, dt = dt, stimulation_time = stimulation_time, baseline = baseline, smooth = smooth)
    
#     ind1 <- as.integer((stimulation_time - baseline) / dt)
#     ind2 <- as.integer(stimulation_time / dt)
#     y2plot <- y - mean(y[ind1:ind2])
    
#     # Draw on the active device (provided by tkrplot)
#     Y <- y2plot[ind1:length(y2plot)]
#     X <- seq(0, dt * (length(Y) - 1), by = dt)
    
#     out <- abline_fun(X, Y, N = N, y_abline = y_abline)
#     A_abline <- out[1]
#     avg_t.abline <- out[2]
#     avg_t.abline <- if (is.na(avg_t.abline)) max(X) else avg_t.abline
    
#     plot(X, Y, col = 'indianred', xlab = 'time (ms)', ylab = ylab, type = 'l', bty = 'l', las = 1, main = '')
#     abline(h = 0, col = 'black', lwd = 1, lty = 1)
    
#     left_axis <- par("usr")[1]
#     bottom_axis <- par("usr")[3]
    
#     lines(c(left_axis, avg_t.abline), c(A_abline, A_abline), col = 'black', lwd = 1, lty = 3)
#     lines(c(avg_t.abline, avg_t.abline), c(A_abline, bottom_axis), col = 'black', lwd = 1, lty = 3)
    
#     ind3 <- as.integer(avg_t.abline / dt)
#     text(x = max(X[ind1:ind3]) * 1.05, y = A_abline * 1.2, 
#          labels = paste0(y_abline * 100, ' %'), pos = 4, cex = 0.6)
#     text(x = max(X[ind1:ind3]) * 1.05, y = bottom_axis * 0.95, 
#          labels = paste0(avg_t.abline, ' ms'), pos = 4, cex = 0.6)
    
#     stim_index <- round(baseline / dt) + 1
#     if (stim_index > length(X)) stim_index <- length(X)
#     points(X[stim_index], Y[stim_index], pch = 8, col = 'darkgray', cex = 1)
#     x_offset <- 0.02 * diff(range(X))
#     text(x = X[stim_index] + x_offset, y = Y[stim_index], labels = "stim", 
#          pos = 4, col = 'darkgray', cex = 0.6)
    
#     x_limit <- avg_t.abline
#   } else {
#     x_limit <- tmax
#   }
  
#   x_limit <- x_limit + stimulation_time - baseline
#   return(x_limit)
# }

# # Modified fit_plot2: Notice that no new device is opened.
# fit_plot2 <- function(traces, func = product2, 
#                      xlab = 'time (ms)', ylab = 'PSC amplitude (pA)', 
#                      xlim = NULL, ylim = NULL, main = '', bl = NULL, 
#                      lwd = 1.2, filter = FALSE, width = 5, height = 5) {
  
#   plot(traces$x, traces$y, col = 'gray', xlab = xlab, ylab = ylab, 
#        xlim = xlim, ylim = ylim, type = 'l', bty = 'l', las = 1, 
#        lwd = lwd, main = main)
  
#   if (filter) {
#     lines(traces$x, traces$yfilter, col = 'black', type = 'l', lwd = lwd)
#   }
  
#   lines(traces$x, traces$yfit, col = 'indianred', lty = 3, lwd = 2 * lwd)
  
#   if (identical(func, product2) || identical(func, product2N)) {
#     lines(traces$x, traces$yfit1, col = '#4C78BC', lty = 3, lwd = 2 * lwd)
#     lines(traces$x, traces$yfit2, col = '#CA92C1', lty = 3, lwd = 2 * lwd)
#   }
  
#   if (identical(func, product3) || identical(func, product3N)) {
#     lines(traces$x, traces$yfit1, col = '#F28E2B', lty = 3, lwd = 2 * lwd)
#     lines(traces$x, traces$yfit2, col = '#4C78BC', lty = 3, lwd = 2 * lwd)
#     lines(traces$x, traces$yfit3, col = '#CA92C1', lty = 3, lwd = 2 * lwd)
#   }
  
#   if (!is.null(bl)) {
#     abline(v = bl, col = 'black', lwd = lwd, lty = 3)
#   }
# }

# # Define drawPlot2 to be used with tkrplot so that the GUI updates with the fit results.
# drawPlot2 <- function(traces, func = product2, xlab = "time (ms)", ylab = "PSC amplitude (pA)",
#                       xlim = NULL, ylim = NULL, main = "", bl = NULL, lwd = 1.2,
#                       filter = FALSE, width = 5, height = 5) {
#   fit_plot2(traces = traces, func = func, xlab = xlab, ylab = ylab,
#             xlim = xlim, ylim = ylim, main = main, bl = bl,
#             lwd = lwd, filter = filter, width = width, height = height)
# }


# PSC_analysis_tk <- function() {
#   # Create the main top-level window
#   tt <- tktoplevel()
#   tkwm.title(tt, "PSC Analysis")
  
#   # Divide the window into a sidebar (for controls) and a main panel (for plot + console)
#   sidebarFrame <- tkframe(tt)
#   mainFrame <- tkframe(tt)
#   tkgrid(sidebarFrame, row = 0, column = 0, sticky = "ns")
#   tkgrid(mainFrame, row = 0, column = 1, sticky = "nsew")
#   tkgrid.rowconfigure(tt, 0, weight = 1)
#   tkgrid.columnconfigure(tt, 1, weight = 1)
  
#   ### Sidebar Controls ###
#   fileLabel <- tklabel(sidebarFrame, text = "Upload CSV or XLSX:")
#   tkgrid(fileLabel, row = 0, column = 0, sticky = "w")
#   filePathVar <- tclVar("")
#   fileEntry <- tkentry(sidebarFrame, textvariable = filePathVar, width = 30)
#   tkgrid(fileEntry, row = 0, column = 1, sticky = "w")
#   browseButton <- tkbutton(sidebarFrame, text = "Browse", 
#                            command = function() {
#                              filePath <- tclvalue(tkgetOpenFile(filetypes = "{{CSV Files} {.csv}} {{Excel Files} {.xlsx .xls}}"))
#                              if (nchar(filePath) > 0) {
#                                tclvalue(filePathVar) <- filePath
#                                ext <- tools::file_ext(filePath)
#                                if (tolower(ext) == "csv") {
#                                  uploaded_data <<- read.csv(filePath)
#                                } else {
#                                  uploaded_data <<- readxl::read_excel(filePath)
#                                }
#                                columns <<- names(uploaded_data)
#                                tkconfigure(columnCombo, values = columns)
#                              }
#                            })
#   tkgrid(browseButton, row = 0, column = 2, padx = 5)
  
#   colLabel <- tklabel(sidebarFrame, text = "Select column:")
#   tkgrid(colLabel, row = 1, column = 0, sticky = "w")
#   columnVar <- tclVar("")
#   columnCombo <- ttkcombobox(sidebarFrame, textvariable = columnVar, values = "", width = 20)
#   tkgrid(columnCombo, row = 1, column = 1, columnspan = 2, sticky = "w")
  
#   nb <- ttknotebook(sidebarFrame)
#   tkgrid(nb, row = 2, column = 0, columnspan = 3, pady = 5, sticky = "nsew")
  
#   mainOptionsFrame   <- tkframe(nb)
#   fitOptionsFrame    <- tkframe(nb)
#   mleSettingsFrame   <- tkframe(nb)
#   advancedFrame      <- tkframe(nb)
#   graphSettingsFrame <- tkframe(nb)
  
#   tkadd(nb, mainOptionsFrame, text = "Main Options")
#   tkadd(nb, fitOptionsFrame, text = "Fit Options")
#   tkadd(nb, mleSettingsFrame, text = "MLE Settings")
#   tkadd(nb, advancedFrame, text = "Advanced")
#   tkadd(nb, graphSettingsFrame, text = "Graph Settings")
  
#   ### Main Options Tab ###
#   dtVar <- tclVar("0.1")
#   stimTimeVar <- tclVar("100")
#   baselineVar <- tclVar("50")
#   nVar <- tclVar("30")
#   yAblineVar <- tclVar("0.1")
#   funcVar <- tclVar("product1N")
  
#   tkgrid(tklabel(mainOptionsFrame, text = "dt (ms):"), row = 0, column = 0, sticky = "w")
#   tkgrid(tkentry(mainOptionsFrame, textvariable = dtVar, width = 10), row = 0, column = 1)
#   tkgrid(tklabel(mainOptionsFrame, text = "Stimulation Time:"), row = 1, column = 0, sticky = "w")
#   tkgrid(tkentry(mainOptionsFrame, textvariable = stimTimeVar, width = 10), row = 1, column = 1)
#   tkgrid(tklabel(mainOptionsFrame, text = "Baseline:"), row = 2, column = 0, sticky = "w")
#   tkgrid(tkentry(mainOptionsFrame, textvariable = baselineVar, width = 10), row = 2, column = 1)
#   tkgrid(tklabel(mainOptionsFrame, text = "n:"), row = 3, column = 0, sticky = "w")
#   tkgrid(tkentry(mainOptionsFrame, textvariable = nVar, width = 10), row = 3, column = 1)
#   tkgrid(tklabel(mainOptionsFrame, text = "Fit cutoff:"), row = 4, column = 0, sticky = "w")
#   tkgrid(tkentry(mainOptionsFrame, textvariable = yAblineVar, width = 10), row = 4, column = 1)
#   tkgrid(tklabel(mainOptionsFrame, text = "Function:"), row = 5, column = 0, sticky = "w")
#   funcChoices <- c("product1N", "product2N", "product3N")
#   funcCombo <- ttkcombobox(mainOptionsFrame, textvariable = funcVar, values = funcChoices, width = 10)
#   tkgrid(funcCombo, row = 5, column = 1)
  
#   ### Fit Options Tab ###
#   NVar <- tclVar("1")
#   IEIVar <- tclVar("50")
#   smoothVar <- tclVar("5")
#   methodVar <- tclVar("BF.LM")
#   weightMethodVar <- tclVar("none")
#   sequentialFitVar <- tclVar("0")
#   intervalMinVar <- tclVar("0.1")
#   intervalMaxVar <- tclVar("0.9")
#   lowerVar <- tclVar("")
#   upperVar <- tclVar("")
#   latencyLimitVar <- tclVar("")
  
#   tkgrid(tklabel(fitOptionsFrame, text = "N:"), row = 0, column = 0, sticky = "w")
#   tkgrid(tkentry(fitOptionsFrame, textvariable = NVar, width = 10), row = 0, column = 1)
#   tkgrid(tklabel(fitOptionsFrame, text = "IEI:"), row = 1, column = 0, sticky = "w")
#   tkgrid(tkentry(fitOptionsFrame, textvariable = IEIVar, width = 10), row = 1, column = 1)
#   tkgrid(tklabel(fitOptionsFrame, text = "Smooth:"), row = 2, column = 0, sticky = "w")
#   tkgrid(tkentry(fitOptionsFrame, textvariable = smoothVar, width = 10), row = 2, column = 1)
#   tkgrid(tklabel(fitOptionsFrame, text = "Method:"), row = 3, column = 0, sticky = "w")
#   methodChoices <- c("BF.LM", "LM", "GN", "port", "robust", "MLE")
#   methodCombo <- ttkcombobox(fitOptionsFrame, textvariable = methodVar, values = methodChoices, width = 10)
#   tkgrid(methodCombo, row = 3, column = 1)
#   tkgrid(tklabel(fitOptionsFrame, text = "Weighting:"), row = 4, column = 0, sticky = "w")
#   weightChoices <- c("none", "~y_sqrt", "~y")
#   weightCombo <- ttkcombobox(fitOptionsFrame, textvariable = weightMethodVar, values = weightChoices, width = 10)
#   tkgrid(weightCombo, row = 4, column = 1)
#   tkgrid(tklabel(fitOptionsFrame, text = "Sequential Fit:"), row = 5, column = 0, sticky = "w")
#   sequentialFitCheck <- tkcheckbutton(fitOptionsFrame, variable = sequentialFitVar)
#   tkgrid(sequentialFitCheck, row = 5, column = 1)
#   tkgrid(tklabel(fitOptionsFrame, text = "Min interval:"), row = 6, column = 0, sticky = "w")
#   tkgrid(tkentry(fitOptionsFrame, textvariable = intervalMinVar, width = 10), row = 6, column = 1)
#   tkgrid(tklabel(fitOptionsFrame, text = "Max interval:"), row = 7, column = 0, sticky = "w")
#   tkgrid(tkentry(fitOptionsFrame, textvariable = intervalMaxVar, width = 10), row = 7, column = 1)
#   tkgrid(tklabel(fitOptionsFrame, text = "Lower bounds:"), row = 8, column = 0, sticky = "w")
#   tkgrid(tkentry(fitOptionsFrame, textvariable = lowerVar, width = 10), row = 8, column = 1)
#   tkgrid(tklabel(fitOptionsFrame, text = "Upper bounds:"), row = 9, column = 0, sticky = "w")
#   tkgrid(tkentry(fitOptionsFrame, textvariable = upperVar, width = 10), row = 9, column = 1)
#   tkgrid(tklabel(fitOptionsFrame, text = "Latency limit:"), row = 10, column = 0, sticky = "w")
#   tkgrid(tkentry(fitOptionsFrame, textvariable = latencyLimitVar, width = 10), row = 10, column = 1)
  
#   ### MLE Settings Tab ###
#   iterVar <- tclVar("1000")
#   metropolisScaleVar <- tclVar("1.5")
#   fitAttemptsVar <- tclVar("10")
#   RWmVar <- tclVar("0")
  
#   tkgrid(tklabel(mleSettingsFrame, text = "MLE Iterations:"), row = 0, column = 0, sticky = "w")
#   tkgrid(tkentry(mleSettingsFrame, textvariable = iterVar, width = 10), row = 0, column = 1)
#   tkgrid(tklabel(mleSettingsFrame, text = "Metropolis Scale:"), row = 1, column = 0, sticky = "w")
#   tkgrid(tkentry(mleSettingsFrame, textvariable = metropolisScaleVar, width = 10), row = 1, column = 1)
#   tkgrid(tklabel(mleSettingsFrame, text = "Fit Attempts:"), row = 2, column = 0, sticky = "w")
#   tkgrid(tkentry(mleSettingsFrame, textvariable = fitAttemptsVar, width = 10), row = 2, column = 1)
#   tkgrid(tklabel(mleSettingsFrame, text = "Random Walk Metropolis:"), row = 3, column = 0, sticky = "w")
#   RWmCheck <- tkcheckbutton(mleSettingsFrame, variable = RWmVar)
#   tkgrid(RWmCheck, row = 3, column = 1)
  
#   ### Advanced Tab ###
#   filterVar <- tclVar("0")
#   fcVar <- tclVar("1000")
#   relDecayFitLimitVar <- tclVar("0.1")
#   halfWidthFitLimitVar <- tclVar("500")
#   seedVar <- tclVar("42")
#   dpVar <- tclVar("3")
#   fastConstraintVar <- tclVar("0")
#   fastConstraintMethodVar <- tclVar("rise")
#   fastDecayLimitVar <- tclVar("")
#   firstDelayConstraintVar <- tclVar("0")
  
#   tkgrid(tklabel(advancedFrame, text = "Filter:"), row = 0, column = 0, sticky = "w")
#   filterCheck <- tkcheckbutton(advancedFrame, variable = filterVar)
#   tkgrid(filterCheck, row = 0, column = 1)
#   tkgrid(tklabel(advancedFrame, text = "Filter cutoff (Hz):"), row = 1, column = 0, sticky = "w")
#   tkgrid(tkentry(advancedFrame, textvariable = fcVar, width = 10), row = 1, column = 1)
#   tkgrid(tklabel(advancedFrame, text = "Relative decay fit limit:"), row = 2, column = 0, sticky = "w")
#   tkgrid(tkentry(advancedFrame, textvariable = relDecayFitLimitVar, width = 10), row = 2, column = 1)
#   tkgrid(tklabel(advancedFrame, text = "Half-width fit limit:"), row = 3, column = 0, sticky = "w")
#   tkgrid(tkentry(advancedFrame, textvariable = halfWidthFitLimitVar, width = 10), row = 3, column = 1)
#   tkgrid(tklabel(advancedFrame, text = "Seed:"), row = 4, column = 0, sticky = "w")
#   tkgrid(tkentry(advancedFrame, textvariable = seedVar, width = 10), row = 4, column = 1)
#   tkgrid(tklabel(advancedFrame, text = "Decimal points:"), row = 5, column = 0, sticky = "w")
#   tkgrid(tkentry(advancedFrame, textvariable = dpVar, width = 10), row = 5, column = 1)
#   tkgrid(tklabel(advancedFrame, text = "Fast constraint:"), row = 6, column = 0, sticky = "w")
#   fastConstraintCheck <- tkcheckbutton(advancedFrame, variable = fastConstraintVar)
#   tkgrid(fastConstraintCheck, row = 6, column = 1)
#   tkgrid(tklabel(advancedFrame, text = "Fast constraint method:"), row = 7, column = 0, sticky = "w")
#   fastConstraintChoices <- c("rise", "peak")
#   fastConstraintCombo <- ttkcombobox(advancedFrame, textvariable = fastConstraintMethodVar, values = fastConstraintChoices, width = 10)
#   tkgrid(fastConstraintCombo, row = 7, column = 1)
#   tkgrid(tklabel(advancedFrame, text = "Fast decay limit(s):"), row = 8, column = 0, sticky = "w")
#   tkgrid(tkentry(advancedFrame, textvariable = fastDecayLimitVar, width = 10), row = 8, column = 1)
#   tkgrid(tklabel(advancedFrame, text = "First delay constraint:"), row = 9, column = 0, sticky = "w")
#   firstDelayCheck <- tkcheckbutton(advancedFrame, variable = firstDelayConstraintVar)
#   tkgrid(firstDelayCheck, row = 9, column = 1)
  
#   ### Graph Settings Tab ###
#   lwdVar <- tclVar("1.2")
#   xlabVar <- tclVar("time (ms)")
#   ylabVar <- tclVar("PSC (pA)")
#   plotHeightVar <- tclVar("5")
#   plotWidthVar <- tclVar("5")
  
#   tkgrid(tklabel(graphSettingsFrame, text = "Line width:"), row = 0, column = 0, sticky = "w")
#   tkgrid(tkentry(graphSettingsFrame, textvariable = lwdVar, width = 10), row = 0, column = 1)
#   tkgrid(tklabel(graphSettingsFrame, text = "x-axis label:"), row = 1, column = 0, sticky = "w")
#   tkgrid(tkentry(graphSettingsFrame, textvariable = xlabVar, width = 10), row = 1, column = 1)
#   tkgrid(tklabel(graphSettingsFrame, text = "y-axis label:"), row = 2, column = 0, sticky = "w")
#   tkgrid(tkentry(graphSettingsFrame, textvariable = ylabVar, width = 10), row = 2, column = 1)
#   tkgrid(tklabel(graphSettingsFrame, text = "Plot height:"), row = 3, column = 0, sticky = "w")
#   tkgrid(tkentry(graphSettingsFrame, textvariable = plotHeightVar, width = 10), row = 3, column = 1)
#   tkgrid(tklabel(graphSettingsFrame, text = "Plot width:"), row = 4, column = 0, sticky = "w")
#   tkgrid(tkentry(graphSettingsFrame, textvariable = plotWidthVar, width = 10), row = 4, column = 1)
  
#   ## Additional sidebar controls for analysis actions:
#   userTmaxVar <- tclVar("")
#   tkgrid(tklabel(sidebarFrame, text = "User Tmax:"), row = 3, column = 0, sticky = "w", pady = 5)
#   tkgrid(tkentry(sidebarFrame, textvariable = userTmaxVar, width = 10), row = 3, column = 1, pady = 5)
  
#   repeatConstraintVar <- tclVar("0")
#   tkgrid(tklabel(sidebarFrame, text = "Repeat with fast constraint:"), row = 4, column = 0, sticky = "w")
#   repeatConstraintCheck <- tkcheckbutton(sidebarFrame, variable = repeatConstraintVar)
#   tkgrid(repeatConstraintCheck, row = 4, column = 1)
  
#   ## Analysis action buttons
  
#   runAnalysisButton <- tkbutton(sidebarFrame, text = "Run Initial Analysis",
#     command = function() {
#       filePath <- tclvalue(filePathVar)
#       if (nchar(filePath) == 0) {
#         tkmessageBox(message = "Please select a file first.")
#         return()
#       }
#       if (nchar(tclvalue(columnVar)) == 0) {
#         tkmessageBox(message = "Please select a column.")
#         return()
#       }
#       ext <- tools::file_ext(filePath)
#       if (tolower(ext) == "csv") {
#         uploaded_data <<- read.csv(filePath)
#       } else {
#         uploaded_data <<- readxl::read_excel(filePath)
#       }
#       response_data <<- uploaded_data[[tclvalue(columnVar)]]
#       tkrreplot(plotWidget, fun = drawPlot1)
#     })
#   tkgrid(runAnalysisButton, row = 5, column = 0, columnspan = 3, pady = 5)
  
#   updatePlotButton <- tkbutton(sidebarFrame, text = "Update Plot",
#     command = function() {
#       tkrreplot(plotWidget, fun = drawPlot1)
#     })
#   tkgrid(updatePlotButton, row = 6, column = 0, columnspan = 3, pady = 5)
  
#   ## Run Main Analysis (inline version of analyse_PSC)
#   runMainAnalysisButton <- tkbutton(sidebarFrame, text = "Run Main Analysis",
#     command = function() {
#       dt                <- as.numeric(tclvalue(dtVar))
#       stimulation_time  <- as.numeric(tclvalue(stimTimeVar))
#       baseline          <- as.numeric(tclvalue(baselineVar))
#       smooth            <- as.numeric(tclvalue(smoothVar))
#       n                 <- as.numeric(tclvalue(nVar))
#       N                 <- as.numeric(tclvalue(NVar))
#       IEI               <- as.numeric(tclvalue(IEIVar))
#       func              <- get(tclvalue(funcVar))
#       method            <- tclvalue(methodVar)
#       weight_method     <- tclvalue(weightMethodVar)
#       sequential.fit    <- as.logical(as.numeric(tclvalue(sequentialFitVar)))
#       fit.limits        <- as.numeric(tclvalue(userTmaxVar))
#       rel.decay.fit.limit <- as.numeric(tclvalue(relDecayFitLimitVar))
#       lwd               <- as.numeric(tclvalue(lwdVar))
#       xlab              <- tclvalue(xlabVar)
#       ylab              <- tclvalue(ylabVar)
#       plotWidth         <- as.numeric(tclvalue(plotWidthVar))
#       plotHeight        <- as.numeric(tclvalue(plotHeightVar))
      
#       fc                <- as.numeric(tclvalue(fcVar))
#       interval          <- c(as.numeric(tclvalue(intervalMinVar)),
#                              as.numeric(tclvalue(intervalMaxVar)))
#       lower             <- if(nchar(tclvalue(lowerVar)) > 0)
#                              as.numeric(unlist(strsplit(tclvalue(lowerVar), ",")))
#                            else NULL
#       upper             <- if(nchar(tclvalue(upperVar)) > 0)
#                              as.numeric(unlist(strsplit(tclvalue(upperVar), ",")))
#                            else NULL
#       iter              <- as.numeric(tclvalue(iterVar))
#       metropolis.scale  <- as.numeric(tclvalue(metropolisScaleVar))
#       fit.attempts      <- as.numeric(tclvalue(fitAttemptsVar))
#       RWm               <- as.logical(as.numeric(tclvalue(RWmVar)))
#       fast.decay.limit  <- if(nchar(tclvalue(fastDecayLimitVar)) > 0)
#                              as.numeric(unlist(strsplit(tclvalue(fastDecayLimitVar), ",")))
#                            else NULL
#       fast.constraint   <- FALSE
#       fast.constraint.method <- tclvalue(fastConstraintMethodVar)
#       first.delay.constraint <- as.logical(as.numeric(tclvalue(firstDelayConstraintVar)))
#       dp                <- as.numeric(tclvalue(dpVar))
#       seed              <- as.numeric(tclvalue(seedVar))
#       filter            <- as.logical(as.numeric(tclvalue(filterVar)))
      
#       y <- response_data
#       if (all(is.na(y[(which(!is.na(y))[length(which(!is.na(y)))] + 1):length(y)]))) {
#         y <- y[!is.na(y)]
#       }
#       x <- seq(0, (length(y) - 1) * dt, by = dt)
      
#       if (!sequential.fit) {
#         tmax <- fit.limits
#         x_limit <- determine_tmax2(
#           y = y,
#           N = N,
#           dt = dt,
#           stimulation_time = stimulation_time,
#           baseline = baseline,
#           smooth = smooth,
#           tmax = tmax,
#           y_abline = rel.decay.fit.limit,
#           ylab = ylab,
#           width = plotWidth,
#           height = plotHeight
#         )
#         adjusted_response <- y[x < x_limit]
        
#         out <- nFIT(
#           response = adjusted_response,
#           n = n,
#           N = N,
#           IEI = IEI,
#           dt = dt,
#           func = func,
#           method = method,
#           weight_method = weight_method,
#           MLEsettings = list(
#             iter = iter,
#             metropolis.scale = metropolis.scale,
#             fit.attempts = fit.attempts,
#             RWm = RWm
#           ),
#           stimulation_time = stimulation_time,
#           baseline = baseline,
#           filter = filter,
#           fc = fc,
#           interval = interval,
#           fast.decay.limit = fast.decay.limit,
#           fast.constraint = fast.constraint,
#           fast.constraint.method = fast.constraint.method,
#           first.delay.constraint = first.delay.constraint,
#           lower = lower,
#           upper = upper,
#           latency.limit = if(nchar(tclvalue(latencyLimitVar)) > 0)
#                             as.numeric(unlist(strsplit(tclvalue(latencyLimitVar), ",")))
#                           else NULL,
#           return.output = TRUE,
#           show.plot = FALSE,
#           half_width_fit_limit = as.numeric(tclvalue(halfWidthFitLimitVar)),
#           dp = dp,
#           height = plotHeight,
#           width = plotWidth,
#           seed = seed
#         )
        
#         out$traces <- traces_fun2(
#           y = y,
#           fits = out$fits,
#           dt = dt,
#           N = N,
#           IEI = IEI,
#           stimulation_time = stimulation_time,
#           baseline = baseline,
#           func = func,
#           filter = filter,
#           fc = fc
#         )
#         # Update the tkrplot using drawPlot2 so that the fit is shown in the GUI
#         tkrreplot(plotWidget, fun = function() {
#           drawPlot2(traces = out$traces, func = func, xlab = xlab, ylab = ylab,
#                     lwd = lwd, filter = filter, width = plotWidth, height = plotHeight)
#         })
#       } else {
#         out <- nFIT_sequential(
#           response = y,
#           n = n,
#           dt = dt,
#           func = func,
#           method = method,
#           weight_method = weight_method,
#           stimulation_time = stimulation_time,
#           baseline = baseline,
#           fit.limits = fit.limits,
#           fast.decay.limit = fast.decay.limit,
#           fast.constraint = as.logical(as.numeric(tclvalue(fastConstraintVar))),
#           fast.constraint.method = fast.constraint.method,
#           first.delay.constraint = first.delay.constraint,
#           latency.limit = if(nchar(tclvalue(latencyLimitVar)) > 0)
#                             as.numeric(unlist(strsplit(tclvalue(latencyLimitVar), ",")))
#                           else NULL,
#           lower = lower,
#           upper = upper,
#           filter = filter,
#           fc = fc,
#           interval = interval,
#           MLEsettings = list(
#             iter = iter,
#             metropolis.scale = metropolis.scale,
#             fit.attempts = fit.attempts,
#             RWm = RWm
#           ),
#           MLE.method = method,
#           half_width_fit_limit = as.numeric(tclvalue(halfWidthFitLimitVar)),
#           dp = dp,
#           lwd = lwd,
#           xlab = xlab,
#           ylab = ylab,
#           width = plotWidth,
#           height = plotHeight,
#           return.output = TRUE,
#           show.output = TRUE,
#           show.plot = TRUE,
#           seed = seed
#         )
#       }
      
#       analysis_output <<- out
#       tkdelete(consoleText, "1.0", "end")
#       tkinsert(consoleText, "end", paste(capture.output(print(out)), collapse = "\n"))
#     }
#   )
#   tkgrid(runMainAnalysisButton, row = 7, column = 0, columnspan = 3, pady = 5)
  
#   ## Run Final Analysis (similar modifications can be made here)
#   runFinalAnalysisButton <- tkbutton(sidebarFrame, text = "Run Final Analysis",
#     command = function() {
#       fast.constraint <- as.logical(as.numeric(tclvalue(repeatConstraintVar)))
      
#       dt                <- as.numeric(tclvalue(dtVar))
#       stimulation_time  <- as.numeric(tclvalue(stimTimeVar))
#       baseline          <- as.numeric(tclvalue(baselineVar))
#       smooth            <- as.numeric(tclvalue(smoothVar))
#       n                 <- as.numeric(tclvalue(nVar))
#       N                 <- as.numeric(tclvalue(NVar))
#       IEI               <- as.numeric(tclvalue(IEIVar))
#       func              <- get(tclvalue(funcVar))
#       method            <- tclvalue(methodVar)
#       weight_method     <- tclvalue(weightMethodVar)
#       sequential.fit    <- as.logical(as.numeric(tclvalue(sequentialFitVar)))
#       fit.limits        <- as.numeric(tclvalue(userTmaxVar))
#       rel.decay.fit.limit <- as.numeric(tclvalue(relDecayFitLimitVar))
#       lwd               <- as.numeric(tclvalue(lwdVar))
#       xlab              <- tclvalue(xlabVar)
#       ylab              <- tclvalue(ylabVar)
#       plotWidth         <- as.numeric(tclvalue(plotWidthVar))
#       plotHeight        <- as.numeric(tclvalue(plotHeightVar))
      
#       fc                <- as.numeric(tclvalue(fcVar))
#       interval          <- c(as.numeric(tclvalue(intervalMinVar)),
#                              as.numeric(tclvalue(intervalMaxVar)))
#       lower             <- if(nchar(tclvalue(lowerVar)) > 0)
#                              as.numeric(unlist(strsplit(tclvalue(lowerVar), ",")))
#                            else NULL
#       upper             <- if(nchar(tclvalue(upperVar)) > 0)
#                              as.numeric(unlist(strsplit(tclvalue(upperVar), ",")))
#                            else NULL
#       iter              <- as.numeric(tclvalue(iterVar))
#       metropolis.scale  <- as.numeric(tclvalue(metropolisScaleVar))
#       fit.attempts      <- as.numeric(tclvalue(fitAttemptsVar))
#       RWm               <- as.logical(as.numeric(tclvalue(RWmVar)))
#       fast.decay.limit  <- if(nchar(tclvalue(fastDecayLimitVar)) > 0)
#                              as.numeric(unlist(strsplit(tclvalue(fastDecayLimitVar), ",")))
#                            else NULL
#       fast.constraint.method <- tclvalue(fastConstraintMethodVar)
#       first.delay.constraint <- as.logical(as.numeric(tclvalue(firstDelayConstraintVar)))
#       dp                <- as.numeric(tclvalue(dpVar))
#       seed              <- as.numeric(tclvalue(seedVar))
#       filter            <- as.logical(as.numeric(tclvalue(filterVar)))
      
#       y <- response_data
#       if (all(is.na(y[(which(!is.na(y))[length(which(!is.na(y)))] + 1):length(y)]))) {
#         y <- y[!is.na(y)]
#       }
#       x <- seq(0, (length(y) - 1) * dt, by = dt)
      
#       if (!sequential.fit) {
#         tmax <- fit.limits
#         x_limit <- determine_tmax2(
#           y = y,
#           N = N,
#           dt = dt,
#           stimulation_time = stimulation_time,
#           baseline = baseline,
#           smooth = smooth,
#           tmax = tmax,
#           y_abline = rel.decay.fit.limit,
#           ylab = ylab,
#           width = plotWidth,
#           height = plotHeight
#         )
#         adjusted_response <- y[x < x_limit]
        
#         out <- nFIT(
#           response = adjusted_response,
#           n = n,
#           N = N,
#           IEI = IEI,
#           dt = dt,
#           func = func,
#           method = method,
#           weight_method = weight_method,
#           MLEsettings = list(
#             iter = iter,
#             metropolis.scale = metropolis.scale,
#             fit.attempts = fit.attempts,
#             RWm = RWm
#           ),
#           stimulation_time = stimulation_time,
#           baseline = baseline,
#           filter = filter,
#           fc = fc,
#           interval = interval,
#           fast.decay.limit = fast.decay.limit,
#           fast.constraint = fast.constraint,
#           fast.constraint.method = fast.constraint.method,
#           first.delay.constraint = first.delay.constraint,
#           lower = lower,
#           upper = upper,
#           latency.limit = if(nchar(tclvalue(latencyLimitVar)) > 0)
#                             as.numeric(unlist(strsplit(tclvalue(latencyLimitVar), ",")))
#                           else NULL,
#           return.output = TRUE,
#           show.plot = FALSE,
#           half_width_fit_limit = as.numeric(tclvalue(halfWidthFitLimitVar)),
#           dp = dp,
#           height = plotHeight,
#           width = plotWidth,
#           seed = seed
#         )
#         out$traces <- traces_fun2(
#           y = y,
#           fits = out$fits,
#           dt = dt,
#           N = N,
#           IEI = IEI,
#           stimulation_time = stimulation_time,
#           baseline = baseline,
#           func = func,
#           filter = filter,
#           fc = fc
#         )
#         tkrreplot(plotWidget, fun = function() {
#           drawPlot2(traces = out$traces, func = func, xlab = xlab, ylab = ylab,
#                     lwd = lwd, filter = filter, width = plotWidth, height = plotHeight)
#         })
#       } else {
#         out <- nFIT_sequential(
#           response = y,
#           n = n,
#           dt = dt,
#           func = func,
#           method = method,
#           weight_method = weight_method,
#           stimulation_time = stimulation_time,
#           baseline = baseline,
#           fit.limits = fit.limits,
#           fast.decay.limit = fast.decay.limit,
#           fast.constraint = fast.constraint,
#           fast.constraint.method = fast.constraint.method,
#           first.delay.constraint = first.delay.constraint,
#           latency.limit = if(nchar(tclvalue(latencyLimitVar)) > 0)
#                             as.numeric(unlist(strsplit(tclvalue(latencyLimitVar), ",")))
#                           else NULL,
#           lower = lower,
#           upper = upper,
#           filter = filter,
#           fc = fc,
#           interval = interval,
#           MLEsettings = list(
#             iter = iter,
#             metropolis.scale = metropolis.scale,
#             fit.attempts = fit.attempts,
#             RWm = RWm
#           ),
#           MLE.method = method,
#           half_width_fit_limit = as.numeric(tclvalue(halfWidthFitLimitVar)),
#           dp = dp,
#           lwd = lwd,
#           xlab = xlab,
#           ylab = ylab,
#           width = plotWidth,
#           height = plotHeight,
#           return.output = TRUE,
#           show.output = TRUE,
#           show.plot = TRUE,
#           seed = seed
#         )
#       }
      
#       analysis_output <<- out
#       tkdelete(consoleText, "1.0", "end")
#       tkinsert(consoleText, "end", paste(capture.output(print(out)), collapse = "\n"))
#     }
#   )
#   tkgrid(runFinalAnalysisButton, row = 8, column = 0, columnspan = 3, pady = 5)

#   clearOutputButton <- tkbutton(sidebarFrame, text = "Clear Output",
#     command = function() {
#       analysis_output <<- NULL
#       tkdelete(consoleText, "1.0", "end")
#       tkrreplot(plotWidget, fun = drawPlot1)
#     }
#   )
#   tkgrid(clearOutputButton, row = 9, column = 0, columnspan = 3, pady = 5)
  
#   ### Main Panel: Plot and Console ###
#   drawPlot1 <- function() {
#     dt <- as.numeric(tclvalue(dtVar))
#     stimTime <- as.numeric(tclvalue(stimTimeVar))
#     baseline <- as.numeric(tclvalue(baselineVar))
#     smooth <- as.numeric(tclvalue(smoothVar))
#     y_abline <- as.numeric(tclvalue(yAblineVar))
    
#     y_val <- if (exists("response_data") && !is.null(response_data)) {
#       response_data
#     } else {
#       rnorm(10000, 0.1)
#     }
    
#     determine_tmax2(
#       y = y_val,
#       N = 1,
#       dt = dt,
#       stimulation_time = stimTime,
#       baseline = baseline,
#       smooth = smooth,
#       tmax = NULL,
#       y_abline = y_abline,
#       ylab = tclvalue(ylabVar),
#       height = 5,
#       width = 5
#     )
#   }
  
#   # Initially, the plot is drawn using drawPlot1.
#   plotWidget <- tkrplot(tt, fun = drawPlot1, hscale = 1.3, vscale = 1.3)
#   tkgrid(plotWidget, row = 0, column = 1, sticky = "nsew")
#   consoleText <- tktext(mainFrame, width = 80, height = 10)
#   tkgrid(consoleText, row = 1, column = 1, sticky = "nsew")
  
#   tkfocus(tt)
# }

# # Launch the PSC Analysis interface
# PSC_analysis_tk()





######################################################################################
# for xquartz to work properly in (some) systems open R from terminal:
# open -n -a R

rm(list = ls(all = TRUE))
graphics.off()

# Load and install necessary packages
load_required_packages <- function(packages) {
  new.packages <- packages[!(packages %in% installed.packages()[, 'Package'])]
  if (length(new.packages)) install.packages(new.packages)
  invisible(lapply(packages, library, character.only = TRUE))
}
required.packages <- c('ARTool', 'robustbase', 'minpack.lm', 'Rcpp', 'signal', 'dbscan')
load_required_packages(required.packages)

# Insert your username here to define the correct path
username <- 'YourUsername'
username <- 'euo9382'

path_repository <- '/Documents/Repositories/Rfits'
# Construct file path
file_path1 <- paste0('/Users/', username, path_repository)

source(paste0(file_path1, '/nNLS functions.R'))

required.packages <- c('tkrplot', 'tcltk', 'readxl')
load_required_packages(required.packages)

determine_tmax2 <- function(y, 
                            N = 1, 
                            dt = 0.1, 
                            stimulation_time = 0, 
                            baseline = 0, 
                            smooth = 5, 
                            tmax = NULL, 
                            y_abline = 0.1, 
                            ylab = NULL, 
                            height = 5, 
                            width = 5) {
  if (is.null(tmax)) {
    # Calculate peak information (assuming peak.fun and abline_fun are defined)
    peak <- peak.fun(y = y, dt = dt, stimulation_time = stimulation_time, baseline = baseline, smooth = smooth)
    
    ind1 <- as.integer((stimulation_time - baseline) / dt)
    ind2 <- as.integer(stimulation_time / dt)
    y2plot <- y - mean(y[ind1:ind2])
    
    # Draw on the active device (provided by tkrplot)
    Y <- y2plot[ind1:length(y2plot)]
    X <- seq(0, dt * (length(Y) - 1), by = dt)
    
    out <- abline_fun(X, Y, N = N, y_abline = y_abline)
    A_abline <- out[1]
    avg_t.abline <- out[2]
    avg_t.abline <- if (is.na(avg_t.abline)) max(X) else avg_t.abline
    
    plot(X, Y, col = 'indianred', xlab = 'time (ms)', ylab = ylab, type = 'l', bty = 'l', las = 1, main = '')
    abline(h = 0, col = 'black', lwd = 1, lty = 1)
    
    left_axis <- par("usr")[1]
    bottom_axis <- par("usr")[3]
    
    lines(c(left_axis, avg_t.abline), c(A_abline, A_abline), col = 'black', lwd = 1, lty = 3)
    lines(c(avg_t.abline, avg_t.abline), c(A_abline, bottom_axis), col = 'black', lwd = 1, lty = 3)
    
    ind3 <- as.integer(avg_t.abline / dt)
    text(x = max(X[ind1:ind3]) * 1.05, y = A_abline * 1.2, 
         labels = paste0(y_abline * 100, ' %'), pos = 4, cex = 0.6)
    text(x = max(X[ind1:ind3]) * 1.05, y = bottom_axis * 0.95, 
         labels = paste0(avg_t.abline, ' ms'), pos = 4, cex = 0.6)
    
    stim_index <- round(baseline / dt) + 1
    if (stim_index > length(X)) stim_index <- length(X)
    points(X[stim_index], Y[stim_index], pch = 8, col = 'darkgray', cex = 1)
    x_offset <- 0.02 * diff(range(X))
    text(x = X[stim_index] + x_offset, y = Y[stim_index], labels = "stim", 
         pos = 4, col = 'darkgray', cex = 0.6)
    
    x_limit <- avg_t.abline
  } else {
    x_limit <- tmax
  }
  
  x_limit <- x_limit + stimulation_time - baseline
  return(x_limit)
}

# Modified fit_plot2: Notice that no new device is opened.
fit_plot2 <- function(traces, func = product2, 
                     xlab = 'time (ms)', ylab = 'PSC amplitude (pA)', 
                     xlim = NULL, ylim = NULL, main = '', bl = NULL, 
                     lwd = 1.2, filter = FALSE, width = 5, height = 5) {
  
  plot(traces$x, traces$y, col = 'gray', xlab = xlab, ylab = ylab, 
       xlim = xlim, ylim = ylim, type = 'l', bty = 'l', las = 1, 
       lwd = lwd, main = main)
  
  if (filter) {
    lines(traces$x, traces$yfilter, col = 'black', type = 'l', lwd = lwd)
  }
  
  lines(traces$x, traces$yfit, col = 'indianred', lty = 3, lwd = 2 * lwd)
  
  if (identical(func, product2) || identical(func, product2N)) {
    lines(traces$x, traces$yfit1, col = '#4C78BC', lty = 3, lwd = 2 * lwd)
    lines(traces$x, traces$yfit2, col = '#CA92C1', lty = 3, lwd = 2 * lwd)
  }
  
  if (identical(func, product3) || identical(func, product3N)) {
    lines(traces$x, traces$yfit1, col = '#F28E2B', lty = 3, lwd = 2 * lwd)
    lines(traces$x, traces$yfit2, col = '#4C78BC', lty = 3, lwd = 2 * lwd)
    lines(traces$x, traces$yfit3, col = '#CA92C1', lty = 3, lwd = 2 * lwd)
  }
  
  if (!is.null(bl)) {
    abline(v = bl, col = 'black', lwd = lwd, lty = 3)
  }
}

# Define drawPlot2 to be used with tkrplot so that the GUI updates with the fit results.
drawPlot2 <- function(traces, func = product2, xlab = "time (ms)", ylab = "PSC amplitude (pA)",
                      xlim = NULL, ylim = NULL, main = "", bl = NULL, lwd = 1.2,
                      filter = FALSE, width = 5, height = 5) {
  fit_plot2(traces = traces, func = func, xlab = xlab, ylab = ylab,
            xlim = xlim, ylim = ylim, main = main, bl = bl,
            lwd = lwd, filter = filter, width = width, height = height)
}

PSC_analysis_tk <- function() {
  # Create the main top-level window
  tt <- tktoplevel()
  tkwm.title(tt, "PSC Analysis")
  
  # Divide the window into a sidebar (for controls) and a main panel (for plot + console)
  sidebarFrame <- tkframe(tt)
  mainFrame <- tkframe(tt)
  tkgrid(sidebarFrame, row = 0, column = 0, sticky = "ns")
  tkgrid(mainFrame, row = 0, column = 1, sticky = "nsew")
  tkgrid.rowconfigure(tt, 0, weight = 1)
  tkgrid.columnconfigure(tt, 1, weight = 1)
  
  ### Sidebar Controls ###
  fileLabel <- tklabel(sidebarFrame, text = "Upload CSV or XLSX:")
  tkgrid(fileLabel, row = 0, column = 0, sticky = "w")
  filePathVar <- tclVar("")
  fileEntry <- tkentry(sidebarFrame, textvariable = filePathVar, width = 30)
  tkgrid(fileEntry, row = 0, column = 1, sticky = "w")
  browseButton <- tkbutton(sidebarFrame, text = "Browse", 
                           command = function() {
                             filePath <- tclvalue(tkgetOpenFile(filetypes = "{{CSV Files} {.csv}} {{Excel Files} {.xlsx .xls}}"))
                             if (nchar(filePath) > 0) {
                               tclvalue(filePathVar) <- filePath
                               ext <- tools::file_ext(filePath)
                               if (tolower(ext) == "csv") {
                                 uploaded_data <<- read.csv(filePath)
                               } else {
                                 uploaded_data <<- readxl::read_excel(filePath)
                               }
                               columns <<- names(uploaded_data)
                               tkconfigure(columnCombo, values = columns)
                             }
                           })
  tkgrid(browseButton, row = 0, column = 2, padx = 5)
  
  colLabel <- tklabel(sidebarFrame, text = "Select column:")
  tkgrid(colLabel, row = 1, column = 0, sticky = "w")
  columnVar <- tclVar("")
  columnCombo <- ttkcombobox(sidebarFrame, textvariable = columnVar, values = "", width = 20)
  tkgrid(columnCombo, row = 1, column = 1, columnspan = 2, sticky = "w")
  
  nb <- ttknotebook(sidebarFrame)
  tkgrid(nb, row = 2, column = 0, columnspan = 3, pady = 5, sticky = "nsew")
  
  mainOptionsFrame   <- tkframe(nb)
  fitOptionsFrame    <- tkframe(nb)
  mleSettingsFrame   <- tkframe(nb)
  advancedFrame      <- tkframe(nb)
  graphSettingsFrame <- tkframe(nb)
  
  tkadd(nb, mainOptionsFrame, text = "Main Options")
  tkadd(nb, fitOptionsFrame, text = "Fit Options")
  tkadd(nb, mleSettingsFrame, text = "MLE Settings")
  tkadd(nb, advancedFrame, text = "Advanced")
  tkadd(nb, graphSettingsFrame, text = "Graph Settings")
  
  ### Main Options Tab ###
  dtVar <- tclVar("0.1")
  stimTimeVar <- tclVar("100")
  baselineVar <- tclVar("50")
  nVar <- tclVar("30")
  yAblineVar <- tclVar("0.1")
  funcVar <- tclVar("product1N")
  
  tkgrid(tklabel(mainOptionsFrame, text = "dt (ms):"), row = 0, column = 0, sticky = "w")
  tkgrid(tkentry(mainOptionsFrame, textvariable = dtVar, width = 10), row = 0, column = 1)
  tkgrid(tklabel(mainOptionsFrame, text = "Stimulation Time:"), row = 1, column = 0, sticky = "w")
  tkgrid(tkentry(mainOptionsFrame, textvariable = stimTimeVar, width = 10), row = 1, column = 1)
  tkgrid(tklabel(mainOptionsFrame, text = "Baseline:"), row = 2, column = 0, sticky = "w")
  tkgrid(tkentry(mainOptionsFrame, textvariable = baselineVar, width = 10), row = 2, column = 1)
  tkgrid(tklabel(mainOptionsFrame, text = "n:"), row = 3, column = 0, sticky = "w")
  tkgrid(tkentry(mainOptionsFrame, textvariable = nVar, width = 10), row = 3, column = 1)
  tkgrid(tklabel(mainOptionsFrame, text = "Fit cutoff:"), row = 4, column = 0, sticky = "w")
  tkgrid(tkentry(mainOptionsFrame, textvariable = yAblineVar, width = 10), row = 4, column = 1)
  tkgrid(tklabel(mainOptionsFrame, text = "Function:"), row = 5, column = 0, sticky = "w")
  funcChoices <- c("product1N", "product2N", "product3N")
  funcCombo <- ttkcombobox(mainOptionsFrame, textvariable = funcVar, values = funcChoices, width = 10)
  tkgrid(funcCombo, row = 5, column = 1)
  
  # --- Added Downsample Factor widget (default = 1) ---
  dsVar <- tclVar("1")
  tkgrid(tklabel(mainOptionsFrame, text = "Downsample Factor:"), row = 6, column = 0, sticky = "w")
  tkgrid(tkentry(mainOptionsFrame, textvariable = dsVar, width = 10), row = 6, column = 1)
  
  ### Fit Options Tab ###
  NVar <- tclVar("1")
  IEIVar <- tclVar("50")
  smoothVar <- tclVar("5")
  methodVar <- tclVar("BF.LM")
  weightMethodVar <- tclVar("none")
  sequentialFitVar <- tclVar("0")
  intervalMinVar <- tclVar("0.1")
  intervalMaxVar <- tclVar("0.9")
  lowerVar <- tclVar("")
  upperVar <- tclVar("")
  latencyLimitVar <- tclVar("")
  
  tkgrid(tklabel(fitOptionsFrame, text = "N:"), row = 0, column = 0, sticky = "w")
  tkgrid(tkentry(fitOptionsFrame, textvariable = NVar, width = 10), row = 0, column = 1)
  tkgrid(tklabel(fitOptionsFrame, text = "IEI:"), row = 1, column = 0, sticky = "w")
  tkgrid(tkentry(fitOptionsFrame, textvariable = IEIVar, width = 10), row = 1, column = 1)
  tkgrid(tklabel(fitOptionsFrame, text = "Smooth:"), row = 2, column = 0, sticky = "w")
  tkgrid(tkentry(fitOptionsFrame, textvariable = smoothVar, width = 10), row = 2, column = 1)
  tkgrid(tklabel(fitOptionsFrame, text = "Method:"), row = 3, column = 0, sticky = "w")
  methodChoices <- c("BF.LM", "LM", "GN", "port", "robust", "MLE")
  methodCombo <- ttkcombobox(fitOptionsFrame, textvariable = methodVar, values = methodChoices, width = 10)
  tkgrid(methodCombo, row = 3, column = 1)
  tkgrid(tklabel(fitOptionsFrame, text = "Weighting:"), row = 4, column = 0, sticky = "w")
  weightChoices <- c("none", "~y_sqrt", "~y")
  weightCombo <- ttkcombobox(fitOptionsFrame, textvariable = weightMethodVar, values = weightChoices, width = 10)
  tkgrid(weightCombo, row = 4, column = 1)
  tkgrid(tklabel(fitOptionsFrame, text = "Sequential Fit:"), row = 5, column = 0, sticky = "w")
  sequentialFitCheck <- tkcheckbutton(fitOptionsFrame, variable = sequentialFitVar)
  tkgrid(sequentialFitCheck, row = 5, column = 1)
  tkgrid(tklabel(fitOptionsFrame, text = "Min interval:"), row = 6, column = 0, sticky = "w")
  tkgrid(tkentry(fitOptionsFrame, textvariable = intervalMinVar, width = 10), row = 6, column = 1)
  tkgrid(tklabel(fitOptionsFrame, text = "Max interval:"), row = 7, column = 0, sticky = "w")
  tkgrid(tkentry(fitOptionsFrame, textvariable = intervalMaxVar, width = 10), row = 7, column = 1)
  tkgrid(tklabel(fitOptionsFrame, text = "Lower bounds:"), row = 8, column = 0, sticky = "w")
  tkgrid(tkentry(fitOptionsFrame, textvariable = lowerVar, width = 10), row = 8, column = 1)
  tkgrid(tklabel(fitOptionsFrame, text = "Upper bounds:"), row = 9, column = 0, sticky = "w")
  tkgrid(tkentry(fitOptionsFrame, textvariable = upperVar, width = 10), row = 9, column = 1)
  tkgrid(tklabel(fitOptionsFrame, text = "Latency limit:"), row = 10, column = 0, sticky = "w")
  tkgrid(tkentry(fitOptionsFrame, textvariable = latencyLimitVar, width = 10), row = 10, column = 1)
  
  ### MLE Settings Tab ###
  iterVar <- tclVar("1000")
  metropolisScaleVar <- tclVar("1.5")
  fitAttemptsVar <- tclVar("10")
  RWmVar <- tclVar("0")
  
  tkgrid(tklabel(mleSettingsFrame, text = "MLE Iterations:"), row = 0, column = 0, sticky = "w")
  tkgrid(tkentry(mleSettingsFrame, textvariable = iterVar, width = 10), row = 0, column = 1)
  tkgrid(tklabel(mleSettingsFrame, text = "Metropolis Scale:"), row = 1, column = 0, sticky = "w")
  tkgrid(tkentry(mleSettingsFrame, textvariable = metropolisScaleVar, width = 10), row = 1, column = 1)
  tkgrid(tklabel(mleSettingsFrame, text = "Fit Attempts:"), row = 2, column = 0, sticky = "w")
  tkgrid(tkentry(mleSettingsFrame, textvariable = fitAttemptsVar, width = 10), row = 2, column = 1)
  tkgrid(tklabel(mleSettingsFrame, text = "Random Walk Metropolis:"), row = 3, column = 0, sticky = "w")
  RWmCheck <- tkcheckbutton(mleSettingsFrame, variable = RWmVar)
  tkgrid(RWmCheck, row = 3, column = 1)
  
  ### Advanced Tab ###
  filterVar <- tclVar("0")
  fcVar <- tclVar("1000")
  relDecayFitLimitVar <- tclVar("0.1")
  halfWidthFitLimitVar <- tclVar("500")
  seedVar <- tclVar("42")
  dpVar <- tclVar("3")
  fastConstraintVar <- tclVar("0")
  fastConstraintMethodVar <- tclVar("rise")
  fastDecayLimitVar <- tclVar("")
  firstDelayConstraintVar <- tclVar("0")
  
  tkgrid(tklabel(advancedFrame, text = "Filter:"), row = 0, column = 0, sticky = "w")
  filterCheck <- tkcheckbutton(advancedFrame, variable = filterVar)
  tkgrid(filterCheck, row = 0, column = 1)
  tkgrid(tklabel(advancedFrame, text = "Filter cutoff (Hz):"), row = 1, column = 0, sticky = "w")
  tkgrid(tkentry(advancedFrame, textvariable = fcVar, width = 10), row = 1, column = 1)
  tkgrid(tklabel(advancedFrame, text = "Relative decay fit limit:"), row = 2, column = 0, sticky = "w")
  tkgrid(tkentry(advancedFrame, textvariable = relDecayFitLimitVar, width = 10), row = 2, column = 1)
  tkgrid(tklabel(advancedFrame, text = "Half-width fit limit:"), row = 3, column = 0, sticky = "w")
  tkgrid(tkentry(advancedFrame, textvariable = halfWidthFitLimitVar, width = 10), row = 3, column = 1)
  tkgrid(tklabel(advancedFrame, text = "Seed:"), row = 4, column = 0, sticky = "w")
  tkgrid(tkentry(advancedFrame, textvariable = seedVar, width = 10), row = 4, column = 1)
  tkgrid(tklabel(advancedFrame, text = "Decimal points:"), row = 5, column = 0, sticky = "w")
  tkgrid(tkentry(advancedFrame, textvariable = dpVar, width = 10), row = 5, column = 1)
  tkgrid(tklabel(advancedFrame, text = "Fast constraint:"), row = 6, column = 0, sticky = "w")
  fastConstraintCheck <- tkcheckbutton(advancedFrame, variable = fastConstraintVar)
  tkgrid(fastConstraintCheck, row = 6, column = 1)
  tkgrid(tklabel(advancedFrame, text = "Fast constraint method:"), row = 7, column = 0, sticky = "w")
  fastConstraintChoices <- c("rise", "peak")
  fastConstraintCombo <- ttkcombobox(advancedFrame, textvariable = fastConstraintMethodVar, values = fastConstraintChoices, width = 10)
  tkgrid(fastConstraintCombo, row = 7, column = 1)
  tkgrid(tklabel(advancedFrame, text = "Fast decay limit(s):"), row = 8, column = 0, sticky = "w")
  tkgrid(tkentry(advancedFrame, textvariable = fastDecayLimitVar, width = 10), row = 8, column = 1)
  tkgrid(tklabel(advancedFrame, text = "First delay constraint:"), row = 9, column = 0, sticky = "w")
  firstDelayCheck <- tkcheckbutton(advancedFrame, variable = firstDelayConstraintVar)
  tkgrid(firstDelayCheck, row = 9, column = 1)
  
  ### Graph Settings Tab ###
  lwdVar <- tclVar("1.2")
  xlabVar <- tclVar("time (ms)")
  ylabVar <- tclVar("PSC (pA)")
  plotHeightVar <- tclVar("5")
  plotWidthVar <- tclVar("5")
  
  tkgrid(tklabel(graphSettingsFrame, text = "Line width:"), row = 0, column = 0, sticky = "w")
  tkgrid(tkentry(graphSettingsFrame, textvariable = lwdVar, width = 10), row = 0, column = 1)
  tkgrid(tklabel(graphSettingsFrame, text = "x-axis label:"), row = 1, column = 0, sticky = "w")
  tkgrid(tkentry(graphSettingsFrame, textvariable = xlabVar, width = 10), row = 1, column = 1)
  tkgrid(tklabel(graphSettingsFrame, text = "y-axis label:"), row = 2, column = 0, sticky = "w")
  tkgrid(tkentry(graphSettingsFrame, textvariable = ylabVar, width = 10), row = 2, column = 1)
  tkgrid(tklabel(graphSettingsFrame, text = "Plot height:"), row = 3, column = 0, sticky = "w")
  tkgrid(tkentry(graphSettingsFrame, textvariable = plotHeightVar, width = 10), row = 3, column = 1)
  tkgrid(tklabel(graphSettingsFrame, text = "Plot width:"), row = 4, column = 0, sticky = "w")
  tkgrid(tkentry(graphSettingsFrame, textvariable = plotWidthVar, width = 10), row = 4, column = 1)
  
  ## Additional sidebar controls for analysis actions:
  userTmaxVar <- tclVar("")
  tkgrid(tklabel(sidebarFrame, text = "User Tmax:"), row = 3, column = 0, sticky = "w", pady = 5)
  tkgrid(tkentry(sidebarFrame, textvariable = userTmaxVar, width = 10), row = 3, column = 1, pady = 5)
  
  repeatConstraintVar <- tclVar("0")
  tkgrid(tklabel(sidebarFrame, text = "Repeat with fast constraint:"), row = 4, column = 0, sticky = "w")
  repeatConstraintCheck <- tkcheckbutton(sidebarFrame, variable = repeatConstraintVar)
  tkgrid(repeatConstraintCheck, row = 4, column = 1)
  
  ## Analysis action buttons
  
  runAnalysisButton <- tkbutton(sidebarFrame, text = "Run Initial Analysis",
    command = function() {
      filePath <- tclvalue(filePathVar)
      if (nchar(filePath) == 0) {
        tkmessageBox(message = "Please select a file first.")
        return()
      }
      if (nchar(tclvalue(columnVar)) == 0) {
        tkmessageBox(message = "Please select a column.")
        return()
      }
      ext <- tools::file_ext(filePath)
      if (tolower(ext) == "csv") {
        uploaded_data <<- read.csv(filePath)
      } else {
        uploaded_data <<- readxl::read_excel(filePath)
      }
      # Load the selected column
      response_data <<- uploaded_data[[tclvalue(columnVar)]]
      
      # --- Apply downsampling using the downsample factor ---
      ds <- as.numeric(tclvalue(dsVar))
      if (ds > 1) {
        response_data <<- response_data[seq(1, length(response_data), by = ds)]
      }
      
      tkrreplot(plotWidget, fun = drawPlot1)
    })
  tkgrid(runAnalysisButton, row = 5, column = 0, columnspan = 3, pady = 5)
    
  ## Run Main Analysis (inline version of analyse_PSC)
  runMainAnalysisButton <- tkbutton(sidebarFrame, text = "Run Main Analysis",
    command = function() {
      fast.constraint <- as.logical(as.numeric(tclvalue(repeatConstraintVar)))
      # --- Adjust dt based on downsample factor ---
      ds <- as.numeric(tclvalue(dsVar))
      dt <- as.numeric(tclvalue(dtVar)) * ds
      stimulation_time  <- as.numeric(tclvalue(stimTimeVar))
      baseline          <- as.numeric(tclvalue(baselineVar))
      smooth            <- as.numeric(tclvalue(smoothVar))
      n                 <- as.numeric(tclvalue(nVar))
      N                 <- as.numeric(tclvalue(NVar))
      IEI               <- as.numeric(tclvalue(IEIVar))
      func              <- get(tclvalue(funcVar))
      method            <- tclvalue(methodVar)
      weight_method     <- tclvalue(weightMethodVar)
      sequential.fit    <- as.logical(as.numeric(tclvalue(sequentialFitVar)))
      fit.limits        <- as.numeric(tclvalue(userTmaxVar))
      rel.decay.fit.limit <- as.numeric(tclvalue(relDecayFitLimitVar))
      lwd               <- as.numeric(tclvalue(lwdVar))
      xlab              <- tclvalue(xlabVar)
      ylab              <- tclvalue(ylabVar)
      plotWidth         <- as.numeric(tclvalue(plotWidthVar))
      plotHeight        <- as.numeric(tclvalue(plotHeightVar))
      
      fc                <- as.numeric(tclvalue(fcVar))
      interval          <- c(as.numeric(tclvalue(intervalMinVar)),
                             as.numeric(tclvalue(intervalMaxVar)))
      lower             <- if(nchar(tclvalue(lowerVar)) > 0)
                             as.numeric(unlist(strsplit(tclvalue(lowerVar), ",")))
                           else NULL
      upper             <- if(nchar(tclvalue(upperVar)) > 0)
                             as.numeric(unlist(strsplit(tclvalue(upperVar), ",")))
                           else NULL
      iter              <- as.numeric(tclvalue(iterVar))
      metropolis.scale  <- as.numeric(tclvalue(metropolisScaleVar))
      fit.attempts      <- as.numeric(tclvalue(fitAttemptsVar))
      RWm               <- as.logical(as.numeric(tclvalue(RWmVar)))
      fast.decay.limit  <- if(nchar(tclvalue(fastDecayLimitVar)) > 0)
                             as.numeric(unlist(strsplit(tclvalue(fastDecayLimitVar), ",")))
                           else NULL
      fast.constraint.method <- tclvalue(fastConstraintMethodVar)
      first.delay.constraint <- as.logical(as.numeric(tclvalue(firstDelayConstraintVar)))
      dp                <- as.numeric(tclvalue(dpVar))
      seed              <- as.numeric(tclvalue(seedVar))
      filter            <- as.logical(as.numeric(tclvalue(filterVar)))
      
      y <- response_data
      if (all(is.na(y[(which(!is.na(y))[length(which(!is.na(y)))] + 1):length(y)]))) {
        y <- y[!is.na(y)]
      }
      # x now uses the adjusted dt value
      x <- seq(0, (length(y) - 1) * dt, by = dt)
      
      if (!sequential.fit) {
        tmax <- fit.limits
        x_limit <- determine_tmax2(
          y = y,
          N = N,
          dt = dt,
          stimulation_time = stimulation_time,
          baseline = baseline,
          smooth = smooth,
          tmax = tmax,
          y_abline = rel.decay.fit.limit,
          ylab = ylab,
          width = plotWidth,
          height = plotHeight
        )
        adjusted_response <- y[x < x_limit]
        
        out <- nFIT(
          response = adjusted_response,
          n = n,
          N = N,
          IEI = IEI,
          dt = dt,
          func = func,
          method = method,
          weight_method = weight_method,
          MLEsettings = list(
            iter = iter,
            metropolis.scale = metropolis.scale,
            fit.attempts = fit.attempts,
            RWm = RWm
          ),
          stimulation_time = stimulation_time,
          baseline = baseline,
          filter = filter,
          fc = fc,
          interval = interval,
          fast.decay.limit = fast.decay.limit,
          fast.constraint = fast.constraint,
          fast.constraint.method = fast.constraint.method,
          first.delay.constraint = first.delay.constraint,
          lower = lower,
          upper = upper,
          latency.limit = if(nchar(tclvalue(latencyLimitVar)) > 0)
                            as.numeric(unlist(strsplit(tclvalue(latencyLimitVar), ",")))
                          else NULL,
          return.output = TRUE,
          show.plot = FALSE,
          half_width_fit_limit = as.numeric(tclvalue(halfWidthFitLimitVar)),
          dp = dp,
          height = plotHeight,
          width = plotWidth,
          seed = seed
        )
        
        out$traces <- traces_fun2(
          y = y,
          fits = out$fits,
          dt = dt,
          N = N,
          IEI = IEI,
          stimulation_time = stimulation_time,
          baseline = baseline,
          func = func,
          filter = filter,
          fc = fc
        )
        # Update the tkrplot using drawPlot2 so that the fit is shown in the GUI
        tkrreplot(plotWidget, fun = function() {
          drawPlot2(traces = out$traces, func = func, xlab = xlab, ylab = ylab,
                    lwd = lwd, filter = filter, width = plotWidth, height = plotHeight)
        })
      } else {
        out <- nFIT_sequential(
          response = y,
          n = n,
          dt = dt,
          func = func,
          method = method,
          weight_method = weight_method,
          stimulation_time = stimulation_time,
          baseline = baseline,
          fit.limits = fit.limits,
          fast.decay.limit = fast.decay.limit,
          fast.constraint = as.logical(as.numeric(tclvalue(fastConstraintVar))),
          fast.constraint.method = fast.constraint.method,
          first.delay.constraint = first.delay.constraint,
          latency.limit = if(nchar(tclvalue(latencyLimitVar)) > 0)
                            as.numeric(unlist(strsplit(tclvalue(latencyLimitVar), ",")))
                          else NULL,
          lower = lower,
          upper = upper,
          filter = filter,
          fc = fc,
          interval = interval,
          MLEsettings = list(
            iter = iter,
            metropolis.scale = metropolis.scale,
            fit.attempts = fit.attempts,
            RWm = RWm
          ),
          MLE.method = method,
          half_width_fit_limit = as.numeric(tclvalue(halfWidthFitLimitVar)),
          dp = dp,
          lwd = lwd,
          xlab = xlab,
          ylab = ylab,
          width = plotWidth,
          height = plotHeight,
          return.output = TRUE,
          show.output = TRUE,
          show.plot = TRUE,
          seed = seed
        )
      }
      
    # Print out$output to the console widget below the graph
    analysis_output <<- out
    tkdelete(consoleText, "1.0", "end")
    # Convert the data frame to a string using capture.output
    output_str <- paste(capture.output(print(out$output)), collapse = "\n")
    tkinsert(consoleText, "end", output_str)
    tcltk::tcl("update", "idletasks")

    }
  )
  tkgrid(runMainAnalysisButton, row = 7, column = 0, columnspan = 3, pady = 5)
  
  clearOutputButton <- tkbutton(sidebarFrame, text = "Clear Output",
    command = function() {
      analysis_output <<- NULL
      tkdelete(consoleText, "1.0", "end")
      tkrreplot(plotWidget, fun = drawPlot1)
    }
  )
  tkgrid(clearOutputButton, row = 9, column = 0, columnspan = 3, pady = 5)
  
  ### Main Panel: Plot and Console ###
  drawPlot1 <- function() {
    # Use downsample factor to adjust dt for plotting
    ds <- as.numeric(tclvalue(dsVar))
    dt <- as.numeric(tclvalue(dtVar)) * ds
    stimTime <- as.numeric(tclvalue(stimTimeVar))
    baseline <- as.numeric(tclvalue(baselineVar))
    smooth <- as.numeric(tclvalue(smoothVar))
    y_abline <- as.numeric(tclvalue(yAblineVar))
    
    y_val <- if (exists("response_data") && !is.null(response_data)) {
      response_data
    } else {
      rnorm(10000, 0.1)
    }
    
    determine_tmax2(
      y = y_val,
      N = 1,
      dt = dt,
      stimulation_time = stimTime,
      baseline = baseline,
      smooth = smooth,
      tmax = NULL,
      y_abline = y_abline,
      ylab = tclvalue(ylabVar),
      height = 5,
      width = 5
    )
  }
  
  # Initially, the plot is drawn using drawPlot1.
  plotWidget <- tkrplot(tt, fun = drawPlot1, hscale = 1.3, vscale = 1.3)
  tkgrid(plotWidget, row = 0, column = 1, sticky = "nsew")
  consoleText <- tktext(mainFrame, width = 80, height = 10)
  tkgrid(consoleText, row = 1, column = 1, sticky = "nsew")
  
  tkfocus(tt)
}

# Launch the PSC Analysis interface
PSC_analysis_tk()




# These steps 
# By manually creating the /tmp/.X11-unix directory with proper permissions, ensuring no conflicting X server processes are running, restarting XQuartz, and setting the DISPLAY variable, you’ve set up the correct environment for your Tcl/Tk applications in R.

#  1. Create the /tmp/.X11-unix directory manually with correct permissions
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






