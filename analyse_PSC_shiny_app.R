# library(shiny)

# ui <- fluidPage(
#   titlePanel("Interactive PSC Analysis with determine_tmax & Fast Constraint"),
#   sidebarLayout(
#     sidebarPanel(
#       fileInput("file", "Choose Response CSV File"),
#       numericInput("dt", "dt (ms):", value = 0.1),
#       numericInput("n", "n:", value = 30),
#       numericInput("N", "N:", value = 1),
#       numericInput("IEI", "IEI:", value = 50),
#       numericInput("stimulation_time", "Stimulation Time:", value = 150),
#       numericInput("baseline", "Baseline:", value = 50),
#       numericInput("smooth", "Smooth:", value = 5),
#       numericInput("y_abline", "y_abline:", value = 0.1),
#       actionButton("run_analysis", "Run Initial Analysis"),
#       conditionalPanel(
#         condition = "output.promptTmax",
#         numericInput("x_limit", "Enter Tmax for fitting:", value = NA),
#         actionButton("run_main_analysis", "Run Main Analysis")
#       ),
#       conditionalPanel(
#         condition = "output.promptFastConstraint",
#         checkboxInput("repeat_constraint", "Repeat with fast constraint?", FALSE),
#         actionButton("run_final_analysis", "Run Final Analysis")
#       ),
#       hr(),
#       downloadButton("download_output", "Download Output")
#     ),
#     mainPanel(
#       plotOutput("plot"),
#       verbatimTextOutput("console")
#     )
#   )
# )

# server <- function(input, output, session) {
#   analysis_output <- reactiveVal(NULL)
#   prompt_tmax <- reactiveVal(FALSE)
#   prompt_fast_constraint <- reactiveVal(FALSE)
#   response_data <- reactiveVal(NULL)

#   observeEvent(input$run_analysis, {
#     req(input$file)
#     data <- read.csv(input$file$datapath, header = FALSE)[,1]
#     response_data(data)

#     determine_tmax(
#       y = data, N = input$N, dt = input$dt,
#       stimulation_time = input$stimulation_time,
#       baseline = input$baseline, smooth = input$smooth,
#       y_abline = input$y_abline
#     )

#     prompt_tmax(TRUE)
#   })

#   observeEvent(input$run_main_analysis, {
#     req(!is.na(input$x_limit))
#     result <- analyse_PSC(
#       response = response_data(), dt = input$dt,
#       n = input$n, N = input$N, IEI = input$IEI,
#       stimulation_time = input$stimulation_time,
#       baseline = input$baseline, smooth = input$smooth,
#       fit.limits = input$x_limit,
#       fast.constraint = FALSE,
#       return.output = TRUE
#     )
#     analysis_output(result)
#     prompt_tmax(FALSE)
#     prompt_fast_constraint(TRUE)
#   })

#   observeEvent(input$run_final_analysis, {
#     result <- analyse_PSC(
#       response = response_data(), dt = input$dt,
#       n = input$n, N = input$N, IEI = input$IEI,
#       stimulation_time = input$stimulation_time,
#       baseline = input$baseline, smooth = input$smooth,
#       fit.limits = input$x_limit,
#       fast.constraint = input$repeat_constraint,
#       return.output = TRUE
#     )
#     analysis_output(result)
#     prompt_fast_constraint(FALSE)
#   })

#   output$plot <- renderPlot({
#     req(analysis_output())
#   })

#   output$console <- renderPrint({
#     req(analysis_output())
#     analysis_output()
#   })

#   output$promptTmax <- reactive({ prompt_tmax() })
#   output$promptFastConstraint <- reactive({ prompt_fast_constraint() })
#   outputOptions(output, "promptTmax", suspendWhenHidden = FALSE)
#   outputOptions(output, "promptFastConstraint", suspendWhenHidden = FALSE)

#   output$download_output <- downloadHandler(
#     filename = function() { paste0("PSC_analysis_", Sys.Date(), ".rds") },
#     content = function(file) {
#       saveRDS(analysis_output(), file)
#     }
#   )
# }

# shinyApp(ui, server)


# 1. use Use tab panels or collapsible sections to neatly group related parameters want to have EVERY input to analyse_PSC

# args(analyse_PSC)
# function (response, dt = 0.1, n = 30, N = 1, IEI = 50, stimulation_time = 150, 
#     baseline = 50, smooth = 5, func = product2N, method = c("BF.LM", 
#         "LM", "GN", "port", "robust", "MLE"), weight_method = c("none", 
#         "~y_sqrt", "~y"), sequential.fit = FALSE, fit.limits = NULL, 
#     MLEsettings = list(iter = 1000, metropolis.scale = 1.5, fit.attempts = 10, 
#         RWm = FALSE), filter = FALSE, fc = 1000, interval = c(0.1, 
#         0.9), lower = NULL, upper = NULL, fast.decay.limit = NULL, 
#     fast.constraint = FALSE, fast.constraint.method = c("rise", 
#         "peak"), first.delay.constraint = FALSE, latency.limit = NULL, 
#     rel.decay.fit.limit = 0.1, half_width_fit_limit = 500, dp = 3, 
#     lwd = 1.2, xlab = "time (ms)", ylab = "PSC (pA)", return.output = TRUE, 
#     height = 5, width = 5, seed = 42) 


# 2. Only dt, stimulation time, baseline, n y_abline on main 

# 3. allow upload of csv OR xlsx

# 4. csv or xlsx of form:


#     time     Control  GABAzine Mecamylamine    NBQX-AP5
# 1  299.0 -0.46735668 0.3154221    1.5275078 -0.39258766
# 2  299.1 -0.79287815 1.6826096    2.9313164  0.33983421
# 3  299.2 -0.38597584 1.7070236    2.3820000  0.66535568
# 4  299.3  0.10230541 0.8037033    1.8937187  1.56053734
# 5  299.4 -0.91494751 0.8769455    0.7950859  1.07225609
# 6  299.5 -0.71149826 1.2675705   -0.3645821  1.03156662
# 7  299.6  0.06161404 1.1455002   -0.4256172  1.15363693
# 8  299.7  0.18368435 1.5605392    1.1612968  0.05500412
# 9  299.8  0.26506519 0.8037033    1.5885429  1.23501778
# 10 299.9  0.71265507 0.8769455    0.6730156  0.74673653

# so must load then allow user to specify which column to analyse


library(shiny)

ui <- fluidPage(
  titlePanel("Interactive PSC Analysis with determine_tmax & Fast Constraint"),
  sidebarLayout(
    sidebarPanel(
      fileInput("file", "Choose Response CSV File"),
      numericInput("dt", "dt (ms):", value = 0.1),
      numericInput("n", "n:", value = 30),
      numericInput("N", "N:", value = 1),
      numericInput("IEI", "IEI:", value = 50),
      numericInput("stimulation_time", "Stimulation Time:", value = 150),
      numericInput("baseline", "Baseline:", value = 50),
      numericInput("smooth", "Smooth:", value = 5),
      numericInput("y_abline", "y_abline:", value = 0.1),
      actionButton("run_analysis", "Run Initial Analysis"),
      conditionalPanel(
        condition = "output.promptTmax",
        numericInput("x_limit", "Enter Tmax for fitting:", value = NA),
        actionButton("run_main_analysis", "Run Main Analysis")
      ),
      conditionalPanel(
        condition = "output.promptFastConstraint",
        checkboxInput("repeat_constraint", "Repeat with fast constraint?", FALSE),
        actionButton("run_final_analysis", "Run Final Analysis")
      ),
      hr(),
      downloadButton("download_output", "Download Output")
    ),
    mainPanel(
      plotOutput("plot"),
      verbatimTextOutput("console")
    )
  )
)

server <- function(input, output, session) {
  analysis_output <- reactiveVal(NULL)
  prompt_tmax <- reactiveVal(FALSE)
  prompt_fast_constraint <- reactiveVal(FALSE)
  response_data <- reactiveVal(NULL)

  observeEvent(input$run_analysis, {
    req(input$file)
    data <- read.csv(input$file$datapath, header = FALSE)[,1]
    response_data(data)

    determine_tmax(
      y = data, N = input$N, dt = input$dt,
      stimulation_time = input$stimulation_time,
      baseline = input$baseline, smooth = input$smooth,
      y_abline = input$y_abline
    )

    prompt_tmax(TRUE)
  })

  observeEvent(input$run_main_analysis, {
    req(!is.na(input$x_limit))
    result <- analyse_PSC(
      response = response_data(), dt = input$dt,
      n = input$n, N = input$N, IEI = input$IEI,
      stimulation_time = input$stimulation_time,
      baseline = input$baseline, smooth = input$smooth,
      fit.limits = input$x_limit,
      fast.constraint = FALSE,
      return.output = TRUE
    )
    analysis_output(result)
    prompt_tmax(FALSE)
    prompt_fast_constraint(TRUE)
  })

  observeEvent(input$run_final_analysis, {
    result <- analyse_PSC(
      response = response_data(), dt = input$dt,
      n = input$n, N = input$N, IEI = input$IEI,
      stimulation_time = input$stimulation_time,
      baseline = input$baseline, smooth = input$smooth,
      fit.limits = input$x_limit,
      fast.constraint = input$repeat_constraint,
      return.output = TRUE
    )
    analysis_output(result)
    prompt_fast_constraint(FALSE)
  })

  output$plot <- renderPlot({
    req(analysis_output())
  })

  output$console <- renderPrint({
    req(analysis_output())
    analysis_output()
  })

  output$promptTmax <- reactive({ prompt_tmax() })
  output$promptFastConstraint <- reactive({ prompt_fast_constraint() })
  outputOptions(output, "promptTmax", suspendWhenHidden = FALSE)
  outputOptions(output, "promptFastConstraint", suspendWhenHidden = FALSE)

  output$download_output <- downloadHandler(
    filename = function() { paste0("PSC_analysis_", Sys.Date(), ".rds") },
    content = function(file) {
      saveRDS(analysis_output(), file)
    }
  )
}

shinyApp(ui, server)


# 1. use Use tab panels or collapsible sections to neatly group related parameters want to have EVERY input to analyse_PSC

# args(analyse_PSC)
# function (response, dt = 0.1, n = 30, N = 1, IEI = 50, stimulation_time = 150, 
#     baseline = 50, smooth = 5, func = product2N, method = c("BF.LM", 
#         "LM", "GN", "port", "robust", "MLE"), weight_method = c("none", 
#         "~y_sqrt", "~y"), sequential.fit = FALSE, fit.limits = NULL, 
#     MLEsettings = list(iter = 1000, metropolis.scale = 1.5, fit.attempts = 10, 
#         RWm = FALSE), filter = FALSE, fc = 1000, interval = c(0.1, 
#         0.9), lower = NULL, upper = NULL, fast.decay.limit = NULL, 
#     fast.constraint = FALSE, fast.constraint.method = c("rise", 
#         "peak"), first.delay.constraint = FALSE, latency.limit = NULL, 
#     rel.decay.fit.limit = 0.1, half_width_fit_limit = 500, dp = 3, 
#     lwd = 1.2, xlab = "time (ms)", ylab = "PSC (pA)", return.output = TRUE, 
#     height = 5, width = 5, seed = 42) 


# 2. Only dt, stimulation time, baseline, n y_abline on main 

# 3. allow upload of csv OR xlsx

# 4. csv or xlsx of form:


#     time     Control  GABAzine Mecamylamine    NBQX-AP5
# 1  299.0 -0.46735668 0.3154221    1.5275078 -0.39258766
# 2  299.1 -0.79287815 1.6826096    2.9313164  0.33983421
# 3  299.2 -0.38597584 1.7070236    2.3820000  0.66535568
# 4  299.3  0.10230541 0.8037033    1.8937187  1.56053734
# 5  299.4 -0.91494751 0.8769455    0.7950859  1.07225609
# 6  299.5 -0.71149826 1.2675705   -0.3645821  1.03156662
# 7  299.6  0.06161404 1.1455002   -0.4256172  1.15363693
# 8  299.7  0.18368435 1.5605392    1.1612968  0.05500412
# 9  299.8  0.26506519 0.8037033    1.5885429  1.23501778
# 10 299.9  0.71265507 0.8769455    0.6730156  0.74673653

# so must load then allow user to specify which column to analyse


# library(shiny)
# library(readxl)

# ui <- fluidPage(
#   titlePanel("Comprehensive PSC Analysis"),
#   sidebarLayout(
#     sidebarPanel(
#       fileInput("file", "Upload CSV or XLSX", accept = c(".csv", ".xlsx")),
#       uiOutput("column_selector"),
      
#       tabsetPanel(
#         tabPanel("Main",
#           numericInput("dt", "dt (ms):", 0.1),
#           numericInput("stimulation_time", "Stimulation Time:", 150),
#           numericInput("baseline", "Baseline:", 50),
#           numericInput("n", "n:", 30),
#           numericInput("y_abline", "y_abline:", 0.1)
#         ),
#         tabPanel("Fit Options",
#           numericInput("N", "N:", 1),
#           numericInput("IEI", "IEI:", 50),
#           numericInput("smooth", "Smooth:", 5),
#           selectInput("method", "Method", c("BF.LM", "LM", "GN", "port", "robust", "MLE")),
#           selectInput("weight_method", "Weight Method", c("none", "~y_sqrt", "~y")),
#           checkboxInput("sequential_fit", "Sequential Fit", FALSE)
#         ),
#         tabPanel("Advanced",
#           checkboxInput("filter", "Filter", FALSE),
#           numericInput("fc", "Filter cutoff (Hz):", 1000),
#           numericInput("rel_decay_fit_limit", "Rel. Decay Fit Limit:", 0.1),
#           numericInput("half_width_fit_limit", "Half-width Fit Limit:", 500),
#           numericInput("seed", "Seed:", 42),
#           numericInput("dp", "Decimal Points:", 3),
#           checkboxInput("fast_constraint", "Fast Constraint", FALSE),
#           selectInput("fast_constraint_method", "Constraint Method", c("rise", "peak"))
#         )
#       ),
      
#       actionButton("run_analysis", "Run Initial Analysis"),
#       conditionalPanel(
#         condition = "output.promptTmax",
#         numericInput("x_limit", "Enter Tmax for fitting:", NA),
#         actionButton("run_main_analysis", "Run Main Analysis")
#       ),
#       conditionalPanel(
#         condition = "output.promptFastConstraint",
#         checkboxInput("repeat_constraint", "Repeat with fast constraint?", FALSE),
#         actionButton("run_final_analysis", "Run Final Analysis")
#       ),
#       hr(),
#       downloadButton("download_output", "Download Output")
#     ),
#     mainPanel(
#       plotOutput("plot"),
#       verbatimTextOutput("console")
#     )
#   )
# )

# server <- function(input, output, session) {
  
#   analysis_output <- reactiveVal(NULL)
#   response_data <- reactiveVal(NULL)
#   prompt_tmax <- reactiveVal(FALSE)
#   suggested_tmax <- reactiveVal(NA)
  
#   data_loaded <- reactive({
#     req(input$file)
#     ext <- tools::file_ext(input$file$name)
#     if(ext == "csv") read.csv(input$file$datapath) else readxl::read_excel(input$file$datapath)
#   })
  
#   output$column_selector <- renderUI({
#     req(data_loaded())
#     selectInput("data_col", "Select column to analyse", choices=names(data_loaded()))
#   })
  
#   observeEvent(input$run_analysis, {
#     req(data_loaded(), input$data_col)
#     response_data(data_loaded()[[input$data_col]])
    
#     output$plot <- renderPlot({
#       suggested <- determine_tmax(
#         y=response_data(), N=input$N, dt=input$dt,
#         stimulation_time=input$stimulation_time,
#         baseline=input$baseline, smooth=input$smooth,
#         y_abline=input$y_abline
#       )
#       suggested_tmax(suggested)
#     })
    
#     prompt_tmax(TRUE)
#   })
  
#   output$tmax_ui <- renderUI({
#     req(prompt_tmax(), suggested_tmax())
#     numericInput("user_tmax", "Enter Tmax for fitting:", value=round(suggested_tmax(),2))
#   })

#   observeEvent(input$run_main_analysis, {
#     req(input$user_tmax)
#     result <- analyse_PSC(
#       response=response_data(), dt=input$dt, n=input$n, N=input$N, IEI=input$IEI,
#       stimulation_time=input$stimulation_time, baseline=input$baseline, smooth=input$smooth,
#       fit.limits=input$user_tmax, fast.constraint=FALSE,
#       method=input$method, weight_method=input$weight_method, filter=input$filter, fc=input$fc,
#       rel.decay.fit.limit=input$rel_decay_fit_limit,
#       half_width_fit_limit=input$half_width_fit_limit, seed=input$seed,
#       sequential.fit=input$sequential_fit, dp=input$dp, return.output=TRUE
#     )
#     analysis_output(result)
#     prompt_tmax(FALSE)
#   })
  
#   output$promptTmax <- reactive({ prompt_tmax() })
#   outputOptions(output, "promptTmax", suspendWhenHidden=FALSE)
  
#   output$console <- renderPrint({ req(analysis_output()) })
  
#   output$download_output <- downloadHandler(
#     filename=function() paste0("PSC_analysis_", Sys.Date(), ".rds"),
#     content=function(file) saveRDS(analysis_output(), file)
#   )
# }

# shinyApp(ui, server)




# # Remove all objects from the environment
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

# # insert your username here to define the correct path
# username <- 'YourUsername'
# username <- 'euo9382'

# path_repository <- '/Documents/Repositories/Rfits'
# # construct file path
# file_path1 <- paste0('/Users/', username, path_repository)

# source(paste0(file_path1, '/nNLS functions.R'))


# library(shiny)
# library(readxl)

# ui <- fluidPage(
#   titlePanel("Comprehensive PSC Analysis"),
#   sidebarLayout(
#     sidebarPanel(
#       fileInput("file", "Upload CSV or XLSX", accept = c(".csv", ".xlsx")),
#       uiOutput("column_selector"),
      
#       tabsetPanel(
#         tabPanel("Main",
#           numericInput("dt", "dt (ms):", 0.1),
#           numericInput("stimulation_time", "Stimulation Time:", 150),
#           numericInput("baseline", "Baseline:", 50),
#           numericInput("n", "n:", 30),
#           numericInput("y_abline", "y_abline:", 0.1)
#         ),
#         tabPanel("Fit Options",
#           numericInput("N", "N:", 1),
#           numericInput("IEI", "IEI:", 50),
#           numericInput("smooth", "Smooth:", 5),
#           selectInput("method", "Method", c("BF.LM", "LM", "GN", "port", "robust", "MLE")),
#           selectInput("weight_method", "Weight Method", c("none", "~y_sqrt", "~y")),
#           checkboxInput("sequential_fit", "Sequential Fit", FALSE)
#         ),
#         tabPanel("Advanced",
#           checkboxInput("filter", "Filter", FALSE),
#           numericInput("fc", "Filter cutoff (Hz):", 1000),
#           numericInput("rel_decay_fit_limit", "Rel. Decay Fit Limit:", 0.1),
#           numericInput("half_width_fit_limit", "Half-width Fit Limit:", 500),
#           numericInput("seed", "Seed:", 42),
#           numericInput("dp", "Decimal Points:", 3),
#           checkboxInput("fast_constraint", "Fast Constraint", FALSE),
#           selectInput("fast_constraint_method", "Constraint Method", c("rise", "peak"))
#         )
#       ),
      
#       actionButton("run_analysis", "Run Initial Analysis"),
#       conditionalPanel(
#         condition = "output.promptTmax",
#         numericInput("user_tmax", "Enter Tmax for fitting:", value=NA),
#         numericInput("stimulation_time_adj", "Adjust Stimulation Time:", value=150),
#         numericInput("baseline_adj", "Adjust Baseline:", value=50),
#         actionButton("update_plot", "Update Plot"),
#         actionButton("run_main_analysis", "Run Main Analysis")
#       ),
#       conditionalPanel(
#         condition = "output.promptFastConstraint",
#         checkboxInput("repeat_constraint", "Repeat with fast constraint?", FALSE),
#         actionButton("run_final_analysis", "Run Final Analysis")
#       ),
#       hr(),
#       downloadButton("download_output", "Download Output")
#     ),
#     mainPanel(
#       plotOutput("plot"),
#       verbatimTextOutput("console")
#     )
#   )
# )

# server <- function(input, output, session) {
  
#   analysis_output <- reactiveVal(NULL)
#   response_data <- reactiveVal(NULL)
#   prompt_tmax <- reactiveVal(FALSE)
#   suggested_tmax <- reactiveVal(NA)

#   data_loaded <- reactive({
#     req(input$file)
#     ext <- tools::file_ext(input$file$name)
#     if(ext == "csv") read.csv(input$file$datapath) else readxl::read_excel(input$file$datapath)
#   })
  
#   output$column_selector <- renderUI({
#     req(data_loaded())
#     selectInput("data_col", "Select column to analyse", choices=names(data_loaded()))
#   })
  
#   current_params <- reactiveValues(
#     stimulation_time = NULL,
#     baseline = NULL
#   )
  
#   observeEvent(input$run_analysis, {
#     req(data_loaded(), input$data_col)
#     response_data(data_loaded()[[input$data_col]])

#     current_params$stimulation_time <- input$stimulation_time
#     current_params$baseline <- input$baseline

#     prompt_tmax(TRUE)
#   })
  
#   observeEvent(input$update_plot, {
#     current_params$stimulation_time <- input$stimulation_time_adj
#     current_params$baseline <- input$baseline_adj
#   })
  
#   output$plot <- renderPlot({
#     req(prompt_tmax())
#     suggested <- determine_tmax(
#       y=response_data(), N=input$N, dt=input$dt,
#       stimulation_time=current_params$stimulation_time,
#       baseline=current_params$baseline, smooth=input$smooth,
#       y_abline=input$y_abline, prompt=FALSE
#     )
#     suggested_tmax(suggested)
#   })

#   output$tmax_ui <- renderUI({
#     req(suggested_tmax())
#     numericInput("user_tmax", "Enter Tmax for fitting:", value=round(suggested_tmax(),2))
#   })

#   observeEvent(input$run_main_analysis, {
#     req(input$user_tmax)
#     result <- analyse_PSC(
#       response=response_data(), dt=input$dt, n=input$n, N=input$N, IEI=input$IEI,
#       stimulation_time=current_params$stimulation_time,
#       baseline=current_params$baseline, smooth=input$smooth,
#       fit.limits=input$user_tmax, fast.constraint=FALSE,
#       method=input$method, weight_method=input$weight_method, filter=input$filter, fc=input$fc,
#       rel.decay.fit.limit=input$rel_decay_fit_limit,
#       half_width_fit_limit=input$half_width_fit_limit, seed=input$seed,
#       sequential.fit=input$sequential_fit, dp=input$dp, return.output=TRUE
#     )
#     analysis_output(result)
#     prompt_tmax(FALSE)
#   })

#   output$promptTmax <- reactive({ prompt_tmax() })
#   outputOptions(output, "promptTmax", suspendWhenHidden=FALSE)
  
#   output$console <- renderPrint({ req(analysis_output()) })
  
#   output$download_output <- downloadHandler(
#     filename=function() paste0("PSC_analysis_", Sys.Date(), ".rds"),
#     content=function(file) saveRDS(analysis_output(), file)
#   )
# }

# shinyApp(ui, server)


# # Remove all objects from the environment
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

# # insert your username here to define the correct path
# username <- 'YourUsername'
# username <- 'euo9382'

# path_repository <- '/Documents/Repositories/Rfits'
# # construct file path
# file_path1 <- paste0('/Users/', username, path_repository)

# source(paste0(file_path1, '/nNLS functions.R'))


# library(shiny)
# library(readxl)

# ui <- fluidPage(
#   titlePanel("Comprehensive PSC Analysis"),
#   sidebarLayout(
#     sidebarPanel(
#       fileInput("file", "Upload CSV or XLSX", accept = c(".csv", ".xlsx")),
#       uiOutput("column_selector"),
      
#       tabsetPanel(
#         tabPanel("Main options",
#           numericInput("dt", "dt (ms):", 0.1),
#           numericInput("stimulation_time", "Stimulation Time:", 100),
#           numericInput("baseline", "Baseline:", 100),
#           numericInput("n", "n:", 30),
#           numericInput("y_abline", "y_abline:", 0.1)
#         ),
#         tabPanel("Fit Options",
#           numericInput("N", "N:", 1),
#           numericInput("IEI", "IEI:", 50),
#           numericInput("smooth", "Smooth:", 5),
#           selectInput("method", "Method", c("BF.LM", "LM", "GN", "port", "robust", "MLE")),
#           selectInput("weight_method", "Weight Method", c("none", "~y_sqrt", "~y")),
#           checkboxInput("sequential_fit", "Sequential Fit", FALSE)
#         ),
#         tabPanel("Advanced",
#           checkboxInput("filter", "Filter", FALSE),
#           numericInput("fc", "Filter cutoff (Hz):", 1000),
#           numericInput("rel_decay_fit_limit", "Rel. Decay Fit Limit:", 0.1),
#           numericInput("half_width_fit_limit", "Half-width Fit Limit:", 500),
#           numericInput("seed", "Seed:", 42),
#           numericInput("dp", "Decimal Points:", 3),
#           checkboxInput("fast_constraint", "Fast Constraint", FALSE),
#           selectInput("fast_constraint_method", "Constraint Method", c("rise", "peak"))
#         )
#       ),
      
#       actionButton("run_analysis", "Run Initial Analysis"),
#       conditionalPanel(
#         condition = "output.promptTmax",
#         numericInput("user_tmax", "Enter Tmax for fitting:", value=NA),
#         numericInput("stimulation_time_adj", "Adjust Stimulation Time:", value=150),
#         numericInput("baseline_adj", "Adjust Baseline:", value=50),
#         actionButton("update_plot", "Update Plot"),
#         actionButton("run_main_analysis", "Run Main Analysis")
#       ),
#       conditionalPanel(
#         condition = "output.promptFastConstraint",
#         checkboxInput("repeat_constraint", "Repeat with fast constraint?", FALSE),
#         actionButton("run_final_analysis", "Run Final Analysis")
#       ),
#       hr(),
#       downloadButton("download_output", "Download Output")
#     ),
#     mainPanel(
#       plotOutput("plot"),
#       verbatimTextOutput("console")
#     )
#   )
# )

# server <- function(input, output, session) {
  
#   analysis_output <- reactiveVal(NULL)
#   response_data <- reactiveVal(NULL)
#   prompt_tmax <- reactiveVal(FALSE)
#   suggested_tmax <- reactiveVal(NA)
#   prompt_fast_constraint <- reactiveVal(FALSE)

#   data_loaded <- reactive({
#     req(input$file)
#     ext <- tools::file_ext(input$file$name)
#     if(ext == "csv") read.csv(input$file$datapath) else readxl::read_excel(input$file$datapath)
#   })
  
#   output$column_selector <- renderUI({
#     req(data_loaded())
#     selectInput("data_col", "Select column to analyse", choices=names(data_loaded()))
#   })
  
#   current_params <- reactiveValues(
#     stimulation_time = NULL,
#     baseline = NULL
#   )
  
#   observeEvent(input$run_analysis, {
#     req(data_loaded(), input$data_col)
#     response_data(data_loaded()[[input$data_col]])

#     current_params$stimulation_time <- input$stimulation_time
#     current_params$baseline <- input$baseline

#     prompt_tmax(TRUE)
#   })
  
#   observeEvent(input$update_plot, {
#     current_params$stimulation_time <- input$stimulation_time_adj
#     current_params$baseline <- input$baseline_adj
#   })
  
#   output$plot <- renderPlot({
#     req(prompt_tmax())
#     suggested <- determine_tmax(
#       y=response_data(), N=input$N, dt=input$dt,
#       stimulation_time=current_params$stimulation_time,
#       baseline=current_params$baseline, smooth=input$smooth,
#       y_abline=input$y_abline, prompt=FALSE
#     )
#     suggested_tmax(suggested)
#   })

#   output$tmax_ui <- renderUI({
#     req(suggested_tmax())
#     numericInput("user_tmax", "Enter Tmax for fitting:", value=round(suggested_tmax(),2))
#   })

#   observeEvent(input$run_main_analysis, {
#     req(input$user_tmax)
#     result <- analyse_PSC(
#       response=response_data(), dt=input$dt, n=input$n, N=input$N, IEI=input$IEI,
#       stimulation_time=current_params$stimulation_time,
#       baseline=current_params$baseline, smooth=input$smooth,
#       fit.limits=input$user_tmax, fast.constraint=FALSE,
#       method=input$method, weight_method=input$weight_method, filter=input$filter, fc=input$fc,
#       rel.decay.fit.limit=input$rel_decay_fit_limit,
#       half_width_fit_limit=input$half_width_fit_limit, seed=input$seed,
#       sequential.fit=input$sequential_fit, dp=input$dp, return.output=TRUE
#     )
#     analysis_output(result)
#     prompt_tmax(FALSE)
#     prompt_fast_constraint(TRUE)
#   })

#   observeEvent(input$run_final_analysis, {
#     result <- analyse_PSC(
#       response=response_data(), dt=input$dt, n=input$n, N=input$N, IEI=input$IEI,
#       stimulation_time=current_params$stimulation_time,
#       baseline=current_params$baseline, smooth=input$smooth,
#       fit.limits=input$user_tmax, fast.constraint=input$repeat_constraint,
#       method=input$method, weight_method=input$weight_method, filter=input$filter, fc=input$fc,
#       rel.decay.fit.limit=input$rel_decay_fit_limit,
#       half_width_fit_limit=input$half_width_fit_limit, seed=input$seed,
#       sequential.fit=input$sequential_fit, dp=input$dp, return.output=TRUE
#     )
#     analysis_output(result)
#     prompt_fast_constraint(FALSE)
#   })

#   output$promptTmax <- reactive({ prompt_tmax() })
#   output$promptFastConstraint <- reactive({ prompt_fast_constraint() })
#   outputOptions(output, "promptTmax", suspendWhenHidden=FALSE)
#   outputOptions(output, "promptFastConstraint", suspendWhenHidden=FALSE)

#   output$console <- renderPrint({ req(analysis_output()) })
  
#   output$download_output <- downloadHandler(
#     filename=function() paste0("PSC_analysis_", Sys.Date(), ".rds"),
#     content=function(file) saveRDS(analysis_output(), file)
#   )
# }

# shinyApp(ui, server)



# Remove all objects from the environment
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

ui <- fluidPage(
  titlePanel("Comprehensive PSC Analysis"),
  sidebarLayout(
    sidebarPanel(
      fileInput("file", "Upload CSV or XLSX", accept = c(".csv", ".xlsx")),
      uiOutput("column_selector"),
      
      tabsetPanel(
        tabPanel("Main Options",
          numericInput("dt", "dt (ms):", 0.1),
          numericInput("stimulation_time", "Stimulation Time:", 150),
          numericInput("baseline", "Baseline:", 50),
          numericInput("n", "n:", 30),
          numericInput("y_abline", "y_abline:", 0.1),
          selectInput("func", "Function:", choices = c("product2N"))
        ),
        tabPanel("Fit Options",
          numericInput("N", "N:", 1),
          numericInput("IEI", "IEI:", 50),
          numericInput("smooth", "Smooth:", 5),
          selectInput("method", "Method", c("BF.LM", "LM", "GN", "port", "robust", "MLE")),
          selectInput("weight_method", "Weight Method", c("none", "~y_sqrt", "~y")),
          checkboxInput("sequential_fit", "Sequential Fit", FALSE),
          numericInput("interval_min", "Interval Min:", 0.1),
          numericInput("interval_max", "Interval Max:", 0.9),
          textInput("lower", "Lower Bounds (comma-separated):", ""),
          textInput("upper", "Upper Bounds (comma-separated):", ""),
          textInput("latency.limit", "latency limit:", "")
        ),
        tabPanel("MLE Settings",
          numericInput("iter", "MLE Iterations:", 1000),
          numericInput("metropolis_scale", "Metropolis Scale:", 1.5),
          numericInput("fit_attempts", "Fit Attempts:", 10),
          checkboxInput("RWm", "Random Walk Metropolis:", FALSE)
        ),
        tabPanel("Advanced",
          checkboxInput("filter", "Filter", FALSE),
          numericInput("fc", "Filter cutoff (Hz):", 1000),
          numericInput("rel_decay_fit_limit", "Rel. Decay Fit Limit:", 0.1),
          numericInput("half_width_fit_limit", "Half-width Fit Limit:", 500),
          numericInput("seed", "Seed:", 42),
          numericInput("dp", "Decimal Points:", 3),
          checkboxInput("fast_constraint", "Fast Constraint", FALSE),
          selectInput("fast_constraint_method", "Constraint Method", c("rise", "peak")),
          textInput("fast_decay_limit", "Fast Decay Limit (comma-separated):", ""),
          checkboxInput("first_delay_constraint", "First Delay Constraint", FALSE)
        ),
        tabPanel("Graph Settings",
          numericInput("lwd", "Line Width:", 1.2),
          textInput("xlab", "X-axis Label:", "time (ms)"),
          textInput("ylab", "Y-axis Label:", "PSC (pA)"),
          numericInput("plot_height", "Plot Height:", 5),
          numericInput("plot_width", "Plot Width:", 5)
        )
      ),
      
      actionButton("run_analysis", "Run Initial Analysis"),
      conditionalPanel(
        condition = "output.promptTmax",
        numericInput("user_tmax", "Enter Tmax for fitting:", value=NA),
        numericInput("stimulation_time_adj", "Adjust Stimulation Time:", value=150),
        numericInput("baseline_adj", "Adjust Baseline:", value=50),
        actionButton("update_plot", "Update Plot"),
        actionButton("run_main_analysis", "Run Main Analysis")
      ),
      conditionalPanel(
        condition = "output.promptFastConstraint",
        checkboxInput("repeat_constraint", "Repeat with fast constraint?", FALSE),
        actionButton("run_final_analysis", "Run Final Analysis")
      ),
      hr(),
      downloadButton("download_output", "Download Output")
    ),
    mainPanel(
      plotOutput("plot"),
      verbatimTextOutput("console")
    )
  )
)

server <- function(input, output, session) {
  
  analysis_output <- reactiveVal(NULL)
  response_data <- reactiveVal(NULL)
  prompt_tmax <- reactiveVal(FALSE)
  suggested_tmax <- reactiveVal(NA)
  prompt_fast_constraint <- reactiveVal(FALSE)

  data_loaded <- reactive({
    req(input$file)
    ext <- tools::file_ext(input$file$name)
    if(ext == "csv") read.csv(input$file$datapath) else readxl::read_excel(input$file$datapath)
  })
  
  output$column_selector <- renderUI({
    req(data_loaded())
    selectInput("data_col", "Select column to analyse", choices=names(data_loaded()))
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
  })
  
  observeEvent(input$update_plot, {
    current_params$stimulation_time <- input$stimulation_time_adj
    current_params$baseline <- input$baseline_adj
  })
  
  output$plot <- renderPlot({
    req(prompt_tmax())
    suggested <- determine_tmax(
      y=response_data(), N=input$N, dt=input$dt,
      stimulation_time=current_params$stimulation_time,
      baseline=current_params$baseline, smooth=input$smooth,
      y_abline=input$y_abline, prompt=FALSE
    )
    suggested_tmax(suggested)
  })

  output$tmax_ui <- renderUI({
    req(suggested_tmax())
    numericInput("user_tmax", "Enter Tmax for fitting:", value=round(suggested_tmax(),2))
  })

  observeEvent(input$run_main_analysis, {
    req(input$user_tmax)

    latency.limit <- if (nzchar(input$latency.limit)) as.numeric(unlist(strsplit(input$latency.limit, ","))) else NULL
    lower <- if (nzchar(input$lower)) as.numeric(unlist(strsplit(input$lower, ","))) else NULL
    upper <- if (nzchar(input$upper)) as.numeric(unlist(strsplit(input$upper, ","))) else NULL
    fast.decay.limit <- if (nzchar(input$fast_decay_limit)) as.numeric(unlist(strsplit(input$fast_decay_limit, ","))) else NULL

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
    latency.limit <- if (nzchar(input$latency.limit)) as.numeric(unlist(strsplit(input$latency.limit, ","))) else NULL
    lower <- if (nzchar(input$lower)) as.numeric(unlist(strsplit(input$lower, ","))) else NULL
    upper <- if (nzchar(input$upper)) as.numeric(unlist(strsplit(input$upper, ","))) else NULL
    fast.decay.limit <- if (nzchar(input$fast_decay_limit)) as.numeric(unlist(strsplit(input$fast_decay_limit, ","))) else NULL

    result <- analyse_PSC(
      response=response_data(), dt=input$dt, n=input$n, N=input$N, IEI=input$IEI,
      stimulation_time=current_params$stimulation_time,
      baseline=current_params$baseline, smooth=input$smooth,
      fit.limits=input$user_tmax, fast.constraint=input$repeat_constraint,
      method=input$method, weight_method=input$weight_method, filter=input$filter, fc=input$fc,
      rel.decay.fit.limit=input$rel_decay_fit_limit,
      half_width_fit_limit=input$half_width_fit_limit, seed=input$seed,
      sequential.fit=input$sequential_fit, dp=input$dp,
      func = get(input$func), interval = c(input$interval_min, input$interval_max),
      lower = lower, upper = upper,
      MLEsettings=list(iter=input$iter, metropolis.scale=input$metropolis_scale, fit.attempts=input$fit_attempts, RWm=input$RWm),
      fast.decay.limit=fast.decay.limit, fast.constraint.method=input$fast_constraint_method,
      first.delay.constraint=input$first_delay_constraint, latency.limit=latency.limit,
      lwd=input$lwd, xlab=input$xlab, ylab=input$ylab, height=input$plot_height, width=input$plot_width,
      return.output=TRUE
    )
    analysis_output(result)
    prompt_fast_constraint(FALSE)
  })

  output$promptTmax <- reactive({ prompt_tmax() })
  output$promptFastConstraint <- reactive({ prompt_fast_constraint() })
  outputOptions(output, "promptTmax", suspendWhenHidden=FALSE)
  outputOptions(output, "promptFastConstraint", suspendWhenHidden=FALSE)

  output$console <- renderPrint({ req(analysis_output()) })
  
  output$download_output <- downloadHandler(
    filename = function() {
      paste0(tools::file_path_sans_ext(basename(input$file$name)), "_", input$data_col, "_PSC_analysis.rds")
    },
    content = function(file) {
      saveRDS(analysis_output(), file)
    }
  )

}

shinyApp(ui, server)


