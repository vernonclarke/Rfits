# R Code for Non-Linear Least Squares Fitting of Signal Waveforms

This repository contains R scripts and functions for performing non-linear least squares fitting of functions that describe signal waveforms. The code was developed using the R GUI and tested on R version 4.4.1 ("Race for Your Life"). It guides you from initial setup through advanced analysis.

---

## Table of Contents

- [Initial Set Up](#initial-set-up)
  - [Setting Up](#setting-up)
- [Initial Guide](#initial-guide)
- [Step-by-Step Guide](#step-by-step-guide)
- [Definitions and Formulae](#definitions-and-formulae)

---

## Initial Set Up

All analyses are executed in R. The code was developed using the R GUI and is known to work with R version 4.4.1. The following software must be installed:

- **R Statistical Software:**  
  [https://www.R-project.org/](https://www.R-project.org/)
- **XQuartz** (for graphical output on macOS):  
  [https://www.xquartz.org/](https://www.xquartz.org/)
- **Sublime Text** (or any text editor; you may also use the default R text editor):  
  [https://www.sublimetext.com/](https://www.sublimetext.com/)

> **Note:** At a minimum, both **R** and **XQuartz** must be installed for this code to work.

### Setting Up

1. **Download the Code:**  
   On GitHub, click the green **<> Code** dropdown and select `Download Zip`.

2. **Create a Working Directory:**  
   Unpack the ZIP file and create a directory (e.g., `/Users/YourUserName/Documents/Rfits`, replacing `YourUserName` with your actual username).

3. **Install Required Packages:**  
   Run the following R code to install (if needed) and load the required packages:
   
   ```r
   # Remove all objects from the environment
   rm(list = ls(all = TRUE))
   
   # Function to load required packages
   load_required_packages <- function(packages) {
      new.packages <- packages[!(packages %in% installed.packages()[, 'Package'])]
      if (length(new.packages)) install.packages(new.packages)
      invisible(lapply(packages, library, character.only = TRUE))
   }
   
   # List of required packages
   required.packages <- c('robustbase', 'minpack.lm', 'Rcpp', 'signal', 'writexl')
   load_required_packages(required.packages)
   ```

4. **Load Custom Functions:**  
   The custom functions are contained in the `nNLS functions.R` file. Adjust the file path as needed:
   
   ```r
   UserName <- 'YourUserName'  # Replace with your actual username
   root_dir <- paste0('/Users/', UserName, '/Documents/Repositories/Rfits')
   path <- file.path(root_dir, 'nNLS functions.R')
   source(path)
   ```

---

## Initial Guide

Follow these steps to run a simple fit example:

1. **Open the R GUI.**

2. **Run the Package Loader Code and Load Custom Functions**  
   (See the "Setting Up" section above.)

3. **Fitting Example:**  
   The example below generates a noisy signal consisting of a train of 3 responses (with an inter-event interval of 50 ms) and performs the fit:
   
   ```r
   dx <- 0.1
   stimulation_time <- 150
   baseline <- 150
   xmax <- 1000
   x <- seq(dx, xmax, dx)
   N <- 3
   IEI <- 50

   # Define parameters for the fit: A1, A2, A3, τrise, τdecay, and delay.
   params1 <- c(-150, -250, -300, 1, 30, 4)
   params1_ <- params1
   params1_[N+3] <- params1_[N+3] + stimulation_time

   std.dev <- 10
   ysignal <- product1N(params = params1_, x = x, N = N, IEI = IEI)
   y <- ysignal + rnorm(length(x), sd = std.dev)

   # (Optional) Quick plot to view the signal
   # plot(x, y, type = 'l')

   # Analyze the signal
   analyse_PSC(
      response = y,
      dt = 0.1,
      n = 30,
      N = 3,
      IEI = 50,
      stimulation_time = 150,
      baseline = 150,
      func = product1N,
      return.output = FALSE
   )
   ```

   When you run the analysis, a graph will appear and you will be prompted in the R console:
   
   ```
   Enter the upper limit for time to use in nFIT (e.g., 400 ms):
   ```
   
   Enter a suitable value (e.g., `330` or `600`) to define the cutoff point for curve fitting. After a short delay, the best fit will be superimposed on the original signal, and a table of parameters (such as amplitudes, time constants, and delay) will be output.

---

## Step-by-Step Guide

This section explains how to simulate a dataset, save it, load it back into R, and analyze multiple responses.

1. **Setup the Environment:**
   
   ```r
   # Remove all objects from the environment
   rm(list = ls(all = TRUE))
   
   # Load required packages
   load_required_packages <- function(packages) {
     new.packages <- packages[!(packages %in% installed.packages()[, 'Package'])]
     if (length(new.packages)) install.packages(new.packages)
     invisible(lapply(packages, library, character.only = TRUE))
   }
   
   required.packages <- c('robustbase', 'minpack.lm', 'Rcpp', 'signal', 'writexl')
   load_required_packages(required.packages)
   
   UserName <- 'YourUserName'  # Replace with your username
   root_dir <- paste0('/Users/', UserName, '/Documents/Repositories/Rfits')
   path <- file.path(root_dir, 'nNLS functions.R')
   source(path)
   ```

2. **Create Data to Analyze:**
   
   The following code simulates 10 responses with added Gaussian noise and saves the dataset as both CSV and XLSX files.
   
   ```r
   # Parameters for modelled response
   dx <- 0.1
   stim_time <- 150
   baseline <- 150
   xmax <- 1000
   x <- seq(dx, xmax, dx)
   
   a1 <- 50; a2 <- 100
   tau1.1 <- 3; tau1.2 <- 30
   tau2.1 <- 10; tau2.2 <- 200
   d1 <- 2; d2 <- 5
   std.dev <- 5
   params <- c(a = a1, b = tau1.1, c = tau1.2, d = d1 + stim_time,
               e = a2, f = tau2.1, g = tau2.2, h = d2 + stim_time)
   
   # Create data
   set.seed(7)
   data <- sapply(1:10, function(ii) {
      ysignal <- product2N(params, x)
      y <- ysignal + rnorm(length(x), sd = std.dev)
      -y  # Invert signal if needed
   })
   
   # Create 'examples' directory if it doesn't exist
   dir.create(file.path(root_dir, 'examples'), showWarnings = FALSE)
   
   # Save data as CSV and XLSX
   csv_file_path <- file.path(root_dir, 'examples', 'data.csv')
   write.csv(data, csv_file_path, row.names = FALSE)
   
   xlsx_file_path <- file.path(root_dir, 'examples', 'data.xlsx')
   write_xlsx(as.data.frame(data), xlsx_file_path)
   ```

3. **View and Save Data:**
   
   ```r
   # View the first 10 rows of data
   data[1:10, ]
   ```

4. **Load Data Using Provided Functions:**
   
   If your data is saved as CSV or XLSX, you can load it with the provided functions:
   
   ```r
   UserName <- 'YourUserName'  # Replace with your username
   wd <- paste0('/Users/', UserName, '/Documents/Repositories/Rfits/examples')
   
   # Load CSV data
   data1 <- load_data(wd = wd, name = 'data')
   
   # Load XLSX data (if there is only one sheet, access it as data2$'Sheet1')
   data2 <- load_data2(wd = wd, name = 'data')
   ```

5. **View Imported Data:**
   
   ```r
   # View first 10 rows of CSV data
   data1[1:10, ]
   
   # View first 10 rows of XLSX data (first sheet)
   data2$'Sheet1'[1:10, ]
   ```

6. **Analyze a Single Column of Data:**
   
   Each column represents a single PSC sampled at 10 KHz (sample interval 0.1 ms).
   
   ```r
   # Analyze the first response (column) from data1
   out1 <- analyse_PSC(
      response = data1[, 1],
      dt = 0.1,
      func = product2N,
      stimulation_time = 150,
      baseline = 50
   )
   ```
   
   If no `fit.limits` is specified, a graph with reference lines will be displayed. For automatic processing, specify `fit.limits`:
   
   ```r
   out1 <- analyse_PSC(
      response = data1[, 1],
      dt = 0.1,
      func = product2N,
      stimulation_time = 150,
      baseline = 50,
      fit.limits = 510
   )
   ```
   
   Retrieve the output table by executing:
   
   ```r
   out1$output
   ```

7. **Analyze an Entire Dataset:**
   
   Loop through each column of data to analyze multiple responses:
   
   ```r
   # Define a vector of fit limits (one per response)
   time.limits <- c(508.4, 507.8, 510.5, 508.6, 508.1, 507.7, 502.7, 498.6, 499.7, 507.0)
   
   # Initialize an empty list for the results
   out_list <- list()
   
   # Loop through columns and analyze each
   for (ii in 1:ncol(data1)) {
      out_list[[ii]] <- analyse_PSC(
         response = data1[, ii],
         dt = 0.1,
         func = product2N,
         stimulation_time = 150,
         baseline = 50,
         fit.limits = time.limits[ii]
      )
   }
   
   # Save the entire environment (optional)
   setwd(wd)
   save.image(file = 'example.RData')
   ```

8. **Retrieve Previously Saved Data:**
   
   ```r
   rm(list = ls(all = TRUE))
   UserName <- 'YourUserName'
   root_dir <- paste0('/Users/', UserName, '/Documents/Repositories/Rfits')
   wd <- paste0(root_dir, '/examples')
   setwd(wd)
   
   load(file.path(wd, 'example.RData'))
   
   # Reload packages if necessary
   load_required_packages(required.packages)
   ```

9. **Examine the Analyzed Data (Summary):**
   
   Create a summary table by extracting outputs from each fit:
   
   ```r
   summary <- t(sapply(1:length(out_list), function(ii) {
      X <- out_list[[ii]]$output
      as.vector(t(X))
   }))
   
   new_colnames <- rep(colnames(out_list[[1]]$output), 2)
   colnames(summary) <- new_colnames
   rownames(summary) <- 1:length(out_list)
   
   summary
   ```

10. **Useful Functions for Analysis:**

    - **`wilcox.test`:**  
      Example for comparing fast and slow amplitudes:
      
      ```r
      wilcox.test(summary[, 1], summary[, 10], paired = TRUE, alternative = 'two.sided', exact = NULL)
      ```
    
    - **`BoxPlot`:**  
      Create a box plot from a formatted data frame:
      
      ```r
      n <- dim(summary)[1]
      A <- data.frame(
         s = rep(1:n, 2),
         x = c(rep(1, n), rep(2, n)),
         y = c(summary[, 1], summary[, 10])
      )
      BoxPlot(A, yrange = c(-120, 0), wid = 0.3, cap = 0.1, y_tick_interval = 20, 
         width = 3, height = 5, tick_length = 0.2, lwd = 1.25, amount = 0.1, p.cex = 0.6)
      ```
    
    - **`ScatterPlot`:**  
      Generate scatter plots for paired data, with optional regression and confidence intervals:
      
      ```r
      ScatterPlot(A, sign = -1, xlim = c(0, 120), ylim = c(0, 120), x_tick_interval = 20, y_tick_interval = 20, 
         xlab = expression(A[fast] * " " * (pA)), ylab = expression(A[slow] * " " * (pA)), 
         lwd = 1.25, p.cex = 0.6, width = 5, height = 5)
      
      ScatterPlot(A, sign = -1, xlim = c(40, 70), ylim = c(90, 110), x_tick_interval = 10, y_tick_interval = 5, 
         xlab = expression(A[fast] * " " * (pA)), ylab = expression(A[slow] * " " * (pA)), 
         lwd = 1.25, p.cex = 0.6, reg = TRUE, plot.CI = TRUE, reg.method = 'Theil-Sen', width = 5, height = 5)
      ```
    
    - **`SingleFitExample`:**  
      Plot fitted responses from a single analysis:
      
      ```r
      # For the 5th fitted response
      SingleFitExample(traces = out_list[[5]]$traces, xlim = c(0, 800), ylim = c(-140, 10), lwd = 1.5, 
         height = 5, width = 5, xbar = 100, ybar = 25, filename = 'egs_control_5.svg', save = FALSE)
      
      # Equivalent semilog plot
      SingleFitExample(traces = out_list[[5]]$traces, xlim = c(50, 400), ylim = NULL, lwd = 1.5, 
         height = 5, width = 5, xbar = 100, filename = 'semilog_egs_5.svg', log_y = TRUE, save = FALSE)
      ```

11. **Output File Structure Explained:**

    The `out_list` object is a list containing outputs for each analyzed trace. For any given trace (e.g., trace 3), you can access various elements:
    
    ```r
    # View available elements for trace 3
    names(out_list[[3]])
    
    # Access the output table for trace 3
    out_list[[3]]$output
    ```
    
    **Key elements include:**
    
    - **output:** Best-fit parameters.
    - **fits:** Parameter fits for the product equation form.
    - **fits.se:** Standard errors for the parameter fits.
    - **gof:** Goodness-of-fit (RMSE) of the model.
    - **AIC / BIC:** Model selection criteria.
    - **model.message:** Termination messages from the fitting algorithm.
    - **sign:** Sign indicator of the fitted response.
    - **traces:** Data traces (original, filtered, and fitted responses).
    - **fit_results:** A list of all individual fits (e.g., when N = 30 fits are performed).

---

## Definitions and Formulae

- The alternative form of the product equation used for fitting is:
  
  $begin:math:display$
  y = A \\, \\left(e^{-t/\\tau_{decay}} - e^{-t/\\tau_{rise}}\\right)
  $end:math:display$

- The version used for fitting is:
  
  $begin:math:display$
  y = A \\, \\left(1 - e^{-t/\\tau_1}\\right) \\, e^{-t/\\tau_2}
  $end:math:display$

- **Goodness-of-Fit (gof):**
  
  $begin:math:display$
  \\text{gof.se} = \\sqrt{\\frac{\\sum (y_i - \\hat{y}_i)^2}{n - k}}
  $end:math:display$

- **Model Selection:**  
  Lower AIC or BIC values indicate a model that best fits the data with minimal information loss.

---

Happy fitting!
````markdown

   
