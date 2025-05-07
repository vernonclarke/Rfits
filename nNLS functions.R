# NLS functions
# Rewrite residFun in C++ and SS.fun in C++


# # Define the C++ code as a string
# cpp_code <- '
# #include <Rcpp.h>
# using namespace Rcpp;

# // Function to calculate residuals
# // [[Rcpp::export]]
# NumericVector residFunCpp(NumericVector params, NumericVector y, NumericVector x, Function func, int N, double IEI) {
#   int n=y.size();
#   NumericVector residuals(n);
  
#   NumericVector fitted_values=as<NumericVector>(func(params, x, N, IEI));
  
#   for (int i=0; i < n; i++) {
#     residuals[i]=y[i] - fitted_values[i];
#   }
  
#   return residuals;
# }

# // Function to calculate sum of squares
# // [[Rcpp::export]]
# double SSfunCpp(NumericVector params, NumericVector x, NumericVector y, Function func, int N, double IEI) {
#   NumericVector residuals=residFunCpp(params, y, x, func, N, IEI);
#   return sum(pow(residuals, 2));
# }
# '

# # Define the C++ code as a string
# cpp_code <- '
# #include <Rcpp.h>
# using namespace Rcpp;

# // Function to calculate residuals
# // [[Rcpp::export]]
# NumericVector residFunCpp(NumericVector params, NumericVector y, NumericVector x, Function func, int N, double IEI) {
#   NumericVector fitted_values = as<NumericVector>(func(params, x, N, IEI));
#   return y - fitted_values;
# }
# // Function to calculate sum of squares
# // [[Rcpp::export]]
# double SSfunCpp(NumericVector params, NumericVector x, NumericVector y, Function func, int N, double IEI) {
#   NumericVector residuals=residFunCpp(params, y, x, func, N, IEI);
#   return sum(pow(residuals, 2));
# }
# '

# Define the C++ code as a string
cpp_code <- '
#include <Rcpp.h>
using namespace Rcpp;

// Function to calculate residuals with weight_method
// [[Rcpp::export]]
NumericVector residFunCpp(NumericVector params, NumericVector y, NumericVector x, Function func, int N, double IEI, Nullable<NumericVector> weights = R_NilValue, std::string weight_method = "none") {
  NumericVector fitted_values = as<NumericVector>(func(params, x, N, IEI));
  NumericVector residuals = y - fitted_values;
  
  // Apply weights based on the method
  if (weight_method == "~y_sqrt") {
    // Use the absolute value of y to avoid NaN issues
    NumericVector y_safe = abs(y);  // Ensure no negative values in y
    residuals = sqrt(y_safe) * residuals;  // Apply sqrt of y as weights
  } else if (weight_method == "~y") {
    residuals = y * residuals;  // 
  } else if (weights.isNotNull()) {
    NumericVector w = weights.get();
    residuals = sqrt(w) * residuals;  // Apply custom weights if provided
  }

  return residuals;
}

// Function to calculate weighted log-likelihood (posterior)
// [[Rcpp::export]]
double logLikPostCpp(NumericVector parameters, NumericVector x, NumericVector y, Function func, Nullable<NumericVector> weights = R_NilValue, std::string weight_method = "none") {
  // Extract sigma (last parameter) and remaining params
  NumericVector params = parameters[seq(0, parameters.size() - 2)];
  double sigma = parameters[parameters.size() - 1];

  // Get the fitted values using the provided function
  NumericVector fitted_values = as<NumericVector>(func(params, x));

  // Initialize log-likelihood
  double loglik = 0.0;

  // Apply weights based on the method or the provided weights vector
  NumericVector w;
  if (weights.isNotNull()) {
    w = weights.get();
  } else {
    w = NumericVector(y.size(), 1.0);  // Default to 1 if no weights are provided
  }

  if (weight_method == "~y_sqrt") {
    w = sqrt(abs(y));  
  } else if (weight_method == "~y") {
    w = abs(y);  // Apply y as weights
  }

  // Calculate the log-likelihood assuming normal distribution with weights
  int n = y.size();
  for (int i = 0; i < n; i++) {
    double diff = y[i] - fitted_values[i];
    loglik += w[i] * (-0.5 * (log(2 * M_PI) + 2 * log(sigma) + pow(diff / sigma, 2)));
  }

  return loglik;
}

// Function to calculate weighted log-likelihood (posterior) with N and IEI
// [[Rcpp::export]]
double logLikPostCpp_N(NumericVector parameters, NumericVector x, NumericVector y, Function func, int N = 1, double IEI = 100, Nullable<NumericVector> weights = R_NilValue, std::string weight_method = "none") {
  // Extract sigma (last parameter) and remaining params
  NumericVector params = parameters[seq(0, parameters.size() - 2)];
  double sigma = parameters[parameters.size() - 1];

  // Get the fitted values using the provided function with N and IEI
  NumericVector fitted_values = as<NumericVector>(func(params, x, N, IEI));

  // Initialize log-likelihood
  double loglik = 0.0;

  // Apply weights based on the method or the provided weights vector
  NumericVector w;
  if (weights.isNotNull()) {
    w = weights.get();
  } else {
    w = NumericVector(y.size(), 1.0);  // Default to 1 if no weights are provided
  }

  if (weight_method == "~y_sqrt") {
    w = sqrt(abs(y));  // Apply sqrt of y as weights
  } else if (weight_method == "~y") {
    w = abs(y);  // Apply y^2 as weights
  }

  // Calculate the log-likelihood assuming normal distribution with weights
  int n = y.size();
  for (int i = 0; i < n; i++) {
    double diff = y[i] - fitted_values[i];
    loglik += w[i] * (-0.5 * (log(2 * M_PI) + 2 * log(sigma) + pow(diff / sigma, 2)));
  }

  return loglik;
}

// Function to calculate sum of squares
// [[Rcpp::export]]
double SSfunCpp(NumericVector params, NumericVector x, NumericVector y, Function func, int N, double IEI, Nullable<NumericVector> weights = R_NilValue) {
  NumericVector residuals = residFunCpp(params, y, x, func, N, IEI, weights);
  return sum(pow(residuals, 2));  // Sum of squared residuals
}
'

# Write the C++ code to a temporary file
cpp_file <- tempfile(fileext=".cpp")
writeLines(cpp_code, cpp_file)

# Source the C++ file
sourceCpp(cpp_file)

FITN <- function(response, dt=0.1, func=product2N, N=1, IEI=50, method=c('BF.LM', 'LM', 'GN', 'port', 'robust', 'MLE'), 
                 weight_method=c('none', '~y_sqrt', '~y'), stimulation_time=0, baseline=0, fast.decay.limit=NULL, latency.limit=NULL, 
                 lower=NULL, upper=NULL, filter=FALSE, fc=1000, interval=c(0.1, 0.9), 
                 MLEsettings=list(iter=1e4, metropolis.scale=1.5, fit.attempts=100, RWm=FALSE), 
                 MLE.method=c('L-BFGS-B', 'Nelder-Mead', 'BFGS','CG', 'SANN', 'Brent'),
                 response_sign_method=c('smooth', 'regression', 'cumsum'), 
                 dp=3, lwd=1.2, xlab='time (ms)', ylab='PSC (pA)', width=5, height=5, 
                 return.output=FALSE, show.output=TRUE, show.plot=TRUE) {

  dx <- dt 
  method <- match.arg(method)
  weight_method <- match.arg(weight_method)
  y <- response
  if (all(is.na(y[(which(!is.na(y))[length(which(!is.na(y)))] + 1):length(y)]))) {
    y <- y[!is.na(y)]
  }

  # sign <- sign_fun(y, direction_method=response_sign_method) 
  # y <- sign * y

  x <- seq(0, (length(y) - 1) * dx, by=dx)

  if (filter) {
    ind=20
    fs=1 / dx * 1000; bf <- butter(2, fc / (fs / 2), type='low')
    yfilter <- signal::filter(bf, y)
  } else {
    ind=1
    yfilter=y
  }

  x.orig <- x
  ind1 <- (stimulation_time - baseline) / dx
  ind2 <- baseline / dx
  
  yorig <- y[ind1:length(y)]
  yfilter <- yfilter[ind1:length(yfilter)]
  xorig <- seq(0, dx * (length(yorig) - 1), by=dx)

  sign <- sign_fun(yfilter - mean(yfilter[1:ind2]), direction_method=response_sign_method) 
  y <- sign * y; yorig <- sign * yorig; yfilter <- sign * yfilter
  
  if (is_product_function(func)) {
    yorig <- yorig - mean(yorig[1:ind2])
    yfilter <- yfilter - mean(yfilter[1:ind2])

    y2fit <- yfilter[ind2:length(yfilter)]
    x2fit <- seq(0, dx * (length(y2fit) - 1), by=dx)

    upper <- adjust_product_bounds_N(bounds=upper, func=func, N=N, upper=TRUE)
    lower <- adjust_product_bounds_N(bounds=lower, func=func, N=N)

    if (is.null(upper) && (!is.null(fast.decay.limit) || !is.null(latency.limit))) {
      if (identical(func, product1N)){
        upper <- c(rep(sign * Inf, N), rep(Inf, 3)) #  N amplitudes plus tau1, tau2, delay
      } else if (identical(func, product2N)) {
        upper <- rep(c(rep(sign * Inf, N), rep(Inf, 3)),2)
      } else if (identical(func, product3N)) {
        upper <- rep(c(rep(sign * Inf, N), rep(Inf, 3)),3)
      }
    
      if (!is.null(fast.decay.limit)){
        upper[N+2] <- fast.decay.limit[1]
        if (identical(func, product2N)){
          upper[2*(N+2)+1] <-  if (length(fast.decay.limit)==1) fast.decay.limit[1] else fast.decay.limit[2]
        }
        if (identical(func, product3N)){
          upper[2*(N+2)+1] <-  if (length(fast.decay.limit)==1) fast.decay.limit[1] else fast.decay.limit[2]
          upper[2*(N+2)+2] <-  if (length(fast.decay.limit)==1) fast.decay.limit[1] else fast.decay.limit[3]
        }
      }

      if (!is.null(latency.limit)){
          if (identical(func, product1N)){
            upper[N+3] <- latency.limit
          } else if (identical(func, product2N)){
            upper[c(N+3, 2*(N+3))] <- latency.limit
          } else if (identical(func, product3N)){
            upper[c(N+3, 2*(N+3), 3*(N+3))] <- latency.limit
          }
        }

    }

    # Define the formulas for each of the product functions
    if (identical(func, product1N)) {
      param_names <- c(paste0('a', 1:N), 'tau1', 'tau2', 'delay')
      N.params <- length(param_names)
 
      if (!is.null(lower) && length(lower) != N.params) stop("number of lower boundaries should be equal to the number of parameters")
      if (!is.null(upper) && length(upper) != N.params) stop("number of upper boundaries should be equal to the number of parameters")
      if (!is.null(lower)) lower[1:N] <- sign * lower[1:N] 
      if (!is.null(upper)) upper[1:N] <- sign * upper[1:N] 
      
    } else if (identical(func, product2N)) {
      param_names <- c(paste0('a', 1:N), 'tau1', 'tau2', 'delay', 
                       paste0('a', 1:N, '2'), 'tau1_2', 'tau2_2', 'delay_2')
      N.params <- length(param_names)
 
      if (!is.null(lower) && length(lower) != N.params) stop("number of lower boundaries should be equal to the number of parameters")
      if (!is.null(upper) && length(upper) != N.params) stop("number of upper boundaries should be equal to the number of parameters")
      if (!is.null(lower)) { lower[1:N] <- sign * lower[1:N]; lower[(N+4):(2*N+3)] <- sign * lower[(N+4):(2*N+3)] }
      if (!is.null(upper)) { upper[1:N] <- sign * upper[1:N]; upper[(N+4):(2*N+3)] <- sign * upper[(N+4):(2*N+3)] }
     
    } else if (identical(func, product3N)) {
      param_names <- c(paste0('a', 1:N), 'tau1', 'tau2', 'delay', 
                       paste0('a', 1:N, '2'), 'tau1_2', 'tau2_2', 'delay_2',
                       paste0('a', 1:N, '3'), 'tau1_3', 'tau2_3', 'delay_3')
      N.params <- length(param_names)

      if (!is.null(lower) && length(lower) != N.params) stop("number of lower boundaries should be equal to the number of parameters")
      if (!is.null(upper) && length(upper) != N.params) stop("number of upper boundaries should be equal to the number of parameters")
      if (!is.null(lower)) { lower[1:N] <- sign * lower[1:N]; lower[(N+4):(2*N+3)] <- sign * lower[(N+4):(2*N+3)]; lower[(2*N+7):(3*N+6)] <- sign * lower[(2*N+7):(3*N+6)] }
      if (!is.null(upper)) { upper[1:N] <- sign * upper[1:N]; upper[(N+4):(2*N+3)] <- sign * upper[(N+4):(2*N+3)]; upper[(2*N+7):(3*N+6)] <- sign * upper[(2*N+7):(3*N+6)] }

    }
  }else{
    fun.2.fit <- func
    y2fit <- yfilter[ind2:length(yfilter)]
    x2fit <- seq(0, dx * (length(y2fit) - 1), by = dx)
    param_names <- generate_param_names(fun.2.fit)
    N.params <- length(param_names)
    form.2.fit <- generate_formula(fun.2.fit, param_names, y_name = 'y', x_name = 'x')
  }


 if (method == 'MLE'){
    
    MLE.method <- match.arg(MLE.method)
    iter <- MLEsettings$iter
    metropolis.scale <- MLEsettings$metropolis.scale
    fit.attempts <- MLEsettings$fit.attempts
    RWm <- MLEsettings$RWm  
    sd.est <- sqrt(mean(yfilter[1:ind2]^2))

    logpost <- if (is_product_function(func)) logLikPostCpp_N else logLikPostCpp

    output <- FIT.MLE(x=x2fit, y=y2fit, func=func, N.params=N.params, N=N, IEI=IEI, sigma=sd.est, iter=iter, 
      metropolis.scale=metropolis.scale, logpost=logpost, params.start=NULL, 
      method=MLE.method, lower=lower, upper=upper, weight_method=weight_method, MLE.fun.attempts=100, fit.attempts=fit.attempts, RWm=RWm) 
  
  } else if (method == 'LM') {
    
    output <- FIT.LM(x=x2fit, y=y2fit, func=func, N.params=N.params, N=N, IEI=IEI, lower=NULL, upper=NULL, 
      weight_method=weight_method, fit.convergence.attempts=10, fit.attempts=10)
  
  } else if (method == 'BF.LM') {
    
    # default method
    bounds <- check_and_set_bounds(x=x2fit, y=y2fit, func=func, N.params=N.params, upper=upper, lower=lower)
    lower <- bounds$lower
    upper <- bounds$upper
    # print(upper)
    
    output <- FIT.LM(x=x2fit, y=y2fit, func=func, N.params=N.params,  N=N, IEI=IEI, lower=lower, upper=upper, 
      weight_method=weight_method, fit.convergence.attempts=10, fit.attempts=10)

  } else if (method == 'GN') {
    
    output <- FIT.NLS(x=x2fit, y=y2fit, func=func, N.params=N.params, N=N, IEI=IEI, lower=NULL, upper=NULL, weight_method=weight_method, algorithm='default')

  } else if (method == 'port') {
    
    bounds <- check_and_set_bounds(x=x, y=y, func=func, N.params=N.params, upper=upper, lower=lower)
    upper <- bounds$upper
    lower <- if (method == 'port') bounds$lower else rep(-Inf, N.params)

    output <- FIT.NLS(x=x2fit, y=y2fit, func=func, N.params=N.params, N=N, IEI=IEI, lower=lower, upper=upper, weight_method=weight_method, algorithm='port')
    
  } else if (method == 'robust') {
    
    output <- FIT.NLSrobust(x=x2fit, y=y2fit, func=func, N.params=N.params, N=N, IEI=IEI, lower=lower, upper=upper, algorithm='port', method = 'M')

  }
    

  if (is_product_function(func)) {
    
    if (identical(func, product1N)){
      
      df_output <- out.fun(params=output$fits[1:(N+3)], interval=interval, dp=dp, sign=sign)
      output$fits[1:N] <- sign * output$fits[1:N]
    
    } else if (identical(func, product2N)){  
      df_output1 <- out.fun(params=output$fits[1:(N+3)], interval=interval, dp=dp, sign=sign)
      df_output2 <- out.fun(params=output$fits[(N+4):(2*N+6)], interval=interval, dp=dp, sign=sign)
            # Create a list of the outputs
      output_list <- list(df_output1, df_output2)
      # Extract the third elements from each output and determine the order
      order_indices <- order(sapply(output_list, function(x) x[[N+2]]))
      # Reorder the output list based on the third element
      output_ordered <- output_list[order_indices]

      # Combine the outputs in the correct order
      df_output <- rbind('fast' = output_ordered[[1]], 'slow' = output_ordered[[2]])

      # Adjust the order of the fits and fits.se based on the sorted order
      output$fits <- c(output$fits[((N+3)*order_indices[1]-(N+2)):((N+3)*order_indices[1])],
                       output$fits[((N+3)*order_indices[2]-(N+2)):((N+3)*order_indices[2])])

      output$fits.se <- c(output$fits.se[((N+3)*order_indices[1]-(N+2)):((N+3)*order_indices[1])],
                       output$fits.se[((N+3)*order_indices[2]-(N+2)):((N+3)*order_indices[2])])

      # adjust signs
      output$fits[1:N] <- sign * output$fits[1:N]
      output$fits[(N+4):(2*N+3)] <- sign * output$fits[(N+4):(2*N+3)]


    } else if (identical(func, product3N)){ 

      df_output1 <- out.fun(params=output$fits[1:(N+3)], interval=interval, dp=dp, sign=sign)
      df_output2 <- out.fun(params=output$fits[(N+4):(2*N+6)], interval=interval, dp=dp, sign=sign)
      df_output3 <- out.fun(params=output$fits[(2*N+7):(3*N+9)], interval=interval, dp=dp, sign=sign)

      # Create a list of the outputs
      output_list <- list(df_output1, df_output2, df_output3)
      # Extract the third elements from each output and determine the order
      order_indices <- order(sapply(output_list, function(x) x[[N+2]]))
      # Reorder the output list based on the third element
      output_ordered <- output_list[order_indices]

      # Combine the outputs in the correct order
      df_output <- rbind('fast' = output_ordered[[1]], 'medium' = output_ordered[[2]], 'slow' = output_ordered[[3]])

      # Adjust the order of the fits and fits.se based on the sorted order
      output$fits <- c(output$fits[((N+3)*order_indices[1]-(N+2)):((N+3)*order_indices[1])],
                       output$fits[((N+3)*order_indices[2]-(N+2)):((N+3)*order_indices[2])],
                       output$fits[((N+3)*order_indices[3]-(N+2)):((N+3)*order_indices[3])])

      output$fits.se <- c(output$fits.se[((N+3)*order_indices[1]-(N+2)):((N+3)*order_indices[1])],
                       output$fits.se[((N+3)*order_indices[2]-(N+2)):((N+3)*order_indices[2])],
                       output$fits.se[((N+3)*order_indices[3]-(N+2)):((N+3)*order_indices[3])])

      # adjust signs
      output$fits[1:N] <- sign * output$fits[1:N]
      output$fits[(N+4):(2*N+3)] <- sign * output$fits[(N+4):(2*N+3)]
      output$fits[(2*N+7):(3*N+6)] <- sign * output$fits[(2*N+7):(3*N+6)]
    }
  
    traces <- data.frame(x = xorig, y = sign * yorig, yfilter = sign * yfilter)
    fits <- output$fits

    # Define the list of functions
    func_list <- list(product1N, product2N, product3N)

    # Find the index of the matching function
    func_index <- sapply(func_list, function(f) identical(f, func))

    # Check if func matches one of the productN functions and apply the corresponding baseline adjustments
    if (any(func_index)) {
      index <- which(func_index)  # Get the index of the matching function
      for (i in 1:index) {
        fits[i * N + (3 * i)] <- fits[i * N + (3 * i)] + baseline
      }
    }

    traces$yfit <- func(params=fits, x=traces$x+dx, N=N, IEI=IEI) 

    if (identical(func, product2N)){
      traces$yfit1 <- product1N(params=fits[1:(N+3)], x=traces$x+dx, N=N, IEI=IEI) 
      traces$yfit2 <- product1N(params=fits[(N+4):(2*N+6)], x=traces$x+dx, N=N, IEI=IEI) 
    } 
    if (identical(func, product3N)){
      traces$yfit1 <- product1N(params=fits[1:(N+3)],  x=traces$x+dx, N=N, IEI=IEI) 
      traces$yfit2 <- product1N(params=fits[(N+4):(2*N+6)],  x=traces$x+dx, N=N, IEI=IEI) 
      traces$yfit3 <- product1N(params=fits[(2*N+7):(3*N+9)], x=traces$x+dx, N=N, IEI=IEI) 
    } 

  }else{
    traces=data.frame(x=xorig, y=sign * yorig, yfilter=sign * yfilter)
    fits <- output$fits
    traces$yfit <- func(fits,x2fit)  
    df_output <- data.frame(fits=fits)
    df_output <- t(df_output)

    df_output <- as.data.frame(df_output)
    colnames(df_output) <- param_names
  }

  if (show.plot) fit_plot(traces=traces, func=func, xlab=xlab, ylab=ylab, lwd=lwd, filter=filter, width=width, height=height)
  
  if (show.output){
    print(df_output)
  }

  if (return.output){
    return(list(output = df_output, fits = output$fits, fits.se = output$fits.se, gof = output$gof, AIC = output$AIC, BIC = output$BIC, model.message = output$model.message, sign=sign, traces=traces))
  }
}


eigen.values <- function(x, only.values = FALSE) {
  z <- .Internal(La_rg(x, only.values))
  ord <- order(Mod(z$values), decreasing = TRUE)
  z$values[ord]
}

residFun <- function(params, y, x, func) {
  y - func(params, x)
}

SS.fun <- function(params, x, y, func) {
  sum((y - func(params, x))^2)
}

# MLE function
MLE.fun <- function(logpost, par, x, y, N=1, IEI=100, gr=NULL, lower=NULL, upper=NULL, weight_method='none', method='Nelder-Mead', func=product1) {
  if (method %in% c('L-BFGS-B', 'Brent')) {
    suppressWarnings(
      MLE.fit <- optim(
        par=par, 
        fn=logpost, 
        gr=gr, 
        method=method, 
        hessian=TRUE, 
        lower=lower, 
        upper=upper, 
        control=list(fnscale=-1, maxit=2000, factr=1e7, pgtol=1e-8, parscale=rep(1, length(par))),
        x=x, 
        y=y, 
        func=func, 
        N=N,         
        IEI=IEI,
        weight_method=weight_method      
      )
    )
  } else {
    suppressWarnings(
      MLE.fit <- optim(
        par=par, 
        fn=logpost, 
        gr=gr, 
        method=method, 
        hessian=TRUE, 
        control=list(fnscale=-1, maxit=2000, reltol=1e-8, parscale=rep(1, length(par))),
        x=x, 
        y=y, 
        func=func, 
        N=N,         
        IEI=IEI,
        weight_method=weight_method      
      )
    )
  }

  fit <- MLE.fit$par
  h <- -solve(MLE.fit$hessian)
  p <- length(fit)
  int <- p / 2 * log(2 * pi) + 0.5 * log(det(h)) + logpost(fit, x, y, func, N, IEI)
  list(fit=fit, fit.se=sqrt(diag(h)), var=h, int=int, converge=MLE.fit$convergence == 0)
}

# Optimized RWmetropolis function
RWmetropolis <- function(logpost, cov.mat, scale, start, m, x, y, func, N=1, IEI=50) {
  pb <- length(start)
  Mpar <- matrix(0, m, pb)
  b <- matrix(t(start))
  lb <- logpost(start, x, y, func)
  a <- chol(cov.mat)
  accept <- 0
  
  random_samples <- scale * t(a) %*% matrix(rnorm(pb * m), nrow = pb)
  
  for (ii in 1:m) {
    bc <- b + random_samples[, ii]
    lbc <- logpost(parameters=t(bc), x=x, y=y, func=func, N=N, IEI=IEI) # lbc <- logpost(t(bc), x, y, func)
    prob <- exp(lbc - lb)
    if (!is.na(prob) && runif(1) < prob) {
      lb <- lbc
      b <- bc
      accept <- accept + 1
    }
    Mpar[ii, ] <- b
  }
  list(par = Mpar, accept = accept / m)
}

# Log-likelihood function
log.lik.post <- function(parameters, x, y, func) {
  params <- parameters[-length(parameters)]
  sigma <- parameters[length(parameters)]
  sum(dnorm(y, mean=func(params, x), sd=sigma, log=TRUE))
}

log.lik.post_N <- function(parameters, x, y, func, N=1, IEI=100) {
  params <- parameters[-length(parameters)]
  sigma <- parameters[length(parameters)]
  sum(dnorm(y, mean=func(params=params, x=x, N=N, IEI=IEI), sd=sigma, log=TRUE))
}

# model.selection.criteria <- function(coeffs, x, y, func, N=1, IEI=100) {
#   res <- residFunCpp(coeffs, y, x, func, N, IEI)
#   k <- length(coeffs)
#   n <- length(res)
#   # definition: log_likelihood <- -(n / 2) * log(2 * pi) - (n / 2) * log(sigma2) - (1 / (2 * sigma2)) * sum(res^2)
#   # log(sigma2) == sum(res^2) / n and (1 / (2 * sigma2)) * sum(res^2) == 0.5 * n
#   loglik <- -0.5 * (n * log(2 * pi) + n * log(sum(res^2) / n) + n)
#   # loglik <- 0.5 * (-n * (log(2 * pi) + 1 - log(n) + log(sum(res^2))))
#   df <- k + 1
#   BIC <- df * log(n) - 2 * loglik
#   AIC <- df * 2 - 2 * loglik
#   c(AIC=AIC, BIC=BIC)
# }

model.selection.criteria <- function(coeffs, x, y, func, N=1, IEI=100) {
  res <- residFunCpp(coeffs, y, x, func, N, IEI)
  k <- length(coeffs)
  n <- length(res)
  # definition: log_likelihood <- -(n / 2) * log(2 * pi) - (n / 2) * log(sigma2) - (1 / (2 * sigma2)) * sum(res^2)
  # log(sigma2) == sum(res^2) / n and (1 / (2 * sigma2)) * sum(res^2) == 0.5 * n
  loglik <- -0.5 * (n * log(2 * pi) + n * log(sum(res^2) / n) + n)
  # loglik <- 0.5 * (-n * (log(2 * pi) + 1 - log(n) + log(sum(res^2))))
  BIC <- k * log(n) - 2 * loglik
  AIC <- k * 2 - 2 * loglik
  c(AIC=AIC, BIC=BIC)
}

# Optimized fit.MLE function
fit.MLE <- function(x, y, func, N.params, N=1, IEI=100, sigma=5, iter=1e4, metropolis.scale=2, logpost=logLikPostCpp, 
  params.start=NULL, method='Nelder-Mead', lower=NULL, upper=NULL, weight_method='none', MLE.fun.attempts=100, RWm=TRUE) {
  
  if (is.null(params.start)) {
     start <- start_optimization(x=x, y=y, func=func, N.params=N.params, SS.fun=SS.fun, lower=rep(0, N.params), upper=upper)
     params.start <- c(start, sigma)
  }
  n <- length(y)
  k <- length(params.start) - 1
  attempts <- 0
  
  l.fit <- list(int=NaN, converge=FALSE)
  while (!(!is.nan(l.fit$int) && l.fit$converge) && (attempts < MLE.fun.attempts)) {
    st1 <- if (is_product_function(func)) params.start * runif(k + 1) else params.start
    l.fit <- suppressWarnings(
      MLE.fun(logpost=logpost, par=st1, x=x, y=y, N=N, IEI=IEI, lower=lower, upper=upper, weight_method=weight_method, method=method, func=func)
    )
    attempts <- attempts + 1
  }
  
  if (!l.fit$converge) stop('MLE.fun does not converge')
  
  e.values <- Re(eigen.values(l.fit$var))
  if (any(e.values <= 0)) stop('MLE.fun covariance matrix not positive definite')

  if (RWm) {
    rw.fit <- RWmetropolis(logpost=logpost, cov.mat=l.fit$var, scale=metropolis.scale, start=l.fit$fit, m=iter, x=x, y=y, func=func, N=N, IEI=IEI)
    
    acceptance.rate <- rw.fit$accept
    parameters <- apply(rw.fit$par, 2, median)
    
    fits <- parameters[1:k]
    fits.se <- apply(rw.fit$par, 2, sd)[1:k]
  } else {
    fits <- l.fit$fit[1:k]
    fits.se <- sqrt(diag(l.fit$var))[1:k]
  }
  
  paramsfit <- fits
  res <- residFunCpp(paramsfit, y, x, func, N, IEI)
  gof.se <- (sum(res^2) / (n - k))^0.5

  msc <- model.selection.criteria(coeffs=paramsfit, x=x, y=y, func, N=N, IEI=IEI)

  if (RWm) {
    list(fits=fits, fits.se=fits.se, acceptance.rate=acceptance.rate, MLE.fun.fit=l.fit$fit, MLE.fun.fit.se=l.fit$fit.se, MLE.fun.var=l.fit$var, MLE.fun.int=l.fit$int, MLE.fun.converge=l.fit$converge, MLE.fun.convergence.attempts=attempts, gof=gof.se, AIC=msc[1], BIC=msc[2])
  } else {
    list(fits=fits, fits.se=fits.se, gof=gof.se, AIC=msc[1], BIC=msc[2], model.info=l.fit$converge, model.message=l.fit$converge)
  }
}

# Wrapper for fit.MLE
FIT.MLE <- function(x, y, func, N.params, N=1, IEI=100, sigma=5, iter=1e4, metropolis.scale=2, logpost=logLikPostCpp_N, params.start=NULL, 
  method='Nelder-Mead', lower=NULL, upper=NULL,  weight_method='none', MLE.fun.attempts=100, fit.attempts=10, RWm=TRUE) {

  bounds <- check_and_set_bounds(x=x, y=y, func=func, N.params=N.params, upper=upper, lower=lower)
  lower <- bounds$lower
  upper <- bounds$upper # some solutions were sitting on this limit so increased 2-fold
  if (method %in% c('L-BFGS-B', 'Brent')) {
    tol_low <- 1e-08
    lower <- if (any(lower <= 0)) lower + tol_low
  }

  output.MLE <- NULL
  attempts <- 0
  while (is.null(output.MLE) && attempts < fit.attempts) {
    tryCatch({
      output.MLE <- fit.MLE(x=x, y=y, func=func, N.params=N.params, N=N, IEI=IEI, sigma=sigma, iter=iter, metropolis.scale=metropolis.scale, logpost=logpost, params.start=params.start, 
        method=method, lower=lower, upper=upper, weight_method=weight_method, MLE.fun.attempts=MLE.fun.attempts, RWm=RWm)
    }, error=function(e) {})
    attempts <- attempts + 1
  }
  c(output.MLE, fit.attempts=attempts)
}

# Function to try optim with N retries
optim.N <- function(N.params, SS.fun, lower, upper, y, x, func, max_attempts=10) {
  attempt <- 1
  success <- FALSE
  est <- NULL

  while (attempt <= max_attempts && !success) {
    tryCatch({
      est <- optim(runif(N.params), SS.fun, method='L-BFGS-B', lower=lower, upper=upper, hessian=TRUE, y=y, x=x, func=func)
      success <- TRUE
    }, error = function(e) {
      if (grepl("L-BFGS-B needs finite values of 'fn'", e$message)) {
        attempt <- attempt + 1
      } else {
        stop(e)
      }
    })
  }

  if (!success) {
    stop("Failed to find a valid solution with 'L-BFGS-B' method after ", max_attempts, " attempts.")
  }

  return(est)
}

# Function to try optim with N retries
start_optimiser <- function(N.params, SS.fun, lower, upper, y, x, func, max_attempts=10, start.method=c('uniform', 'lognormal'), cv=0.4){
  start.method <- match.arg(start.method)
  if (identical(func, product3) || identical(func, product3N)) {
    N <- (N.params - 9)/3
    ests <- ests_fun(x=x, y=y, N=N, showplot=FALSE)
    pars <- c(rep(ests[1], N), ests[2],  ests[3]) 
    pars1 <- c(0.25 * pars, 15)
    if (!is.null(upper)) pars1 <- pmin(pars1, upper[1:(N+3)])
    pars2 <- c(pars, 15)
    if (!is.null(upper)) pars2 <- pmin(pars2, upper[(N+4):(2*(N+3))]) 
    pars3 <- c(pars, 15)
    if (!is.null(upper)) pars3 <- pmin(pars3, upper[(2*(N+3)+1):(3*(N+3))]) 
    if (start.method=='uniform'){
      par <- (c(pars1, pars2, pars3) * runif(N.params))
    }
    else if (start.method=='lognormal'){
      par <- generate_lognormal_samples(means=c(pars1, pars2, pars3), cv=cv, n=1)
    }

  } else if (identical(func, product2) || identical(func, product2N)) {
    N <- (N.params - 6)/2
    ests <- ests_fun(x=x, y=y, N=N, showplot=FALSE)
    pars <- c(rep(ests[1], N), ests[2],  ests[3]) 
    pars1 <- c(0.25 * pars, 15)
    if (!is.null(upper)) pars1 <- pmin(pars1, upper[1:(N+3)])
    pars2 <- c(pars, 15)
    if (!is.null(upper)) pars2 <- pmin(pars2, upper[(N+4):(2*(N+3))]) 
    if (start.method=='uniform'){
      par <- (c(pars1, pars2) * runif(N.params))
    }
    else if (start.method=='lognormal'){
      par <- generate_lognormal_samples(means=c(pars1, pars2), cv=cv, n=1)
    }

  } else if (identical(func, product1) || identical(func, product1N)) {
    N <- (N.params - 3)
    ests <- ests_fun(x=x, y=y, N=N, showplot=FALSE)
    pars1 <- c(rep(ests[1], N), ests[2],  ests[3], 15) 
    if (!is.null(upper)) pars1 <- pmin(pars1, upper[1:(N+3)])
    if (start.method=='uniform'){
      par <- pars1 * runif(N.params)
    }
    else if (start.method=='lognormal'){
      par <- generate_lognormal_samples(means=pars1, cv=cv, n=1)
    }

  }
  est <- list(par=par)
  return(est)
}

# function to check if a given function is product1 or product2
is_product_function <- function(func) {
  valid_funcs <- list(product1, product2, product3, product1N, product2N, product3N)
  any(vapply(valid_funcs, identical, logical(1), func))
}

check_and_set_bounds <- function(x, y, func, N.params, upper=NULL, lower=NULL) {
  if (is_product_function(func)) {
    if (is.null(upper)) {
      if (identical(func, product2) || identical(func, product2N)) {
        N <- (N.params - 6)/2
        ests <- ests_fun(x=x, y=y, N=N)
        ests <- rep(c(rep(ests[1], N), Inf,  ests[3], Inf) ,2) 
        upper <- sapply(6 * ests, round_up)
      } else if (identical(func, product1) || identical(func, product1N)) {
        N <- N.params - 3
        ests <- ests_fun(x=x, y=y, N=N)
        ests <- c(rep(ests[1], N), Inf,  ests[3], Inf)
        upper <- sapply(6 * ests, round_up)
      } else if (identical(func, product3) || identical(func, product3N)){
        N <- (N.params - 9)/3
        ests <- ests_fun(x=x, y=y, N=N)
        ests <- rep(c(rep(ests[1], N), Inf,  ests[3], Inf) ,3) 
        upper <- sapply(6 * ests, round_up)
      }  
    }
    if (is.null(lower)) lower <- rep(0, N.params)
  } else {
    upper <- rep(Inf, N.params)
    lower <- rep(-Inf, N.params)
  }
  return(list(upper=upper, lower=lower))
}

start_optimization <- function(x, y, func, N.params, SS.fun, lower, upper) {
  if (is_product_function(func)) {
    est <- start_optimiser(N.params = N.params, SS.fun = SS.fun, lower = lower, upper = upper, y = y, x = x, func = func, max_attempts = 10)
    start <- est$par
  } else {
    est <- optim(par = runif(N.params), SS.fun, method = "L-BFGS-B", lower= lower, upper = upper, hessian = TRUE, y = y, x = x, func = func) 
    start <- est$par
  }
  return(start)
}

fit.LM <- function(x, y, func,  N.params, N=1, IEI=100, lower=NULL, upper=NULL, weight_method='none') {
  # est <- optim(runif(N.params), SS.fun, method='L-BFGS-B', lower=lower, upper=upper, hessian=TRUE, y=y, x=x, func=func)
 
  start <- start_optimization(x=x, y=y, func=func, N.params=N.params, SS.fun=SS.fun, lower=lower, upper=upper)
  model <- minpack.lm::nls.lm(
    par=start, 
    fn=residFunCpp, 
    y=y, 
    x=x, 
    lower=lower, 
    upper=upper, 
    control=list(maxiter=1024, tol=1e-05, warnOnly=TRUE),
    func=func,  
    N=N,             
    IEI=IEI,
    weight_method=weight_method     
  )
  e.values <- Re(eigen.values(model$hessian))
  if (any(e.values <= 0)) stop('Hessian not positive definite')
  fits <- summary(model)$coefficients[, 'Estimate']
  fits.se <- summary(model)$coefficients[, 'Std. Error']
  gof.se <- summary(model)$sigma
  msc <- model.selection.criteria(coeffs=fits, x=x, y=y, func=func, N=N, IEI=IEI)
  list(fits=fits, fits.se=fits.se, gof=gof.se, AIC=msc[1], BIC=msc[2], model.info=model$info, model.message=model$message)
}

round_up <- function(x, factor=10) factor * ceiling( x / factor)
  
# Wrapper for fit.LM
FIT.LM <- function(x, y, func, N.params, N=1, IEI=100, lower=NULL, upper=NULL, weight_method='none', fit.convergence.attempts=10, fit.attempts=10) {
  convergence.attempts <- 0
  modelinfo <- NULL
  
  while ((is.null(modelinfo) || !any(modelinfo == c(1, 2, 3, 4))) && convergence.attempts < fit.convergence.attempts) {
    output.LM <- NULL
    attempts <- 0
    while (is.null(output.LM) && attempts < fit.attempts) {
      tryCatch({
        output.LM <- fit.LM(x=x, y=y, func=func, N.params=N.params, N=N, IEI=IEI, lower=lower, upper=upper, weight_method=weight_method)
      }, error=function(e) {})
      attempts <- attempts + 1
    }
    modelinfo <- output.LM$model.info
    convergence.attempts <- convergence.attempts + 1
  }
  c(output.LM, convergence.attempts=convergence.attempts)
}

# FIT.bf.LM <- function(x, y, func, N.params, N=1, IEI=100, lower=NULL, upper=NULL, fit.convergence.attempts = 10, fit.attempts = 10) {
#   convergence.attempts <- 0
#   modelinfo <- NULL

#   # bounds <- check_and_set_bounds(x=x, y=y, func=func, N.params=N.params, upper=upper, lower=lower)
#   # lower <- bounds$lower
#   # upper <- bounds$upper

#   while ((is.null(modelinfo) || !any(modelinfo == c(1, 2, 3, 4))) && convergence.attempts < fit.convergence.attempts) {
#     output.LM <- NULL
#     attempts <- 0
#     while (is.null(output.LM) && attempts < fit.attempts) {
#       tryCatch({
#         output.LM <- fit.bf.LM(x=x, y=y, func=func, N.params=N.params, N=N, IEI=IEI, lower=lower, upper=upper)
#       }, error = function(e) {})
#       attempts <- attempts + 1
#     }
#     modelinfo <- output.LM$model.info
#     convergence.attempts <- convergence.attempts + 1
#   }
#   c(output.LM, convergence.attempts = convergence.attempts)
# }

# # Wrapper for FIT.bf.LM
# FIT.BF.LM <- function(x, y, func, N.params, lower = NULL, upper = NULL, fit.convergence.attempts = 10, fit.attempts = 10, N.repeat = 25) {
#   fits <- matrix(NA, N.repeat, N.params)
#   fits.se <- matrix(NA, N.repeat, N.params)
#   gof.se <- rep(NA, N.repeat)
#   BIC <- rep(NA, N.repeat)
#   AIC <- rep(NA, N.repeat)
#   fit.convergence.attempts <- rep(NA, N.repeat)
#   model.info <- rep(NA, N.repeat)
#   model.message <- rep(NA, N.repeat)
#   convergence.attempts <- rep(NA, N.repeat)
  
#   for (iii in 1:N.repeat) {
#     output1 <- FIT.bf.LM(x=x, y=y, func=func, N.params=N.params, lower=lower, upper=upper, fit.convergence.attempts=fit.convergence.attempts, fit.attempts=fit.attempts)
#     if (!is.null(output1)) {
#       fits[iii, ] <- output1$fits
#       fits.se[iii, ] <- output1$fits.se
#       gof.se[iii] <- output1$gof
#       AIC[iii] <- output1$AIC
#       BIC[iii] <- output1$BIC
#       fit.convergence.attempts[iii] <- output1$convergence.attempts
#       model.info[iii] <- output1$model.info
#       model.message[iii] <- output1$model.message
#       convergence.attempts[iii] <- output1$convergence.attempts
#     }
#   }
  
#   ind <- which.min(gof.se)
  
#   list(fits = fits[ind, ], fits.se = fits.se[ind, ], gof = gof.se[ind], AIC = AIC[ind], BIC = BIC[ind], model.info = model.info[ind], model.message = model.message[ind], convergence.attempts = convergence.attempts[ind])
# }

# # Algorithm = c('default', 'port')
# fit.NLS <- function(x, y, func,  N.params, N=1, IEI=100, lower=NULL, upper=NULL, algorithm='default') {

#   if (is_product_function(func)){
#     if (algorithm == 'port'){ 
#       start <- start_optimization(x=x, y=y, func=func, N.params=N.params, SS.fun=SS.fun, lower=lower, upper=upper)
#     }else{          # default is Gauss-Newton
#       est <- start_optimiser(N.params, SS.fun=SS.fun, lower=rep(0, N.params), upper=upper, y=y, x=x, func=func, max_attempts=10)
#       start <- est$par
#     }
#   }else{
#     start <- start_optimization(x=x, y=y, func=func, N.params=N.params, SS.fun=SS.fun, lower=lower, upper=upper)   
#   }
  
#   # names(start) <- param.names[!param.names == 'x']
  
#   # must declare form.2.fit formula in local environment
  
#   form.2.fit <- as.formula(y ~ func(params, x, N, IEI))
  
#   model <- suppressWarnings(
#           nls(
#             form.2.fit, 
#             start=list(params = start),
#             lower=lower,
#             upper=upper,
#             algorithm=algorithm
#             ) 
#           )

#   e.values <- Re(eigen.values(summary(model)$cov.unscaled))
#   if (any(e.values <= 0)) stop('Hessian not positive definite')
  
#   fits <- summary(model)$coefficients[, 'Estimate']
#   fits.se <- summary(model)$coefficients[, 'Std. Error']
#   gof.se <- summary(model)$sigma
#   msc <- model.selection.criteria(coeffs=fits, x=x, y=y, func, N=N, IEI=IEI)
#   list(fits=unname(fits), fits.se=unname(fits.se), gof=gof.se, AIC=msc[1], BIC=msc[2], model.info=model$convInfo$stopCode, model.message=model$convInfo$stopMessage)
# }

# Algorithm = c('default', 'port')
fit.NLS <- function(x, y, func,  N.params, N=1, IEI=100, lower=NULL, upper=NULL, weight_method='none', algorithm='default') {

  if (is_product_function(func)){
    if (algorithm == 'port'){ 
      start <- start_optimization(x=x, y=y, func=func, N.params=N.params, SS.fun=SS.fun, lower=lower, upper=upper)
    }else{          # default is Gauss-Newton
      est <- start_optimiser(N.params, SS.fun=SS.fun, lower=rep(0, N.params), upper=upper, y=y, x=x, func=func, max_attempts=10)
      start <- est$par
    }
  }else{
    start <- start_optimization(x=x, y=y, func=func, N.params=N.params, SS.fun=SS.fun, lower=lower, upper=upper)   
  }
  
  # weights
  if (weight_method == 'none') {
    weights <- rep(1, length(y)) 
  } else if (weight_method == '~y_sqrt') {
    weights <- abs(y)  
  } else if (weight_method == '~y') {
    weights <- y^2
  }
  
  # must declare form.2.fit formula in local environment
  
  form.2.fit <- as.formula(y ~ func(params, x, N, IEI))
  
  model <- suppressWarnings(
          nls(
            form.2.fit, 
            start=list(params = start),
            lower=lower,
            upper=upper,
            weights=weights,
            algorithm=algorithm
            ) 
          )

  e.values <- Re(eigen.values(summary(model)$cov.unscaled))
  if (any(e.values <= 0)) stop('Hessian not positive definite')
  
  fits <- summary(model)$coefficients[, 'Estimate']
  fits.se <- summary(model)$coefficients[, 'Std. Error']
  gof.se <- summary(model)$sigma
  msc <- model.selection.criteria(coeffs=fits, x=x, y=y, func, N=N, IEI=IEI)
  list(fits=unname(fits), fits.se=unname(fits.se), gof=gof.se, AIC=msc[1], BIC=msc[2], model.info=model$convInfo$stopCode, model.message=model$convInfo$stopMessage)
}

# Wrapper for fit.NLS
FIT.NLS <- function(x, y, func,  N.params, N=1, IEI=100, lower=NULL, upper=NULL, weight_method='none', algorithm='default', fit.convergence.attempts=10, fit.attempts=10) {
  convergence.attempts <- 0
  modelinfo <- NULL
  
  while ((is.null(modelinfo) || !any(modelinfo == c(0, 1, 2, 3, 4))) && convergence.attempts < fit.convergence.attempts) {
    output.NLS <- NULL
    attempts <- 0
    while (is.null(output.NLS) && attempts < fit.attempts) {
      tryCatch({
        output.NLS <- fit.NLS(x=x, y=y, func=func, N.params=N.params, N=N, IEI=IEI, lower=lower, upper=upper, weight_method=weight_method, algorithm=algorithm)
      }, error = function(e) {})
      attempts <- attempts + 1
    }
    modelinfo <- output.NLS$model.info
    convergence.attempts <- convergence.attempts + 1
  }
  c(output.NLS, convergence.attempts = convergence.attempts)
}


fit.NLSrobust <- function(x, y, func,  N.params, N=1, IEI=100, lower=NULL, upper=NULL, algorithm='default', method = 'M') {

  if (method == 'M') {
    if (is_product_function(func)){
      start <- start_optimization(x=x, y=y, func=func, N.params=N.params, SS.fun=SS.fun, lower=lower, upper=upper)

      # Define the formula, but only in terms of params and x
      func_wrapper <- function(params, x) {
        func(params, x, N, IEI)  # Pass N and IEI directly here
      }
      form.2.fit <- as.formula(y ~ func_wrapper(params, x))

    }else{   
      est <- start_optimiser(N.params, SS.fun=SS.fun, lower=rep(0, N.params), upper=upper, y=y, x=x, func=func, max_attempts=10)
      start <- est$par
      form.2.fit <- as.formula(y ~ func_wrapper(params, x))
    }


    control=robustbase:::nlrob.control(method)

    model <- suppressWarnings(
          nlrob(
            formula=form.2.fit, 
            start = list(params = start),
            data=data.frame(x=x, y=y),
            algorithm=algorithm, 
            method=method, 
            control=control, 
            lower=lower, 
            upper=upper,
            )
          )  

  }else{
    start <- start_optimization(x=x, y=y, func=func, N.params=N.params, SS.fun=SS.fun, lower=lower, upper=upper)   

    model <- suppressWarnings(
          robustbase::nlrob(
            formula=form.2.fit, 
            data=data.frame(x=x, y=y),
            algorithm=algorithm, 
            control=control, 
            lower=lower, 
            upper=upper
            )
          )  

    }

  e.values <- Re(eigen.values(vcov(model)))
  if (any(e.values <= 0)) stop('Hessian not positive definite')
  
  fits <- summary(model)$coefficients[, 'Estimate']
  fits.se <- summary(model)$coefficients[, 'Std. Error']
  gof.se <- summary(model)$Scale
  msc <- model.selection.criteria(coeffs=fits, x=x, y=y, func=func, N=N, IEI=IEI)
  list(fits=unname(fits), fits.se=unname(fits.se), gof=gof.se, AIC=msc[1], BIC=msc[2], model.message=summary(model)$status)
}

# Wrapper for fit.NLSrobust
FIT.NLSrobust <- function(x, y, func, N.params, N=1, IEI=100, lower=NULL, upper=NULL, algorithm = 'default', method = 'M', fit.convergence.attempts = 10, fit.attempts = 10) {
  
  convergence.attempts <- 0
  modelinfo <- NULL

  while ((is.null(modelinfo) || modelinfo != 'converged') && convergence.attempts < fit.convergence.attempts) {
    output.NLS <- NULL
    attempts <- 0
    while (is.null(output.NLS) && attempts < fit.attempts) {
      {
        output.NLS <- fit.NLSrobust(x=x, y=y, func=func, N.params=N.params, N=N, IEI=IEI, lower=lower, upper=upper, algorithm=algorithm, method=method)
      }
      attempts <- attempts + 1
    }
    modelinfo <- output.NLS$model.message
    convergence.attempts <- convergence.attempts + 1
  }
  c(output.NLS, convergence.attempts = convergence.attempts)
}

sort2.fun <- function(x){
  index = seq(1:length(x))
  if (any(is.na(x))) {
    return(x)
  }
  if (x[4] > x[2]) {
    x <- x[c(3,4,1,2,index[5:length(x)])]
  }
  x
}

sort2 <- function(mat, col.names=NULL){
  mat <- t(apply(mat, 1, function(row) {
    if (any(is.na(row))) {
      return(row)
    } else {
      return(sort2.fun(row))
    }
  }))
  colnames(mat) <- col.names
  mat
}

sort3.fun <- function(x){
  index=seq(1:length(x))
  if ( (x[2] > x[4]) && (x[4] > x[6]) ) {
    x <- x[c(1,2,3,4,5,6,7)]
  } else if ( (x[2] > x[6]) && (x[6] > x[4]) ) {
    x <- x[c(1,2,5,6,3,4,7)]
  } else if ( (x[4] > x[2]) && (x[2] > x[6]) ) {
    x <- x[c(3,4,1,2,5,6,7)]
  } else if ( (x[4] > x[6]) && (x[6] > x[2]) ) {
    x <- x[c(3,4,5,6,1,2,7)]
  } else if ( (x[6] > x[2]) && (x[2] > x[4]) ) {
    x <- x[c(5,6,1,2,3,4,7)]
  } else if ( (x[6] > x[4]) && (x[4] > x[2]) ) {
    x <- x[c(5,6,3,4,1,2,7)]
  }
  x
}

sort3 <- function(mat, col.names=NULL){
  mat <- t(apply(mat,1,sort3.fun))
  colnames(mat) <- col.names
  mat
  }

generate_formula <- function(fun, param_names, y_name = 'yfilter', x_name = 'x') {
  # Extract the body of the function as a character string
  fun_body <- deparse(body(fun))
  
  # Replace 'params' with the corresponding parameter names
  for (i in seq_along(param_names)) {
    fun_body <- gsub(paste0('params\\[', i, '\\]'), param_names[i], fun_body)
  }
  
  # Remove the function declaration and braces
  fun_body <- gsub('\\{', '', fun_body)
  fun_body <- gsub('\\}', '', fun_body)
  fun_body <- trimws(fun_body)
  fun_body <- paste(fun_body, collapse = ' ')
  
  # Create the formula
  formula_str <- paste(y_name, '~', fun_body)
  formula <- as.formula(formula_str)
  
  return(formula)
}

# generate parameter names
generate_param_names <- function(fun) {
  # Get the function body
  fun_body <- body(fun)
  
  # Convert the function body to a character vector
  fun_body_char <- deparse(fun_body)
  
  # Find all the params indices used in the function
  param_indices <- regmatches(fun_body_char, gregexpr("params\\[[0-9]+\\]", fun_body_char))
  param_indices <- unlist(param_indices)
  param_indices <- as.numeric(gsub("params\\[([0-9]+)\\]", "\\1", param_indices))
  
  # Handle any NAs or duplicates
  param_indices <- unique(param_indices[!is.na(param_indices)])
  
  # Generate names for the parameters
  param_names <- letters[seq_along(param_indices)]
  
  # Return a named vector of parameter names
  setNames(param_names, paste0("params[", param_indices, "]"))
}

tau_rise <- function(tau1, tau2) tau1 * tau2 / ( tau1 + tau2 )

tau1_fun <- function(tau_rise, tau_decay) tau_rise * tau_decay / ( tau_decay - tau_rise )

# Define the function to solve
find_tau1 <- function(tpeak, tau2) {
  # Define the equation to solve
  equation <- function(tau1) {
    tau1 * (exp(tpeak / tau1) - 1) - tau2
  }
  
  # Use uniroot to solve for tau1, give an interval of reasonable guesses
  result <- uniroot(equation, interval = c(1e-6, 100), tol = 1e-9)
  
  # Return the root found for tau1
  return(result$root)
}

product_area <- function(Apeak, tau1, tau2){
  f <- ((tau1 / (tau1 + tau2)) ^ (tau1 / tau2)) * tau2 / (tau1 + tau2) 
  abs(Apeak / f * tau2^2 / (tau1 + tau2))
}

product_tpeak <- function(tau1, tau2) tau1 * log( ( tau1 + tau2 ) / tau1 )

# calculate interval 10-90% rise and 90-10% decay based on default interval <- c(0.1, 0.9) from tau1 and tau2
product_rise_and_decay_percent <- function(tau1, tau2, interval=c(0.1, 0.9), showplot=FALSE, verbose=FALSE) {
  # Calculate tpeak
  tpeak <- product_tpeak(tau1, tau2)

  # Calculate f
  f <- ((tau1 / (tau1 + tau2)) ^ (tau1 / tau2)) * tau2 / (tau1 + tau2)

  # Define the function to find the root of
  target_function <- function(x, p) {
    exp(-x / tau2) - exp(-x * (tau1 + tau2) / (tau1 * tau2)) - f * p
  }

  # Function to find roots
  find_roots <- function(p, end_time) {
    # Define the intervals based on the plot
    interval1 <- c(0, tpeak)  # Interval before the peak
    interval2 <- c(tpeak, end_time)  # Interval after the peak
    
    # Find the roots
    root1 <- tryCatch(uniroot(function(x) target_function(x, p), interval1)$root, 
                      error = function(e) NA)
    root2 <- tryCatch(uniroot(function(x) target_function(x, p), interval2)$root, 
                      error = function(e) NA)
    
    if (is.na(root1) | is.na(root2)) {
      if (verbose) {
        cat("Failed to find roots for p =", p, "\n")
        cat("Interval1:", interval1, "\n")
        cat("Interval2:", interval2, "\n")
      }
    }
    
    return(c(root1, root2))
  }

  # Determine the end time of the signal
  end_time <- 10 * max(tau1, tau2) # This is a heuristic; adjust based on your signal's characteristics

  # Find roots for both p values and plot if required
  if (showplot) {
    x <- seq(0, end_time, length.out = 1000)
    yplot <- exp(-x / tau2) - exp(-x * (tau1 + tau2) / (tau1 * tau2))

    plot(x, yplot, type = 'l', xlab = 'x', ylab = 'f(x)', bty='l', las=1, main = '')
    abline(h = interval[1] * f, col = 'Indianred', lty = 3)
    abline(h = interval[2] * f, col = 'Indianred', lty = 3)

    # Add text labels
    text(x = max(x) * 0.95, y = interval[1] * f, labels = 'p=0.1', pos = 3, col = 'Indianred')
    text(x = max(x) * 0.95, y = interval[2] * f, labels = 'p=0.9', pos = 3, col = 'Indianred')
  }

  # Find roots for both p values
  roots_lower <- find_roots(interval[1], end_time)
  roots_upper <- find_roots(interval[2], end_time)

  # Compute the absolute differences
  diffs <- abs(roots_lower - roots_upper)

  return(diffs)
}

# product1 <- function(params, x) {
#   a1_max <- params[1]
#   tau_rise <- params[2]
#   tau2 <- params[3]
#   delay <- params[4]
  
#   tau1 <- tau1_fun(tau_rise, tau2) # tau2 = tau_decay

#   x_adjusted <- pmax(0, x - delay)
#   f <- ((tau1 / (tau1 + tau2)) ^ (tau1 / tau2)) * tau2 / (tau1 + tau2)
#   a1 <- a1_max / f
#   a1 * (1 - exp(-x_adjusted / tau1)) * exp(-x_adjusted / tau2)
# }

product1 <- function(params, x) {
  a1_max <- params[1]
  tau1 <- params[2]
  tau2 <- params[3]
  delay <- params[4]
  
  x_adjusted <- pmax(0, x - delay)
  f <- ((tau1 / (tau1 + tau2)) ^ (tau1 / tau2)) * tau2 / (tau1 + tau2)
  a1 <- a1_max / f
  a1 * (1 - exp(-x_adjusted / tau1)) * exp(-x_adjusted / tau2)
}

product2 <- function(params, x) {
  product1(params[1:4], x) + product1(params[5:8], x)
}

# product1N <- function(params, x, N=1, IEI=50) {
#   # Extract parameters
#   amplitudes <- params[1:N]  
#   shared_params <- params[(N + 1):length(params)]  
#   out <- sapply(1:N, function(ii){
#     shifted_x <- x - (ii - 1) * IEI 
#     current_params <- c(amplitudes[ii], shared_params)  
#     product1(current_params, shifted_x)  
#     })
#   response <- apply(out, 1, sum) 
#   return(response)
# }

# product2N <- function(params, x, N=1, IEI=50) {
  
#   response1 <- product1N(params[1:(N+3)], x=x, N=N, IEI=IEI)
#   response2 <- product1N(params[(N+4):(2*N+6)], x=x, N=N, IEI=IEI)
#   response <- response1 + response2
  
#   return(response) 

# }

product1N <- function(params, x, N=1, IEI=50) {
  # Extract parameters
  amplitudes <- params[1:N]  
  shared_params <- params[(N + 1):length(params)]  
  
  response <- numeric(length(x))  
  
  # Loop over the number of events
  for (i in seq_len(N)) {
    shifted_x <- x - (i - 1) * IEI 
    current_params <- c(amplitudes[i], shared_params)  
    current_response <- product1(current_params, shifted_x)  
    response <- response + current_response  
  }
  
  return(response)
}

product2N <- function(params, x, N=1, IEI=50) {
  # parameters for the first product1
  amplitudes1 <- params[1:N]                
  tau1_1 <- params[N + 1]                   
  tau2_1 <- params[N + 2]                   
  delay1 <- params[N + 3]                   
  
  # parameters for the second product1
  amplitudes2 <- params[(N + 4):(2 * N + 3)]  
  tau1_2 <- params[2 * N + 4]                
  tau2_2 <- params[2 * N + 5]                
  delay2 <- params[2 * N + 6]                 
  
  response <- numeric(length(x)) 
  
  # loop over the number of events
  for (i in seq_len(N)) {
    shifted_x <- x - (i - 1) * IEI 
    
    # first product1
    current_params1 <- c(amplitudes1[i], tau1_1, tau2_1, delay1)
    current_response1 <- product1(current_params1, shifted_x)
    
    # second product1
    current_params2 <- c(amplitudes2[i], tau1_2, tau2_2, delay2)
    current_response2 <- product1(current_params2, shifted_x)
    
    response <- response + current_response1 + current_response2
  }
  
  return(response) 
}

product3 <- function(params, x) {
  product1(params[1:4], x) + product1(params[5:8], x) + product1(params[9:12], x)
}

product3N <- function(params, x, N=1, IEI=50) {
  # parameters for the first product1
  amplitudes1 <- params[1:N]   
  tau1_1 <- params[N + 1]
  tau2_1 <- params[N + 2]
  delay1 <- params[N + 3]
  
  # parameters for the second product1
  amplitudes2 <- params[(N + 4):(2 * N + 3)] 
  tau1_2 <- params[2 * N + 4]                
  tau2_2 <- params[2 * N + 5]                
  delay2 <- params[2 * N + 6]                
  
  # parameters for the third product1
  amplitudes3 <- params[(2 * N + 7):(3 * N + 6)] 
  tau1_3 <- params[3 * N + 7]                    
  tau2_3 <- params[3 * N + 8]                    
  delay3 <- params[3 * N + 9]                    
  
  response <- numeric(length(x))  
  
  # loop over the number of events
  for (i in seq_len(N)) {
    shifted_x <- x - (i - 1) * IEI 
    
    # first product1
    current_params1 <- c(amplitudes1[i], tau1_1, tau2_1, delay1)
    current_response1 <- product1(current_params1, shifted_x)
    
    # second product1
    current_params2 <- c(amplitudes2[i], tau1_2, tau2_2, delay2)
    current_response2 <- product1(current_params2, shifted_x)
    
    # third product1
    current_params3 <- c(amplitudes3[i], tau1_3, tau2_3, delay3)
    current_response3 <- product1(current_params3, shifted_x)
    
    response <- response + current_response1 + current_response2 + current_response3
  }
  
  return(response)  # Return the total response
}

# sign_fun <- function(y, direction_method = c('smooth', 'regression', 'cumsum'), k=5) {
#   method <- match.arg(direction_method)
  
#   # Check if the length of y is sufficient for processing
#   if (length(y) < 10) {
#     stop("The length of y must be at least 10.")
#   }
  
#   n <- length(y)
#   peak_value <- NA
  
#   if (method == 'smooth') {
#     # Calculate the smoothed signal using a simple moving average
#     smoothed_signal <- rep(NA, n)
#     for (i in 1:(n - k + 1)) {
#       smoothed_signal[i + floor(k/2)] <- mean(y[i:(i + k - 1)])
#     }
#     # Find the peak value
#     peak_value <- max(abs(smoothed_signal), na.rm = TRUE)
#     peak_value <- smoothed_signal[which.max(abs(smoothed_signal))]
    
#   } else if (method == 'diff') {
#     # Calculate the differences of the signal
#     diff_signal <- diff(y)
#     # Identify both the maximum and minimum differences
#     max_diff <- max(diff_signal, na.rm = TRUE)
#     min_diff <- min(diff_signal, na.rm = TRUE)
    
#     # Determine the peak value considering the direction
#     peak_value <- ifelse(abs(max_diff) > abs(min_diff), max_diff, min_diff)
     
#   } else if (method == 'regression') {
#     # Fit a quadratic regression to the signal
#     x <- 1:n
#     model <- lm(y ~ I(x^2) + x)
    
#     # Use the fitted values to find the peak
#     fitted_values <- predict(model)
#     peak_value <- max(abs(fitted_values), na.rm = TRUE)
#     peak_value <- fitted_values[which.max(abs(fitted_values))]
    
#   } else if (method == 'cumsum') {
#     # Calculate the cumulative sum of the signal
#     cumsum_signal <- cumsum(y)
#     # Find the peak value of the cumulative sum
#     peak_value <- max(abs(cumsum_signal), na.rm = TRUE)
#     peak_value <- cumsum_signal[which.max(abs(cumsum_signal))]
#   }
  
#   # Determine if the peak is positive or negative
#   peak_direction <- ifelse(peak_value > 0, 1, -1)
  
#   # Return the sign of the peak direction
#   return(peak_direction[[1]])
# }


sign_fun <- function(y, direction_method = c('smooth', 'regression', 'cumsum'), k = 5) {
  
  method <- match.arg(direction_method)
  
  # Check if the length of y is sufficient for processing
  if (length(y) < 10) {
    stop("The length of y must be at least 10.")
  }
  
  n <- length(y)
  peak_value <- NA
  
  if (method == 'smooth') {
    # Calculate the smoothed signal using a simple moving average
    smoothed_signal <- rep(NA, n)
    for (i in 1:(n - k + 1)) {
      smoothed_signal[i + floor(k/2)] <- mean(y[i:(i + k - 1)])
    }
    # Remove NA values before finding the indices of the maximum and minimum values
    valid_indices <- which(!is.na(smoothed_signal))
    max_idx <- valid_indices[which.max(smoothed_signal[valid_indices])]
    min_idx <- valid_indices[which.min(smoothed_signal[valid_indices])]
    
    # Determine the peak value based on the magnitude comparison
    if (abs(smoothed_signal[max_idx]) > abs(smoothed_signal[min_idx])) {
      peak_value <- smoothed_signal[max_idx]
    } else {
      peak_value <- smoothed_signal[min_idx]
    }
    
  } else if (method == 'diff') {
    # Calculate the differences of the signal
    diff_signal <- diff(y)
    # Identify both the maximum and minimum differences
    max_diff <- max(diff_signal, na.rm = TRUE)
    min_diff <- min(diff_signal, na.rm = TRUE)
    
    # Determine the peak value considering the direction
    peak_value <- ifelse(abs(max_diff) > abs(min_diff), max_diff, min_diff)
     
  } else if (method == 'regression') {
    # Fit a quadratic regression to the signal
    x <- 1:n
    model <- lm(y ~ I(x^2) + x)
    
    # Use the fitted values to find the peak
    fitted_values <- predict(model)
    max_idx <- which.max(fitted_values)
    min_idx <- which.min(fitted_values)
    
    # Determine the peak value based on the magnitude comparison
    if (abs(fitted_values[max_idx]) > abs(fitted_values[min_idx])) {
      peak_value <- fitted_values[max_idx]
    } else {
      peak_value <- fitted_values[min_idx]
    }
    
  } else if (method == 'cumsum') {
    # Calculate the cumulative sum of the signal
    cumsum_signal <- cumsum(y)
    # Find the indices of the maximum and minimum values
    max_idx <- which.max(cumsum_signal)
    min_idx <- which.min(cumsum_signal)
    
    # Determine the peak value based on the magnitude comparison
    if (abs(cumsum_signal[max_idx]) > abs(cumsum_signal[min_idx])) {
      peak_value <- cumsum_signal[max_idx]
    } else {
      peak_value <- cumsum_signal[min_idx]
    }
  }
  
  # Determine if the peak is positive or negative
  peak_direction <- ifelse(peak_value > 0, 1, -1)
  
  # Return the sign of the peak direction
  return(peak_direction)
}

out.fun <- function(params, interval = c(0.1, 0.9), dp = 3, sign=1) {
  
  # Calculate derived metrics
  N <- length(params)-3

  A <- sign*params[1:N]
  tau.rise <- tau_rise(params[N+1], params[N+2])
  tau.decay <- params[N+2]
  tpeak <- product_tpeak(params[N+1], params[N+2])
  
  diffs <- product_rise_and_decay_percent(tau1 = params[N+1], tau2 = params[N+2], interval = interval, showplot = FALSE, verbose=FALSE)
  trise.percent <- diffs[1]
  tdecay.percent <- diffs[2]
  
  area <- sapply(1:N, function(ii) product_area(params[ii], params[N+1], params[N+2]))
  delay <- params[[N+3]]
  
  # Calculate rise and decay interval labels
  r_percent_start <- interval[1] * 100
  r_percent_end <- interval[2] * 100
  r_label <- paste0('r', r_percent_start, '_', r_percent_end)
  
  d_label <- paste0('d', interval[2] * 100, '_', interval[1] * 100)
  
  # Create data frame with calculated values
  output <- c(round(A, dp), round(tau.rise, dp), round(tau.decay, dp), round(tpeak, dp), round(trise.percent, dp), round(tdecay.percent, dp), round(delay, dp), round(area, dp))
  col.names <- c(paste0('A', 1:N), 'τrise', 'τdecay', 'tpeak', r_label, d_label, 'delay', paste0('area', 1:N))
  names(output) <- col.names
  df_output <- data.frame(t(output))
  
  return(df_output)
}
  
adjust_product_bounds <- function(bounds, func, upper=FALSE) {
  tau1_fun <- function(tau_rise, tau_decay) {
    return(tau_rise * tau_decay / (tau_decay - tau_rise))
  }
  adjust_tau <- function(value1, value2, upper=FALSE) {
    tau <- tau1_fun(value1, value2)
    if (is.nan(tau)) {
      tau <- if (upper) Inf else 0
    } else if (!upper && tau == Inf) {
      tau <- 0
    }
    return(tau)
  }
  
  if (!is.null(bounds)) {
    if (identical(func, product1)) {
      bounds[2]  <- adjust_tau(bounds[2], bounds[3], upper=upper)
    } else if (identical(func, product2)) {
      bounds[2]  <- adjust_tau(bounds[2], bounds[3], upper=upper)
      bounds[6]  <- adjust_tau(bounds[6], bounds[7], upper=upper)
    } else if (identical(func, product3)) {
      bounds[2]  <- adjust_tau(bounds[2], bounds[3], upper=upper)
      bounds[6]  <- adjust_tau(bounds[6], bounds[7], upper=upper)
      bounds[10] <- adjust_tau(bounds[10], bounds[11], upper=upper)
    }
  }
  return(bounds)
}

# adjust_product_bounds_N <- function(bounds, func, N, upper=FALSE) {
#   # Function to calculate tau1 based on tau_rise and tau_decay
#   tau1_fun <- function(tau_rise, tau_decay) {
#     return(tau_rise * tau_decay / (tau_decay - tau_rise))
#   }
  
#   # Function to adjust tau1 and handle upper or lower bounds
#   adjust_tau <- function(value1, value2, upper=FALSE) {
#     tau <- tau1_fun(value1, value2)
#     if (is.nan(tau)) {
#       tau <- if (upper) Inf else 0
#     } else if (!upper && tau == Inf) {
#       tau <- 0
#     }
#     return(tau)
#   }
  
#   # Adjust bounds based on the product function type
#   if (!is.null(bounds)) {
#     if (identical(func, product1N)) {
#       for (i in 1:length(bounds[1:N])) {
#         bounds[2 + (i - 1) * 3] <- adjust_tau(bounds[2 + (i - 1) * 3], bounds[3 + (i - 1) * 3], upper=upper)
#       }
#     } else if (identical(func, product2N)) {
#       for (i in 1:N) {
#         bounds[2 + (i - 1) * 3] <- adjust_tau(bounds[2 + (i - 1) * 3], bounds[3 + (i - 1) * 3], upper=upper)
#         bounds[(N + 5) + (i - 1) * 3] <- adjust_tau(bounds[(N + 5) + (i - 1) * 3], bounds[(N + 6) + (i - 1) * 3], upper=upper)
#       }
#     } else if (identical(func, product3N)) {
#       for (i in 1:N) {
#         bounds[2 + (i - 1) * 3] <- adjust_tau(bounds[2 + (i - 1) * 3], bounds[3 + (i - 1) * 3], upper=upper)
#         bounds[(N + 5) + (i - 1) * 3] <- adjust_tau(bounds[(N + 5) + (i - 1) * 3], bounds[(N + 6) + (i - 1) * 3], upper=upper)
#         bounds[(2 * N + 7) + (i - 1) * 3] <- adjust_tau(bounds[(2 * N + 7) + (i - 1) * 3], bounds[(2 * N + 8) + (i - 1) * 3], upper=upper)
#       }
#     }
#   }
  
#   return(bounds)
# }


adjust_product_bounds_N <- function(bounds, func, N, upper = FALSE) {
  # Function to calculate tau1 based on tau_rise and tau_decay
  tau1_fun <- function(tau_rise, tau_decay) {
    tau_rise * tau_decay / (tau_decay - tau_rise)
  }

  # Function to adjust tau1 and handle upper or lower bounds
  adjust_tau <- function(tau_rise, tau_decay, upper = FALSE) {
    tau <- tau1_fun(tau_rise, tau_decay)
    if (is.nan(tau) || tau < 0) {
      tau <- if (upper) Inf else 0
    }
    tau
  }

  # Adjust bounds based on the product function type
  if (!is.null(bounds)) {
    if (identical(func, product1N)) {
      # Adjust tau1
      bounds[N + 1] <- adjust_tau(bounds[N + 1], bounds[N + 2], upper = upper)
    } else if (identical(func, product2N)) {
      # Adjust tau1 for first component
      bounds[N + 1] <- adjust_tau(bounds[N + 1], bounds[N + 2], upper = upper)
      # Adjust tau1 for second component
      bounds[2 * N + 4] <- adjust_tau(bounds[2 * N + 4], bounds[2 * N + 5], upper = upper)
    } else if (identical(func, product3N)) {
      # Adjust tau1 for first component
      bounds[N + 1] <- adjust_tau(bounds[N + 1], bounds[N + 2], upper = upper)
      # Adjust tau1 for second component
      bounds[2 * N + 4] <- adjust_tau(bounds[2 * N + 4], bounds[2 * N + 5], upper = upper)
      # Adjust tau1 for third component
      bounds[3 * N + 7] <- adjust_tau(bounds[3 * N + 7], bounds[3 * N + 8], upper = upper)
    }
  }

  return(bounds)
}

# FIT <- function(response, dt=0.1, func=product2, method= c('LM', 'BF.LM', 'GN', 'port', 'robust', 'MLE'), 
#   stimulation_time=0, baseline=0, fast.decay.limit=NULL, latency.limit=NULL, lower=NULL, upper=NULL, filter=FALSE, fc=1000, interval=c(0.1, 0.9), 
#   MLEsettings=list(iter=1e4, metropolis.scale=1.5, fit.attempts=100, RWm=FALSE), 
#   MLE.method=c('L-BFGS-B', 'Nelder-Mead', 'BFGS','CG', 'SANN', 'Brent'),
#   response_sign_method = c('smooth', 'regression', 'cumsum'), 
#   dp=3, lwd=1.2, xlab='time (ms)', ylab='PSC (pA)', width=5, height=5, 
#   return.output=FALSE, show.output=TRUE, show.plot=TRUE){

#   dx <- dt 
#   method <- match.arg(method)
#   y <- response
#   if (all(is.na(y[(which(!is.na(y))[length(which(!is.na(y)))] + 1):length(y)]))) {
#     y <- y[!is.na(y)]
#   }

#   sign <- sign_fun(y, direction_method=response_sign_method) 
#   y <- sign * y

#   x <- seq(0, (length(y) - 1) * dx, by = dx)

#   if (filter){
#     ind = 20
#     fc = fc; fs = 1/dx*1000; bf <- butter(2, fc/(fs/2), type='low')
#     yfilter <- signal::filter(bf, y)
#   } else {
#     ind=1
#     yfilter=y
#   }
#   x.orig <- x

#   ind1 <- (stimulation_time - baseline)/dx
#   ind2 <- baseline/dx
  
#   yorig <- y[ind1:length(y)]
#   yfilter <- yfilter[ind1:length(yfilter)]
#   xorig <- seq(0, dx * (length(yorig) - 1), by = dx)

#   if (is_product_function(func)) {
#     yorig <- yorig - mean(yorig[1:ind2])
#     yfilter <- yfilter - mean(yfilter[1:ind2])

#     y2fit <- yfilter[ind2:length(yfilter)]
#     x2fit <- seq(0, dx * (length(y2fit) - 1), by = dx)

#     # check if both upper and fast.decay.limit are specified
#     if ( !is.null(upper) && (!is.null(fast.decay.limit) || !is.null(latency.limit)) ) {
#       warning("both 'upper' boundaries and at either 'fast decay limit' and/or 'latency limit' are specified:\ninputs 'fast.decay.limit' and/or 'latency.limit' will be ignored")
#       fast.decay.limit <- NULL
#       latency.limit <- NULL
#     }

#     # adjusts bounds for tau_rise to correct form involving tau1
#     upper <- adjust_product_bounds(bounds=upper, func=func, upper=TRUE)
#     lower <- adjust_product_bounds(bounds=lower, func=func)
      
#     if (is.null(upper) && (!is.null(fast.decay.limit) || !is.null(latency.limit))){

#       if (identical(func, product1)){
#         upper <- c(sign * Inf, Inf, Inf, Inf) 
#       } else if (identical(func, product2)) {
#         upper <- c(sign * Inf, Inf, Inf, Inf, sign * Inf, Inf, Inf, Inf)
#       } else if (identical(func, product3)) {
#         upper <- c(sign * Inf, Inf, Inf, Inf, sign * Inf, Inf, Inf, Inf, sign * Inf, Inf, Inf, Inf)
#       }

#       if (!is.null(fast.decay.limit)){
#         upper[3] <- fast.decay.limit[1]
#         if (identical(func, product3)){
#           upper[7] <-  if (length(fast.decay.limit)==1) fast.decay.limit[1] else fast.decay.limit[2]
#         }
#       }

#       if (!is.null(latency.limit)){
#         if (identical(func, product1)){
#           upper[4] <- latency.limit
#         } else if (identical(func, product2)){
#           upper[c(4,8)] <- latency.limit
#         } else if (identical(func, product3)){
#           upper[c(4,8,12)] <- latency.limit
#         }
#       }

#     }
  
#     if (identical(func, product1)){
#       param_names <- c('a', 'b', 'c', 'd')
#       N.params <- length(param_names)
#       form.2.fit <- y ~ (a / (((b / (b + c)) ^ (b / c)) * c / (b + c))) * (1 - exp(-(x - d) / b)) * exp(-(x - d) / c) * (x >= d)

#       if (!is.null(lower) && length(lower) != N.params) stop("number of lower boundaries should be equal to the number of parameters")
#       if (!is.null(upper) && length(upper) != N.params) stop("number of upper boundaries should be equal to the number of parameters")
#       if (!is.null(lower)) lower[1] <- sign * lower[1] 
#       if (!is.null(upper)) upper[1] <- sign * upper[1] 

#     } else if (identical(func, product2)){
#       param_names <- c('a', 'b', 'c', 'd', 'e', 'f', 'g', 'h')
#       N.params <- length(param_names)
#       form.2.fit <- y ~ (a / (((b / (b + c)) ^ (b / c)) * c / (b + c))) * (1 - exp(-(x - d) / b)) * exp(-(x - d) / c) * (x >= d) + 
#                         (e / (((f / (f + g)) ^ (f / g)) * g / (f + g))) * (1 - exp(-(x - h) / f)) * exp(-(x - h) / g) * (x >= h)
#       if (!is.null(lower) && length(lower) != N.params) stop("number of lower boundaries should be equal to the number of parameters")
#       if (!is.null(upper) && length(upper) != N.params) stop("number of upper boundaries should be equal to the number of parameters")
#       if (!is.null(lower)) { lower[1] <- sign * lower[1]; lower[5] <- sign * lower[5] }
#       if (!is.null(upper)) { upper[1] <- sign * upper[1]; upper[5] <- sign * upper[5] }
    
#     } else if (identical(func, product3)){
#       param_names <- c('a', 'b', 'c', 'd', 'e', 'f', 'g', 'h', 'i', 'j', 'k', 'l')
#       N.params <- length(param_names)
#       form.2.fit <- y ~ (a / (((b / (b + c)) ^ (b / c)) * c / (b + c))) * (1 - exp(-(x - d) / b)) * exp(-(x - d) / c) * (x >= d) + 
#                         (e / (((f / (f + g)) ^ (f / g)) * g / (f + g))) * (1 - exp(-(x - h) / f)) * exp(-(x - h) / g) * (x >= h) +
#                         (i / (((j / (j + k)) ^ (j / k)) * g / (j + k))) * (1 - exp(-(x - l) / j)) * exp(-(x - l) / k) * (x >= l)
#       if (!is.null(lower) && length(lower) != N.params) stop("number of lower boundaries should be equal to the number of parameters")
#       if (!is.null(upper) && length(upper) != N.params) stop("number of upper boundaries should be equal to the number of parameters")
#       if (!is.null(lower)) { lower[1] <- sign * lower[1]; lower[5] <- sign * lower[5]; lower[9] <- sign * lower[9] }
#       if (!is.null(upper)) { upper[1] <- sign * upper[1]; upper[5] <- sign * upper[5]; upper[9] <- sign * upper[9] }
#     }
#   }else{
#     fun.2.fit <- func
#     y2fit <- yfilter[ind2:length(yfilter)]
#     x2fit <- seq(0, dx * (length(y2fit) - 1), by = dx)
#     param_names <- generate_param_names(fun.2.fit)
#     N.params <- length(param_names)
#     form.2.fit <- generate_formula(fun.2.fit, param_names, y_name = 'y', x_name = 'x')
#   }

#  if (method == 'MLE'){
    
#     MLE.method <- match.arg(MLE.method)
#     iter <- MLEsettings$iter
#     metropolis.scale <- MLEsettings$metropolis.scale
#     fit.attempts <- MLEsettings$fit.attempts
#     RWm <- MLEsettings$RWm  
#     sd.est <- sqrt(mean(yfilter[1:ind2]^2))

#     output <- FIT.MLE(x=x2fit, y=y2fit, func=func, N.params=N.params, sigma=sd.est, iter=iter, metropolis.scale=metropolis.scale, 
#       logpost=log.lik.post, params.start=NULL, method=MLE.method, lower=lower, upper=upper, MLE.fun.attempts=100, fit.attempts=fit.attempts, RWm=RWm) 

#   } else if (method == 'LM'){
#      output <- FIT.LM(x=x2fit, y=y2fit, func=func, N.params=N.params, lower=lower, upper=upper, fit.convergence.attempts=10, fit.attempts=10)

#   } else if (method == 'BF.LM'){
#     output <- FIT.bf.LM(x=x2fit, y=y2fit, func=func, N.params=N.params, lower=lower, upper=upper, fit.convergence.attempts=10, fit.attempts=10)

#   } else if (method == 'GN'){
#     output <- FIT.NLS(x=x2fit, y=y2fit, func=func, form=form.2.fit, N.params=N.params, lower=lower, upper=upper)

#   } else if (method == 'port'){
#     output <- FIT.NLS(x=x2fit, y=y2fit, func=func, form=form.2.fit, N.params=N.params, lower=lower, upper=upper, algorithm='port')

#   } else if (method == 'robust'){
#     algorithm <- if (!is.null(lower) || !is.null(upper)) 'port' else 'default'
#     # currently only 'M' working
#     algorithm <-  'port'
#     output <- FIT.NLSrobust(x=x2fit, y=y2fit, func=func, form=form.2.fit, N.params=N.params, lower=lower, upper=upper, algorithm=algorithm, method='M')
#   }
    

#   if (is_product_function(func)) {
#     if (identical(func, product1)){
#       df_output <- out.fun(params=output$fits[1:4], interval=interval, dp=dp, sign=sign)
#       output$fits[1] <- sign * output$fits[1]
#     } else if (identical(func, product2)){  
#       df_output1 <- out.fun(params=output$fits[1:4], interval=interval, dp=dp, sign=sign)
#       df_output2 <- out.fun(params=output$fits[5:8], interval=interval, dp=dp, sign=sign)
#       df_output <- if (df_output1[3] < df_output2[3]) rbind('fast' = df_output1, 'slow' = df_output2) else rbind('fast' = df_output2, 'slow' = df_output1)
#       output$fits[1] <- sign * output$fits[1]; output$fits[5] <- sign * output$fits[5]
#       if (output$fits[3] > output$fits[7]){
#         output$fits <- c(output$fits[5:8], output$fits[1:4])
#         output$fits.se <- c(output$fits.se[5:8], output$fits.se[1:4])
#       }
#     } else if (identical(func, product3)){  
#       # Compute the output data frames
#       df_output1 <- out.fun(params=output$fits[1:4],  interval=interval, dp=dp, sign=sign)
#       df_output2 <- out.fun(params=output$fits[5:8],  interval=interval, dp=dp, sign=sign)
#       df_output3 <- out.fun(params=output$fits[9:12], interval=interval, dp=dp, sign=sign)

#       # Create a list of the outputs
#       output_list <- list(df_output1, df_output2, df_output3)

#       # Extract the third elements from each output and determine the order
#       order_indices <- order(sapply(output_list, function(x) x[[3]]))

#       # Reorder the output list based on the third element
#       output_ordered <- output_list[order_indices]

#       # Combine the outputs in the correct order
#       df_output <- rbind('fast' = output_ordered[[1]], 'medium' = output_ordered[[2]], 'slow' = output_ordered[[3]])

#       # Adjust the order of the fits and fits.se based on the sorted order
#       output$fits <- c(output$fits[(4*order_indices[1]-3):(4*order_indices[1])],
#                        output$fits[(4*order_indices[2]-3):(4*order_indices[2])],
#                        output$fits[(4*order_indices[3]-3):(4*order_indices[3])])

#       output$fits.se <- c(output$fits.se[(4*order_indices[1]-3):(4*order_indices[1])],
#                           output$fits.se[(4*order_indices[2]-3):(4*order_indices[2])],
#                           output$fits.se[(4*order_indices[3]-3):(4*order_indices[3])])

#       # adjust signs
#       output$fits[1] <- sign * output$fits[1]
#       output$fits[5] <- sign * output$fits[5]
#       output$fits[9] <- sign * output$fits[9]
#     }

#     traces=data.frame(x=xorig, y= sign * yorig, yfilter= sign * yfilter)
#     fits <- output$fits
#     if (identical(func, product1)){
#       fits[4] <- fits[4] + baseline
#     } else if (identical(func, product2)){
#       fits[4] <- fits[4] + baseline; fits[8] <- fits[8] + baseline
#     } else if (identical(func, product3)){
#       fits[4] <- fits[4] + baseline; fits[8] <- fits[8] + baseline; fits[12] <- fits[12] + baseline
#     }    
#     traces$yfit <- func(fits,traces$x+dx)  
#     if (identical(func, product2)){
#       traces$yfit1 <- product1(fits[1:4],traces$x+dx) 
#       traces$yfit2 <- product1(fits[5:8],traces$x+dx) 
#     } 
#     if (identical(func, product3)){
#       traces$yfit1 <- product1(fits[1:4],traces$x+dx) 
#       traces$yfit2 <- product1(fits[5:8],traces$x+dx) 
#       traces$yfit3 <- product1(fits[9:12],traces$x+dx) 
#     } 

#   }else{
#     traces=data.frame(x=xorig, y=sign * yorig, yfilter=sign * yfilter)
#     fits <- output$fits
#     traces$yfit <- func(fits,x2fit)  
#     df_output <- data.frame(fits=fits)
#     df_output <- t(df_output)

#     df_output <- as.data.frame(df_output)
#     colnames(df_output) <- param_names
#   }

#   if (show.plot) fit_plot(traces=traces, func=func, xlab=xlab, ylab=ylab, lwd=lwd, filter=filter, width=width, height=height)

#   if (show.output){
#     print(df_output)
#   }

#   if (return.output){
#     return(list(output = df_output, fits = output$fits, fits.se = output$fits.se, gof = output$gof, AIC = output$AIC, BIC = output$BIC, model.message = output$model.message, sign=sign, traces=traces))
#   }
# }


fit_plot <- function(traces, func=product2, xlab='time (ms)', ylab='PSC amplitude (pA)', xlim=NULL, ylim=NULL, main='', bl=NULL, lwd=1.2, filter=FALSE, width=5, height=5, bg='transparent', filename='trace.svg', save=FALSE) {
  
  if (save) {
    # Open SVG device
    svg(file=filename, width=width, height=height, bg=bg)
  } else {
    dev.new(width=width, height=height, noRStudioGD=TRUE)
  }

  plot(traces$x, traces$y, col='gray', xlab=xlab, ylab=ylab, xlim=xlim, ylim=ylim, type='l', bty='l', las=1, lwd=lwd, main=main)
  
  if (filter) {
    lines(traces$x, traces$yfilter, col='black', type='l', lwd=lwd)
  }
  
  lines(traces$x, traces$yfit, col='indianred', lty=3, lwd=2 * lwd)
  
  if (identical(func, product2) || identical(func, product2N)) {
    lines(traces$x, traces$yfit1, col='#4C78BC', lty=3, lwd=2 * lwd)
    lines(traces$x, traces$yfit2, col='#CA92C1', lty=3, lwd=2 * lwd)
  }

  if (identical(func, product3) || identical(func, product3N)) {
    lines(traces$x, traces$yfit1, col='#F28E2B', lty=3, lwd=2 * lwd)
    lines(traces$x, traces$yfit2, col='#4C78BC', lty=3, lwd=2 * lwd)
    lines(traces$x, traces$yfit3, col='#CA92C1', lty=3, lwd=2 * lwd)
  }

  if (!is.null(bl)) abline(v=bl, col='black', lwd=lwd, lty=3)

  if (save) {
    dev.off()
  }
}

nFIT <- function(response, n=30, N=1, IEI=50, dt=0.1, func=product2N, method= c('BF.LM', 'LM', 'GN', 'port', 'robust', 'MLE'), weight_method=c('none', '~y_sqrt', '~y'),
  stimulation_time=0, baseline=0, fast.decay.limit=NULL, fast.constraint=FALSE, fast.constraint.method=c('rise', 'peak'), first.delay.constraint=FALSE,
  latency.limit=NULL, lower=NULL, upper=NULL, filter=FALSE, fc=1000, interval=c(0.1, 0.9), MLEsettings=list(iter=1e4, metropolis.scale=1.5, fit.attempts=100, RWm=FALSE),  
  MLE.method=c('L-BFGS-B', 'Nelder-Mead', 'BFGS','CG', 'SANN', 'Brent'), response_sign_method = c('smooth', 'regression', 'cumsum'), half_width_fit_limit=500, dp=3, lwd=1.2, 
  xlab='time (ms)', ylab='PSC (pA)', width=5, height=5, return.output=FALSE, show.output=TRUE, show.plot=TRUE, seed=42) {
  
  set.seed(seed)
  
  output <- NULL
  gof <- Inf
  fit_results <- list()
  df_output <- NULL
  sign <- NULL
  traces <- NULL

  fast.constraint.method <- match.arg(fast.constraint.method)

  for (i in 1:n) {
    fit_result <- tryCatch({
      
      FITN(response=response, dt=dt, func=func, N=N, IEI=IEI, method=method, weight_method=weight_method, stimulation_time=stimulation_time, baseline=baseline, fast.decay.limit=fast.decay.limit, 
                      latency.limit=latency.limit, lower=lower, upper=upper, filter=filter, fc=fc, interval=interval, MLEsettings=MLEsettings, 
                      MLE.method=MLE.method, response_sign_method=response_sign_method, dp=10, lwd=lwd, xlab=xlab, ylab=ylab, width=width, height=height, 
                      return.output=TRUE, show.output=FALSE, show.plot=FALSE)
    }, error = function(e) {
      NULL
    })
    
    if (!is.null(fit_result)) {
      fit_results[[i]] <- fit_result[!(names(fit_result) %in% c('sign', 'traces'))]
      
      if (first.delay.constraint){
        d1 <- fit_result$output[1, 'delay'] 
        d2 <- fit_result$output[2, 'delay'] 
      }
      
      if (identical(func, product2N)){

        if (fast.constraint){

          r_percent_start <- interval[1] * 100
          r_percent_end   <- interval[2] * 100
          r_label <- paste0('r', r_percent_start, '_', r_percent_end)

          if (fast.constraint.method == 'rise'){
            fast_ <- fit_result$output[1, r_label]
            slow_ <- fit_result$output[2, r_label]
          }else if (fast.constraint.method == 'peak'){
            fast_ <- fit_result$output[1, 'tpeak'] 
            slow_ <- fit_result$output[2, 'tpeak'] 
          }
          
          if (first.delay.constraint){
            if (!(fast_ < slow_ && d1 < d2)) {
              next  # Skip updating if fast r10_90 is greater than slow r10_90
            }
          }else{
            if (fast_ > slow_) {
              next  # Skip updating if fast r10_90 is greater than slow r10_90
            }
          }

        }

      }else if (identical(func, product3N)){
      
        if (fast.constraint){
          r_percent_start <- interval[1] * 100
          r_percent_end   <- interval[2] * 100
          r_label <- paste0('r', r_percent_start, '_', r_percent_end)

          if (fast.constraint.method == 'rise'){
            fast_  <- fit_result$output[1, r_label]
            medium_ <- fit_result$output[2, r_label]
            slow_   <- fit_result$output[3, r_label]
          }else if (fast.constraint.method == 'peak'){
            fast_   <- fit_result$output[1, 'tpeak'] 
            medium_ <- fit_result$output[2, 'tpeak'] 
            slow_   <- fit_result$output[3, 'tpeak']
          }
          
          if (first.delay.constraint){
            d3 <- fit_result$output[3, 'delay'] 
            if (!(fast_ < medium_ && medium_ < slow_ && d1 < d2 && d1 < d3)) {
              next  # Skip updating if fast r10_90 is greater than slow r10_90
            }
          }else{
            if (!(fast_ < medium_ && medium_ < slow_)) {
              next  # Skip updating if fast r10_90 is greater than slow r10_90
            }
          }
         
        }   
      }
  
      if (fit_result$gof < gof) {
        df_output <- fit_result$output
        output <- fit_result
        gof <- fit_result$gof
        sign <- fit_result$sign
        traces <- fit_result$traces      
      }
    }
  }
  
  if (is.null(output)) {
    stop("All fit attempts failed")
  }
  
  # apply half_width function and extract its value
  df_output$half_width <- mapply(
    function(A, tau1, tau2) half_width(A, tau1, tau2, limit = half_width_fit_limit)[['half_width']],
    A = sign * df_output$A1, 
    tau1 = df_output$τrise,
    tau2 = df_output$τdecay
  )

  # reorder
  # df_output <- df_output[, c(names(df_output)[1:6], 'half_width', 'delay', 'area1')]
  # df_output <- round(df_output, digits=dp)

  cols <- names(df_output)
  area_idx <- which(grepl("^area", cols))[1]
  new_order <- append(cols[-which(cols == "half_width")], "half_width", after = area_idx - 1)
  df_output <- df_output[, new_order]
  df_output <- round(df_output, digits=dp)

  if (show.plot) fit_plot(traces=traces, func=func, xlab=xlab, ylab=ylab, lwd=lwd, filter=filter, width=width, height=height)
  
  if (show.output) print(df_output)
  
  if (return.output) {
    return(list(output = df_output, fits = output$fits, fits.se = output$fits.se, gof = output$gof, AIC = output$AIC, BIC = output$BIC, model.message = output$model.message, sign=sign, traces=traces, fit_results=fit_results))
  }
}

first.peak <- function(y, threshold=0.1) {
  
  sign <- sign_fun(y)
  y <- sign * y
  
  # Define a simple threshold to ignore tiny fluctuations (for example, 10% of max)
  threshold_value <- threshold * max(abs(y), na.rm=TRUE)

  # Initialize peak_index to NA, which will be returned if no peak is found
  peak_index <- NA
  
  # Loop through the data and find the first peak (positive value followed by a negative one)
  for (i in 2:(length(y) - 1)) {
    if (!is.na(y[i]) && y[i] > y[i - 1] && y[i] > y[i + 1] && y[i] > threshold_value) {
      peak_index <- i
      break  # Exit the loop once the first peak is found
    }
  }

  # Return the index of the first peak or NA if none is found
  return(peak_index)
}



last.peak <- function(y, threshold=0.9) {
  
  sign <- sign_fun(y)
  y <- sign * y

  # Define a simple threshold to ignore tiny fluctuations (for example, 10% of max)
  threshold_value <- threshold * max(abs(y), na.rm=TRUE)

  # Reverse the data to find the "first peak" in the reversed vector, which corresponds to the last peak
  y_rev <- rev(y)

  peak_index_rev <- NA
  for (i in 2:(length(y_rev) - 1)) {
    if (!is.na(y_rev[i]) && y_rev[i] > y_rev[i - 1] && y_rev[i] > y_rev[i + 1] && y_rev[i] > threshold_value) {
      peak_index_rev <- i
      break
    }
  }

  # Convert the index from the reversed vector back to the original vector's index
  if (!is.na(peak_index_rev)) {
    peak_index <- length(y) - peak_index_rev + 1
  } else {
    peak_index <- NA
  }

  return(peak_index)
}

ests_fun <- function(x, y, N=1, showplot=FALSE) {
  n <- length(y)
  
  dx=x[2]-x[1]
  k <- round(ceiling(1/dx))

  ysmooth <- rep(NA, n)
  for (i in 1:(n - k + 1)) {
    ysmooth[i + floor(k/2)] <- mean(y[i:(i + k - 1)])
  }

  idx.peak <- which.max(abs(ysmooth))
  t.peak <- x[idx.peak]
  A <- ysmooth[idx.peak]

  # Find the index for peak value; if N=1 then indices are the same else indices for first and last peak
  idx.1st <- ifelse(N == 1 || is.na(first.peak(ysmooth)), which.max(abs(ysmooth)), first.peak(ysmooth))
  idx.Nth <- ifelse(N == 1 || is.na(last.peak(ysmooth)), idx.1st, last.peak(ysmooth))

  t.peak1 <- x[idx.1st]
  A1 <- ysmooth[idx.1st]
  A.Nth <- ysmooth[idx.Nth]
  
  # Calculate 20% and 80% of the 1st and Nth peak amplitude
  A1.20 <- 0.20 * A1; A.Nth.20 <- 0.20 * A.Nth 
  A1.80 <- 0.80 * A1; A.Nth.80 <- 0.80 * A.Nth 
  
  # Find all indices where y crosses 20% and 80% amplitude levels during rising phase
  rise_idx.20 <- which(diff(sign(ysmooth[1:idx.1st] - A1.20)) != 0)
  rise_idx.80 <-  which(diff(sign(ysmooth[1:idx.1st] - A1.80)) != 0)
  
  # Find all indices where y crosses 20% and 80% amplitude levels during decaying phase
  decay_idx.80 <- which(diff(sign(ysmooth[idx.Nth:length(y)] - A.Nth.80)) != 0) + idx.Nth
  decay_idx.20 <- which(diff(sign(ysmooth[idx.Nth:length(y)] - A.Nth.20)) != 0) + idx.Nth
  
  # Handle the case where decay_idx.20 does not exist
  if (length(decay_idx.20) == 0) {
    # Perform linear extrapolation
    if (length(decay_idx.80) > 1) {
      # Get the last two points near the 80% level
      x1 <- x[decay_idx.80[length(decay_idx.80) - 1]]
      x2 <- x[decay_idx.80[length(decay_idx.80)]]
      y1 <- y[decay_idx.80[length(decay_idx.80) - 1]]
      y2 <- y[decay_idx.80[length(decay_idx.80)]]
      
      # Calculate the slope of the line
      slope <- (y2 - y1) / (x2 - x1)
      
      # Extrapolate to the 20% level
      t.decay.20 <- x2 + (A1.20 - y2) / slope
    } else {
      # If we don't have enough points to extrapolate, use the last time point
      t.decay.20 <- x[length(x)]
    }
  } else {
    t.decay.20 <- x[decay_idx.20]
  }
  
  # Calculate times for these points
  t.rise.20 <- x[rise_idx.20]
  t.rise.80 <- x[rise_idx.80]
  t.decay.80 <- x[decay_idx.80]
  
  # Average the times if there are multiple points
  avg_t.rise.20 <- median(t.rise.20)
  avg_t.rise.80 <- median(t.rise.80)
  avg_t.decay.80 <- median(t.decay.80)
  avg_t.decay.20 <- median(t.decay.20)
  
  # Calculate 20-80% rise time and 80-20% decay time
  rise.time <- avg_t.rise.80 - avg_t.rise.20

  if (rise.time==0) rise.time <- 0.5 * t.peak
  decay.time <- avg_t.decay.20 - avg_t.decay.80

  # Plot the original and smoothed signals
  if (showplot) {
    plot(x, y, col='indianred', xlab='x', type='l', bty='l', las=1, main='')
    lines(x, ysmooth, col='lightgray')
    # abline(h=A1.20, col='black', lty=3)
    # abline(h=A1.80, col='black', lty=3)
    abline(h= A.Nth.20, col='black', lty=3)
    abline(h= A.Nth.80, col='black', lty=3)
  }

  return(c(A, rise.time, decay.time))
}


# Use generate_lognormal_samples to produce random starting points balsed on estimated starting values
generate_lognormal_samples <- function(means, cv=0.4, n=1) {
  # Calculate parameters for the lognormal distribution
  calculate_lognormal_params <- function(mean, cv) {
    mu <- log(mean) - 0.5 * log(1 + cv^2)
    sigma <- sqrt(log(1 + cv^2))
    list(mu = mu, sigma = sigma)
  }
  
  # Generate samples for each mean
  samples <- sapply(means, function(mean) {
    params <- calculate_lognormal_params(mean, cv)
    rlnorm(n, meanlog = params$mu, sdlog = params$sigma)
  })
  return(samples)
}

adjust_upper1 <- function(upper, ests, factor) {
  changed <- FALSE
  for (i in 1:length(upper)) {
    if (upper[i] < ests[i]) {
      upper[i] <- factor * ceiling(5 * ests[i] / factor)
      changed <- TRUE
    }
  }
  return(list(upper = upper, changed = changed))
}

bounds_check <- function(ests, model='product2', lower=NULL, upper=NULL, fast.decay.limit=NULL, latency.limit=NULL){
  if (model == 'product') {
    result <- adjust_upper1(upper, ests, 10)
    N.params <- 4
    upper <- result$upper
  } else if (model == 'product2') {
    result <- adjust_upper1(upper, ests, 10)
    N.params <- 8
    upper <- result$upper
    if (upper[1] < upper[5]) upper[1] <- upper[5] else if (upper[5] < upper[1]) upper[5] <- upper[1]
  }
  if (result$changed) {
    message("upper bounds have been changed to: [", paste(upper, collapse = ", "), "]")
  }

  if (is.null(lower)) lower <- rep(0, N.params) 
  return(list(lower=lower, upper=upper))
  }

adjust_upper <- function(upper, ests, indices, factor) {
  changed <- FALSE
  for (i in indices) {
    est_index <- ifelse(i %in% c(5, 6, 7), (i - 4), (i - 1) %% length(ests) + 1)
    if (upper[i] < ests[est_index]) {
      upper[i] <- factor * round(5 * ests[est_index] / factor)
      changed <- TRUE
    }
  }
  return(list(upper = upper, changed = changed))
}

start_fun <- function(x, y, cv=0.4, showplot=FALSE, model='product2', lower=NULL, upper=NULL, fast.decay.limit=NULL, latency.limit=NULL){

  ests <- ests_fun(x, y)

  if (!is.null(upper)) {
    if (model == 'product') {
      result <- adjust_upper(upper, ests, c(1, 2, 3), 10)
    } else if (model == 'product2') {
      result <- adjust_upper(upper, ests, c(1, 2, 3, 5, 6, 7), 10)
    }
    
    upper <- result$upper
    if (result$changed) {
      message("upper bounds have been changed to: [", paste(upper, collapse = ", "), "]")
    }

    if (model == 'product'){
      success<- FALSE
      for (i in 1:1e4) {
        st <- generate_lognormal_samples(means=ests, cv=cv, n=1)
        if (all(st[1:3] < upper[1:3])){
          success <- TRUE
          break
        }
      }
      if (!success) {
        stop("failed to generate initial starting values (1e4 attempts)")
      }

      st <- c(st, runif(1)*upper[4])

    }  else if (model == 'product2'){
        success <- FALSE
        ests1 <- ests; ests1 <- ests1/2
        ests2 <- ests; ests2[1] <- ests2[1]/2

        for (i in 1:1e4) {
        st <- generate_lognormal_samples(means = c(ests1, ests2), cv = cv, n = 1)
        if ((all(st[1:3] < upper[1:3])) && all(st[4:6] < upper[5:7])) {
          success <- TRUE
          break
        }
      }
      if (!success) {
        stop("failed to generate initial starting values (1e4 attempts)")
      }
      st <- c(st[1:3], runif(1)*upper[4], st[4:6], runif(1)*upper[4])
    } 

  } else {
  
    if (model == 'product'){
        success <- FALSE
        if (!is.null(fast.decay.limit)) {
          for (i in 1:1e4) {
            st <- generate_lognormal_samples(means=ests, cv=cv, n=1)
            if (st[3] < fast.decay.limit) {
              success <- TRUE
              break
            }
          }
          if (!success) {
            stop("failed to generate initial starting values where τdecay is less than the fast decay limit (1e4 attempts)")
          }
        }else{
          st <- generate_lognormal_samples(means=ests, cv=cv, n=1)
        }

        st <- if (!is.null(latency.limit)) c(st, runif(1)*latency.limit) else c(st, runif(1)*5)

        if (is.null(lower)) lower= c(0, 0, 0, 0) # rep(0, N.params)
        if (is.null(upper)) upper = c(ests * 10, Inf)
        if (!is.null(fast.decay.limit)) upper[3] <- fast.decay.limit
        if (!is.null(latency.limit)) upper[4] <- latency.limit

      } else if (model == 'product2'){
        success <- FALSE
        # 
        ests1 <- ests; ests1 <- ests1/2
        ests2 <- ests; ests2[1] <- ests2[1]/2
        
        if (!is.null(fast.decay.limit)) {
          for (i in 1:1e4) {
            st <- generate_lognormal_samples(means = c(ests1, ests2), cv = cv, n = 1)
            if (st[3] < fast.decay.limit) {
              success <- TRUE
              break
            }
          }
          if (!success) {
            stop("failed to generate initial starting values where τdecay is less than the fast decay limit (1e4 attempts)")
          }
        } else {
          st <- generate_lognormal_samples(means = c(ests1, ests2), cv = cv, n = 1)
        }

        st <- if (!is.null(latency.limit)) c(st[1:3], runif(1)*latency.limit, st[4:6], runif(1)*latency.limit) else c(st[1:3], runif(1)*5, st[4:6], runif(1)*5)

        if (is.null(lower)) lower <- c(0, 0, 0, 0, 0, 0, 0, 0) # rep(0, N.params)
        if (is.null(upper)) upper <- rep(c(ests * 10, Inf),2)
        if (!is.null(fast.decay.limit)) upper[3] <- fast.decay.limit
        if (!is.null(latency.limit)) upper[c(4,8)] <- latency.limit

      }
    }
    return(list(start=st, lower=lower, upper=upper))
  }

# analysis functions
load_data <- function(wd, name) {
    
  # Create the file path
  file_path <- file.path(wd, paste0(name, '.', 'csv'))
  
  # Load the csv data
  data <- read.csv(file_path, stringsAsFactors = FALSE)
  
  return(data)
}

moving_avg <- function(y, n = 5) {
  sign <- sign_fun(y)
  y <- y * sign
  y_length <- length(y)
  result <- rep(NA, y_length)
  
  for (i in 1:y_length) {
    # Determine the start and end indices for the window
    start_idx <- max(1, i - floor(n / 2))
    end_idx <- min(y_length, i + floor(n / 2))
    
    # Calculate the mean for the current window
    result[i] <- mean(y[start_idx:end_idx], na.rm = TRUE)
  }
  
  return(max(result) * sign)
}

# load_data2 <- function(wd, name) {
#   # Create the file path
#   file_path <- file.path(wd, paste0(name, '.', 'xlsx'))
  
#   # Load the Excel file
#   workbook <- openxlsx::loadWorkbook(file_path)
  
#   # Get the sheet names
#   sheet_names <- openxlsx::getSheetNames(file_path)
  
#   # Initialize an empty list to store each sheet's data
#   data_list <- list()
  
#   # Loop through each sheet and read the data into the list
#   for (sheet in sheet_names) {
#     data_list[[sheet]] <- openxlsx::read.xlsx(file_path, sheet = sheet)
#   }
  
#   return(data_list)
# }


load_data2 <- function(wd, name, header = TRUE) {
  # Create the file path
  file_path <- file.path(wd, paste0(name, '.', 'xlsx'))
  
  # Load the Excel file
  workbook <- openxlsx::loadWorkbook(file_path)
  
  # Get the sheet names
  sheet_names <- openxlsx::getSheetNames(file_path)
  
  # Initialize an empty list to store each sheet's data
  data_list <- list()
  
  # Loop through each sheet and read the data into the list
  for (sheet in sheet_names) {
    data_list[[sheet]] <- openxlsx::read.xlsx(file_path, sheet = sheet, colNames = header)
  }
  
  return(data_list)
}

# save list to excel spreadsheet
list2excel <- function(data_list, file_name, wd = getwd(), center_align = TRUE) {
  # Load the openxlsx library
  library(openxlsx)
  
  # Create a new workbook
  workbook <- createWorkbook()
  
  # Define styles: bold and centered for headers, and centered for data
  header_style <- createStyle(textDecoration = "bold", halign = "center", valign = "center")
  center_style <- createStyle(halign = "center", valign = "center")
  
  # Loop over each element in the list
  for (i in seq_along(data_list)) {
    # Use the name of the list element as the sheet name
    sheet_name <- names(data_list)[i]
    
    # Default to "Sheet1", "Sheet2", etc., if the name is missing
    if (is.null(sheet_name) || sheet_name == "") {
      sheet_name <- paste0("Sheet", i)
    }
    
    # Add a new sheet with the specified name to the workbook
    addWorksheet(workbook, sheet_name)
    
    # Write the data to the sheet
    writeData(workbook, sheet_name, data_list[[i]])
    
    # Apply header style to make headers bold and centered
    addStyle(workbook, sheet_name, style = header_style, rows = 1, cols = 1:ncol(data_list[[i]]), gridExpand = TRUE)
    
    # Apply center alignment to all cells if center_align is TRUE
    if (center_align) {
      addStyle(workbook, sheet_name, style = center_style, rows = 1:(nrow(data_list[[i]]) + 1), 
               cols = 1:ncol(data_list[[i]]), gridExpand = TRUE)
    }
  }
  
  # Create the full file path
  file_path <- file.path(wd, file_name)
  
  # Save the workbook
  saveWorkbook(workbook, file_path, overwrite = TRUE)
}


peak.fun <- function(y, dt, stimulation_time, baseline, smooth=5){
  
  idx1 <- (stimulation_time - baseline) / dt
  idx2 <- baseline / dt

  y1 <- y[idx1:length(y)]

  window_idx <- smooth/dt
  y1 <- y1 - mean(y1[0:idx2])
  return(moving_avg(y1, n = smooth)) 
}

raw_plot <- function(response, dt=0.1, stimulation_time=0, baseline=0, smooth=5, y_abline=0.1, height=5, width=5){

    y <- response

    if (all(is.na(y[(which(!is.na(y))[length(which(!is.na(y)))] + 1):length(y)]))) {
      y <- y[!is.na(y)]
    }
  
    peak <- peak.fun(y=y, dt=dt, stimulation_time=stimulation_time, baseline=baseline, smooth=smooth)

    ind1 <- (stimulation_time - baseline)/dt
    ind2 <- stimulation_time/dt
    y2plot <- y - mean(y[ind1:ind2])

    dev.new(width=width, height=height, noRStudioGD=TRUE)
    Y <- y2plot[ind1:length(y2plot)]
    X <- seq(0, dt * (length(Y) - 1), by = dt)
    plot(X, Y, col='indianred', xlab='time (ms)', type='l', bty='l', las=1, main='')
    abline(h = 0, col = 'black', lwd = 1, lty=1)
    abline(h = peak * y_abline, col = 'black', lwd = 1, lty=3)

    # Add a label to the abline
    text(x=max(X[ind1:length(X)]) * 0.95, y=peak * (y_abline - 0.025), labels=y_abline, pos=4)

}

# determine_tmax <- function(y, dt=0.1, stimulation_time=0, baseline=0, smooth=5, tmax=NULL, y_abline=0.1, height=5, width=5){ 
#   if (is.null(tmax)){
#     # Plot the data
#     # plot(x, y, col='indianred', xlab='x', type='l', bty='l', las=1, main=paste('trace', ii))
#     peak <- peak.fun(y=y, dt=dt, stimulation_time=stimulation_time, baseline=baseline, smooth=smooth)

#     # Get one fit and assess
#     ind1 <- (stimulation_time - baseline)/dt
#     ind2 <- stimulation_time/dt
#     y2plot <- y - mean(y[ind1:ind2])

#     dev.new(width=width, height=height, noRStudioGD=TRUE)
#     Y <- y2plot[ind1:length(y2plot)]
#     X <- seq(0, dt * (length(Y) - 1), by = dt)
#     plot(X, Y, col='indianred', xlab='time (ms)', type='l', bty='l', las=1, main='')
#     abline(h = 0, col = 'black', lwd = 1, lty=1)
#     abline(h = peak * y_abline, col = 'black', lwd = 1, lty=3)

#     # Add a label to the abline
#     text(x=max(X[ind1:length(X)]) * 0.95, y=peak * (y_abline - 0.025), labels=y_abline, pos=4)

#     # Prompt user for the range of x to use for nFIT
#     x_limit <- NA
#     while (is.na(x_limit)) {
#       cat('\nEnter the upper limit for time to use in nFIT (e.g., 400 ms): ')
#       x_limit <- as.numeric(readLines(n = 1))
#       if (is.na(x_limit)) {
#         cat('\nInvalid input. Please enter a numeric value.\n')
#       }
#     }
#     dev.off()
#   }else{
#     x_limit <- tmax
#   }

#   x_limit <- x_limit + stimulation_time - baseline
#   return(x_limit)
# }

# determine_tmax <- function(y, N=1, dt=0.1, stimulation_time=0, baseline=0, smooth=5, tmax=NULL, y_abline=0.1, height=5, width=5, prompt=TRUE) { 
#   if (is.null(tmax)){
#     peak <- peak.fun(y=y, dt=dt, stimulation_time=stimulation_time, baseline=baseline, smooth=smooth)

#     ind1 <- as.integer((stimulation_time - baseline)/dt)
#     ind2 <- as.integer(stimulation_time/dt)
#     y2plot <- y - mean(y[ind1:ind2])

#     dev.new(width=width, height=height, noRStudioGD=TRUE)
#     Y <- y2plot[ind1:length(y2plot)]
#     X <- seq(0, dt * (length(Y) - 1), by = dt)

#     out <- abline_fun(X, Y, N=N, y_abline=y_abline) 
#     A_abline <- out[1]
#     avg_t.abline <- out[2]
#     avg_t.abline <- if (is.na(avg_t.abline)) max(X) else avg_t.abline

#     plot(X, Y, col='indianred', xlab='time (ms)', type='l', bty='l', las=1, main='')
#     abline(h = 0, col = 'black', lwd = 1, lty=1)
    
#     # Get the left and bottom of the plot
#     left_axis <- par("usr")[1]
#     bottom_axis <- par("usr")[3]
    
#     # Horizontal dotted line from X=0 to X=avg_t.abline at height A_abline
#     lines(c(left_axis, avg_t.abline), c(A_abline, A_abline), col = 'black', lwd = 1, lty = 3)
    
#     # Vertical dotted line down to the bottom of the plot
#     lines(c(avg_t.abline, avg_t.abline), c(A_abline, bottom_axis), col = 'black', lwd = 1, lty = 3)

#     # Add a label to the abline
#     ind3 <- as.integer(avg_t.abline/dt)
#     text(x=max(X[ind1:ind3])*1.05, y=A_abline * 1.2, labels=paste0(y_abline*100, ' %'), pos=4, cex=0.6)

#     text(x=max(X[ind1:ind3])*1.05, y=bottom_axis*0.95, labels=paste0(avg_t.abline, ' ms'), pos=4, cex=0.6)

#     if (prompt) {
#       x_limit <- NA
#       while (is.na(x_limit)) {
#         cat('\nEnter the upper limit for time to use in nFIT (e.g., 400 ms): ')
#         x_limit <- as.numeric(readLines(n = 1))
#         if (is.na(x_limit)) {
#           cat('\nInvalid input. Please enter a numeric value.\n')
#         }
#       }
#       dev.off()
#     } else {
#       x_limit <- avg_t.abline
#     }
#   } else {
#     x_limit <- tmax
#   }

#   x_limit <- x_limit + stimulation_time - baseline
#   return(x_limit)
# }


determine_tmax <- function(y, N=1, dt=0.1, stimulation_time=0, baseline=0, smooth=5, tmax=NULL, y_abline=0.1, ylab=NULL, height=5, width=5, prompt=TRUE) { 
  if (is.null(tmax)){
    peak <- peak.fun(y=y, dt=dt, stimulation_time=stimulation_time, baseline=baseline, smooth=smooth)
    
    ind1 <- as.integer((stimulation_time - baseline)/dt)
    ind2 <- as.integer(stimulation_time/dt)
    y2plot <- y - mean(y[ind1:ind2])
    
    dev.new(width=width, height=height, noRStudioGD=TRUE)
    Y <- y2plot[ind1:length(y2plot)]
    X <- seq(0, dt * (length(Y) - 1), by = dt)
    
    out <- abline_fun(X, Y, N=N, y_abline=y_abline) 
    A_abline <- out[1]
    avg_t.abline <- out[2]
    avg_t.abline <- if (is.na(avg_t.abline)) max(X) else avg_t.abline
    
    plot(X, Y, col='indianred', xlab='time (ms)', ylab=ylab, type='l', bty='l', las=1, main='')
    abline(h = 0, col = 'black', lwd = 1, lty=1)
    
    # Get the left and bottom of the plot
    left_axis <- par("usr")[1]
    bottom_axis <- par("usr")[3]
    
    # Horizontal dotted line from X=0 to X=avg_t.abline at height A_abline
    lines(c(left_axis, avg_t.abline), c(A_abline, A_abline), col = 'black', lwd = 1, lty = 3)
    
    # Vertical dotted line down to the bottom of the plot for the abline
    lines(c(avg_t.abline, avg_t.abline), c(A_abline, bottom_axis), col = 'black', lwd = 1, lty = 3)
    
    # Add a label to the abline
    ind3 <- as.integer(avg_t.abline/dt)
    text(x=max(X[ind1:ind3])*1.05, y=A_abline * 1.2, labels=paste0(y_abline*100, ' %'), pos=4, cex=0.6)
    text(x=max(X[ind1:ind3])*1.05, y=bottom_axis*0.95, labels=paste0(avg_t.abline, ' ms'), pos=4, cex=0.6)
    
    # asterisk and label stim
    stim_index <- round(baseline/dt) + 1
    if (stim_index > length(X)) stim_index <- length(X)
    points(X[stim_index], Y[stim_index], pch=8, col='darkgray', cex=1)  # pch=8 is a star symbol.

    # Place the label "stim" on the same y level as the star, slightly to the right.
    x_offset <- 0.02 * diff(range(X))  # horizontal offset based on the range of X
    text(x = X[stim_index] + x_offset, y = Y[stim_index], labels = "stim", pos = 4, col = 'darkgray', cex = 0.6)    

    if (prompt) {
      x_limit <- NA
      while (is.na(x_limit)) {
        cat('\nEnter the upper limit for time to use in nFIT (e.g., 400 ms): ')
        x_limit <- as.numeric(readLines(n = 1))
        if (is.na(x_limit)) {
          cat('\nInvalid input. Please enter a numeric value.\n')
        }
      }
      dev.off()
    } else {
      x_limit <- avg_t.abline
    }
  } else {
    x_limit <- tmax
  }
  
  x_limit <- x_limit + stimulation_time - baseline
  return(x_limit)
}


abline_fun <- function(x, y, N = 1, y_abline = 0.1) {
  sign <- sign_fun(y)
  y <- y * sign
  n <- length(y)
  
  if (n < 2) return(NA)  # Ensure input length is valid
  
  dx <- x[2] - x[1]
  k <- round(ceiling(1 / dx))

  # Ensure we have enough points to perform the moving average
  if (n < k) return(NA)
  
  # Vectorized smoothing using a moving average filter
  ysmooth <- stats::filter(y, rep(1/k, k), sides = 2) # Use two-sided filtering for better alignment
  
  # Ensure the smoothed signal is of the correct length
  ysmooth <- as.numeric(ysmooth) # Convert from 'ts' object to numeric
  ysmooth[is.na(ysmooth)] <- 0 # Replace NAs from edges with 0 (or you can opt for another value)
  
  # Identify the first and Nth peaks
  idx.1st <- ifelse(N == 1 || is.na(first.peak(ysmooth)), which.max(abs(ysmooth)), first.peak(ysmooth))
  idx.Nth <- ifelse(N == 1 || is.na(last.peak(ysmooth)), idx.1st, last.peak(ysmooth))

  # Check if peaks were found
  if (is.na(idx.1st) || is.na(idx.Nth)) return(NA)

  A.Nth <- ysmooth[idx.Nth]
  
  # Calculate threshold amplitude level based on y_abline
  A_abline <- y_abline * A.Nth 
    
  # Find indices where the smoothed signal crosses the threshold during the decaying phase
  abline_idx <- which(diff(sign(ysmooth[idx.Nth:length(y)] - A_abline)) != 0) + idx.Nth
  
  # If no crossing points are found, return NA
  if (length(abline_idx) == 0) return(NA)

  # Find corresponding times and average if there are multiple crossings
  t.abline <- x[abline_idx]
  avg_t.abline <- median(t.abline)
  A_abline <- sign * A_abline

  return(c(A_abline, avg_t.abline))
}

sequential_fit <- function(response, n=30, dt=0.1, func=product2N, N=1, IEI=50, method= c('LM', 'BF.LM', 'GN', 'port', 'robust', 'MLE'), weight_method=c('none', '~y_sqrt', '~y'),
    stimulation_time=0, baseline=0, tmax=NULL, y_abline=0.1, fast.decay.limit=NULL, fast.constraint=FALSE, fast.constraint.method=c('rise', 'peak'), first.delay.constraint=FALSE,
    latency.limit=NULL, lower=NULL, upper=NULL,  filter=FALSE, fc=1000, interval=c(0.1, 0.9), MLEsettings=list(iter=1e4, metropolis.scale=1.5, fit.attempts=100, RWm=FALSE), 
    MLE.method=c('L-BFGS-B', 'Nelder-Mead', 'BFGS','CG', 'SANN', 'Brent'), smooth=5, response_sign_method = c('smooth', 'regression', 'cumsum'), dp=3, 
    lwd=1.2, xlab='time (ms)', ylab='PSC (pA)', width=5, height=5, return.output=FALSE, show.output=TRUE, show.plot=TRUE, seed=42){
  
  y <- response
  x_limit <- determine_tmax(y=y, N=N, dt=dt, stimulation_time=stimulation_time, baseline=baseline, smooth=smooth, tmax=tmax, y_abline=y_abline, width=width, height=height)

  x <- seq(0, (length(y) - 1) * dt, by = dt)

  # Adjust response based on user input
  y2fit <- y[x < x_limit]
  x2fit <- seq(0, dt * (length(y2fit)-1), by = dt)

  ind1 <- (stimulation_time - baseline)/dt
  ind2 <- stimulation_time/dt
  ind3 <- length(y2fit)
  # ind4 <- length(y)
   
  # bl <- mean(y[ind1:ind2])

  # xfit <- seq(0, dt * (ind3 - ind2), by = dt)
  out <- nFIT(response=y2fit, n=n, dt=dt, N=N, IEI=IEI, func=func, method=method, weight_method=weight_method, stimulation_time=stimulation_time, baseline=baseline, 
    fast.decay.limit=fast.decay.limit, fast.constraint=fast.constraint, fast.constraint.method=fast.constraint.method, 
    first.delay.constraint=first.delay.constraint, latency.limit=latency.limit, lower=lower, upper=upper, filter=filter, 
    fc=fc, interval=interval, MLEsettings=MLEsettings, MLE.method=MLE.method, response_sign_method=response_sign_method, 
    dp=dp, lwd=lwd, xlab=xlab, ylab=ylab, width=width, height=height, return.output=return.output, show.output=FALSE, show.plot=FALSE, seed=seed)

  if (show.plot) fit_plot(traces=out$traces, func=func, xlab=xlab, ylab=ylab, lwd=lwd, filter=filter, width=width, height=height)
  if (show.output) print(out$output)

  if (!fast.constraint && (!identical(func, product1))){
      
      # Define the specific message based on the fast.constraint.method
      message_part <- switch(fast.constraint.method,
                       rise = "the fastest rise",
                       peak = "the fastest time to peak")

      # Prompt user if they want to repeat with fast.rise.constraint=TRUE
      cat('\nDo you want to repeat with "fast constraint" turned on?',
          '\nThis constraint ensures the response with the fastest decay also has', message_part, "(y/n): ")
      
      repeat_with_constraint <- tolower(readLines(n = 1))
      
      if (repeat_with_constraint == 'y') {
      dev.off()
      out <- nFIT(response=y2fit, n=n, dt=dt, N=N, IEI=IEI, func=func, method=method, weight_method=weight_method, stimulation_time=stimulation_time, baseline=baseline, 
        fast.decay.limit=fast.decay.limit, fast.constraint=TRUE, fast.constraint.method=fast.constraint.method, 
        first.delay.constraint=first.delay.constraint,latency.limit=latency.limit, lower=lower, upper=upper, filter=filter, 
        fc=fc, interval=interval, MLEsettings=MLEsettings, MLE.method=MLE.method, response_sign_method=response_sign_method, 
        dp=dp, lwd=lwd, xlab=xlab, ylab=ylab, width=width, height=height, return.output=return.output, show.output=FALSE, show.plot=FALSE, seed=seed)

      if (show.plot) fit_plot(traces=out$traces, func=func, xlab=xlab, ylab=ylab, lwd=lwd, filter=filter, width=width, height=height)
      if (show.output) print(out$output)
    }
  }
 
  N <- length(out$fits)
  params <- out$fits[(N-3):N]

  xfit <- seq(0, dt * (ind3 - ind2), by = dt)
  yfit <- product1(params, xfit)

  ynew <- c(y2fit[1:(ind2-1)], (y2fit[ind2:ind3] - yfit))

  traces <- data.frame(x=x2fit, y=y2fit, yreduced=ynew)
    
  # plot(x2fit,y2fit,type='l')
  # lines(x2fit, ynew, type='l')

  return(list(params=params, traces=traces, tmax=x_limit))

}


nFIT_sequential <- function(response, n=30, N=1, IEI=50, dt=0.1, func=product2N, method= c('LM', 'BF.LM', 'GN', 'port', 'robust', 'MLE'), 
    weight_method = c('none', '~y_sqrt', '~y'), stimulation_time=0, baseline=0, fit.limits=NULL, fast.decay.limit=NULL, fast.constraint=FALSE, 
    fast.constraint.method=c('rise', 'peak'), first.delay.constraint=FALSE, latency.limit=NULL, lower=NULL, upper=NULL, filter=FALSE, fc=1000, interval=c(0.1, 0.9), 
    MLEsettings=list(iter=1e4, metropolis.scale=1.5, fit.attempts=100, RWm=FALSE), MLE.method=c('L-BFGS-B', 'Nelder-Mead', 'BFGS','CG', 'SANN', 'Brent'), 
    response_sign_method = c('smooth', 'regression', 'cumsum'), dp=3, lwd=1.2, xlab='time (ms)', ylab='PSC (pA)', width=5, height=5, return.output=FALSE, 
    show.output=TRUE, show.plot=TRUE, seed=42){

  fast.constraint.method <- match.arg(fast.constraint.method)

  y <- response
  if (identical(func, product1N)){
    functions <- list(product1N)
  }else if (identical(func, product2N)){
    functions <- list(product2N, product1N)
  }else if (identical(func, product3N)){
    functions <- list(product3N, product2N, product1N)
  }

  if (!is.null(fit.limits) && length(fit.limits) != length(functions)){
    warning("fit.limits must contain same number of elements as fitted product equations; fit.limits will be ignored")
    fit.limits=NULL
  }
  outputs <- vector("list", length(functions))

  for (ii in 1:length(functions)) {
    if (ii == 1) {
        y2fit <- y
    } else {
        y2fit <- outputs[[ii - 1]]$traces$yreduced
        if (!is.null(upper)) upper <- upper[1:(length(upper) - 4)]
        if (!is.null(lower)) lower <- lower[1:(length(lower) - 4)]
    }
    
    # fit.limits is used to reproduce the fits rapidly; do not show plots - stimulation_time + baseline
    tmax <- if (is.null(fit.limits)) NULL else fit.limits[ii] #  + stimulation_time - baseline
    output_logic <- if (is.null(tmax)) TRUE else FALSE 
    outputs[[ii]] <- sequential_fit(response=y2fit, n=n, N=N, IEI=IEI, dt=dt, func=functions[[ii]], method=method, weight_method=weight_method, stimulation_time=stimulation_time, baseline=baseline, tmax=tmax, 
                    fast.decay.limit=fast.decay.limit, fast.constraint=fast.constraint, fast.constraint.method=fast.constraint.method, first.delay.constraint=first.delay.constraint,
                    latency.limit=latency.limit, lower=lower, upper=upper, filter=filter, fc=fc, interval=interval, MLEsettings=MLEsettings, MLE.method=MLE.method, 
                    response_sign_method=response_sign_method, dp=dp, lwd=lwd, xlab=xlab, ylab=ylab, width=width, height=height, return.output=TRUE, show.output=output_logic, show.plot=output_logic, seed=seed)

    if (is.null(tmax)) { # Prompt user if they want to repeat
      cat('Do you want to repeat fit with new time base? (y/n): ')
      repeat_fit <- tolower(readLines(n = 1))
      while (repeat_fit == 'y') {
          dev.off()
          outputs[[ii]] <- sequential_fit(response=y2fit, n=n, N=N, IEI=IEI, dt=dt, func=functions[[ii]], method=method, stimulation_time=stimulation_time, baseline=baseline, tmax=tmax, 
                        fast.decay.limit=fast.decay.limit, fast.constraint=fast.constraint, fast.constraint.method=fast.constraint.method, first.delay.constraint=first.delay.constraint,
                        latency.limit=latency.limit, lower=lower, upper=upper, filter=filter, fc=fc, interval=interval, MLEsettings=MLEsettings, MLE.method=MLE.method, 
                        response_sign_method=response_sign_method, dp=dp, lwd=lwd, xlab=xlab, ylab=ylab, width=width, height=height, return.output=TRUE, show.output=output_logic, show.plot=output_logic, seed=seed)

          cat('Do you want to repeat fit with new time base? (y/n): ')
        repeat_fit <- tolower(readLines(n = 1))
        }
        dev.off()
    }
  }
  
  fits_list <- lapply(outputs, function(out) out$params)
  # Sort the list by the increasing decay (3rd element)
  fits_list_sorted <- fits_list[order(sapply(fits_list, function(x) x[3]))]
  fits <- unlist(fits_list_sorted)

  t_limits <- sapply(outputs, function(out) out$tmax)

  # if (identical(func, product1N)){
  #   df_output <- out.fun(params=fits[1:4], interval=interval, dp=dp, sign=1)
  # } else if (identical(func, product2)){  
  #   df_output1 <- out.fun(params=fits[1:4], interval=interval, dp=dp, sign=1)
  #   df_output2 <- out.fun(params=fits[5:8], interval=interval, dp=dp, sign=1)
  #   df_output <- rbind('fast' = df_output1, 'slow' = df_output2)
  # } else if (identical(func, product3N)){  
  #   # Compute the output data frames
  #   df_output1 <- out.fun(params=fits[1:4],  interval=interval, dp=dp, sign=1)
  #   df_output2 <- out.fun(params=fits[5:8],  interval=interval, dp=dp, sign=1)
  #   df_output3 <- out.fun(params=fits[9:12], interval=interval, dp=dp, sign=1)
  #   df_output <- rbind('fast' = df_output1, 'medium' =  df_output2, 'slow' =  df_output3)
  # }
  
  if (identical(func, product1N)){
    
    df_output <- out.fun(params=fits[1:(N+3)], interval=interval, dp=dp, sign=1)
  
  } else if (identical(func, product2N)){  
    df_output1 <- out.fun(params=fits[1:(N+3)], interval=interval, dp=dp, sign=1)
    df_output2 <- out.fun(params=fits[(N+4):(2*N+6)], interval=interval, dp=dp, sign=1)
      # Create a list of the outputs
    output_list <- list(df_output1, df_output2)
    # Extract the third elements from each output and determine the order
    order_indices <- order(sapply(output_list, function(x) x[[N+2]]))
    # Reorder the output list based on the third element
    output_ordered <- output_list[order_indices]

    # Combine the outputs in the correct order
    df_output <- rbind('fast' = output_ordered[[1]], 'slow' = output_ordered[[2]])



  } else if (identical(func, product3N)){ 

    df_output1 <- out.fun(params=fits[1:(N+3)], interval=interval, dp=dp, sign=1)
    df_output2 <- out.fun(params=fits[(N+4):(2*N+6)], interval=interval, dp=dp, sign=1)
    df_output3 <- out.fun(params=fits[(2*N+7):(3*N+9)], interval=interval, dp=dp, sign=1)

    # Create a list of the outputs
    output_list <- list(df_output1, df_output2, df_output3)
    # Extract the third elements from each output and determine the order
    order_indices <- order(sapply(output_list, function(x) x[[N+2]]))
    # Reorder the output list based on the third element
    output_ordered <- output_list[order_indices]

    # Combine the outputs in the correct order
    df_output <- rbind('fast' = output_ordered[[1]], 'medium' = output_ordered[[2]], 'slow' = output_ordered[[3]])
  }

  traces <- traces_fun2(y=y, fits=fits, dt=dt, N=N, IEI=IEI, stimulation_time=stimulation_time, baseline=baseline, func=func, filter=filter, fc=fc)

  if (show.plot) fit_plot(traces=traces, func=func, xlab=xlab, ylab=ylab, lwd=lwd, filter=filter, width=width, height=height)
  if (show.output) print(df_output)
    
  idx1 <- baseline/dt
  idx2 <- max(t_limits)/dt

  k = length(fits)
  gof.se <- sqrt(sum((traces$y[idx1:idx2] - traces$yfit[idx1:idx2])^2) / (length(traces$y[idx1:idx2])-k))
  msc <- model.selection.criteria(coeffs=fits, x=traces$x[idx1:idx2]-traces$x[idx1], y=traces$y[idx1:idx2], func=func, N=N, IEI=IEI)

  if (return.output) {
    out <- list(output=df_output, fits=fits, gof=gof.se, AIC=msc[1], BIC=msc[2], traces=traces, fit.limits=t_limits - stimulation_time + baseline)
    return(out)
  }
}

analyse_PSC <- function(response, dt=0.1, n=30, N=1, IEI=50, stimulation_time=150, baseline=50, smooth=5, func=product2N,  method=c("BF.LM", "LM", "GN", "port", "robust", "MLE"), 
  weight_method=c('none', '~y_sqrt', '~y'), sequential.fit=FALSE, fit.limits=NULL, MLEsettings=list(iter=1e3, metropolis.scale=1.5, fit.attempts=10, RWm=FALSE), 
  filter=FALSE, fc=1000, interval=c(0.1, 0.9), lower=NULL, upper=NULL,  fast.decay.limit=NULL, fast.constraint=FALSE, fast.constraint.method=c('rise', 'peak'), 
  first.delay.constraint=FALSE, latency.limit=NULL, rel.decay.fit.limit=0.1, half_width_fit_limit=500, dp=3, lwd=1.2, xlab='time (ms)', ylab='PSC (pA)', return.output=TRUE, height=5, width=5, seed=42) {
  
  y <- response
  if (all(is.na(y[(which(!is.na(y))[length(which(!is.na(y)))] + 1):length(y)]))) {
    y <- y[!is.na(y)]
  }

  x <- seq(0, (length(y) - 1) * dt, by = dt)

  if (!sequential.fit){
    
    tmax <- fit.limits
    x_limit <- determine_tmax(y=y, N=N, dt=dt, stimulation_time=stimulation_time, baseline=baseline, smooth=smooth, tmax=tmax, y_abline=rel.decay.fit.limit, ylab=ylab, width=width, height=height)   
   
    adjusted_response <- y[x < x_limit]
    
    # Execute nFIT
    out <- nFIT(response=adjusted_response, n=n, N=N, IEI=IEI, dt=dt, func=func, method=method, weight_method=weight_method, MLEsettings=MLEsettings, stimulation_time=stimulation_time, baseline=baseline, 
      filter=filter, fc=fc, interval=interval, fast.decay.limit=fast.decay.limit, fast.constraint=fast.constraint, fast.constraint.method=fast.constraint.method,
      first.delay.constraint=first.delay.constraint, lower=lower, upper=upper, latency.limit=latency.limit, return.output=TRUE, show.plot=FALSE, half_width_fit_limit=half_width_fit_limit, dp=dp, height=height, width=width, seed=seed)

    out$traces <- traces_fun2(y=y, fits=out$fits, dt=dt, N=N, IEI=IEI, stimulation_time=stimulation_time, baseline=baseline, func=func, filter=filter, fc=fc)
    fit_plot(traces=out$traces, func=func, xlab=xlab, ylab=ylab, lwd=lwd, filter=filter, width=width, height=height)
    
    if (!fast.constraint && (!identical(func, product1)) && is.null(fit.limits)){
      
      fast.constraint.method <- match.arg(fast.constraint.method)
      # Define the specific message based on the fast.constraint.method
      message_part <- switch(fast.constraint.method,
                       rise = "the fastest rise",
                       peak = "the fastest time to peak")

      # Prompt user if they want to repeat with fast.rise.constraint=TRUE
      cat('\nDo you want to repeat with "fast constraint" turned on?',
          '\nThis constraint ensures the response with the fastest decay also has', message_part, "(y/n): ")

      repeat_with_constraint <- tolower(readLines(n = 1))
      
      if (repeat_with_constraint == 'y') {
        dev.off()
        out <- nFIT(response=adjusted_response, n=n,  N=N, IEI=IEI,dt=dt, func=func, method=method, weight_method=weight_method, MLEsettings=MLEsettings, stimulation_time=stimulation_time, baseline=baseline, 
          filter=filter, fc=fc, interval=interval, fast.decay.limit=fast.decay.limit, fast.constraint=TRUE, fast.constraint.method=fast.constraint.method, half_width_fit_limit=half_width_fit_limit,
          first.delay.constraint=first.delay.constraint, lower=lower, upper=upper, latency.limit=latency.limit, return.output=TRUE, show.plot=FALSE, dp=dp, height=height, width=width, seed=seed)

        out$traces <- traces_fun2(y=y, fits=out$fits, dt=dt, N=N, IEI=IEI, stimulation_time=stimulation_time, baseline=baseline, func=func, filter=filter, fc=fc)
        fit_plot(traces=out$traces, func=func, xlab=xlab, ylab=ylab, lwd=lwd, filter=filter, width=width, height=height)
    
      }
    }
    
  }else{
    
    out <- nFIT_sequential(response=y, n=n, dt=dt, func=func, method=method, weight_method=weight_method, stimulation_time=stimulation_time, baseline=baseline, fit.limits=fit.limits, 
      fast.decay.limit=fast.decay.limit, fast.constraint=fast.constraint, fast.constraint.method=fast.constraint.method, first.delay.constraint=first.delay.constraint,
      latency.limit=latency.limit, lower=lower, upper=upper, filter=filter, fc=fc, interval=interval, MLEsettings=MLEsettings, MLE.method=MLE.method,  half_width_fit_limit=half_width_fit_limit, 
      dp=dp, lwd=lwd, xlab=xlab, ylab=ylab, width=width, height=height, return.output=TRUE, show.output=TRUE, show.plot=TRUE, seed=seed)
  }

  if (return.output) return(out)

}

traces_fun <- function(y, fits, dt=0.1,  stimulation_time=150, baseline=50, func=product2, filter=FALSE, fc=1000, height=5, width=5){
  
  dx <- dt  
  x <- seq(0, (length(y) - 1) * dx, by = dx)

  if (filter){
    ind = 20
    fc = fc; fs = 1/dx*1000; bf <- butter(2, fc/(fs/2), type='low')
    yfilter <- signal::filter(bf, y)
  } else {
    ind=1
    yfilter=y
  }

  ind1 <- (stimulation_time - baseline)/dx
  ind2 <- baseline/dx
  
  yorig <- y[ind1:length(y)]
  yfilter <- yfilter[ind1:length(yfilter)]
  xorig <- seq(0, dx * (length(yorig) - 1), by = dx)


  yorig <- yorig - mean(yorig[1:ind2])
  yfilter <- yfilter - mean(yfilter[1:ind2])

  traces=data.frame(x=xorig, y=yorig, yfilter=yfilter)
  if (identical(func, product1)){
    fits[4] <- fits[4] + baseline
  } else if (identical(func, product2)){
    fits[4] <- fits[4] + baseline; fits[8] <- fits[8] + baseline
  } else if (identical(func, product3)){
    fits[4] <- fits[4] + baseline; fits[8] <- fits[8] + baseline; fits[12] <- fits[12] + baseline
  }    
  traces$yfit <- func(fits,traces$x+dx)  
  if (identical(func, product2)){
    traces$yfit1 <- product1(fits[1:4],traces$x+dx) 
    traces$yfit2 <- product1(fits[5:8],traces$x+dx) 
  } 
  if (identical(func, product3)){
    traces$yfit1 <- product1(fits[1:4],traces$x+dx) 
    traces$yfit2 <- product1(fits[5:8],traces$x+dx) 
    traces$yfit3 <- product1(fits[9:12],traces$x+dx) 
  } 

  return(traces)
}

traces_fun2 <- function(y, fits, dt=0.1,  N=1, IEI=50, stimulation_time=150, baseline=50, func=product2N, filter=FALSE, fc=1000){
  
  dx <- dt  
  x <- seq(0, (length(y) - 1) * dx, by = dx)

  if (filter){
    ind = 20
    fc = fc; fs = 1/dx*1000; bf <- butter(2, fc/(fs/2), type='low')
    yfilter <- signal::filter(bf, y)
  } else {
    ind=1
    yfilter=y
  }

  ind1 <- (stimulation_time - baseline)/dx
  ind2 <- baseline/dx
  
  yorig <- y[ind1:length(y)]
  yfilter <- yfilter[ind1:length(yfilter)]
  xorig <- seq(0, dx * (length(yorig) - 1), by = dx)

  yorig <- yorig - mean(yorig[1:ind2])
  yfilter <- yfilter - mean(yfilter[1:ind2])

  traces <- data.frame(x = xorig, y = yorig, yfilter = yfilter)
  
  # Define the list of functions
  func_list <- list(product1N, product2N, product3N)

  # Find the index of the matching function
  func_index <- sapply(func_list, function(f) identical(f, func))

  # Check if func matches one of the productN functions and apply the corresponding baseline adjustments
  if (any(func_index)) {
    index <- which(func_index)  # Get the index of the matching function
    for (i in 1:index) {
      fits[i * N + (3 * i)] <- fits[i * N + (3 * i)] + baseline
    }
  }

  traces$yfit <- func(params=fits, x=traces$x+dx, N=N, IEI=IEI) 

  if (identical(func, product2N)){
    traces$yfit1 <- product1N(params=fits[1:(N+3)], x=traces$x+dx, N=N, IEI=IEI) 
    traces$yfit2 <- product1N(params=fits[(N+4):(2*N+6)], x=traces$x+dx, N=N, IEI=IEI) 
  } 
  if (identical(func, product3N)){
    traces$yfit1 <- product1N(params=fits[1:(N+3)],  x=traces$x+dx, N=N, IEI=IEI) 
    traces$yfit2 <- product1N(params=fits[(N+4):(2*N+6)],  x=traces$x+dx, N=N, IEI=IEI) 
    traces$yfit3 <- product1N(params=fits[(2*N+7):(3*N+9)], x=traces$x+dx, N=N, IEI=IEI) 
  } 

  return(traces)
}

# Loop through each dataset
datasets2list <- function(datasets, id='_fits'){
  out <- list()
  for (dataset in datasets) {
    # Construct the name of the fits object
    name <- paste0(dataset, id)
    
    # Check if the fits object exists
    if (exists(name)) {
      # Add the fits object to the list
      out[[dataset]] <- get(name)
    } else {
      cat(paste("Object", name, "does not exist.\n"))
    }
  }
  return(out)
}


# function to extract from list of summaries
extract_variable <- function(data_list, id = 'area', rename_datasets = TRUE) {
  # Initialize empty vectors to store the combined data
  dataset_names <- c()
  values1 <- c()
  values2 <- c()
  subjects <- c()
  
  # Create a mapping for dataset renaming based on the order in the data_list
  dataset_mapping <- setNames(seq_along(names(data_list)), names(data_list))
  
  # Loop through each item in data_list
  for (name in names(data_list)) {
    # Extract the dataframe
    df <- data_list[[name]]
    
    # Extract the columns with the specified id
    variable_columns <- df[, grep(paste0("^", id, "$"), colnames(df))]
    
    # Ensure we have exactly two columns with the specified id
    if (ncol(variable_columns) == 2) {
      # Determine the dataset name
      dataset_name <- if (rename_datasets) {
        dataset_mapping[[name]]
      } else {
        name
      }
      
      # Append data to vectors
      dataset_names <- c(dataset_names, rep(dataset_name, nrow(df)))
      values1 <- c(values1, variable_columns[, 1])
      values2 <- c(values2, variable_columns[, 2])
    } else {
      cat(paste("Dataframe", name, "does not contain exactly two '", id, "' columns.\n"))
    }
  }

  # Create the combined dataframe
  combined_df <- data.frame(
    dataset = dataset_names,
    s = 1:length(dataset_names),
    value1 = values1,
    value2 = values2
  )
  
  # Rename the value1 and value2 columns to id1 and id2
  colnames(combined_df)[3:4] <- paste0(id, 1:2)

  return(combined_df)
}


boxplot_calculator <- function(data, type = 6, na.rm=FALSE) {
  unique_x <- unique(data$x)
  result <- data.frame(x = numeric(), Q1 = numeric(), Q3 = numeric(), Median = numeric(), Min = numeric(), Max = numeric(), MAD = numeric())
  
  for (i in 1:length(unique_x)) {
    current_x <- unique_x[i]
    d <- data$y[data$x == current_x]
    
    q1 <- quantile(d, probs = 0.25, type = type, na.rm = na.rm)
    q3 <- quantile(d, probs = 0.75, type = type, na.rm = na.rm)
    iqr <- q3 - q1  # Calculate IQR
    
    lower_bound <- q1 - 1.5 * iqr  # Lower bound for outliers
    upper_bound <- q3 + 1.5 * iqr  # Upper bound for outliers
    
    # Exclude outliers
    d_filtered <- d[d >= lower_bound & d <= upper_bound]
    
    median_val <- median(d, na.rm = na.rm)
    min_val <- min(d_filtered, na.rm = na.rm)
    max_val <- max(d_filtered, na.rm = na.rm)
    
    # Calculate MAD
    mad <- median(abs(d - median_val), na.rm = na.rm)
    
    result <- rbind(result, data.frame(x = current_x, Q1 = q1, Q3 = q3, Median = median_val, Min = min_val, Max = max_val, MAD = mad))
  }
  
  rownames(result) <- NULL  # Remove row names
  return(result)
}


WBplot <- function(data, wid = 0.2, cap = 0.05, xlab = '', ylab = 'PSP amplitude (mV)', 
                   xrange = c(0.75, 2.25), yrange = c(0, 400), main = '', tick_length = 0.02, 
                   x_tick_interval = NULL, y_tick_interval = 100, lwd = 0.8, type = 6, na.rm=FALSE) {
  
  boxplot_values <- boxplot_calculator(data=data, type=type, na.rm=na.rm)
  
  if (is.null(x_tick_interval)) {
    x_ticks <- unique(data$x)
  } else {
    x_ticks <- seq(xrange[1], xrange[2], by = x_tick_interval)
  }
  xrange <- xrange + c(-wid, wid)
  
  # Ensure background is off and plot area is clear
  par(bg = NA)
  plot(1, type = 'n', ylim = yrange, xlim = xrange, xlab = xlab, ylab = ylab, 
       main = main, xaxt = 'n', yaxt = 'n', bty = 'n', lwd = lwd)
  
  for (i in 1:nrow(boxplot_values)) {
    # Convert current_x to numeric for arithmetic operations
    current_x <- as.numeric(boxplot_values$x[i])
    
    rect(current_x - wid, boxplot_values$Q1[i], current_x + wid, boxplot_values$Q3[i], col = 'white', lwd = lwd)
    segments(current_x, boxplot_values$Q1[i], current_x, boxplot_values$Min[i], lwd = lwd)
    segments(current_x, boxplot_values$Q3[i], current_x, boxplot_values$Max[i], lwd = lwd)
    segments(current_x - cap, boxplot_values$Min[i], current_x + cap, boxplot_values$Min[i], lwd = lwd)
    segments(current_x - cap, boxplot_values$Max[i], current_x + cap, boxplot_values$Max[i], lwd = lwd)
    segments(current_x - wid * 1.1, boxplot_values$Median[i], current_x + wid * 1.1, boxplot_values$Median[i], col = 'black', lwd = 3 * lwd)
  }
  
  # Set the x-axis ticks
  axis(1, at = x_ticks, labels = x_ticks, tcl = -tick_length, lwd = lwd)
  
  # Set the y-axis ticks
  y_ticks <- seq(yrange[1], yrange[2], by = y_tick_interval)
  axis(2, at = y_ticks, tcl = -tick_length, las = 1, lwd = lwd)
}


BoxPlot <- function(data, wid=0.2, cap=0.05, xlab='', ylab='PSC amplitude (pA)', main='', 
                    xrange=c(0.75,2.25), yrange=c(-400, 0), tick_length=0.2, 
                    x_tick_interval = NULL, y_tick_interval=100, lwd=1, 
                    type=6, amount=0.05, p.cex=0.5, filename='boxplot.svg', 
                    height=2.5, width=4, bg='transparent', alpha=0.6, na.rm=FALSE, save=FALSE){
  
  if (save) {
    # Open SVG device
    svg(file=filename, width=width, height=height, bg=bg)
  } else {
    dev.new(width=width, height=height, noRStudioGD=TRUE)
  }

  # data1 <- if ("s" %in% colnames(data)) data[, !colnames(data) %in% "s"] else data

  WBplot(data=data, wid=wid, cap=cap, xlab=xlab, ylab=ylab, main=main, xrange=xrange, yrange=yrange, 
         tick_length=tick_length, x_tick_interval=x_tick_interval, y_tick_interval=y_tick_interval, 
         lwd=lwd, type=type, na.rm=na.rm)
  
  set.seed(42)
  data$x_jitter <- jitter(data$x, amount=amount)
  
  # Set color with alpha transparency for points
  point_color <- rgb(169/255, 169/255, 169/255, alpha=alpha)  # darkgray with alpha=0.4
  points(data$x_jitter, data$y, pch=19, bg='transparent', col=point_color, lwd=lwd/3, cex=p.cex)

  if ("s" %in% colnames(data)) {

    # Connect data points within subjects with gray dotted lines
    line <- TRUE
    if (line) {
      subjects <- unique(data$s)
      for (subj in subjects) {
        subset_data <- data[data$s == subj, ]
        lines(subset_data$x_jitter, subset_data$y, col='darkgray', lwd=lwd, lty=3)  # lty=3 for dotted line
      }
    }
  }

  if (save) {
    dev.off()
  }
}



# BoxPlot2 <- function(formula, data, wid = 0.2, cap = 0.05, xlab = '', ylab = 'PSC amplitude (pA)', 
#                     main = '', xrange = NULL, yrange = c(-400, 0), tick_length = 0.2, 
#                     x_tick_interval = NULL, y_tick_interval = 100, lwd = 1, 
#                     type = 6, amount = 0.05, p.cex = 0.5, filename = 'boxplot.svg', 
#                     height = 2.5, width = 4, bg = 'transparent', alpha = 0.6, na.rm=FALSE, save = FALSE) {
  
#   # Parse the formula to extract response and predictors
#   response <- as.character(formula[[2]])
#   predictors <- all.vars(formula[[3]]) # Get the predictor variables
  
#   # Check if the specified columns exist in the data
#   if (!all(c(response, predictors) %in% colnames(data))) {
#     stop("The specified response or predictor variables are not found in the data.")
#   }
  
#   # Handle grouping if interaction is specified or random effects are included
#   if (any(grepl("\\|", predictors))) {
#     # Mixed effects formula with random effects
#     fixed_effects <- sub(" \\+ \\(1\\|.*\\)", "", predictors)
#     group_vars <- strsplit(fixed_effects, " \\* | \\+ ")[[1]]
#     subject_var <- gsub(".*\\|", "", predictors)
#     data$s <- as.factor(data[[subject_var]])
#   } else {
#     # Standard formula without random effects
#     group_vars <- predictors
#   }
  
#   # Determine the number of grouping factors
#   if (length(group_vars) == 1) {
#     # Single grouping variable
#     data$x <- as.factor(data[[group_vars[1]]])
#   } else {
#     # Interaction of two grouping variables
#     data$x <- interaction(data[[group_vars[1]]], data[[group_vars[2]]], sep = " : ")
#   }
  
#   # Set the response variable
#   data$y <- data[[response]]

#   # Set x range based on the unique levels of x
#   if (is.null(xrange)) {
#     xrange <- range(as.numeric(data$x)) + c(-wid, wid)
#   }
  
#   # Handle saving the plot
#   if (save) {
#     svg(file = filename, width = width, height = height, bg = bg)
#   } else {
#     dev.new(width = width, height = height, noRStudioGD = TRUE)
#   }
  
#   # Create the box plot using the WBplot function
#   WBplot(data = data, wid = wid, cap = cap, xlab = xlab, ylab = ylab, main = main, 
#          xrange = xrange, yrange = yrange, tick_length = tick_length, x_tick_interval = x_tick_interval, 
#          y_tick_interval = y_tick_interval, lwd = lwd, type = type, na.rm=na.rm)
  
#   # Jitter x-values for plotting individual points
#   set.seed(42)
#   data$x_jitter <- jitter(as.numeric(data$x), amount = amount)
  
#   # Set the color with alpha transparency for the points
#   point_color <- rgb(169 / 255, 169 / 255, 169 / 255, alpha = alpha)  # darkgray with alpha transparency
#   points(data$x_jitter, data$y, pch = 19, col = point_color, lwd = lwd / 3, cex = p.cex)
  
#   # Connect data points for repeated measures (if subject information is provided)
#   if ("s" %in% colnames(data)) {
#     subjects <- unique(data$s)
#     for (subj in subjects) {
#       subset_data <- data[data$s == subj, ]
#       lines(subset_data$x_jitter, subset_data$y, col = 'darkgray', lwd = lwd, lty = 3)  # lty=3 for dotted line
#     }
#   }
  
#   # Close the SVG device if saving
#   if (save) {
#     dev.off()
#   }
# }

BoxPlot2 <- function(formula, data, wid = 0.2, cap = 0.05,
                     xlab = '', ylab = 'PSC amplitude (pA)', main = '',
                     xrange = NULL, yrange = c(-400, 0), tick_length = 0.2,
                     x_tick_interval = NULL, y_tick_interval = 100,
                     xlabel_angle = NULL,    # NEW: angle in degrees, or NULL for horizontal
                     lwd = 1, type = 6, amount = 0.05, p.cex = 0.5,
                     filename = 'boxplot.svg', height = 2.5, width = 4,
                     bg = 'transparent', alpha = 0.6, na.rm = FALSE, save = FALSE) {

  # formula
  response   <- as.character(formula[[2]])
  predictors <- all.vars(formula[[3]])
  if (!all(c(response, predictors) %in% names(data))) {
    stop('Response or predictor not found in data.')
  }

  if (length(predictors) == 1) {
    data$x <- factor(data[[predictors[1]]])
  } else {
    data$x <- interaction(data[[predictors[1]]], data[[predictors[2]]], sep=' : ')
  }
  data$y <- data[[response]]

  if (is.null(xrange)) {
    xrange <- range(as.numeric(data$x), na.rm=TRUE) + c(-wid, wid)
  }

  if (save) {
    svg(filename, width=width, height=height, bg=bg)
    on.exit(dev.off(), add=TRUE)
  } else {
    dev.new(width=width, height=height, noRStudioGD=TRUE)
  }

  # suppress default x‐labels
  orig_axis <- graphics::axis
  assign('axis',
    function(side, at, labels = TRUE, tcl = NA, ...) {
      if (side == 1) {
        orig_axis(side, at = at, labels = FALSE, tcl = -tick_length, ...)
      } else {
        orig_axis(side, at = at, labels = labels, tcl = -tick_length, ...)
      }
    },
    envir = .GlobalEnv
  )
  on.exit(assign('axis', orig_axis, envir = .GlobalEnv), add = TRUE)

  # draw the boxplot (WBplot calls axis(1) & axis(2) internally)
  WBplot(data = data, wid = wid, cap = cap, xlab = xlab, ylab = ylab, main = main,
         xrange = xrange, yrange = yrange, tick_length = tick_length,
         x_tick_interval = x_tick_interval, y_tick_interval = y_tick_interval,
         lwd = lwd, type = type, na.rm = na.rm)

  # restore axis()
  assign('axis', orig_axis, envir = .GlobalEnv)


  # draw custom x‐labels
  x_labels <- levels(data$x)

  # strip off everything from the first ' :'
  labs <- sub('\\s*:\\s*.*$', '', x_labels)

  # draw axis
  at <- seq_along(labs)
  if (is.null(xlabel_angle)) {
    axis(1, at = at, labels = labs, tcl = -tick_length, lwd = lwd)
  } else {
    axis(1, at = at, labels = FALSE, tcl = -tick_length, lwd = lwd)
    usr <- par('usr')
    y0  <- usr[3] - 0.1 * diff(usr[3:4])
    text(x = at, y = y0, labels = labs,
         srt = xlabel_angle, adj = 1, xpd = TRUE)
  }

  # —————— jitter & paired/unpaired block ——————
  set.seed(42)
  data$x_jitter <- jitter(as.numeric(data$x), amount = amount)
  
  # identify paired subjects (>=2 non-NA y’s)
  if ('s' %in% names(data)) {
    counts      <- ave(!is.na(data$y), data$s, FUN = sum)
    data$paired <- counts >= 2
  } else {
    data$paired <- rep(FALSE, nrow(data))
  }
  
  # draw lines for repeated measures (>=2)
  if ('s' %in% names(data)) {
    for (subj in unique(data$s[data$paired])) {
      sd <- subset(data, s == subj & !is.na(y))
      sd <- sd[order(as.numeric(sd$x)), ]
      lines(sd$x_jitter, sd$y,
            col = 'darkgray', lwd = lwd, lty = 3)
    }
  }
  
  # draw all paired points (or, if no 's', all points) in gray
  gray_col <- rgb(169/255,169/255,169/255, alpha = alpha)
  points(data$x_jitter[data$paired | !('s' %in% names(data))],
         data$y    [data$paired | !('s' %in% names(data))],
         pch = 19, col = gray_col, cex = p.cex, lwd = lwd/3)
  
  # draw the strictly unpaired points in red—but only if you really had 's'
  if ('s' %in% names(data)) {
    unpaired <- !data$paired & !is.na(data$y)
    red_col  <- rgb(205/255,92/255,92/255, alpha = alpha)
    points(data$x_jitter[unpaired],
           data$y    [unpaired],
           pch = 19, col = red_col, cex = p.cex, lwd = lwd/3)
  }
  
  if (save) dev.off()
}

BoxPlot3 <- function(formula, data, wid = 0.2, cap = 0.05, xlab = '', ylab = 'PSC amplitude (pA)', main = '',
                     xrange = NULL, yrange = c(-400, 0), xlabel_angle = NULL, tick_length = 0.2, 
                     x_tick_interval = NULL, y_tick_interval = 100, lwd = 1, type = 6, amount = 0.05, p.cex = 0.5,
                     height = 2.5, width = 4, bg = 'transparent', alpha = 0.6, na_rm_subjects = FALSE,
                     test_result, alpha_level = 0.05, group_names = NULL, sig_offset = NULL) {

  
  f_str <- deparse(formula)
  has_error <- grepl('Error', f_str)

  if (has_error) {
    err_part <- sub('.*Error\\((.*)\\).*', '\\1', f_str)
    subject_var <- strsplit(err_part, '/')[[1]][1]
    subject_var <- gsub('[[:space:]]', '', subject_var)
    main_formula_str <- sub('\\+\\s*Error\\(.*\\)', '', f_str)
    main_formula <- as.formula(main_formula_str)
  } else {
    subject_var <- NULL
    main_formula <- formula
  }

  response_var <- all.vars(formula(main_formula))[1]
  predictors <- all.vars(formula(main_formula))[-1]

  if (na_rm_subjects && !is.null(subject_var)) {
    data <- df[ !ave(is.na(data[[response_var]]), data[[subject_var]], FUN = any), ]
  }






  BoxPlot2(formula = formula, data = data, wid = wid, cap = cap, xlab = xlab, ylab = ylab, xlabel_angle = xlabel_angle, main = main,
           xrange = xrange, yrange = yrange, tick_length = tick_length, y_tick_interval = y_tick_interval,
           lwd = lwd, type = type, amount = amount, p.cex = p.cex, height = height, width = width,
           bg = bg, na.rm = TRUE)

  has_error <- grepl("Error", deparse(formula))
  if (missing(test_result) || is.null(test_result) || nrow(test_result) == 0) {
    return()
  }

  # Prepare data$x & data$y
  response   <- response_var
  # predictors <- all.vars(formula[[3]])
  data$y     <- data[[response]]

  if (length(predictors) == 1) {
    data$x        <- factor(data[[predictors[1]]])
    single_factor <- TRUE
  } else {
    data$x        <- interaction(
      data[[predictors[1]]],
      data[[predictors[2]]],
      sep = ' : '
    )
    single_factor <- FALSE
  }

  x_labels    <- levels(data$x)
  x_positions <- setNames(seq_along(x_labels), x_labels)

  # compute base offset and tick
  y_span <- diff(range(data$y, na.rm = TRUE))
  offset <- if (is.null(sig_offset)) 0.05 * y_span else sig_offset
  tick   <- 0.25 * offset

  # Loop over tests
  for (i in seq_len(nrow(test_result))) {
    p_adj <- test_result$`p adjusted`[i]
    if (is.na(p_adj) || p_adj >= alpha_level) next

    parts <- strsplit(as.character(test_result$contrast[i]), ' vs ')[[1]]
    if (length(parts) != 2) next
    lev1 <- trimws(parts[1])
    lev2 <- trimws(parts[2])

    # Build the two labels that match interaction()
    if (single_factor) {
      label1 <- lev1
      label2 <- lev2
    } else {
      comp      <- as.character(test_result$comparison[i])
      outer_lev <- sub('.*? ([^ ]+) \\(.*', '\\1', comp)
      if (grepl('\\(paired\\)', comp)) {
        label1 <- paste0(lev1, ' : ', outer_lev)
        label2 <- paste0(lev2, ' : ', outer_lev)
      } else {
        label1 <- paste0(outer_lev, ' : ', lev1)
        label2 <- paste0(outer_lev, ' : ', lev2)
      }
    }

    if (!(label1 %in% x_labels) || !(label2 %in% x_labels)) next
    x1 <- x_positions[label1]
    x2 <- x_positions[label2]

    yvals <- c(data$y[data$x == label1], data$y[data$x == label2])
    yvals <- yvals[!is.na(yvals)]
    if (length(yvals) == 0) next
    y_max <- max(yvals)
    y_min <- min(yvals)

    # shift unpaired sig bars by 6*tick
    if (single_factor) {
      shift_amt <- 0
    } else {
      paired    <- grepl('\\(paired\\)', comp)
      shift_amt <- if (!paired) 6 * tick else 0
    }

    if (y_max > 0) {
      y_line <- y_max + offset + shift_amt
      segments(x1, y_line, x1, y_line - tick, lwd = lwd)
      segments(x2, y_line, x2, y_line - tick, lwd = lwd)
      text_y <- y_line + tick
    } else {
      y_line <- y_min - offset - shift_amt
      segments(x1, y_line, x1, y_line + tick, lwd = lwd)
      segments(x2, y_line, x2, y_line + tick, lwd = lwd)
      text_y <- y_line - tick
    }

    segments(x1, y_line, x2, y_line, lwd = lwd)
    text((x1 + x2) / 2, text_y, labels = '*', cex = 1.2)
  }
}


# BoxPlot3 <- function(formula, data, wid=0.2, cap=0.05, xlab='', ylab='PSC amplitude (pA)', main='', xrange=NULL, 
#   yrange=c(-400, 0), tick_length=0.2, x_tick_interval=NULL, y_tick_interval=100, lwd=1, type=6, amount=0.05, 
#     p.cex=0.5, filename='boxplot.svg', height=2.5, width=4, bg='transparent', alpha=0.6, na.rm=FALSE, 
#     test_results=NULL, alpha_level=0.05, group_names=NULL, sig_offset=NULL, save=FALSE) {
  
#   BoxPlot2(formula=formula, data=data, wid=wid, cap=cap, xlab=xlab, 
#            ylab=ylab, xrange=xrange, yrange=yrange, tick_length=tick_length, 
#            y_tick_interval=y_tick_interval, lwd=lwd, type=type, amount=amount, 
#            p.cex=p.cex, filename=filename, height=height, width=width, 
#            na.rm=na.rm, save=save)
  
#   if (!is.null(test_results)){
#     response <- as.character(formula[[2]])
#     data$y <- data[[response]]
    
#     # grouping variable (assumed to be the first predictor)
#     predictors <- all.vars(formula[[3]])
#     group_col <- predictors[1]
    
#     # if group_names is provided, recode else use the numeric levels
#     orig_levels <- sort(unique(data[[group_col]]))
#     if (!is.null(group_names)) {
#       if (length(group_names) != length(orig_levels)) {
#         stop("number of group_names should match the number of unique groups in the group column")
#       }
#       data$x <- factor(data[[group_col]], levels=orig_levels, labels=group_names)
#       groups <- group_names
#     } else {
#       groups <- as.character(orig_levels)
#       data$x <- factor(data[[group_col]], levels=groups)
#     }
    
#     # map groups to x positions
#     x_positions <- setNames(seq_along(groups), groups)
#     y_range <- diff(range(data$y, na.rm=TRUE))
    
#     # offset sig bars
#     offset <- if (is.null(sig_offset)) 0.05 * y_range else sig_offset
#     tick <- 0.25 * offset  # fixed tick height
    
#     for (i in 1:nrow(test_results)) {
#       p_val <- as.numeric(test_results[i, "p adjusted"])
#       if (p_val < alpha_level) {
#         # Determine groups to compare.
#         if (is.null(group_names)) {
#           # if group_names are provided, assumes numeric order corresponds to testing order
#           if (i >= length(groups)) {
#             warning("Test result index exceeds number of available group pairs.")
#             next
#           }
#           group1 <- groups[i]
#           group2 <- groups[i + 1]
#         } else {
#           contrast_str <- as.character(test_results[i, "contrast"])
#           groups_in_contrast <- strsplit(contrast_str, " vs ")[[1]]
#           if (length(groups_in_contrast) != 2) next
#           group1 <- groups_in_contrast[1]
#           group2 <- groups_in_contrast[2]
#         }
        
#         # get x positions
#         if (!(group1 %in% names(x_positions)) || !(group2 %in% names(x_positions))) {
#           warning(paste("One of the groups in contrast", group1, "vs", group2, "was not found."))
#           next
#         }
#         x1 <- x_positions[group1]
#         x2 <- x_positions[group2]
        
#         # Retrieve y-values for these groups.
#         y_vals_group1 <- data$y[as.character(data$x) == group1]
#         y_vals_group2 <- data$y[as.character(data$x) == group2]
#         y_vals <- c(y_vals_group1, y_vals_group2)
#         y_max <- max(y_vals, na.rm=TRUE)
#         y_min <- min(y_vals, na.rm=TRUE)
        
#         if (is.infinite(y_max)) {
#           warning(paste("No valid y-values found for contrast", group1, "vs", group2))
#           next
#         }
        
#         # significance bar above the max if positive
#         if (y_max > 0) {
#           y_line <- y_max + offset
#           segments(x0=x1, y0=y_line, x1=x1, y1=y_line - tick, lwd=lwd)
#           segments(x0=x2, y0=y_line, x1=x2, y1=y_line - tick, lwd=lwd)
#           text_y <- y_line + tick 
#         } else {
#           # significance bar below the min if negative
#           y_line <- y_min - offset
#           segments(x0=x1, y0=y_line, x1=x1, y1=y_line + tick, lwd=lwd)
#           segments(x0=x2, y0=y_line, x1=x2, y1=y_line + tick, lwd=lwd)
#           text_y <- y_line - tick 
#         }
        
#         # horizontal line connecting boxes
#         segments(x0=x1, y0=y_line, x1=x2, y1=y_line, lwd=lwd)
        
#         # 
#         star_label <- "*"
        
#         # show sig *
#         text(x=(x1 + x2) / 2, y=text_y, labels=star_label)
#       }
#     }
#   }
# }

scatter_plot <- function(scatter, xlim=c(0, 400), ylim=c(0, 400), x_tick_interval=100, y_tick_interval=100, height=4, width=4, main='',
                         colors=c("black", "indianred"), open_symbols=FALSE, lwd=1, p.cex=0.5, filename='scatter.svg', save=FALSE) {
  # Create a basic scatter plot
  if (save) {
    # Open SVG device
    svg(file=filename, width=width, height=height, bg='transparent')
  } else {
    dev.new(width=width, height=height, noRStudioGD=TRUE)
  }

  # Map levels to 1 and 2 alternately
  unique_levels <- unique(scatter$level)
  n <- length(unique_levels)
  scatter$level <- as.numeric(factor(scatter$level, levels = unique_levels, labels = rep(1:n, length.out = length(unique_levels))))

  # Determine plot symbols (only open circles if open_symbols is TRUE)
  pch <- if (open_symbols) 1 else 19

  # Plot without the box and customize axes
  xlab = expression(A[fast] * " " * (pA))
  ylab = expression(A[slow] * " " * (pA))

  cols <- hex_palette(n=n, color1=colors[1], color2=colors[2], reverse = FALSE)


  plot(scatter$x, scatter$y, col = cols[scatter$level], 
       pch = pch, cex = p.cex, xlim=xlim, ylim=ylim, xlab = xlab, ylab = ylab, 
       main = main, xaxt='n', yaxt='n', bty='n')

  # Add abline y = x as light gray; use segments as abline function extends beyond axes
  segments(min(xlim), min(ylim), max(xlim), max(ylim), col = "lightgray")

  # Define tick intervals and lengths
  x_ticks <- seq(min(xlim), max(xlim), by=x_tick_interval)
  y_ticks <- seq(min(ylim), max(ylim), by=y_tick_interval)
  tick_length <- -0.2

  # Customize x-axis
  axis(1, at=x_ticks, labels=x_ticks, tcl=tick_length, lwd=lwd)

  # Customize y-axis with horizontal labels
  axis(2, at=y_ticks, labels=y_ticks, tcl=tick_length, las=1, lwd=lwd)

  if (save) {
    dev.off()
  }
}

# Custom Theil-Sen function
theil_sen <- function(x, y) {
  if (length(x) != length(y)) {
    stop('x and y must have the same length')
  }
  
  n <- length(x)
  slopes <- c()
  
  for (i in 1:(n-1)) {
    for (j in (i+1):n) {
      if (x[j] != x[i]) {
        slopes <- c(slopes, (y[j] - y[i]) / (x[j] - x[i]))
      }
    }
  }
  
  slope <- median(slopes, na.rm=TRUE)
  intercept <- median(y - slope * x)
  
  return(list(slope=slope, intercept=intercept))
}

# Custom function to compute z-scores and p-values
compute_regression_stats <- function(estimate, bootstrap_estimates, conf_level=0.95) {
  # Standard error
  std_err <- sd(bootstrap_estimates, na.rm=TRUE)
  
  # Z-score
  z_value <- estimate / std_err
  
  # P-value (two-tailed)
  p_value <- 2 * (1 - pnorm(abs(z_value)))
  
  # Confidence intervals
  lower_ci <- quantile(bootstrap_estimates, (1 - conf_level) / 2, na.rm=TRUE)
  upper_ci <- quantile(bootstrap_estimates, 1 - (1 - conf_level) / 2, na.rm=TRUE)
  
  return(list(
    estimate=estimate,
    std_err=std_err,
    z_value=z_value,
    p_value=p_value,
    lower_ci=lower_ci,
    upper_ci=upper_ci
  ))
}

# Function to format p-values, showing as 0.00000 when below threshold
format_p_value <- function(p, dp=5) {
  if (p < 10^(-dp)) {
    return(formatC(0, format='f', digits=dp))
  } else {
    return(formatC(p, format='f', digits=dp))
  }
}

# Custom Siegel Estimator function
siegel_sen <- function(x, y) {
  if (length(x) != length(y)) {
    stop('x and y must have the same length')
  }
  
  n <- length(x)
  slopes_per_point <- numeric(n)  # Vector to store median slopes for each point
  
  for (i in 1:n) {
    slopes <- c()
    for (j in 1:n) {
      if (i != j && x[j] != x[i]) {
        slopes <- c(slopes, (y[j] - y[i]) / (x[j] - x[i]))
      }
    }
    slopes_per_point[i] <- median(slopes, na.rm=TRUE)  # Store the median slope for point i
  }
  
  slope <- median(slopes_per_point, na.rm=TRUE)  # Final slope is the median of individual medians
  intercept <- median(y - slope * x, na.rm=TRUE)
  
  return(list(slope=slope, intercept=intercept))
}

theil_sen_with_ci <- function(x, y, n_bootstrap=1e4, conf_level=0.95, seed=42, dp=5, n_points=1000) {
  if (length(x) != length(y)) {
    stop('x and y must have the same length')
  }

  # Create a finer sequence for smoother CI plotting
  x1 <- seq(min(x), max(x), length.out=n_points)

  # Calculate Theil-Sen slope and intercept
  n <- length(x)
  slopes <- c()
  for (i in 1:(n-1)) {
    for (j in (i+1):n) {
      if (x[j] != x[i]) {
        slopes <- c(slopes, (y[j] - y[i]) / (x[j] - x[i]))
      }
    }
  }

  slope <- median(slopes, na.rm=TRUE)
  intercept <- median(y - slope * x)

  # Bootstrap to estimate confidence intervals
  set.seed(seed)  # For reproducibility
  boot_slopes <- numeric(n_bootstrap)
  boot_intercepts <- numeric(n_bootstrap)
  preds <- matrix(NA, nrow=n_bootstrap, ncol=length(x1))

  for (b in 1:n_bootstrap) {
    sample_indices <- sample(1:n, replace=TRUE)
    x_boot <- x[sample_indices]
    y_boot <- y[sample_indices]
    
    slopes_boot <- c()
    if (length(unique(x_boot)) > 1) {
      for (i in 1:(n-1)) {
        for (j in (i+1):n) {
          if (x_boot[j] != x_boot[i]) {
            slopes_boot <- c(slopes_boot, (y_boot[j] - y_boot[i]) / (x_boot[j] - x_boot[i]))
          }
        }
      }
      
      boot_slope <- median(slopes_boot, na.rm=TRUE)
      boot_intercept <- median(y_boot - boot_slope * x_boot)
      
      # Store bootstrap results
      boot_slopes[b] <- boot_slope
      boot_intercepts[b] <- boot_intercept
      
      # Predict y values for finer sequence x1
      preds[b, ] <- boot_intercept + boot_slope * x1
    }
  }

  # Remove NAs from bootstrap results
  boot_slopes <- boot_slopes[!is.na(boot_slopes)]
  boot_intercepts <- boot_intercepts[!is.na(boot_intercepts)]

  # Compute confidence intervals for predicted y values at x1
  ci <- apply(preds, 2, quantile, probs=c((1 - conf_level) / 2, 1 - (1 - conf_level) / 2), na.rm=TRUE)

  # Compute regression stats for slope and intercept
  slope_stats <- compute_regression_stats(slope, boot_slopes, conf_level)
  intercept_stats <- compute_regression_stats(intercept, boot_intercepts, conf_level)

  # Dynamic confidence interval labels based on conf_level
  lower_ci_label <- paste0(round((1 - conf_level) / 2 * 100, 1), '%')
  upper_ci_label <- paste0(round((1 + conf_level) / 2 * 100, 1), '%')

  # Format results, rounding to the specified decimal places and applying p-value formatting
  summary_table <- data.frame(
    est=round(c(intercept_stats$estimate, slope_stats$estimate), dp),
    `se`=round(c(intercept_stats$std_err, slope_stats$std_err), dp),
    `z`=round(c(intercept_stats$z_value, slope_stats$z_value), dp),
    `P(>|z|)`=c(format_p_value(intercept_stats$p_value, dp), format_p_value(slope_stats$p_value, dp)),
    check.names=FALSE
  )
  
  # Add dynamically labeled confidence intervals to the table
  summary_table[lower_ci_label] <- round(c(intercept_stats$lower_ci, slope_stats$lower_ci), dp)
  summary_table[upper_ci_label] <- round(c(intercept_stats$upper_ci, slope_stats$upper_ci), dp)
  
  # Set row names
  rownames(summary_table) <- c('(intercept)', 'slope')
  
  return(list(
    summary_table=summary_table,
    lower_ci=ci[1, ],
    upper_ci=ci[2, ],
    predicted_y=intercept + slope * x1,
    x1=x1  # Return x1 for plotting purposes
  ))
}

# Custom Siegel Estimator with confidence intervals and regression stats
siegel_sen_with_ci <- function(x, y, n_bootstrap=1e4, conf_level=0.95, seed=42, dp=5, n_points=1000) {
  if (length(x) != length(y)) {
    stop('x and y must have the same length')
  }

  # Create a finer sequence for smoother CI plotting
  x1 <- seq(min(x), max(x), length.out=n_points)

  n <- length(x)
  slopes_per_point <- numeric(n)
  
  # Original Siegel estimate
  for (i in 1:n) {
    slopes <- c()
    for (j in 1:n) {
      if (i != j && x[j] != x[i]) {
        slopes <- c(slopes, (y[j] - y[i]) / (x[j] - x[i]))
      }
    }
    slopes_per_point[i] <- median(slopes, na.rm=TRUE)
  }
  
  slope <- median(slopes_per_point, na.rm=TRUE)
  intercept <- median(y - slope * x)

  # Bootstrap to estimate confidence intervals
  set.seed(seed)  # For reproducibility
  boot_slopes <- numeric(n_bootstrap)
  boot_intercepts <- numeric(n_bootstrap)
  preds <- matrix(NA, nrow=n_bootstrap, ncol=length(x1))

  for (b in 1:n_bootstrap) {
    sample_indices <- sample(1:n, replace=TRUE)
    x_boot <- x[sample_indices]
    y_boot <- y[sample_indices]
    
    slopes_per_point <- numeric(n)
    if (length(unique(x_boot)) > 1) {
      for (i in 1:n) {
        slopes <- c()
        for (j in 1:n) {
          if (i != j && x_boot[j] != x_boot[i]) {
            slopes <- c(slopes, (y_boot[j] - y_boot[i]) / (x_boot[j] - x_boot[i]))
          }
        }
        slopes_per_point[i] <- median(slopes, na.rm=TRUE)
      }
      
      boot_slope <- median(slopes_per_point, na.rm=TRUE)
      boot_intercept <- median(y_boot - boot_slope * x_boot)
      
      # Store bootstrap results
      boot_slopes[b] <- boot_slope
      boot_intercepts[b] <- boot_intercept
      
      # Predict y values for finer sequence x1
      preds[b, ] <- boot_intercept + boot_slope * x1
    }
  }

  # Remove NAs from bootstrap results
  boot_slopes <- boot_slopes[!is.na(boot_slopes)]
  boot_intercepts <- boot_intercepts[!is.na(boot_intercepts)]

  # Compute confidence intervals for predicted y values at x1
  ci <- apply(preds, 2, quantile, probs=c((1 - conf_level) / 2, 1 - (1 - conf_level) / 2), na.rm=TRUE)

  # Compute regression stats for slope and intercept
  slope_stats <- compute_regression_stats(slope, boot_slopes, conf_level)
  intercept_stats <- compute_regression_stats(intercept, boot_intercepts, conf_level)

  # Dynamic confidence interval labels based on conf_level
  lower_ci_label <- paste0(round((1 - conf_level) / 2 * 100, 1), '%')
  upper_ci_label <- paste0(round((1 + conf_level) / 2 * 100, 1), '%')

  # Format results, rounding to the specified decimal places and applying p-value formatting
  summary_table <- data.frame(
    est=round(c(intercept_stats$estimate, slope_stats$estimate), dp),
    `se`=round(c(intercept_stats$std_err, slope_stats$std_err), dp),
    `z`=round(c(intercept_stats$z_value, slope_stats$z_value), dp),
    `P(>|z|)`=c(format_p_value(intercept_stats$p_value, dp), format_p_value(slope_stats$p_value, dp)),
    check.names=FALSE
  )
  
  # Add dynamically labeled confidence intervals to the table
  summary_table[lower_ci_label] <- round(c(intercept_stats$lower_ci, slope_stats$lower_ci), dp)
  summary_table[upper_ci_label] <- round(c(intercept_stats$upper_ci, slope_stats$upper_ci), dp)
  
  # Set row names
  rownames(summary_table) <- c('(intercept)', 'slope')
  
  # Return the summary table and the confidence intervals for predictions
  return(list(
    summary_table=summary_table,
    lower_ci=ci[1, ],
    upper_ci=ci[2, ],
    predicted_y=intercept + slope * x1,
    x1=x1  # Return x1 for plotting purposes
  ))
}

convert_A_to_scatter <- function(A, sign=1) {
  # Get the unique 's' values
  s_vals <- unique(A$s)
  
  # Initialize empty vectors for x and y
  x_vals <- numeric(length(s_vals))
  y_vals <- numeric(length(s_vals))
  
  # Loop through each 's' and assign the correct values to x and y
  for (i in seq_along(s_vals)) {
    x_vals[i] <- A$y[A$s == s_vals[i] & A$x == 1]
    y_vals[i] <- A$y[A$s == s_vals[i] & A$x == 2]
  }
  
  # Create the scatter data frame
  scatter <- data.frame(s=s_vals, x=sign*x_vals, y=sign*y_vals)
  
  return(scatter)
}

convert_to_scatter <- function(A, sign=1) {
  # Split data frame A by 's' and 'x' values to get unique levels
  split_data <- split(A, A$s)
  
  # Create scatter data frame
  scatter <- do.call(rbind, lapply(split_data, function(df) {
    n <- nrow(df) / 2  # Calculate the number of levels
    data.frame(
      s = df$s[1:n],
      level = df$x[1:n],
      x = sign * (df$y[1:n]),
      y = sign *(df$y[(n + 1):(2 * n)])
    )
  }))
  
  return(scatter)
}


# ScatterPlot <- function(A, sign=1, xlim=c(0, 400), ylim=c(0, 400), x_tick_interval=100, y_tick_interval=100, tick_length=0.2, 
#   height=4, width=4, xlab='', ylab='', main='', colors=c('black', 'indianred'), open_symbols=FALSE, lwd=1, p.cex=0.5, 
#   filename='scatter.svg', reg=FALSE, plot.CI=FALSE, reg.points=1e3, reg.color='darkgray', reg.method=c('Siegel', 'Theil-Sen'), reg.CI.settings=list(nboot=1e4, conf_level=0.95), 
#   save=FALSE, bg='transparent', dp=3, return.output=FALSE) {

#   # scatter <- convert_A_to_scatter(A=A, sign=sign)
#   scatter <- convert_to_scatter(A=A, sign=sign)

#   # Create the scatter plot
#   if (save) {
#     svg(file=filename, width=width, height=height, bg=bg)
#   } else {
#     dev.new(width=width, height=height, noRStudioGD=TRUE)
#   }

#   # Check if 'level' column exists and map levels to 1 and 2 alternately, if not set n=1
#   if ('level' %in% colnames(scatter)) {
#     unique_levels <- unique(scatter$level)
#     n <- length(unique_levels)
#     scatter$level <- as.numeric(factor(scatter$level, levels=unique_levels, labels=rep(1:n, length.out=length(unique_levels))))
#   } else {
#     n <- 1  # If 'level' column does not exist
#     scatter$level <- rep(1, dim(scatter)[1])
#   }

#   # Determine plot symbols (only open circles if open_symbols is TRUE)
#   pch <- if (open_symbols) 1 else 19

#   cols <- hex_palette(n=n, color1=colors[1], color2=colors[2], reverse=FALSE)

#   plot(scatter$x, scatter$y, col=cols[scatter$level], pch=pch, cex=p.cex, xlim=xlim, ylim=ylim, 
#        xlab=xlab, ylab=ylab, main=main, xaxt='n', yaxt='n', bty='n')

#   # Add regression line or non-parametric line based on 'reg' parameter
#   if (!reg) {
#     segments(min(xlim), min(ylim), max(xlim), max(ylim), lwd=lwd, col=reg.color, lty=3)  # lty=3 for dotted line
#   } else {
#     reg.method <- match.arg(reg.method)

#     if (reg.method == 'Siegel') {
#       reg_func <- if (plot.CI || return.output) siegel_sen_with_ci else siegel_sen
#     } else if (reg.method == 'Theil-Sen') {
#       reg_func <- if (plot.CI || return.output) theil_sen_with_ci else theil_sen
#     }

#     # Perform regression and extract results
#     if (plot.CI || return.output) {
#       reg_results <- reg_func(scatter$x, scatter$y, n_bootstrap=reg.CI.settings$nboot, 
#                               conf_level=reg.CI.settings$conf_level, dp=dp, n_points=reg.points)
      
#       summary_table <- reg_results$summary_table
#       intercept <- summary_table["(intercept)", "est"]
#       slope <- summary_table["slope", "est"]
      
#       # Confidence intervals for predictions
#       lower_ci <- reg_results$lower_ci
#       upper_ci <- reg_results$upper_ci
#       predicted_y <- reg_results$predicted_y
#       x1 <- reg_results$x1
#     } else {
#       out <- reg_func(scatter$x, scatter$y)
#       intercept <- out$intercept
#       slope <- out$slope
#       x1 <- seq(min(scatter$x), max(scatter$x), length.out=1000)  # Fallback for plotting
#       predicted_y <- intercept + slope * x1
#     }

#     # Plot the regression line using x1
#     lines(x1, predicted_y, lwd=lwd, col=reg.color, lty=3)

#     # Plot confidence intervals if requested
#     if (plot.CI) {
#       lines(x1, lower_ci, col=reg.color, lty=2)
#       lines(x1, upper_ci, col=reg.color, lty=2)
#     }
#   }

#   # Define tick intervals and lengths
#   x_ticks <- seq(min(xlim), max(xlim), by=x_tick_interval)
#   y_ticks <- seq(min(ylim), max(ylim), by=y_tick_interval)
  
#   # Customize x-axis
#   axis(1, at=x_ticks, labels=x_ticks, tcl=-tick_length, lwd=lwd)

#   # Customize y-axis with horizontal labels
#   axis(2, at=y_ticks, labels=y_ticks, tcl=-tick_length, las=1, lwd=lwd)

#   if (save) {
#     dev.off()
#   }

#   if (reg && return.output) {
#     return(summary_table)
#   }
# }


# ScatterPlot <- function(A, sign=1, xlim=c(0, 400), ylim=c(0, 400), x_tick_interval=100, y_tick_interval=100, tick_length=0.2, 
#                         height=4, width=4, xlab='', ylab='', main='', colors=c('black', 'indianred'), open_symbols=FALSE, 
#                         lwd=1, p.cex=0.5, filename='scatter.svg', reg=FALSE, plot.CI=FALSE, reg.points=1e3, reg.color='darkgray', 
#                         reg.method=c('Siegel', 'Theil-Sen'), reg.CI.settings=list(nboot=1e4, conf_level=0.95), save=FALSE, 
#                         bg='transparent', dp=3, return.output=FALSE) {
  
#   # Convert A to scatter
#   scatter <- convert_to_scatter(A=A, sign=sign)

#   # Create the scatter plot
#   if (save) {
#     svg(file=filename, width=width, height=height, bg=bg)
#   } else {
#     dev.new(width=width, height=height, noRStudioGD=TRUE)
#   }

#   # Check if 'level' column exists and map levels to 1 and 2 alternately, if not set n=1
#   if ('level' %in% colnames(scatter)) {
#     unique_levels <- unique(scatter$level)
#     n <- length(unique_levels)
#     scatter$level <- as.numeric(factor(scatter$level, levels=unique_levels, labels=rep(1:n, length.out=length(unique_levels))))
#   } else {
#     n <- 1  # If 'level' column does not exist
#     scatter$level <- rep(1, dim(scatter)[1])
#   }

#   # Determine plot symbols
#   pch <- if (open_symbols) 1 else 19
#   cols <- hex_palette(n=2, color1=colors[1], color2=colors[2], reverse=FALSE)
  
#   # Plot scatter points
#   plot(scatter$x, scatter$y, col=cols[scatter$level], pch=pch, cex=p.cex, xlim=xlim, ylim=ylim, 
#        xlab=xlab, ylab=ylab, main=main, xaxt='n', yaxt='n', bty='n')

#   # Initialize list to store summary tables for each level
#   summary_tables <- list()

#   # Add regression line or non-parametric line based on 'reg' parameter
#   if (!reg) {
#     segments(min(xlim), min(ylim), max(xlim), max(ylim), lwd=lwd, col=reg.color, lty=3)  # lty=3 for dotted line
#   } else {
#     reg.method <- match.arg(reg.method)
#     levels <- unique(scatter$level)

#     reg_func <- switch(reg.method,
#                        'Siegel' = if (plot.CI || return.output) siegel_sen_with_ci else siegel_sen,
#                        'Theil-Sen' = if (plot.CI || return.output) theil_sen_with_ci else theil_sen)

#     for (level in levels) {
#       # Filter data for the current level
#       level_data <- scatter[scatter$level==level,]
      
#       if (plot.CI || return.output) {
#         # Perform regression with confidence intervals for the current level
#         reg_results <- reg_func(level_data$x, level_data$y, n_bootstrap=reg.CI.settings$nboot, 
#                                 conf_level=reg.CI.settings$conf_level, dp=dp, n_points=reg.points)
        
#         # Store the summary table for this level
#         summary_tables[[as.character(level)]] <- reg_results$summary_table
        
#         intercept <- reg_results$summary_table["(intercept)", "est"]
#         slope <- reg_results$summary_table["slope", "est"]
        
#         # Confidence intervals for predictions
#         lower_ci <- reg_results$lower_ci
#         upper_ci <- reg_results$upper_ci
#         predicted_y <- reg_results$predicted_y
#         x1 <- reg_results$x1
#       } else {
#         # Perform standard regression without confidence intervals for the current level
#         out <- reg_func(level_data$x, level_data$y)
#         intercept <- out$intercept
#         slope <- out$slope
#         x1 <- seq(min(level_data$x), max(level_data$x), length.out=reg.points)
#         predicted_y <- intercept + slope * x1
#       }
      
#       # Plot the regression line for the current level
#       lines(x1, predicted_y, lwd=lwd, col=cols[as.integer(level)], lty=3)
      
#       # Plot confidence intervals if requested
#       if (plot.CI) {
#         lines(x1, lower_ci, col=cols[as.integer(level)], lty=2)
#         lines(x1, upper_ci, col=cols[as.integer(level)], lty=2)
#       }
#     }
#   }

#   # Customize axes
#   axis(1, at=seq(min(xlim), max(xlim), by=x_tick_interval), tcl=-tick_length, lwd=lwd)
#   axis(2, at=seq(min(ylim), max(ylim), by=y_tick_interval), las=1, tcl=-tick_length, lwd=lwd)

#   if (save) {
#     dev.off()
#   }

#   # Return the list of summary tables if requested
#   if (reg && return.output) {
#     return(summary_tables)
#   }
# }


ScatterPlot <- function(A, sign=1, xlim=c(0, 400), ylim=c(0, 400), x_tick_interval=100, y_tick_interval=100, tick_length=0.2, 
                        height=4, width=4, xlab='', ylab='', main='', colors=c('black', 'indianred'), open_symbols=FALSE, 
                        lwd=1, p.cex=0.5, filename='scatter.svg', reg=FALSE, plot.CI=FALSE, reg.points=1e3, reg.color='darkgray', 
                        reg.method=c('Siegel', 'Theil-Sen'), reg.CI.settings=list(nboot=1e4, conf_level=0.95), save=FALSE, 
                        show_pairs=FALSE, bg='transparent', dp=3, return.output=FALSE) {
  
  # Convert A to scatter
  scatter <- convert_to_scatter(A=A, sign=sign)

  # Create the scatter plot
  if (save) {
    svg(file=filename, width=width, height=height, bg=bg)
  } else {
    dev.new(width=width, height=height, noRStudioGD=TRUE)
  }

  # Check if 'level' column exists and map levels to 1 and 2 alternately, if not set n=1
  if ('level' %in% colnames(scatter)) {
    unique_levels <- unique(scatter$level)
    n <- length(unique_levels)
    scatter$level <- as.numeric(factor(scatter$level, levels=unique_levels, labels=rep(1:n, length.out=length(unique_levels))))
  } else {
    n <- 1  # If 'level' column does not exist
    scatter$level <- rep(1, dim(scatter)[1])
  }

  # Determine plot symbols
  pch <- if (open_symbols) 1 else 19
  cols <- hex_palette(n=2, color1=colors[1], color2=colors[2], reverse=FALSE)
  
  # Plot scatter points
  plot(scatter$x, scatter$y, col=cols[scatter$level], pch=pch, cex=p.cex, xlim=xlim, ylim=ylim, 
       xlab=xlab, ylab=ylab, main=main, xaxt='n', yaxt='n', bty='n')

  # Initialize list to store summary tables for each level
  summary_tables <- list()

  # Add regression line or non-parametric line based on 'reg' parameter
  if (!reg) {
    segments(min(xlim), min(ylim), max(xlim), max(ylim), lwd=lwd, col=reg.color, lty=3)  # lty=3 for dotted line
  } else {
    reg.method <- match.arg(reg.method)
    levels <- unique(scatter$level)

    reg_func <- switch(reg.method,
                       'Siegel' = if (plot.CI || return.output) siegel_sen_with_ci else siegel_sen,
                       'Theil-Sen' = if (plot.CI || return.output) theil_sen_with_ci else theil_sen)

    for (level in levels) {
      # Filter data for the current level
      level_data <- scatter[scatter$level==level,]
      
      if (plot.CI || return.output) {
        # Perform regression with confidence intervals for the current level
        reg_results <- reg_func(level_data$x, level_data$y, n_bootstrap=reg.CI.settings$nboot, 
                                conf_level=reg.CI.settings$conf_level, dp=dp, n_points=reg.points)
        
        # Store the summary table for this level
        summary_tables[[as.character(level)]] <- reg_results$summary_table
        
        intercept <- reg_results$summary_table["(intercept)", "est"]
        slope <- reg_results$summary_table["slope", "est"]
        
        # Confidence intervals for predictions
        lower_ci <- reg_results$lower_ci
        upper_ci <- reg_results$upper_ci
        predicted_y <- reg_results$predicted_y
        x1 <- reg_results$x1
      } else {
        # Perform standard regression without confidence intervals for the current level
        out <- reg_func(level_data$x, level_data$y)
        intercept <- out$intercept
        slope <- out$slope
        x1 <- seq(min(level_data$x), max(level_data$x), length.out=reg.points)
        predicted_y <- intercept + slope * x1
      }
      
      # Plot the regression line for the current level
      lines(x1, predicted_y, lwd=lwd, col=cols[as.integer(level)], lty=3)
      
      # Plot confidence intervals if requested
      if (plot.CI) {
        lines(x1, lower_ci, col=cols[as.integer(level)], lty=2)
        lines(x1, upper_ci, col=cols[as.integer(level)], lty=2)
      }
    }
  }

  if (show_pairs && "s" %in% colnames(scatter)) {
    subjects <- unique(scatter$s)
    for (subj in subjects) {
      subj_data <- scatter[scatter$s == subj, ]
      # Optionally, order the data by x or some grouping variable:
      # subj_data <- subj_data[order(subj_data$x), ]
      # Connect them with lines
      if (nrow(subj_data) > 1) {
        lines(subj_data$x, subj_data$y, col='darkgray', lty=3, lwd=lwd)
      }
    }
  }

  # Customize axes
  axis(1, at=seq(min(xlim), max(xlim), by=x_tick_interval), tcl=-tick_length, lwd=lwd)
  axis(2, at=seq(min(ylim), max(ylim), by=y_tick_interval), las=1, tcl=-tick_length, lwd=lwd)

  if (save) {
    dev.off()
  }

  # Return the list of summary tables if requested
  if (reg && return.output) {
    return(summary_tables)
  }
}
# Define the start and end colors of your palette Slate Blue to Indian Red
hex_palette <- function(n, color1='#6A5ACD', color2='#CD5C5C', reverse = FALSE) {
  
  # Create a sequence of colors
  colors <- colorRampPalette(c(color1, color2))(n)
  
  # reverse colors if reverse=TRUE
  if (reverse) {
    colors <- rev(colors)
  }
  
  return(colors)
}

traces2plot <- function(V, traces, offsets=NULL, color1='#6A5ACD', color2='#CD5C5C', xlim=NULL, ylim=NULL, lwd=1, 
  xbar=100, ybar=50, reverse=TRUE, show_text=FALSE, normalise=FALSE, height=3, width=3, 
  filename='traces.svg', dt=0.1, save=FALSE){
  
  if (save) {
    # Open SVG device
    svg(file=filename, width=width, height=height, bg='transparent')
  } else {
    dev.new(width=width, height=height, noRStudioGD=TRUE)
  }

  n <- length(traces)
  colors <- hex_palette(n=n, color1=color1, color2=color2, reverse=reverse)
  if (reverse) traces = rev(traces)
  if (is.null(offsets)) offsets <- rep(0, n)

  x = 0:dt:(dim(V)[1]-1)*dt

  if (is.null(xlim)) xlim <- c(min(x), max(x))
  if (is.null(ylim)) ylim <- c(0, max(apply(V[, traces],2,max)))

  idx1 <- which.min(abs(x - xlim[1]))
  idx2 <- which.min(abs(x - xlim[2]))

  y <- V[, traces[1]]
  
  # # Design a low-pass filter
  # fs <- 10  # Sampling frequency (e.g., 1 Hz if data points are 1 second apart)
  # cutoff <- 0.1  # Cutoff frequency (e.g., 0.1 Hz)
  # butter_order <- 5  # Order of the Butterworth filter
  # bf <- butter(butter_order, cutoff / (0.5 * fs), type = "low")

  # # Apply the low-pass filter
  # y_smooth <- filtfilt(bf, y)


  plot(x[idx1:idx2]+offsets[1], y[idx1:idx2], type='l', col=colors[1], xlim=xlim, ylim=ylim, bty='n', lwd=lwd, lty=1, axes=FALSE, frame=FALSE, xlab = '', ylab = '')

  # Loop through remaining traces and add them to the plot
  for (i in 2:n) {
    y <- V[, traces[i]]
    lines(x[idx1:idx2]+offsets[i], y[idx1:idx2], col=colors[i], lwd=lwd, lty=1)
  }

  # Define scale bar lengths and ybar position
  # ybar_start <- (max(ylim) - min(ylim)) / 10
  ybar_start <- min(ylim) + (max(ylim) - min(ylim)) / 20

  # Add scale bars at the bottom right
  x_start <- max(xlim) - xbar - 50
  y_start <- ybar_start
  x_end <- x_start + xbar
  y_end <- y_start + ybar

  # Draw the scale bars
  segments(x_start, y_start, x_end, y_start, lwd=lwd, col='black') # Horizontal scale bar
  if (!normalise){
    segments(x_start, y_start, x_start, y_end, lwd=lwd, col='black') # Vertical scale bar
  }

  # Add labels to the scale bars
  if (show_text){
    text(x = (x_start + x_end) / 2, y = y_start - ybar / 20, labels = paste(xbar, 'ms'), adj = c(0.5, 1))
    if (!normalise) text(x = x_start -xbar/4, y = (y_start + y_end) / 2, labels = paste(ybar, 'mV'), adj = c(0.5, 0.5), srt = 90)
  }

  if (save) {
    dev.off()
  }

}


fun_single_example <- function(rawdata, fits, start_time=50, baseline=50, idx=1, model='product2', response_sign_method = c('smooth', 'regression', 'cumsum')){
  
  ind1 <- (start_time - baseline)/dx
  ind2 <- start_time/dx

  y <- rawdata[, idx+1]
  y <- na.omit(y)
  y <- y[ind1:length(y)]
  y <- y - mean(y[ind1:ind2])
  x <- seq(0, (length(y) - 1) * dx, by = dx)

  # sign <- sign_fun(y, direction_method=response_sign_method) 
  # single_example <- data.frame('x'=x, 'y'= sign * y)
  single_example <- data.frame('x'=x, 'y'= y)

  if (model == 'product') {
    fun.2.fit <- product1
  } else if (model == 'product2') {
    fun.2.fit <- product2
  }

  pars <- fits[[idx]]$fits
  if (model == 'product') {
    pars[4] <- pars[4] + baseline
  } else if (model == 'product2') {
    pars[4] <- pars[4] + baseline
    pars[8] <- pars[8] + baseline
  }    
  yfit <- fun.2.fit(pars, x) 
  single_example$y_fit <- yfit

  if (model == 'product2'){
    yfit1 <- product1(pars[1:4], x) 
    yfit2 <- product1(pars[5:8], x) 
    single_example$y_fit1 <- yfit1
    single_example$y_fit2 <- yfit2
  }
  return(single_example)
}

# single_fit_egs <- function(traces, xlim=NULL, ylim=NULL, lwd=1, show_text=FALSE, normalise=FALSE, func=product2N, height=4, width=2.5, xbar=100, ybar=50, log_y=FALSE, colors=c('#4C77BB', '#CA92C1', '#F28E2B'), filename='plot.svg', bg='transparent', save=FALSE) {
  
#   if (save) {
#     # Open SVG device
#     svg(file=filename, width=width, height=height, bg=bg)
#   } else {
#     dev.new(width=width, height=height, noRStudioGD=TRUE)
#   }
  
#   x <- traces$x
#   y <- traces$y
  
#   fit1 <- traces$yfit1
#   if (identical(func, product2N)){
#     fit2 <- traces$yfit2
#   }
  
#   if (identical(func, product3N)){
#     fit2 <- traces$yfit2
#     fit3 <- traces$yfit3
#   }

#   if (is.null(xlim)) xlim <- c(min(x), max(x))
  
#   if (is.null(ylim)) {
#     if (log_y) {
#       y <- -y
#       fit1 <- -fit1
#       if (identical(func, product2N)){
#         fit2 <- -fit2
#       }
#       if (identical(func, product3N)){
#         fit2 <- -fit2
#         fit3 <- -fit3
#       }

#       # Define custom major tick positions
#       y_ticks <- c(1, 10, 100, 1000)  # Example custom major ticks
#       log_y_ticks <- log(y_ticks)
#       valid_ticks <- log_y_ticks[log_y_ticks <= log(max(y[y > 0], na.rm=TRUE))]  # Get ticks up to the maximum y

#       # Set ylim based on valid ticks, using the lowest major tick for the lower bound
#       ylim <- c(min(valid_ticks), log(max(y[y > 0], na.rm=TRUE)))
#     } else {
#       ylim <- c(-max(y, na.rm=TRUE), 0)
#     }
#   }

#   if (log_y) {
#     y <- ifelse(y > 0, log(pmax(y, .Machine$double.eps)), NA)
#     fit1 <- ifelse(fit1 > 0, log(fit1), NA)
#     if (identical(func, product2N)){
#       fit2 <- ifelse(fit2 > 0, log(fit2), NA)
#     }
#     if (identical(func, product3N)){
#       fit2 <- ifelse(fit2 > 0, log(fit2), NA)
#       fit3 <- ifelse(fit3 > 0, log(fit3), NA)
#     }
#   } 

#   idx1 <- which.min(abs(x - xlim[1]))
#   idx2 <- which.min(abs(x - xlim[2]))

#   plot(x[idx1:idx2], y[idx1:idx2], type='l', col='#A6A8AA', xlim=xlim, ylim=ylim, bty='n', lwd=lwd, lty=1, axes=FALSE, frame=FALSE, xlab='', ylab='')

#   if (identical(func, product1N)){
#     fits <- cbind(fit)
#   }else if (identical(func, product2N)){
#     fits <- cbind(fit1, fit2)
#   }else if (identical(func, product3N)){
#     fits <- cbind(fit1, fit2, fit3)
#   }

  
#   # Loop through remaining traces and add them to the plot
#   for (i in 1:dim(fits)[2]) {
#     y_fit <- fits[, i]
#     lines(x[idx1:idx2], y_fit[idx1:idx2], col=colors[i], lwd=lwd, lty=1)
#   }
  
#   # Define scale bar lengths and ybar position
#   ybar <- ifelse(log_y, exp(1), ybar)
#   ybar_start <- ifelse(log_y, log(1) + (log(max(exp(ylim))) - log(1)) / 20, min(ylim) + (max(ylim) - min(ylim)) / 20)
  
#   # Add scale bars at the bottom right
#   x_start <- max(xlim) - xbar - 50
#   y_start <- ybar_start
#   x_end <- x_start + xbar
#   y_end <- ifelse(log_y, y_start + log(ybar), y_start + ybar)
  
#   # Draw the scale bars
#   segments(x_start, y_start, x_end, y_start, lwd=lwd, col='black') # Horizontal scale bar
#   if (!normalise) {
#     segments(x_start, y_start, x_start, y_end, lwd=lwd, col='black') # Vertical scale bar
#   }
  
#   # Add labels to the scale bars
#   if (show_text) {
#     text(x = (x_start + x_end) / 2, y = y_start - ybar / 20, labels = paste(xbar, 'ms'), adj = c(0.5, 1))
#     if (!normalise) text(x = x_start - xbar / 4, y = (y_start + y_end) / 2, labels = ifelse(log_y, "e-fold change", paste(ybar, 'pA')), adj = c(0.5, 0.5), srt = 90)
#   }
  
#   # Add the y-axis only if log_y is TRUE
#   if (log_y) {
#     tick_length <- -0.2
#     minor_tick_length <- -0.1
    
#     # Major tick positions and labels for the log scale
#     y_ticks <- c(1, 10, 100, 1000)  # Example custom major ticks
#     log_y_ticks <- log(y_ticks)
#     valid_ticks <- log_y_ticks[log_y_ticks >= ylim[1] & log_y_ticks <= ylim[2]]  # Filter major ticks within the plot range

#     # Minor tick positions for the log scale
#     minor_y_ticks <- c(2, 3, 4, 5, 6, 7, 8, 9, 
#                        20, 30, 40, 50, 60, 70, 80, 90, 
#                        200, 300, 400, 500, 600, 700, 800, 900)  # Example custom minor ticks
#     log_minor_y_ticks <- log(minor_y_ticks)
#     valid_minor_ticks <- log_minor_y_ticks[log_minor_y_ticks >= ylim[1] & log_minor_y_ticks <= ylim[2]]  # Filter minor ticks within the plot range

#     # Add major ticks
#     axis(2, at=valid_ticks, labels=y_ticks[log_y_ticks >= ylim[1] & log_y_ticks <= ylim[2]], tcl=tick_length, las=1)

#     # Add minor ticks
#     axis(2, at=valid_minor_ticks, labels=NA, tcl=minor_tick_length, las=1)
    
#     # Add y-axis label
#     mtext('PSC amplitude (pA)', side=2, line=2.5)
#   }
  
#   if (save) {
#     dev.off()
#   }
# }

single_fit_egs <- function(traces, sign=-1, xlim=NULL, ylim=NULL, lwd=1, show_text=FALSE, normalise=FALSE, func=product2N, height=4, width=2.5, 
  xbar=100, ybar=50, ybar_units = 'pA', log_y=FALSE, colors=c('#4C77BB', '#CA92C1', '#F28E2B'), filename='plot.svg', bg='transparent', save=FALSE) {
  
  if (save) {
    # Open SVG device
    svg(file=filename, width=width, height=height, bg=bg)
  } else {
    dev.new(width=width, height=height, noRStudioGD=TRUE)
  }

  x <- traces$x
  y <- traces$y
  
  fit1 <- if (identical(func, product1N)) traces$yfit else traces$yfit1
  
  if (identical(func, product2N)){
    fit2 <- traces$yfit2
  }
  
  if (identical(func, product3N)){
    fit2 <- traces$yfit2
    fit3 <- traces$yfit3
  }

  if (is.null(xlim)) xlim <- c(min(x), max(x))
  
  if (is.null(ylim)) {
    if (log_y) {
      y <- sign * y
      fit1 <- sign * fit1
      if (identical(func, product2N)){
        fit2 <- sign * fit2
      }
      if (identical(func, product3N)){
        fit2 <- sign * fit2
        fit3 <- sign * fit3
      }

      # Define custom major tick positions
      y_ticks <- c(1, 10, 100, 1000)  # Example custom major ticks
      log_y_ticks <- log(y_ticks)
      valid_ticks <- log_y_ticks[log_y_ticks <= log(max(y[y > 0], na.rm=TRUE))]  # Get ticks up to the maximum y

      # Set ylim based on valid ticks, using the lowest major tick for the lower bound
      ylim <- c(min(valid_ticks), log(max(y[y > 0], na.rm=TRUE)))
    } else {
      ylim <- c(-max(y, na.rm=TRUE), 0)
    }
  }

  if (log_y) {
    y <- ifelse(y > 0, log(pmax(y, .Machine$double.eps)), NA)
    fit1 <- ifelse(fit1 > 0, log(pmax(fit1, .Machine$double.eps)), NA)
    if (identical(func, product2N)){
      fit2 <- ifelse(fit2 > 0, log(pmax(fit2, .Machine$double.eps)), NA)
    }
    if (identical(func, product3N)){
      fit2 <- ifelse(fit2 > 0, log(pmax(fit2, .Machine$double.eps)), NA)
      fit3 <- ifelse(fit3 > 0, log(pmax(fit3, .Machine$double.eps)), NA)
    }
    y[ y < 0 ] <- NA
    y[which.min(abs(x - xlim[1]))] <- 0
  } 

  idx1 <- which.min(abs(x - xlim[1]))
  idx2 <- which.min(abs(x - xlim[2]))

  plot(x[idx1:idx2], y[idx1:idx2], type='l', col='#A6A8AA', xlim=xlim, ylim=ylim, bty='n', lwd=lwd, lty=1, axes=FALSE, frame=FALSE, xlab='', ylab='')

  if (identical(func, product1N)){
    fits <- cbind(fit1)
  }else if (identical(func, product2N)){
    fits <- cbind(fit1, fit2)
  }else if (identical(func, product3N)){
    fits <- cbind(fit1, fit2, fit3)
  }

  if (log_y) {
    fits[ fits < 0 ] <- NA
    # Step 2: For each column, set the last NA before the first numeric to 1, if it exists.
    for (col_idx in seq_len(ncol(fits))) {
      col_data <- fits[, col_idx]
      
      # Find the first non-NA value's index
      f <- which(!is.na(col_data))[1]
      
      # If there's a first numeric value and it's not in row 1...
      if (!is.na(f) && f > 0) {
        # Check if row f-1 is NA
        if (is.na(col_data[f - 1])) {
          col_data[f - 1] <- 0
        }
      }
      
      # Write the updated column back
      fits[, col_idx] <- col_data
    }
  }

  # Loop through remaining traces and add them to the plot
  for (i in 1:dim(fits)[2]) {
    y_fit <- fits[, i]
    lines(x[idx1:idx2], y_fit[idx1:idx2], col=colors[i], lwd=lwd, lty=1)
  }
  
  # Define scale bar lengths and ybar position
  ybar <- ifelse(log_y, exp(1), ybar)
  ybar_start <- ifelse(log_y, log(1) + (log(max(exp(ylim))) - log(1)) / 20, min(ylim) + (max(ylim) - min(ylim)) / 20)
  
  # Add scale bars at the bottom right
  x_start <- max(xlim) - xbar - 50
  y_start <- ybar_start
  x_end <- x_start + xbar
  y_end <- ifelse(log_y, y_start + log(ybar), y_start + ybar)
  
  # Draw the scale bars
  segments(x_start, y_start, x_end, y_start, lwd=lwd, col='black') # Horizontal scale bar
  if (!normalise) {
    segments(x_start, y_start, x_start, y_end, lwd=lwd, col='black') # Vertical scale bar
  }
  
  # Add labels to the scale bars
  if (show_text) {
    text(x = (x_start + x_end) / 2, y = y_start - ybar / 20, labels = paste(xbar, 'ms'), adj = c(0.5, 1))
    if (!normalise) text(x = x_start - xbar / 4, y = (y_start + y_end) / 2, labels = ifelse(log_y, "e-fold change", paste(ybar, ybar_units)), adj = c(0.5, 0.5), srt = 90)
  }
  
  # Add the y-axis only if log_y is TRUE
  if (log_y) {
    tick_length <- -0.2
    minor_tick_length <- -0.1
    
    # Major tick positions and labels for the log scale
    y_ticks <- c(1, 10, 100, 1000)  # Example custom major ticks
    log_y_ticks <- log(y_ticks)
    valid_ticks <- log_y_ticks[log_y_ticks >= ylim[1] & log_y_ticks <= ylim[2]]  # Filter major ticks within the plot range

    # Minor tick positions for the log scale
    minor_y_ticks <- c(2, 3, 4, 5, 6, 7, 8, 9, 
                       20, 30, 40, 50, 60, 70, 80, 90, 
                       200, 300, 400, 500, 600, 700, 800, 900)  # Example custom minor ticks
    log_minor_y_ticks <- log(minor_y_ticks)
    valid_minor_ticks <- log_minor_y_ticks[log_minor_y_ticks >= ylim[1] & log_minor_y_ticks <= ylim[2]]  # Filter minor ticks within the plot range

    # Add major ticks
    axis(2, at=valid_ticks, labels=y_ticks[log_y_ticks >= ylim[1] & log_y_ticks <= ylim[2]], tcl=tick_length, las=1)

    # Add minor ticks
    axis(2, at=valid_minor_ticks, labels=NA, tcl=minor_tick_length, las=1)
    
    # Add y-axis label
    if (grepl("A$", ybar_units)) {
      mtext(bquote(PSC~amplitude~"(" * plain(.(ybar_units)) * ")"), side = 2, line = 2.5)
    } else if (grepl("V$", ybar_units)) {
      mtext(bquote(PSP~amplitude~"(" * plain(.(ybar_units)) * ")"), side = 2, line = 2.5)
    }

  }
  
  if (save) {
    dev.off()
  }
}

SingleFitExample <- function(traces, xlim=NULL, ylim=NULL, ylab='PSC amplitude (pA)', tick_length=0.2, lwd=1, show_text=FALSE, normalise=FALSE, func=product2N, 
  height=4, width=2.5, xbar=100, ybar=50, log_y=FALSE, colors=c('#4C77BB', '#CA92C1', '#F28E2B'), filename='plot.svg', bg='transparent', save=FALSE) {
  
  if (save) {
    # Open SVG device
    svg(file=filename, width=width, height=height, bg=bg)
  } else {
    dev.new(width=width, height=height, noRStudioGD=TRUE)
  }
  
  x <- traces$x
  y <- traces$y
  
  fit1 <- traces$yfit1
  if (identical(func, product2N)){
    fit2 <- traces$yfit2
  }
  
  if (identical(func, product3N)){
    fit2 <- traces$yfit2
    fit3 <- traces$yfit3
  }

  if (is.null(xlim)) xlim <- c(min(x), max(x))
  
  if (is.null(ylim)) {
    if (log_y) {
      y <- -y
      fit1 <- -fit1
      if (identical(func, product2N)){
        fit2 <- -fit2
      }
      if (identical(func, product3N)){
        fit2 <- -fit2
        fit3 <- -fit3
      }

      # Define custom major tick positions
      y_ticks <- c(1, 10, 100, 1000)  # Example custom major ticks
      log_y_ticks <- log(y_ticks)
      valid_ticks <- log_y_ticks[log_y_ticks <= log(max(y[y > 0], na.rm=TRUE))]  # Get ticks up to the maximum y

      # Set ylim based on valid ticks, using the lowest major tick for the lower bound
      ylim <- c(min(valid_ticks), log(max(y[y > 0], na.rm=TRUE)))
    } else {
      ylim <- c(-max(y, na.rm=TRUE), 0)
    }
  }

  if (log_y) {
    y <- ifelse(y > 0, log(pmax(y, .Machine$double.eps)), NA)
    fit1 <- ifelse(fit1 > 0, log(fit1), NA)
    if (identical(func, product2N)){
      fit2 <- ifelse(fit2 > 0, log(fit2), NA)
    }
    if (identical(func, product3N)){
      fit2 <- ifelse(fit2 > 0, log(fit2), NA)
      fit3 <- ifelse(fit3 > 0, log(fit3), NA)
    }
  } 

  idx1 <- which.min(abs(x - xlim[1]))
  idx2 <- which.min(abs(x - xlim[2]))

  plot(x[idx1:idx2], y[idx1:idx2], type='l', col='#A6A8AA', xlim=xlim, ylim=ylim, bty='n', lwd=lwd, lty=1, axes=FALSE, frame=FALSE, xlab='', ylab='')

  if (identical(func, product1N)){
    fits <- cbind(fit)
  }else if (identical(func, product2N)){
    fits <- cbind(fit1, fit2)
  }else if (identical(func, product3N)){
    fits <- cbind(fit1, fit2, fit3)
  }

  
  # Loop through remaining traces and add them to the plot
  for (i in 1:dim(fits)[2]) {
    y_fit <- fits[, i]
    lines(x[idx1:idx2], y_fit[idx1:idx2], col=colors[i], lwd=lwd, lty=1)
  }
  
  # Define scale bar lengths and ybar position
  ybar <- ifelse(log_y, exp(1), ybar)
  ybar_start <- ifelse(log_y, log(1) + (log(max(exp(ylim))) - log(1)) / 20, min(ylim) + (max(ylim) - min(ylim)) / 20)
  
  # Add scale bars at the bottom right
  x_start <- max(xlim) - xbar - 50
  y_start <- ybar_start
  x_end <- x_start + xbar
  y_end <- ifelse(log_y, y_start + log(ybar), y_start + ybar)
  
  # Draw the scale bars
  segments(x_start, y_start, x_end, y_start, lwd=lwd, col='black') # Horizontal scale bar
  if (!normalise) {
    segments(x_start, y_start, x_start, y_end, lwd=lwd, col='black') # Vertical scale bar
  }
  
  # Add labels to the scale bars
  if (show_text) {
    text(x = (x_start + x_end) / 2, y = y_start - ybar / 20, labels = paste(xbar, 'ms'), adj = c(0.5, 1))
    if (!normalise) text(x = x_start - xbar / 4, y = (y_start + y_end) / 2, labels = ifelse(log_y, "e-fold change", paste(ybar, 'pA')), adj = c(0.5, 0.5), srt = 90)
  }
  
  # Add the y-axis only if log_y is TRUE
  if (log_y) {
    tick_length <- -tick_length
    minor_tick_length <- -tick_length/2
    
    # Major tick positions and labels for the log scale
    y_ticks <- c(1, 10, 100, 1000)  # Example custom major ticks
    log_y_ticks <- log(y_ticks)
    valid_ticks <- log_y_ticks[log_y_ticks >= ylim[1] & log_y_ticks <= ylim[2]]  # Filter major ticks within the plot range

    # Minor tick positions for the log scale
    minor_y_ticks <- c(2, 3, 4, 5, 6, 7, 8, 9, 
                       20, 30, 40, 50, 60, 70, 80, 90, 
                       200, 300, 400, 500, 600, 700, 800, 900)  # Example custom minor ticks
    log_minor_y_ticks <- log(minor_y_ticks)
    valid_minor_ticks <- log_minor_y_ticks[log_minor_y_ticks >= ylim[1] & log_minor_y_ticks <= ylim[2]]  # Filter minor ticks within the plot range

    # Add major ticks
    axis(2, at=valid_ticks, labels=y_ticks[log_y_ticks >= ylim[1] & log_y_ticks <= ylim[2]], tcl=tick_length, las=1)

    # Add minor ticks
    axis(2, at=valid_minor_ticks, labels=NA, tcl=minor_tick_length, las=1)
    
    # Add y-axis label
    mtext(ylab, side=2, line=2.5)
  }
  
  if (save) {
    dev.off()
  }
}


trace_extract <- function(rawdata, start_time=50, baseline=50, idx=1, response_sign_method = c('smooth', 'regression', 'cumsum')){
  ind1 <- (start_time - baseline)/dx
  ind2 <- start_time/dx

  y <- rawdata[, idx+1]
  y <- na.omit(y)
  y <- y[ind1:length(y)]
  y <- y - mean(y[ind1:ind2])
  # x <- 0:dx:(length(y)-1)*dx

  sign <- sign_fun(y, direction_method=response_sign_method) 
  return(sign * y)
}

# Function to process each sheet with user-defined base name
process_sheet <- function(sheet_name, summary, data, dt, stimulation_time, baseline, smooth, base_name, func=product2) {
  # Get the summary list
  summary_list <- get(paste0(base_name, '_summary'), envir = .GlobalEnv)
  summary_list[[sheet_name]] <- summary
  assign(paste0(base_name, '_summary'), summary_list, envir = .GlobalEnv)
  
  # Process fits
  fits <- t(sapply(1:length(summary), function(ii){
    X <- summary[[ii]]$output
    as.vector(t(X))
  }))
  
  # Create new column names by appending 1 and 2 to the original names
  if (identical(func, product1)){
    new_colnames <- rep(colnames(summary[[1]]$output), 2)
    new_colnames[new_colnames == 'amp'] <- c('A1', 'A2')   
  }else if (identical(func, product2)){
    new_colnames <- rep(colnames(summary[[1]]$output), 2)
    new_colnames[new_colnames == 'amp'] <- c('A1', 'A2')
  }else if (identical(func, product3)){
    new_colnames <- rep(colnames(summary[[1]]$output), 3)
    new_colnames[new_colnames == 'amp'] <- c('A1', 'A2', 'A3')
  }
  
  colnames(fits) <- new_colnames
  rownames(fits) <- 1:dim(fits)[1]
  
  fits_list <- get(paste0(base_name, '_fits'), envir = .GlobalEnv)
  fits_list[[sheet_name]] <- fits
  assign(paste0(base_name, '_fits'), fits_list, envir = .GlobalEnv)
  
  # Process peaks
  peaks <- sapply(2:dim(data[[sheet_name]])[2], function(ii){
    x <- data[[sheet_name]][,1]
    y <- data[[sheet_name]][,ii]
    
    x <- x[!is.na(y)]
    y <- y[!is.na(y)]
    peak.fun(y=y, dt=dt, stimulation_time=stimulation_time, baseline=baseline, smooth=smooth)
  })
  
  peaks_list <- get(paste0(base_name, '_peaks'), envir = .GlobalEnv)
  peaks_list[[sheet_name]] <- peaks
  assign(paste0(base_name, '_peaks'), peaks_list, envir = .GlobalEnv)
}



# # for analysis of noise baseline to determoine cut-off frequencies

# FFT_plot <- function(response, response0=NULL, dx=0.1, xlim=c(0, 3000), ylim=c(10^-5, 10^2), xlab='frequency (Hz)', ylab='amplitude (pA)', col='gray', lwd=1, width=6, height=6, filename="plot.svg", save=FALSE, freq4baseline=25, spar=0.85) {
  
#   y <- response

#   N <- length(y)
#   fft_result <- fft(y)
#   fft_magnitude <- Mod(fft_result) / N

#   # frequency axis
#   sampling_rate <- 1e3 / dx
#   frequencies <- seq(0, sampling_rate / 2, length.out = N / 2 + 1)

#   # Plot the frequency spectrum
#   X <- frequencies
#   Y <- fft_magnitude[1:(N / 2 + 1)]

#   # Filter frequencies within the desired range
#   in_range <- X >= xlim[1] & X <= xlim[2]
#   X <- X[in_range]
#   Y <- Y[in_range]

#   # Compute gain in dB if response0 is provided
#   if (!is.null(response0)) {
#     N0 <- length(response0)
#     fft_result0 <- fft(response0)
#     fft_magnitude0 <- Mod(fft_result0) / N0
#     Y0 <- fft_magnitude0[1:(N0 / 2 + 1)][in_range]
#     gain_db <- 20 * log10(Y / Y0)
#     ylab <- 'gain (dB)'
#     Y <- gain_db
#     if (is.null(ylim)) ylim <- range(c(10 * floor(min(Y) / 10), 1, gain_db))  # Adjust ylim for gain in dB
#     log_scale <- ""
#   } else {
#     log_scale <- "y"
#   }

#   if (save) {
#     # Open SVG device
#     svg(file = filename, width = width, height = height, bg = 'transparent')
#   } else {
#     dev.new(width = width, height = height, noRStudioGD = TRUE)
#   }

#   # Plot the FFT or gain
#   plot(X, Y, col = col, type = 'l', bty = 'l', las = 1, lwd = lwd, 
#        main = '', xlim = xlim, ylim = ylim, xlab = xlab, ylab = ylab, log = log_scale, axes = FALSE, frame = FALSE)

#   # Adding custom ticks and grid
#   axis(1, at = seq(xlim[1], xlim[2], by = 1000), tcl = -0.2)  # x-axis with custom tick length

#   if (is.null(response0)) {
#     # Major ticks (decades) for amplitude scale
#     major_ticks <- 10^(log10(ylim[1]):log10(ylim[2]))
#     # Minor ticks (in-between decades)
#     minor_ticks <- c(2:9 %o% 10^(log10(ylim[1]):log10(ylim[2])))

#     # Custom y-axis labels using expressions
#     major_labels <- sapply(major_ticks, function(tick) as.expression(bquote(10^.(log10(tick)))))

#     # Add major ticks and labels
#     axis(2, at = major_ticks, labels = major_labels, tcl = -0.2, las = 1)

#     # Add minor ticks without labels
#     axis(2, at = minor_ticks, labels = NA, tcl = -0.1)
#   } else {
#     # Major ticks for dB scale
#     major_ticks <- seq(10 * floor(min(Y) / 10), 0, by = 10)
#     axis(2, at = major_ticks, tcl = -0.2, las = 1)
#   }

#   # Adding grid
#   grid()

#   # Add a smooth spline interpolation line for the periodogram
#   if (is.null(response0)) {
#     smooth_fit <- smooth.spline(X, Y, spar = spar)
#     lines(smooth_fit, col = 'slateblue', lty = 3, lwd = 2)

#     # Calculate drop point at 1/sqrt(2) based on spline fit
#     smooth_values <- predict(smooth_fit, X)$y
#     drop_point <- 1 / sqrt(2) * mean(smooth_values[X<freq4baseline])
#     crossings <- which(diff(sign(smooth_values - drop_point)) != 0)  # Find all crossing points
    
#     # Linear interpolation to find more accurate crossing points
#     if (length(crossings) > 0) {
#       intersection_freqs <- numeric()
#       for (i in crossings) {
#         x1 <- X[i]
#         x2 <- X[i + 1]
#         y1 <- smooth_values[i]
#         y2 <- smooth_values[i + 1]
#         intersection_freqs <- c(intersection_freqs, x1 + (x2 - x1) * (drop_point - y1) / (y2 - y1))
#       }

#       if (length(intersection_freqs) > 0) {
#         median_freq <- median(intersection_freqs)  # Find the median frequency
#         print(paste0('cutoff frequency is ', round(median_freq,2)))
#         # Draw lines from axes to the median point
#         segments(x0 = xlim[1], y0 = drop_point, x1 = median_freq, y1 = drop_point, col = 'indianred', lty = 2)
#         segments(x0 = median_freq, y0 = ylim[1], x1 = median_freq, y1 = drop_point, col = 'indianred', lty = 2)
#       }
#     }
#   }

#   # Add -3 dB line if response0 is provided
#   if (!is.null(response0)) {
#     # Calculate -3 dB point
#     minus_3db <- -3
#     crossings <- which(diff(sign(Y - minus_3db)) != 0)  # Find all crossing points
    
#     # Linear interpolation to find more accurate crossing points
#     intersection_freqs <- numeric()
#     for (i in crossings) {
#       x1 <- X[i]
#       x2 <- X[i + 1]
#       y1 <- Y[i]
#       y2 <- Y[i + 1]
#       intersection_freqs <- c(intersection_freqs, x1 + (x2 - x1) * (minus_3db - y1) / (y2 - y1))
#     }

#     if (length(intersection_freqs) > 0) {
#       median_freq <- median(intersection_freqs)  # Find the median frequency
#       print(paste0('cutoff frequency is ', round(median_freq,2)))
#       # Draw lines from axes to the median point
#       segments(x0 = xlim[1], y0 = minus_3db, x1 = median_freq, y1 = minus_3db, col = 'indianred', lty = 2)
#       segments(x0 = median_freq, y0 = ylim[1], x1 = median_freq, y1 = minus_3db, col = 'indianred', lty = 2)
#     }
#   }

#   if (save) {
#     dev.off()
#   }
# }


FFT_plot <- function(response, response0 = NULL, dx = 0.1, xlim = c(0, 3000), ylim = c(10^-5, 10^2), 
                     xlab = 'frequency (Hz)', ylab = 'amplitude (pA)', col = 'gray', lwd = 1, 
                     width = 6, height = 6, filename = "plot.svg", save = FALSE, 
                     freq4baseline = 25, spar = 0.85, logx = FALSE) {
  
  y <- response
  N <- length(y)
  fft_result <- fft(y)
  fft_magnitude <- Mod(fft_result) / N

  # frequency axis
  sampling_rate <- 1e3 / dx
  frequencies <- seq(0, sampling_rate / 2, length.out = N / 2 + 1)

  # Plot the frequency spectrum
  X <- frequencies
  Y <- fft_magnitude[1:(N / 2 + 1)]

  # Filter frequencies within the desired range
  if (logx) xlim <- c(max(1, xlim[1]), xlim[2])
  in_range <- X >= xlim[1] & X <= xlim[2]
  X <- X[in_range]
  Y <- Y[in_range]

  # Compute gain in dB if response0 is provided
  if (!is.null(response0)) {
    N0 <- length(response0)
    fft_result0 <- fft(response0)
    fft_magnitude0 <- Mod(fft_result0) / N0
    Y0 <- fft_magnitude0[1:(N0 / 2 + 1)][in_range]
    gain_db <- 20 * log10(Y / Y0)
    ylab <- 'gain (dB)'
    Y <- gain_db
    if (is.null(ylim)) ylim <- range(c(10 * floor(min(Y) / 10), 1, gain_db))  # Adjust ylim for gain in dB
    log_scale <- ""
  } else {
    log_scale <- "y"
  }

  if (logx) {
    # Filter out zero values to avoid log(0)
    nonzero_idx <- X > 0
    X <- X[nonzero_idx]
    Y <- Y[nonzero_idx]
    # X <- log10(X)
    xlim <- c(max(1, xlim[1]), xlim[2])
    log_scale <- if (log_scale == "y") "xy" else "x"
  }

  if (save) {
    # Open SVG device
    svg(file = filename, width = width, height = height, bg = 'transparent')
  } else {
    dev.new(width = width, height = height, noRStudioGD = TRUE)
  }

  # Plot the FFT or gain
  plot(X, Y, col = col, type = 'l', bty = 'l', las = 1, lwd = lwd, 
       main = '', xlim = xlim, ylim = ylim, xlab = xlab, ylab = ylab, log = log_scale, axes = FALSE, frame = FALSE)

  # Adding custom ticks and grid
  if (logx) {
    major_ticks <- 10^(log10(xlim[1]):log10(xlim[2]))

    # Minor ticks (in-between decades)
    minor_ticks <- c(2:9 %o% 10^(log10(xlim[1]):log10(xlim[2])))
    # Custom y-axis labels using expressions
    major_labels <- sapply(major_ticks, function(tick) as.expression(bquote(10^.(log10(tick)))))
    # Add major ticks and labels
    axis(1, at = major_ticks, labels = major_labels, tcl = -0.2)
    # Add minor ticks without labels
    axis(1, at = minor_ticks, labels = NA, tcl = -0.1)
  } else {
    axis(1, at = seq(xlim[1], xlim[2], by = 1000), tcl = -0.2)  # x-axis with custom tick length
  }

  if (is.null(response0)) {
    # Major ticks (decades) for amplitude scale
    major_ticks <- 10^(log10(ylim[1]):log10(ylim[2]))
    # Minor ticks (in-between decades)
    minor_ticks <- c(2:9 %o% 10^(log10(ylim[1]):log10(ylim[2])))
    # Custom y-axis labels using expressions
    major_labels <- sapply(major_ticks, function(tick) as.expression(bquote(10^.(log10(tick)))))
    # Add major ticks and labels
    axis(2, at = major_ticks, labels = major_labels, tcl = -0.2, las = 1)
   # Add minor ticks without labels
    axis(2, at = minor_ticks, labels = NA, tcl = -0.1)
  } else {
    # Major ticks for dB scale
    major_ticks <- seq(10 * floor(min(Y) / 10), 0, by = 10)
    axis(2, at = major_ticks, tcl = -0.2, las = 1)
  }

  # Adding grid
  grid()

  # Add a smooth spline interpolation line for the periodogram
  if (is.null(response0)) {
    smooth_fit <- smooth.spline(X, Y, spar = spar)
    lines(smooth_fit, col = 'slateblue', lty = 3, lwd = 2)

    # Calculate drop point at 1/sqrt(2) based on spline fit
    smooth_values <- predict(smooth_fit, X)$y
    drop_point <- if (logx)  1 / sqrt(2) * mean(smooth_values[X < log10(freq4baseline)]) else 1 / sqrt(2) * mean(smooth_values[X<freq4baseline])
    crossings <- which(diff(sign(smooth_values - drop_point)) != 0)  # Find all crossing points
    
    # Linear interpolation to find more accurate crossing points
    if (length(crossings) > 0) {
      intersection_freqs <- numeric()
      for (i in crossings) {
        x1 <- X[i]
        x2 <- X[i + 1]
        y1 <- Y[i]     # smooth_values[i]
        y2 <- Y[i + 1] # smooth_values[i + 1]
        intersection_freqs <- c(intersection_freqs, x1 + (x2 - x1) * (drop_point - y1) / (y2 - y1))
      }

      if (length(intersection_freqs) > 0) {
        median_freq <- median(intersection_freqs)  # Find the median frequency
        # Draw lines from axes to the median point
        segments(x0 = xlim[1], y0 = drop_point, x1 = median_freq, y1 = drop_point, col = 'indianred', lty = 2)
        segments(x0 = median_freq, y0 = ylim[1], x1 = median_freq, y1 = drop_point, col = 'indianred', lty = 2)
        print(paste0('cutoff frequency is ', round(median_freq, 2)))       
      }
    }
  }

  # Add -3 dB line if response0 is provided
  if (!is.null(response0)) {
    # Calculate -3 dB point
    drop_point <- if (logx)  mean(Y[X < log10(freq4baseline)]) - 3 else 1 / sqrt(2) * mean(Y[X<freq4baseline])
    # if logx and input tends to output then drop_point is same as minus_3db i.e. -3
    crossings <- which(diff(sign(Y - drop_point)) != 0)  # Find all crossing points
    
    # Linear interpolation to find more accurate crossing points
    intersection_freqs <- numeric()
    for (i in crossings) {
      x1 <- X[i]
      x2 <- X[i + 1]
      y1 <- Y[i]
      y2 <- Y[i + 1]
      intersection_freqs <- c(intersection_freqs, x1 + (x2 - x1) * (drop_point - y1) / (y2 - y1))
    }

    if (length(intersection_freqs) > 0) {
      median_freq <- median(intersection_freqs)  # Find the median frequency
      # Draw lines from axes to the median point
      segments(x0 = xlim[1], y0 = drop_point, x1 = median_freq, y1 = drop_point, col = 'indianred', lty = 2)
      segments(x0 = median_freq, y0 = ylim[1], x1 = median_freq, y1 = drop_point, col = 'indianred', lty = 2)
      print(paste0('cutoff frequency is ', round(median_freq, 2)))      
    }
  }

  if (save) {
    dev.off()
  }
}


# for power gain_db <- 10 * log10(power / original_power)
periodogram_plot <- function(response, response0=NULL, dx=0.1, xlim=c(0, 3000), ylim=c(10^-10, 10^2), xlab='frequency (Hz)', ylab='power (pA^2/Hz)', col='gray', lwd=1, width=6, height=6, filename='plot.svg', save=FALSE, spar=0.85, freq4baseline=25) {
  
  y <- response
  fs <- 1e3 / dx  # Compute sampling frequency

  # Compute periodogram for the response
  periodogram <- spectrum(y, plot = FALSE)
  freqs <- periodogram$freq * fs
  power <- periodogram$spec

  # Filter frequencies within the desired range
  in_range <- freqs >= xlim[1] & freqs <= xlim[2]
  freqs <- freqs[in_range]
  power <- power[in_range]

  # Compute gain in dB if response0 is provided
  if (!is.null(response0)) {
    original_periodogram <- spectrum(response0, plot = FALSE)
    original_power <- original_periodogram$spec[in_range]
    gain_db <- 10 * log10(power / original_power)
    ylab <- 'gain (dB)'
    power <- gain_db
    if (is.null(ylim)) ylim <- range(c(10*floor(min(power)/10), 1, gain_db))  # Adjust ylim for gain in dB
    log_scale <- ""
  } else {
    log_scale <- "y"
  }

  if (save) {
    # Open SVG device
    svg(file=filename, width=width, height=height, bg='transparent')
  } else {
    dev.new(width=width, height=height, noRStudioGD=TRUE)
  }

  # Plot the periodogram or gain
  plot(freqs, power, col=col, type = 'l', bty='l', las=1, lwd=lwd, 
       main='', xlim=xlim, ylim=ylim, xlab=xlab, ylab=ylab, log=log_scale, axes=FALSE, frame=FALSE)

  # Adding custom ticks and grid
  axis(1, at = seq(xlim[1], xlim[2], by = 500), tcl = -0.2)  # x-axis with custom tick length

  if (is.null(response0)) {
    # Major ticks (decades) for linear scale
    major_ticks <- 10^(floor(log10(ylim[1])):ceiling(log10(ylim[2])))
    # Minor ticks (in-between decades)
    minor_ticks <- c(outer(2:9, 10^(floor(log10(ylim[1])):ceiling(log10(ylim[2])))))

    # Custom y-axis labels using expressions
    major_labels <- sapply(major_ticks, function(tick) as.expression(bquote(10^.(log10(tick)))))

    # Add major ticks and labels
    axis(2, at = major_ticks, labels = major_labels, tcl = -0.2, las=1)

    # Add minor ticks without labels
    axis(2, at = minor_ticks, labels = NA, tcl = -0.1)
  } else {
    # Major ticks for dB scale
    major_ticks <- seq(10*floor(min(power)/10), 0, by=10)
    axis(2, at = major_ticks, tcl = -0.2, las=1)
  }

  # Adding grid
  grid()

  # Add a smooth spline interpolation line for the periodogram
  if (is.null(response0)){
    smooth_fit <- smooth.spline(freqs, power, spar = spar)
    lines(smooth_fit, col='slateblue', lty=3, lwd=2)

    # Calculate drop point at 1/sqrt(2) based on spline fit
    smooth_values <- predict(smooth_fit, freqs)$y
    drop_point <- 1 / 2 * mean(smooth_values[freqs<freq4baseline])
    crossings <- which(diff(sign(smooth_values - drop_point)) != 0)  # Find all crossing points
    # Linear interpolation to find more accurate crossing points
    if (length(crossings) > 0) {
      intersection_freqs <- numeric()
      for (i in crossings) {
        x1 <- freqs[i]
        x2 <- freqs[i + 1]
        y1 <- smooth_values[i]
        y2 <- smooth_values[i + 1]
        intersection_freqs <- c(intersection_freqs, x1 + (x2 - x1) * (drop_point - y1) / (y2 - y1))
      }
    if (length(intersection_freqs) > 0) {
        median_freq <- median(intersection_freqs)  # Find the median frequency
        print(paste0('cutoff frequency is ', round(median_freq,2)))
        # Draw lines from axes to the median point
        segments(x0 = xlim[1], y0 = drop_point, x1 = median_freq, y1 = drop_point, col = 'indianred', lty = 2)
        segments(x0 = median_freq, y0 = ylim[1], x1 = median_freq, y1 = drop_point, col = 'indianred', lty = 2)
      }
    }
  }
  # Add -3 dB line if response0 is provided
  if (!is.null(response0)) {
    # Calculate -3 dB point
    minus_3db <- -3
    crossings <- which(diff(sign(power - minus_3db)) != 0)  # Find all crossing points
    
    # Linear interpolation to find more accurate crossing points
    intersection_freqs <- numeric()
    for (i in crossings) {
      x1 <- freqs[i]
      x2 <- freqs[i + 1]
      y1 <- power[i]
      y2 <- power[i + 1]
      intersection_freqs <- c(intersection_freqs, x1 + (x2 - x1) * (minus_3db - y1) / (y2 - y1))
    }

    if (length(intersection_freqs) > 0) {
      median_freq <- median(intersection_freqs)  # Find the median frequency
      print(paste0('cutoff frequency is ', round(median_freq,2)))
      # Draw lines from axes to the median point
      segments(x0 = xlim[1], y0 = minus_3db, x1 = median_freq, y1 = minus_3db, col = 'indianred', lty = 2)
      segments(x0 = median_freq, y0 = ylim[1], x1 = median_freq, y1 = minus_3db, col = 'indianred', lty = 2)
    }
  }

  if (save) {
    dev.off()
  }
}

processed_signals_plot <- function(input, output, width=6, height=6, main=''){
  dev.new(width=width, height=height, noRStudioGD=TRUE)
  plot(t * 1e3, input, type = 'l', col = 'darkgray', xlab = '', ylab = '', main = main, axes = FALSE, frame = FALSE, ylim = range(c(input, output)))
  lines(t * 1e3, output, col = 'lightgray')
  lines(c(max(t) * 1e3 - 1000, max(t) * 1e3), c(-3, -3), col = 'black', lwd = 2)
  lines(c(max(t) * 1e3 - 1000, max(t) * 1e3 - 1000), c(-3, -1), col = 'black', lwd = 2)
}


FFT_plot2 <- function(output, input=NULL, t, fs, width=12, height=6, xlab = 'frequency [Hz]',  ylab = 'Vm/Vc', spar = 0.05) {
  compute_fft <- function(signal, N) {
    fft_result <- fft(signal)
    freq <- (0:(length(t) - 1)) * (fs / length(t))
    valid_indices <- freq < 3e3
    list(
      fft = fft_result[valid_indices],
      freq = freq[valid_indices]
    )
  }
  
  add_epsilon <- function(magnitude, N, epsilon = 1e-6) {
    abs(magnitude) / N + epsilon
  }
  
  find_cutoff_frequency <- function(positive_freq, positive_vm_vc_ratio) {
    spline_fit <- smooth.spline(positive_freq, positive_vm_vc_ratio, spar = spar)
    fitted_spline <- predict(spline_fit, positive_freq)$y
    max_response <- max(positive_vm_vc_ratio)
    cutoff_level <- max_response / sqrt(2)
    cutoff_index <- which.min(abs(fitted_spline - cutoff_level))
    positive_freq[cutoff_index]
  }
  
  plot_signals <- function(t, input, output) {
    plot(t * 1e3, input, type = "l", col = "darkgray", xlab = "", ylab = "", main = "", axes = FALSE, frame = FALSE, ylim = range(c(input, output)))
    lines(t * 1e3, output, col = "lightgray")
    lines(c(max(t) * 1e3 - 1000, max(t) * 1e3), c(-3, -3), col = "black", lwd = 2)
    lines(c(max(t) * 1e3 - 1000, max(t) * 1e3 - 1000), c(-3, -1), col = "black", lwd = 2)
  }
  
  plot_frequency_response <- function(freq, vm_vc_ratio, fitted_spline, cutoff_frequency, xlab, ylab) {
    plot(freq, vm_vc_ratio, type = "l", col = "slateblue", log = "x", xlab = xlab, ylab = ylab, main = "", frame = FALSE)
    lines(freq, fitted_spline, col = "indianred", lty = "dotted")
    abline(v = cutoff_frequency, col = "darkgray", lty = "dotted")
  }
  
  if (!is.null(input)) {
    N <- length(input)
    input_fft_res <- compute_fft(input, N)
    output_fft_res <- compute_fft(output, N)
    
    input_fft_magnitude <- add_epsilon(input_fft_res$fft, N)
    output_fft_magnitude <- add_epsilon(output_fft_res$fft, N)
    
    vm_vc_ratio <- output_fft_magnitude / input_fft_magnitude
    positive_freq_indices <- input_fft_res$freq > 0
    positive_freq <- input_fft_res$freq[positive_freq_indices]
    positive_vm_vc_ratio <- vm_vc_ratio[positive_freq_indices]
    
    cutoff_frequency <- find_cutoff_frequency(positive_freq, positive_vm_vc_ratio)
    
    dev.new(width=width, height=height, noRStudioGD=TRUE)
    par(mfrow = c(1, 2))  
    plot_signals(t, input, output)
    plot_frequency_response(positive_freq, positive_vm_vc_ratio, predict(smooth.spline(positive_freq, positive_vm_vc_ratio, spar = spar), positive_freq)$y, cutoff_frequency, xlab, ylab)
    
    print(sprintf("-3 dB cutoff frequency is approximately %.2f Hz", cutoff_frequency))
  } else {
    N <- length(output)
    output_fft_res <- compute_fft(output, N)
    
    output_fft_magnitude <- add_epsilon(output_fft_res$fft, N)
    positive_freq_indices <- output_fft_res$freq > 0
    positive_freq <- output_fft_res$freq[positive_freq_indices]
    positive_output_fft_magnitude <- output_fft_magnitude[positive_freq_indices]
    
    cutoff_frequency <- find_cutoff_frequency(positive_freq, positive_output_fft_magnitude)
    
    par(mfrow = c(1, 2))
    plot_signals(t, output, output)
    plot_frequency_response(positive_freq, positive_output_fft_magnitude, predict(smooth.spline(positive_freq, positive_output_fft_magnitude, spar = spar), positive_freq)$y, cutoff_frequency, xlab, ylab)
    
    print(sprintf("-3 dB cutoff frequency is approximately %.2f Hz", cutoff_frequency))
  }
}


# functions for voltage step analysis

# function to plot step for visualisation
step_plot <- function(x, y, tstep=c(50, 250), xlab='time (ms)', ylab='PSC amplitude (pA)', ylim=NULL, lwd=1.2, filter=FALSE, width=5, height=5, save=FALSE) {
  
  xfit <- x[x<sum(tstep)]
  yfit <- y[x<sum(tstep)]
  
  xlim <- c(10*floor(min(xfit)/10), 10*ceiling(max(xfit)/10))

  dev.new(width=width, height=height, noRStudioGD=TRUE)
  plot(xfit, yfit, col='gray', xlab=xlab, ylab=ylab, ylim=ylim, xlim=xlim, type='l', bty='l', las=1, lwd=lwd, main='')
  abline(v = tstep[1], col='darkgray', lwd=lwd, lty=3) 
  abline(v = tstep[2], col='darkgray', lwd=lwd, lty=3) 

}

vc_step_fit <- function(response, dt=0.1, tstep=c(50, 250), dV=5, t_interval=25, method= c('LM', 'BF.LM', 'GN', 'port', 'robust', 'MLE'), 
  xlab='time (ms)', ylab='PSC amplitude (pA)', ylim=NULL, response_sign_method='smooth', lwd=1.2, width=8, height=4, 
  return.output=FALSE, show.output=TRUE, show.plot=TRUE) {

  method <- match.arg(method)
  fun.2.fit <- function(params,x){
    -params[1]*exp(-params[2]*x) -params[3]*exp(-params[4]*x) + params[1] + params[3] + params[5]
  }
  form.2.fit <- y ~ -a*exp(-b*x) -c*exp(-d*x) + a + c + e

  y <- response
  x <- seq(0, length(y)-1) * dt

  X <- x[x<(tstep[2]-dt)]
  Y <- y[x<(tstep[2]-dt)]

  sign <- sign_fun(Y, direction_method=response_sign_method) 
  Y <- sign * Y

  idx1 <- which.max(Y)
  A <- Y[idx1]

  A.fit <- 0.9 * A
  idx2 <-which.min(abs(Y[idx1:length(Y)] - A.fit)) + idx1 - 1
  xfit <- X[idx2:length(X)]
  yfit <- Y[idx2:length(X)]

  output <- nFIT(response=yfit, dt=dt, func=fun.2.fit, method='LM', lower=NULL, upper=NULL, filter=FALSE, return.output=TRUE, show.output=FALSE, show.plot=FALSE)

  fits <- output$fits
  traces <- output$traces

  n <- idx2-idx1
  x2fit <- c(rev(-seq(1:n)) * dt, traces$x)
  y2fit <- fun.2.fit(fits, x2fit)  

  traces=data.frame(x=x2fit - x2fit[1] + tstep[1], y=sign * Y[idx1:length(X)], yfit=sign * y2fit)
    
  dI <- fits[1] + fits[3] + fits[5]
  Rm <- dV/dI * 1e3
  Rs <- dV/y2fit[1] * 1e3
  Cm <- abs(fits[1])/fits[2]/dV + abs(fits[3])/fits[4]/dV # area under given exponetial of form a * exp(-bx) is a/b pAms ie fC then divide by 5 mV gives pC
  
  # to check area use trapz in pracma package:
  Cm_est <- ( pracma:::trapz(x2fit, y2fit) - dI * (max(x2fit) - min(x2fit)) ) / dV
  
  df_output <- data.frame(Cm=Cm, Cm_est=Cm_est, Rm=Rm, Rs=Rs)
  row.names(df_output) <- ''

  # Display the output
  df_output

  if (show.plot){
    dev.new(width=width, height=height, noRStudioGD=TRUE)
    par(mfrow = c(1, 2))
    
    step_plot2(x=x[x<sum(tstep)], y=y[x<sum(tstep)], x2fit=traces$x, y2fit=traces$yfit, tstep=tstep)

    ylim <- c(5*floor(min(traces$y)/5), 0)
    step_plot2(x=traces$x[traces$x > tstep[1] & traces$x < tstep[1] + t_interval], y=traces$y[traces$x > tstep[1] & traces$x < tstep[1] + t_interval], 
      x2fit=traces$x[traces$x > tstep[1] & traces$x < tstep[1] + t_interval], y2fit=traces$yfit[traces$x > tstep[1] & traces$x < tstep[1] + t_interval], tstep=tstep, ylim=ylim, round2=5)
  }

  if (show.output){
    print(df_output)
  }
  
  if (return.output){
    return(list(pars=df_output, fits=fits))
  }

}

step_plot2 <- function(x, y, x2fit, y2fit, tstep=c(50, 250), xlab='time (ms)', ylab='PSC amplitude (pA)', ylim=NULL, lwd=1.2, filter=FALSE, width=5, height=5, tcl = -0.2, round2=10, save=FALSE) {
  xlim <- c(round2*floor(min(x)/round2), round2*ceiling(max(x)/round2))
  ylim <- if (is.null(ylim)) c(round2*floor(min(y2fit)/round2), round2*ceiling(max(y)/round2)) else ylim
  
  plot(x, y, col='darkgray', xlab=xlab, ylab=ylab, ylim=ylim, xlim=xlim, type='l', bty='l', las=1, lwd=lwd, main='',
       tcl=tcl, xaxs='i', yaxs='i')
  
  lines(x2fit, y2fit, col='indianred', lty=3, lwd=lwd)
}



# Function to wait for a mouse click
wait_for_click <- function() {
  cat("Click on the plot to continue\n")
  locator(1)
}


# Function to generate initial start values
generate_start_values <- function(y, x) {
  xmin <- min(x) * runif(1, 0.9, 1.1)
  Imax <- max(y) * runif(1, 0.9, 1.1)
  tau <- runif(1, 0, 1)
  list(xmin = xmin, Imax = Imax, tau=tau)
}

# Function to fit the model with a fallback to fixing n = 1 if fitting fails
fit_fun <- function(y, x, start_values, lower_bounds, upper_bounds) {
  fit <- tryCatch({
    nlsLM(y ~ exp_fun(x, xmin, Imax, tau), 
          start = start_values, 
          lower = lower_bounds, 
          upper = upper_bounds, 
          control = nls.lm.control(maxiter = 500))
  }, error = function(e) {
    cat("Initial fit failed with error:", e$message, "\n")
    NULL
  })
  
  
  return(fit)
}

fits_fun <- function(mat, ylab = "amplitude (pA)", ids = NULL, attempts = 10, seed = 7) {
  # Run the plotting and fitting in lapply with mouse click pause
  set.seed(seed)
  N <- dim(mat)[2]
  
  # Open a graphical device
  dev.new(width = 4, height = 4, noRStudioGD = TRUE)
  
  fits <- lapply(1:N, function(ii) {
    y <- c(rev(mat[, ii]))
    attempt <- 0
    fit <- NULL

    while (is.null(fit) && attempt < attempts) {
      start_values <- generate_start_values(y, x)
      lower_bounds <- c(0, 0, 0)
      upper_bounds <- c(1, 10 * max(y), 1)
      
      fit <- fit_fun(y, x, start_values, lower_bounds, upper_bounds)
      attempt <- attempt + 1
    }

    if (!is.null(fit)) {
      plot(x, y, main = "", xlab = "input (intensity)", ylab = ylab, 
           xlim = c(0, 1), ylim = c(0, max(y) * 1.1), bty = 'n')
      curve(exp_fun(x, coef(fit)[1], coef(fit)[2], coef(fit)[3]), add = TRUE, 
            col = "slateblue", lty = 3, lwd = 2)
    } else {
      plot(x, y, main = "fit fails", xlab = "input (intensity)", ylab = ylab, 
           xlim = c(0, 1), ylim = c(0, max(y) * 1.1), bty = 'n')
    }

    wait_for_click()

    fit
  })
  
  # Clear the final plot
  dev.off()
  names(fits) <- ids

  # Print fit summaries
  lapply(fits, function(fit) {
    if (!is.null(fit)) {
      summary(fit)
    } else {
      cat("fit fails\n")
    }
  })

  # Extract coefficients
  coefficients <- do.call(rbind, lapply(fits, function(fit) {
    if (!is.null(fit)) {
      as.data.frame(t(coef(fit)))
    } else {
      data.frame(xmin = NA, Imax = NA, tau = NA)
    }
  }))

  return(coefficients)
}

exp_fun <- function(x, xmin, Imax, tau)
  Imax * (1 - exp(-(x-xmin)/tau))



IO_plot <- function(boxplot_values, main='', xlab='relative LED intensity', ylab='relative amplitude', col='indianred', 
  xlim=c(0, 1.05), ylim=c(0, 1.05), x_tick_interval=0.2, y_tick_interval=0.2, tick_length=-0.2, xround=1, yround=1, 
  seed=42, fit.attempts=10, maxiter=500, pch=20, lwd=1.2, height=4, width=4, filename='IO.svg', save=FALSE){
  
  y <- abs(boxplot_values$Median)
  x <- boxplot_values$x

  if (save) {
    # Open SVG device
    svg(file=filename, width=width, height=height, bg='transparent')
  } else {
    dev.new(width=width, height=height, noRStudioGD=TRUE)
  }

  plot(x, y, main=main, xlab=xlab, ylab=ylab, col=col, pch=pch, xlim=xlim, ylim=ylim, lwd=lwd, xaxt='n', yaxt='n', bty='n')

  # Fit the exponential to median values
  lower_bounds <- c(0, 0, 0)  
  upper_bounds <- c(1, 2*max(y), 10) 

  set.seed(seed)
  
  attempt <- 0
  fit <- NULL
  while (is.null(fit) && attempt < fit.attempts) {
    start_values <- generate_start_values(y, x)
    fit <- nlsLM(y ~ exp_fun(x, xmin, Imax, tau), 
                 start = start_values, 
                 lower = lower_bounds, 
                 upper = upper_bounds, 
                 control = nls.lm.control(maxiter=maxiter))
    attempt <- attempt + 1
  }

  # Fit model to the Q1 values
  y_q1 <- abs(boxplot_values$Q1)
  upper_bounds <- c(1, 2*max(y_q1), 10) 
  
  attempt <- 0
  fit1 <- NULL
  while (is.null(fit1) && attempt < fit.attempts) {
    start_values <- generate_start_values(y_q1, x)
    fit1 <- nlsLM(y_q1 ~ exp_fun(x, xmin, Imax, tau), 
                 start = start_values, 
                 lower = lower_bounds, 
                 upper = upper_bounds, 
                 control = nls.lm.control(maxiter=maxiter))
    attempt <- attempt + 1
  }

  # Fit model to the Q3 values
  y_q3 <- abs(boxplot_values$Q3)
  upper_bounds <- c(1, 2*max(y_q3), 10) 
  
  attempt <- 0
  fit3 <- NULL
  while (is.null(fit3) && attempt < fit.attempts) {
    start_values <- generate_start_values(y_q3, x)
    fit3 <- nlsLM(y_q3 ~ exp_fun(x, xmin, Imax, tau), 
                 start = start_values, 
                 lower = lower_bounds, 
                 upper = upper_bounds, 
                 control = nls.lm.control(maxiter=maxiter))
    attempt <- attempt + 1
  }

  # Filter x-values for which the y-values of the fitted curve are positive
  x_curve <- seq(0, max(x), length.out = 200)
  
  y_fit <- exp_fun(x_curve, coef(fit)[1], coef(fit)[2], coef(fit)[3])
  x_curve_pos <- x_curve[y_fit > 0]
  y_fit_pos <- y_fit[y_fit > 0]
  lines(x_curve_pos, y_fit_pos, col = col, lty = 1, lwd = lwd)
  
  y_fit1 <- exp_fun(x_curve, coef(fit1)[1], coef(fit1)[2], coef(fit1)[3])
  x_curve_pos1 <- x_curve[y_fit1 > 0]
  y_fit_pos1 <- y_fit1[y_fit1 > 0]
  lines(x_curve_pos1, y_fit_pos1, col = "gray", lty = 3, lwd = lwd)

  y_fit3 <- exp_fun(x_curve, coef(fit3)[1], coef(fit3)[2], coef(fit3)[3])
  x_curve_pos3 <- x_curve[y_fit3 > 0]
  y_fit_pos3 <- y_fit3[y_fit3 > 0]
  lines(x_curve_pos3, y_fit_pos3, col = "gray", lty = 3, lwd = lwd)

  # Define tick intervals and lengths
  x_ticks <- seq(xround*floor(xlim[1]/xround), xround*ceiling(xlim[2]/xround), by=x_tick_interval)
  y_ticks <- seq(yround*floor(ylim[1]/yround), yround*ceiling(ylim[2]/yround), by=y_tick_interval)
  
  # Customize x-axis
  axis(1, at=x_ticks, labels=x_ticks, tcl=tick_length, lwd=lwd)

  # Customize y-axis with horizontal labels
  axis(2, at=y_ticks, labels=y_ticks, tcl=tick_length, las=1, lwd=lwd)

  if (save) {
    dev.off()
  }

}


IO_plot2 <- function(boxplot_values, main='', xlab='relative LED intensity', ylab='relative amplitude', col='indianred', 
  xlim=c(0, 1.05), ylim=c(0, 1.05), x_tick_interval=0.2, y_tick_interval=0.2, tick_length=-0.2, xround=1, yround=1, 
  seed=42, fit.attempts=10, maxiter=500, pch=20, lwd=1.2, height=4, width=4, spar=0.5, filename='IO.svg', save=FALSE){

  y <- abs(boxplot_values$Median)
  x <- boxplot_values$x

  if (save) {
    # Open SVG device
    svg(file=filename, width=width, height=height, bg='transparent')
  } else {
    dev.new(width=width, height=height, noRStudioGD=TRUE)
  }

  plot(x, y, main=main, xlab=xlab, ylab=ylab, col=col, pch=pch, xlim=xlim, ylim=ylim, lwd=lwd, xaxt='n', yaxt='n', bty='n')

  # Fit the exponential to median values
  lower_bounds <- c(0, 0, 0)  
  upper_bounds <- c(1, 2*max(y), 1) 

  set.seed(seed)
  
  attempt <- 0
  fit <- NULL
  while (is.null(fit) && attempt < fit.attempts) {
    start_values <- generate_start_values(y, x)
    fit <- nlsLM(y ~ exp_fun(x, xmin, Imax, tau), 
                 start = start_values, 
                 lower = lower_bounds, 
                 upper = upper_bounds, 
                 control = nls.lm.control(maxiter=maxiter))
    attempt <- attempt + 1
  }

  # Fit spline to Q1 values
  y_q1 <- abs(boxplot_values$Q1)
  spline_fit1 <- smooth.spline(x, y_q1, spar=spar)

  # Fit spline to Q3 values
  y_q3 <- abs(boxplot_values$Q3)
  spline_fit3 <- smooth.spline(x, y_q3, spar=spar)

  # Filter x-values for which the y-values of the fitted exponential curve are positive
  x_curve <- seq(0, max(x), length.out = 200)
  
  y_fit <- exp_fun(x_curve, coef(fit)[1], coef(fit)[2], coef(fit)[3])
  x_curve_pos <- x_curve[y_fit > 0]
  y_fit_pos <- y_fit[y_fit > 0]
  lines(x_curve_pos, y_fit_pos, col = col, lty = 1, lwd = lwd)

  # Filter and plot the spline for Q1 where the y-values are positive
  spline_fit1$y <- spline_fit1$y[spline_fit1$y > 0]
  spline_fit1$x <- spline_fit1$x[1:length(spline_fit1$y)]
  lines(spline_fit1, col="gray", lty=3, lwd=lwd)

  # Filter and plot the spline for Q3 where the y-values are positive
  spline_fit3$y <- spline_fit3$y[spline_fit3$y > 0]
  spline_fit3$x <- spline_fit3$x[1:length(spline_fit3$y)]
  lines(spline_fit3, col="gray", lty=3, lwd=lwd)

  # Define tick intervals and lengths
  x_ticks <- seq(xround*floor(xlim[1]/xround), xround*ceiling(xlim[2]/xround), by=x_tick_interval)
  y_ticks <- seq(yround*floor(ylim[1]/yround), yround*ceiling(ylim[2]/yround), by=y_tick_interval)
  
  # Customize x-axis
  axis(1, at=x_ticks, labels=x_ticks, tcl=tick_length, lwd=lwd)

  # Customize y-axis with horizontal labels
  axis(2, at=y_ticks, labels=y_ticks, tcl=tick_length, las=1, lwd=lwd)

  if (save) {
    dev.off()
  }
}

datasets2list_2 <- function(names, idxs, id='_fits'){
  all_objects <- ls(envir = .GlobalEnv)
  objects <- all_objects[grep(paste0(id, "$"), all_objects)]
  out_list <- list()
  # Loop through each name and id to group the objects
  for (name in names) {
    for (idx in idxs) {
      # Construct the object name pattern
      pattern <- paste0(name, idx, id)
      
      # Find the object that matches the pattern
      matched_object <- objects[objects == pattern]
      
      # If a match is found, add it to the grouped_fits list
      if (length(matched_object) > 0) {
        if (is.null(out_list[[name]])) {
          out_list[[name]] <- list()
        }
        out_list[[name]][[idx]] <- get(matched_object)
      }
    }
  }
  return(out_list)
}



traces_smoothfits <- function(y, fits, dt=0.1, N=1, IEI=50, stimulation_time=150, baseline=50, func=product1N, filter=FALSE, fc=1000, upsample.fit = c(upsample=TRUE, factor=100)){
  
  if (all(is.na(y[(which(!is.na(y))[length(which(!is.na(y)))] + 1):length(y)]))) {
    y <- y[!is.na(y)]
  }

  upsample <- upsample.fit[['upsample']]
  factor <- upsample.fit[['factor']]

  dx <- dt  
  x <- seq(0, (length(y) - 1) * dx, by = dx)

  if (filter){
    ind = 20
    fc = fc; fs = 1/dx*1000; bf <- butter(2, fc/(fs/2), type='low')
    yfilter <- signal::filter(bf, y)
  } else {
    ind=1
    yfilter=y
  }

  ind1 <- (stimulation_time - baseline)/dx
  ind2 <- baseline/dx
  
  yorig <- y[ind1:length(y)]
  yfilter <- yfilter[ind1:length(yfilter)]
  xorig <- seq(0, dx * (length(yorig) - 1), by = dx)

  xfit <- if (upsample) seq(0, dx * (length(yorig) - 1), by = dx/factor) else xorig


  yorig <- yorig - mean(yorig[1:ind2])
  yfilter <- yfilter - mean(yfilter[1:ind2])

  traces1 <- data.frame(x=xorig, y=yorig, yfilter=yfilter)
  traces2 <- data.frame(xfit=xfit)


  if (identical(func, product1N)){
    fits[N+3] <- fits[N+3] + baseline
  } else if (identical(func, product2N)){
    fits[N+3] <- fits[N+3] + baseline; fits[2*N+6] <- fits[2*N+6] + baseline
  } else if (identical(func, product3N)){
    fits[N+3] <- fits[N+3] + baseline; fits[2*N+6] <- fits[2*N+6] + baseline; fits[3*N+9] <- fits[3*N+9] + baseline
  }    
  traces2$yfit <- func(fits,traces2$x+dx, N=N, IEI=IEI)  
  if (identical(func, product2N)){
    traces2$yfit1 <- product1N(fits[1:(N+3)],traces2$x+dx, N=N, IEI=IEI) 
    traces2$yfit2 <- product1N(fits[(N+4):(2*N+6)],traces2$x+dx, N=N, IEI=IEI) 
  } 
  if (identical(func, product3N)){
    traces2$yfit1 <- product1N(fits[1:(N+3)],traces2$x+dx, N=N, IEI=IEI) 
    traces2$yfit2 <- product1N(fits[(N+4):(2*N+6)],traces2$x+dx, N=N, IEI=IEI)
    traces2$yfit3 <- product1N(fits[(2*N+7):(3*N+9)],traces2$x+dx, N=N, IEI=IEI) 
  } 

  return(list(original=traces1, fit=traces2))
}



fit_plot2 <- function(traces, func=product1N, xlab='time (ms)', ylab='PSC amplitude (pA)', xlim=NULL, ylim=NULL, bl=NULL, lwd=1.2, filter=FALSE, width=4, height=4, bg='transparent', filename='trace.svg', save=FALSE) {
  
  if (save) {
    # Open SVG device
    svg(file=filename, width=width, height=height, bg=bg)
  } else {
    dev.new(width=width, height=height, noRStudioGD=TRUE)
  }

  plot(traces$original$x, traces$original$y, col='gray', xlab=xlab, ylab=ylab, xlim=xlim, ylim=ylim, type='l', bty='l', las=1, lwd=lwd, main='')
  
  if (filter) {
    lines(traces$original$x, traces$original$yfilter, col='black', type='l', lwd=lwd)
  }
  
  lines(traces$fit$xfit, traces$fit$yfit, col='#CD5C5C', lty=3, lwd=2 * lwd)
  
  if (identical(func, product2N)) {
    lines(traces$fit$x, traces$fit$yfit1, col='#4C78BC', lty=3, lwd=2 * lwd)
    lines(traces$fit$x, traces$fit$yfit2, col='#CA92C1', lty=3, lwd=2 * lwd)
  }

  if (identical(func, product3N)) {
    lines(traces$fit$x, traces$fit$yfit1, col='#F28E2B', lty=3, lwd=2 * lwd)
    lines(traces$fit$x, traces$fit$yfit2, col='#4C78BC', lty=3, lwd=2 * lwd)
    lines(traces$fit$x, traces$fit$yfit3, col='#CA92C1', lty=3, lwd=2 * lwd)
  }

  if (!is.null(bl)) abline(v=bl, col='black', lwd=lwd, lty=3)

  if (save) {
    dev.off()
  }
}



smooth.plots <- function(y, fits, N=1, IEI=50, dt=0.1,  stimulation_time=150, baseline=50, func=product1N, filter=FALSE, fc=1000, upsample.fit = c(upsample=FALSE, factor=100),
  xlab='time (ms)', ylab='', xlim=NULL, ylim=NULL, lwd=1.2, width=5, height=5, bg='transparent', filename='trace.svg', save=FALSE){

  traces <- traces_smoothfits(y=y, fits=fits, N=N, IEI=IEI, dt=dt,  stimulation_time=stimulation_time, baseline=baseline, func=func, filter=filter, fc=fc, upsample.fit = upsample.fit)

  fit_plot2(traces=traces, func=func, xlab=xlab, ylab=ylab, xlim=xlim, ylim=ylim, lwd=lwd, filter=filter, width=width, height=height, bg=bg, filename=filename, save=save) 
  
}


# DBSCAN (Density-Based Spatial Clustering of Applications with Noise)
# DBSCAN is an algorithm that discovers clusters by identifying regions of high density, which are separated by areas of lower density. 
# It is effective for finding clusters of varying shapes and sizes in large datasets, even when noise or outliers are present.
# The algorithm operates using two key parameters:
# minPts: The minimum number of points required for a region to be considered a cluster.
# epsilon (eps): Defines the neighborhood radius around a point for identifying neighboring points, often measured using Euclidean distance.

# k-NN distance: In general, for a given value of k, the k-NN distance is the distance from a point to its k-th nearest neighbor

# The kNN distance plot helps identify the “elbow” point, which corresponds to a good eps value
# Use k = minPts - 1; minPts <- 4  is default for 2D data
# Small minPts values (e.g., 3 or 4) work well for data with:
#   Small clusters
#   High density
# Large minPts values (e.g., 10 or 15) work well for:
#   Large datasets
#   Clusters with high variability in density

# dbscan::dbscan performs the following operations
#   1.  Core Points Identification:
#   • For each point in the dataset:
#   • Compute the number of points within a radius of eps (including the point itself).
#   • If this count is ≥ minPts, the point is classified as a core point.
#   2.  Cluster Formation:
#   • Starting from a core point, expand the cluster by adding all points that are directly reachable from the core point:
#   • A point is directly reachable if it lies within the eps radius of the core point.
#   • If this point is also a core point, repeat the process to find more directly reachable points (density reachability).
#   3.  Noise Points Identification:
#   • Any point that:
#   • Is not a core point.
#   • Is not directly reachable from any core point.
#   • These points are classified as noise (or outliers) and assigned a cluster label of 0.


kNNdistplot2 <- function(x, k, minPts,  bty="n", lwd=1, lty=1, axes=FALSE, frame=FALSE, xtick=1, ...) {
    
    if (missing(k) && missing(minPts)) 
        stop("k or minPts need to be specified.")
    if (missing(k)) 
        k <- minPts - 1
    
    # Sort the k-NN distances
    kNNdist <- sort(dbscan::kNNdist(x, k, ...))
    factor <- 10^floor(log10(max(kNNdist)))  # Dynamically determine rounding factor
    ylim <- c(0, factor * ceiling(max(kNNdist) / factor))
    # Dynamically calculate tick intervals based on the unified limit
    yticks <- pretty(ylim, n = 5)
    xlim <- c(1, length(kNNdist))
    # Plot with additional parameters passed using ...
    plot(kNNdist, type="l", ylab=paste0(k, "-NN distance"), 
         xlab="points sorted by distance", xlim=xlim, ylim=ylim, 
         bty=bty, lwd=lwd, lty=lty, axes=axes, frame=frame, ...)
 
  if (!axes){
    axis(1, at=seq(xlim[1], xlim[2], by=xtick), tcl=-0.2, las = 1)  
    axis(2, at=yticks, tcl=-0.2, las = 1)  
  }

}

DBSCAN_analyse <- function(data, minPts = 4, k = NULL, height = 5, width = 10, bg = 'transparent', 
                           eps = NA, filename = 'DBscan.svg', save = FALSE) {
  # Ensure k is derived from minPts if k is not provided
  if (is.null(k)) {
    k <- minPts - 1
  }
  
  # Open a wider plotting window
  dev.new(width = width, height = height, noRStudioGD = TRUE)
  
  # Split the plotting region into two panels
  par(mfrow = c(1, 2))  # 1 row, 2 columns
  
  # Panel 1: kNN plot
  kNNdistplot2(data, k = k)
  
  proceed <- if (is.na(eps)) FALSE else TRUE
  
  if (proceed) {
    # Draw a horizontal line at the specified eps
    abline(h = eps, col = 'indianred', lty = 3)
  } else {
    while (!proceed) {
      # Prompt for eps value
      while (is.na(eps)) {
        cat('\nEnter eps: ')
        eps <- as.numeric(readLines(n = 1))
        if (is.na(eps)) {
          cat('\nInvalid input. Please enter a numeric value.\n')
        }
      }

      # Redraw the kNN plot in the left panel
      par(mfg = c(1, 1))  # Focus back on the first panel
      kNNdistplot2(data, k = k)
      abline(h = eps, col = 'indianred', lty = 3)

      # Ask user if they're happy with the eps
      cat('\nAre you happy with the position of the line (y/n)? ')
      response <- tolower(readLines(n = 1))

      if (response == 'y') {
        proceed <- TRUE
      } else {
        eps <- NA
        cat('\nTry again...\n')
      }
    }
  }

  # Apply DBSCAN on the dataset
  dbscan_result <- dbscan::dbscan(x=data, eps=eps, minPts=minPts)
  
  # Panel 2: DBSCAN plot
  par(mfg = c(1, 2))  # Switch to the second panel
  DBSCAN_plot(data, dbscan_result)

  # Save the plot to an SVG file if save = TRUE
  if (save) {
    # Open SVG device
    svg(file = filename, width = width, height = height, bg = bg)
    # Recreate the plots
    par(mfrow = c(1, 2))
    # Replot kNN plot
    kNNdistplot2(data, k = k)
    abline(h = eps, col = 'indianred', lty = 3)
    # Replot DBSCAN plot
    par(mfg = c(1, 2))
    DBSCAN_plot(data, dbscan_result)
    dev.off()  # Close SVG device
  }
}

DBSCAN_plot <- function(data, dbscan_result) {
  # Ensure `data` is a matrix for proper subsetting
  if (!is.matrix(data)) {
    data <- as.matrix(data)
  }
  
  # Extract column names for axis labels
  xlab <- colnames(data)[1]
  ylab <- colnames(data)[2]
  
  # Determine the maximum value across both axes
  # max_limit <- 100 * ceiling(max(data) / 100)
  factor <- 10^floor(log10(max(data))) 
  max_limit <- factor * ceiling(max(data) / factor)

  lim <- c(0, max_limit)
  
  # Dynamically calculate tick intervals based on the unified limit
  ticks <- pretty(lim, n = 5)
  
  # Plot all data points
  plot(data, col = 'black', pch = 19, cex = 0.75, main = 'DBSCAN clustering results',
       xlim = lim, ylim = lim, bty = 'n', lwd = 1, lty = 1, axes = FALSE, frame = FALSE,
       xlab = xlab, ylab = ylab)
  
  # Highlight noise points (cluster == 0)
  noise_indices <- which(dbscan_result$cluster == 0)
  if (length(noise_indices) > 0) {
    noise_points <- data[noise_indices, , drop = FALSE]
    points(noise_points, col = 'indianred', pch = 19)
    text(noise_points, labels = noise_indices, pos = 4, col = 'indianred', cex = 0.75)
  }
  
  # Add axes with adaptive ticks
  axis(1, at = ticks, tcl = -0.2, las = 1)  # x-axis
  axis(2, at = ticks, tcl = -0.2, las = 1)  # y-axis
}

# DBSCAN_plot <- function(data, dbscan_result, height = height, width = width) {
#   dev.new(width = width, height = height, noRStudioGD = TRUE)
  
#   # Ensure `data` is a matrix for proper subsetting
#   if (!is.matrix(data)) {
#     data <- as.matrix(data)
#   }
  
#   # Extract column names for axis labels
#   xlab <- colnames(data)[1]
#   ylab <- colnames(data)[2]
  
#   # Determine the maximum value across both axes
#   max_limit <- 100 * ceiling(max(data) / 100)
#   lim <- c(0, max_limit)
  
#   # Dynamically calculate tick intervals based on the unified limit
#   ticks <- pretty(lim, n = 5)
  
#   # Plot all data points
#   plot(data, col = 'black', pch = 19, cex = 0.75, main = '', xlim = lim, ylim = lim, 
#        bty = 'n', lwd = 1.2, lty = 1, axes = FALSE, frame = FALSE, 
#        xlab = xlab, ylab = ylab)
  
#   # Identify noise points (cluster == 0)
#   noise_indices <- which(dbscan_result$cluster == 0)
#   if (length(noise_indices) > 0) {
#     noise_points <- data[noise_indices, , drop = FALSE]
    
#     # Highlight noise points
#     points(noise_points, col = "indianred", pch = 19)
    
#     # Add labels to the right of noise points using row identifiers
#     text(noise_points, labels = noise_indices, pos = 4, col = "indianred", cex = 0.75)
#   }
  
#   # Add axes with adaptive ticks
#   axis(1, at = ticks, tcl = -0.2)  # x-axis
#   axis(2, at = ticks, tcl = -0.2)  # y-axis
# }

# DBSCAN_analyse <- function(data, height=5, width=5) {

#   # Use kNNdistplot to select eps
#   dev.new(width=width, height=height, noRStudioGD=TRUE)
#   kNNdistplot2(data, k=2)

#   eps <- NA
#   proceed <- FALSE

#   # Loop until the user is happy with the abline position
#   while (!proceed) {
#     # Prompt for eps value
#     while (is.na(eps)) {
#       cat('\nEnter eps: ')
#       eps <- as.numeric(readLines(n=1))
#       if (is.na(eps)) {
#         cat('\nInvalid input. Please enter a numeric value.\n')
#       }
#     }

#     # Get current x-axis limits
#     x_limits <- par("usr")[1:2]  # This gets the x-axis limits from the plot (xmin and xmax)

#     # Draw the abline with the current eps value, restricted to x-axis limits
#     segments(x0=x_limits[1], y0=eps, x1=x_limits[2], y1=eps, col="indianred", lty=3)

#     # Ask user if they're happy with the line
#     cat('\nAre you happy with the position of the line (y/n)? ')
#     response <- tolower(readLines(n=1))

#     # Check if user is happy
#     if (response == 'y') {
#       proceed <- TRUE
#     } else {
#       # Reset eps and prompt again
#       eps <- NA
#       cat('\ntry again...\n')
#     }
#   }

#   dev.off()

#   # Apply DBSCAN on the dataset
#   dbscan_result <- dbscan::dbscan(data, eps=eps, minPts=3)

#   DBSCAN_plot(data, dbscan_result, height=height, width=width)
# }

create_test_output <- function(parameter, test_result) {
  # Create an empty data frame for mixed types
  output_matrix <- data.frame(
    parameter = parameter,  # Add the parameter name
    test = as.character(test_result$method),
    alternative = as.character(test_result$alternative),
    W = as.numeric(test_result$statistic),
    p.value = as.numeric(test_result$p.value),
    stringsAsFactors = FALSE  # Avoid factors for character columns
  )
  
  return(output_matrix)
}

create_art_output <- function(formula, data, parameter=NULL, dp=5) {
  
  model <- NULL  # Clear any old model

  temp_model <- tryCatch(
    {
      # Use `try(..., silent = TRUE)` to suppress messages
      res <- try(art(formula = formula, data = data), silent = TRUE)
      
      # 2) Check if `res` is a `try-error`; if so, raise an error for `tryCatch` to handle
      if (inherits(res, "try-error")) {
        # Extract the error message from the `try-error` object
        error_msg <- attr(res, "condition")$message
        stop(error_msg)
      }
      
      # If successful, `res` is the fitted 'art' model
      res
    },
    error = function(e) {
      e
    }
  )

  if (!inherits(temp_model, "error")) {
    # Assign the successful model to 'model'
    model <- temp_model
    anova_res <- anova(model)
  }else{
    # Return a single-row data frame with the error
    return(data.frame(
      parameter = parameter,
      test      = 'Analysis of Variance of Aligned Rank Transformed Data', 
      factor    = '',
      df        = '',
      dfe       = '',
      `F value` = '',
      `Pr(>F)`  = '',
      Note      = gsub("\\.$", "", error_msg),   # Store the error text here
      stringsAsFactors = FALSE,
      check.names = FALSE
    ))
  }
  
  anova_res <- anova(model)
  
  # Figure out if it's an anova.art object or a regular anova
  if ('anova.art' %in% class(anova_res)) {
    header_text <- 'Analysis of Variance of Aligned Rank Transformed Data'
  } else {
    header_text <- 'ANOVA'
  }
  
  # Convert the ANOVA result to a data frame
  anova_df <- as.data.frame(anova_res)
  factor_names <- rownames(anova_df)
  
  # Build the output data frame with multiple rows for each factor
  # and only the first row containing the parameter & test name
  output_matrix <- data.frame(
    parameter = c(parameter, rep("", nrow(anova_df) - 1)),  # Only first row has 'parameter'
    test      = c(header_text, rep("", nrow(anova_df) - 1)), # Only first row has 'test'
    factor    = factor_names,                                # e.g., celltype, condition, etc.
    df        = round(anova_df$Df, digits=dp),               # Degrees of freedom
    dfe       = round(anova_df$Df.res, digits=dp),           # Residual degrees of freedom
    `F value` = round(anova_df$`F value`, digits=dp),        # F statistic
    `Pr(>F)`  = round(anova_df$`Pr(>F)`, digits=dp),         # p-value
    Note      = "",                                          # Blank note for successful runs
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
  
  return(output_matrix)
}

# if ms and pA then output would be fC so 1e3 corrects to pC
trap_fun <- function(x, y) {
  sum(diff(x) * (head(y, -1) + tail(y, -1)) / 2) / 1e3
}

fit_plot2 <- function(traces, func=product2, xlab='time (ms)', ylab='PSC amplitude (pA)', xlim=NULL, ylim=NULL, main='', bl=NULL, lwd=1.2, filter=FALSE, width=5, height=5, bg='transparent', filename='trace.svg', save=FALSE) {
  
  if (save) {
    # Open SVG device
    svg(file=filename, width=width, height=height, bg=bg)
  } else {
    dev.new(width=width, height=height, noRStudioGD=TRUE)
  }

  plot(traces$x, traces$y, col='gray', xlab=xlab, ylab=ylab, xlim=xlim, ylim=ylim, type='l', bty='l', las=1, lwd=lwd, main=main)
  lines(traces$x, traces$yfilter, col='black', type='l', lwd=lwd)
  if (!is.null(bl)) abline(v=bl, col='black', lwd=lwd, lty=3)

  if (save) {
    dev.off()
  }
}

single_egs <- function(x, y, sign=-1, xlim=NULL, ylim=NULL, lwd=1, show_text=FALSE, height=4, width=2.5, xbar=100, ybar=50,  color='#4C77BB', filename='trace1.svg', bg='transparent', save=FALSE) {
  
  if (save) {
    svg(file=filename, width=width, height=height, bg=bg)
  } else {
    dev.new(width=width, height=height, noRStudioGD=TRUE)
  }

  if (is.null(ylim)) ylim <- if (sign == 1) c(0, max(y)) else c(-max(-y), 0)

  if (is.null(xlim)) xlim <- c(min(x), max(x))
  idx1 <- which.min(abs(x - xlim[1]))
  idx2 <- which.min(abs(x - xlim[2]))

  plot(x[idx1:idx2], y[idx1:idx2], type='l', col=color, xlim=xlim, ylim=ylim, bty='n', lwd=lwd, lty=1, axes=FALSE, frame=FALSE, xlab='', ylab='')

  #  scale bar lengths and ybar position
  ybar_start <- min(ylim) + (max(ylim) - min(ylim)) / 20
  
  # Add scale bars at the bottom right
  x_start <- max(xlim) - xbar - 50
  y_start <- ybar_start
  x_end <- x_start + xbar
  y_end <- y_start + ybar
  
  # Draw the scale bars
  segments(x_start, y_start, x_end, y_start, lwd=lwd, col='black')
  segments(x_start, y_start, x_start, y_end, lwd=lwd, col='black')
  
  # Add labels to the scale bars
  if (show_text) {
    text(x = (x_start + x_end) / 2, y = y_start - ybar / 20, labels = paste(xbar, 'ms'), adj = c(0.5, 1))
    text(x = x_start - xbar / 4, y = (y_start + y_end) / 2,  labels = paste(ybar, 'pA'), adj = c(0.5, 0.5), srt = 90)
  }
    
  if (save) {
    dev.off()
  }
}
charge_transfer_fun <- function(x, y, fc=300, dt=0.1, baseline=40, filter=TRUE, width=5, height=5, bg='transparent', filename='trace.svg', showplot=FALSE, save=FALSE){

  idx <- baseline / dt
  y <- y - mean(y[1:idx])

  if (filter) {
    fs=1 / dt * 1000; bf <- butter(2, fc / (fs / 2), type='low')
    yfilter <- signal::filter(bf, y)
  } else {
    ind=1
    yfilter=y
  }

  if (showplot){
    traces <- data.frame(x=x, y=y, yfilter=yfilter)
    fit_plot2(traces, width=width, height=height, bg=bg, filename=filename, save=save)
  }

  trap_fun(x, yfilter)

}

charge_fun <- function(data_list, condition='Control', fc=300, dt=0.1, baseline=10, filter=TRUE, showplot=TRUE) {
  sapply(data_list, function(df) {
    if (condition %in% colnames(df)) {
      x <- df$time
      y <- df[[condition]]  # Select column dynamically
      auc <- charge_transfer_fun(x, y, fc=fc, dt=dt, baseline=baseline, filter=filter, showplot=showplot)  # Compute charge transfer
      return(auc)
    } else {
      return(NA)  # Return NA if column is missing
    }
  })
}

# Function to calculate the half-width of y(t) = A * (exp(-t / tau2) - exp(-t / tau1))
half_width <- function(A, tau1, tau2, limit=100) {
  # Define the response function y(t)
  y <- function(t) {
    A * (exp(-t / tau2) - exp(-t / tau1))
  }
  
  # Find the peak value of y(t) and the corresponding time (t_peak)
  opt <- optimize(y, interval = c(0, 100), maximum = TRUE)
  t_peak <- opt$maximum
  y_max <- opt$objective
  
  # Define the target half-maximum value
  y_half_max <- y_max / 2
  
  # Define a function for the difference from half-maximum
  half_max_eq <- function(t) {
    y(t) - y_half_max
  }
  
  # Solve for t1 (before the peak) where y(t) = y_max / 2
  t1 <- uniroot(half_max_eq, interval = c(0, t_peak))$root
  
  # Solve for t2 (after the peak) where y(t) = y_max / 2
  t2 <- uniroot(half_max_eq, interval = c(t_peak, limit))$root
  
  # Calculate the half-width
  half_width <- t2 - t1
  
  # Return results as a numeric vector
  c(t1 = t1, t2 = t2, half_width = half_width)
}


amplifier_gain <- function(dataset = NULL, headstage_gain = 0.5, additional_gain = NULL,
                           AD_range = c(-10, 10), AD_bits = 16,
                           dp = 3, tol = 1e-3, VClamp = TRUE) {
  
  if (!is.null(dataset) && !is.null(additional_gain)) {
    amplifier_gain3(dataset = dataset,
                    headstage_gain = headstage_gain,
                    additional_gain = additional_gain,
                    AD_range = AD_range,
                    AD_bits = AD_bits,
                    tol = tol,
                    dp = dp,
                    VClamp = VClamp)
    
  } else if (!is.null(dataset)) {
    amplifier_gain1(dataset = dataset,
                    headstage_gain = headstage_gain,
                    AD_range = AD_range,
                    AD_bits = AD_bits,
                    dp = dp,
                    tol = tol,
                    VClamp = VClamp)
    
  } else {
    amplifier_gain2(headstage_gain = headstage_gain,
                    additional_gain = additional_gain,
                    AD_range = AD_range,
                    AD_bits = AD_bits,
                    VClamp = VClamp)
  }
}

dpA_fun <- function(dataset, tol=1e-2){
  sapply(1:dim(dataset)[2], function(ii){
    vec <- diff(sort(unique(dataset[,ii])))
    vec <- vec[vec>tol]
    min(vec)
    }
  )
}

amplifier_gain1 <- function(dataset, headstage_gain=0.5, AD_range=c(-10, 10), AD_bits=16, dp=3, tol=1e-3, VClamp=TRUE) {
  
  digitiser_range <- abs(diff(AD_range))
  min_A_D <- rep(AD_range[1], ncol(dataset))
  max_A_D <- rep(AD_range[2], ncol(dataset))
  
  if (VClamp) {
    dpA <- dpA_fun(dataset=dataset, tol=tol)
    recording_range <- dpA * 2^AD_bits
    min_recording <- -recording_range / 2
    max_recording <-  recording_range / 2
    final_gain <- digitiser_range * 1e3 / recording_range
    additional_gain <- final_gain / headstage_gain
    
    output <- data.frame(
      'R GOhms' = rep(headstage_gain, ncol(dataset)),
      'gain mV/pA' = rep(headstage_gain, ncol(dataset)),
      'additional gain' = round(additional_gain, dp),
      'final gain mV/pA' = round(final_gain, dp),
      'min A-D board V' = min_A_D,
      'max A-D board V' = max_A_D,
      'A-D board range V' = rep(digitiser_range, ncol(dataset)),
      'A-D bits' = rep(AD_bits, ncol(dataset)),
      'min recording pA' = min_recording,
      'max recording pA' = max_recording,
      'recording range pA' = recording_range,
      'digitisation pA/bit' = dpA,
      check.names = FALSE
    )
    
  } else {
    dV <- sapply(1:dim(dataset)[2], function(ii) min(diff(sort(unique(dataset[, ii])))) )
    recording_range <- dV * 2^AD_bits
    min_recording <- -recording_range / 2
    max_recording <-  recording_range / 2
    final_gain <- digitiser_range * 1e3 / recording_range
    additional_gain <- final_gain / headstage_gain
    
    output <- data.frame(
      'R GOhms' = rep(headstage_gain, ncol(dataset)),
      'gain V/V' = rep(headstage_gain, ncol(dataset)),
      'additional gain' = round(additional_gain, dp),
      'final gain V/V' = round(final_gain, dp),
      'min A-D board V' = min_A_D,
      'max A-D board V' = max_A_D,
      'A-D board range V' = rep(digitiser_range, ncol(dataset)),
      'A-D bits' = rep(AD_bits, ncol(dataset)),
      'min recording mV' = min_recording,
      'max recording mV' = max_recording,
      'recording range mV' = recording_range,
      'digitisation mV/bit' = dV,
      check.names = FALSE
    )
  }
  
  return(output)
}

amplifier_gain2 <- function(headstage_gain=0.5, additional_gain=20, AD_range=c(-10, 10), AD_bits=16, VClamp=TRUE) {
  
  if (length(additional_gain) == 1 && length(headstage_gain) > 1) {
    additional_gain <- rep(additional_gain, length(headstage_gain))
  }
  
  if (length(headstage_gain) != length(additional_gain)) {
    stop("if 'additional_gain' is a vector, it must be the same length as 'headstage_gain'")
  }
  
  min_A_D <- AD_range[1]
  max_A_D <- AD_range[2]
  digitiser_range <- abs(diff(AD_range))
  final_gain <- headstage_gain * additional_gain
  recording_range <- digitiser_range * 1e3 / final_gain
  min_recording <- -recording_range / 2
  max_recording <-  recording_range / 2
  dUnit <- recording_range / 2^AD_bits
  
  if (VClamp) {
    output <- data.frame(
      'R GOhms' = headstage_gain,
      'gain mV/pA' = headstage_gain,
      'additional gain' = additional_gain,
      'final gain mV/pA' = final_gain,
      'min A-D board V' = rep(min_A_D, length(headstage_gain)),
      'max A-D board V' = rep(max_A_D, length(headstage_gain)),
      'A-D board range V' = rep(digitiser_range, length(headstage_gain)),
      'A-D bits' = rep(AD_bits, length(headstage_gain)),
      'min recording pA' = min_recording,
      'max recording pA' = max_recording,
      'recording range pA' = recording_range,
      'digitisation pA/bit' = dUnit,
      check.names = FALSE
    )
  } else {
    output <- data.frame(
      'R GOhms' = headstage_gain,
      'gain V/V' = headstage_gain,
      'additional gain' = additional_gain,
      'final gain V/V' = final_gain,
      'min A-D board V' = rep(min_A_D, length(headstage_gain)),
      'max A-D board V' = rep(max_A_D, length(headstage_gain)),
      'A-D board range V' = rep(digitiser_range, length(headstage_gain)),
      'A-D bits' = rep(AD_bits, length(headstage_gain)),
      'min recording mV' = min_recording,
      'max recording mV' = max_recording,
      'recording range mV' = recording_range,
      'digitisation mV/bit' = dUnit,
      check.names = FALSE
    )
  }
  return(output)
}

amplifier_gain3 <- function(dataset, headstage_gain = 0.5, additional_gain = 20,
                            AD_range = c(-10, 10), AD_bits = 16,
                            tol = 1e-3, dp=3, VClamp = TRUE) {
  
  digitiser_range <- abs(diff(AD_range))
  min_A_D <- rep(AD_range[1], ncol(dataset))
  max_A_D <- rep(AD_range[2], ncol(dataset))
  final_gain <- headstage_gain * additional_gain
  recording_range <- digitiser_range * 1e3 / final_gain
  min_recording <- -recording_range / 2
  max_recording <-  recording_range / 2
  dUnit <- recording_range / 2^AD_bits  # actual theoretical digitisation
  
  if (VClamp) {
    dpA_actual <- dpA_fun(dataset = dataset, tol = tol)
    n <- dUnit / dpA_actual 
    
    output <- data.frame(
      'R GOhms' = rep(headstage_gain, ncol(dataset)),
      'gain mV/pA' = rep(headstage_gain, ncol(dataset)),
      'additional gain' = rep(additional_gain, ncol(dataset)),
      'final gain mV/pA' = final_gain,
      'min A-D board V' = min_A_D,
      'max A-D board V' = max_A_D,
      'A-D board range V' = rep(digitiser_range, ncol(dataset)),
      'A-D bits' = rep(AD_bits, ncol(dataset)),
      'min recording pA' = min_recording,
      'max recording pA' = max_recording,
      'recording range pA' = recording_range,
      'digitisation pA/bit' = dUnit,
      'n' = round(n, dp),
      check.names = FALSE
    )
    
  } else {
    dV_actual <- sapply(1:dim(dataset)[2], function(ii) min(diff(sort(unique(dataset[, ii])))) )
    n <- dUnit / dV_actual
    
    output <- data.frame(
      'R GOhms' = rep(headstage_gain, ncol(dataset)),
      'gain V/V' = rep(headstage_gain, ncol(dataset)),
      'additional gain' = rep(additional_gain, ncol(dataset)),
      'final gain V/V' = final_gain,
      'min A-D board V' = min_A_D,
      'max A-D board V' = max_A_D,
      'A-D board range V' = rep(digitiser_range, ncol(dataset)),
      'A-D bits' = rep(AD_bits, ncol(dataset)),
      'min recording mV' = min_recording,
      'max recording mV' = max_recording,
      'recording range mV' = recording_range,
      'digitisation mV/bit' = dUnit,
      'n' = round(n, dp),
      check.names = FALSE
    )
  }
  
  return(output)
}


fit_limit <- function(y, N=1, dt=0.1, stimulation_time=0, baseline=0, smooth=5, 
                      y_abline=0.1, height=4, width=4, show_plot=FALSE) { 
  # Calculate peak (unused in the remainder but may be important elsewhere)
  peak <- peak.fun(y=y, dt=dt, stimulation_time=stimulation_time, 
                   baseline=baseline, smooth=smooth)
  
  ind1 <- as.integer((stimulation_time - baseline)/dt)
  ind2 <- as.integer(stimulation_time/dt)
  y2plot <- y - mean(y[ind1:ind2])
  
  Y <- y2plot[ind1:length(y2plot)]
  X <- seq(0, dt * (length(Y) - 1), by = dt)

  out <- abline_fun(X, Y, N=N, y_abline=y_abline) 
  A_abline <- out[1]
  avg_t.abline <- out[2]
  avg_t.abline <- if (is.na(avg_t.abline)) max(X) else avg_t.abline

  if (show_plot) {
    dev.new(width=width, height=height, noRStudioGD=TRUE)
    
    plot(X, Y, col='indianred', xlab='time (ms)', type='l', bty='l', las=1, main='')
    abline(h = 0, col = 'black', lwd = 1, lty = 3)
    
    # Get the left and bottom boundaries of the plot
    left_axis <- par("usr")[1]
    bottom_axis <- par("usr")[3]
    
    # Draw horizontal dotted line from the left to avg_t.abline at height A_abline
    lines(c(left_axis, avg_t.abline), c(A_abline, A_abline), col = 'black', lwd = 1, lty = 3)
    
    # Draw vertical dotted line from avg_t.abline down to the bottom of the plot
    lines(c(avg_t.abline, avg_t.abline), c(A_abline, bottom_axis), col = 'black', lwd = 1, lty = 3)
    
    # Determine index for labeling and add labels
    ind3 <- as.integer(avg_t.abline/dt)
    text(x = max(X[ind1:ind3]) * 1.05, y = A_abline * 1.2, 
         labels = paste0(y_abline * 100, " %"), pos = 4, cex = 0.6)
    text(x = max(X[ind1:ind3]) * 1.05, y = bottom_axis * 0.95, 
         labels = paste0(avg_t.abline, " ms"), pos = 4, cex = 0.6)
    
    # Prompt user for the range of x to use for nFIT
    x_limit <- NA
    while (is.na(x_limit)) {
      cat('\nEnter the upper limit for time to use in nFIT (e.g., 400 ms): ')
      x_limit <- as.numeric(readLines(n = 1))
      if (is.na(x_limit)) {
        cat('\nInvalid input. Please enter a numeric value.\n')
      }
    }
    dev.off()
  } else {
    x_limit <- avg_t.abline
  }
  
  return(x_limit)
}

smooth_moving_avg <- function(y, n = 5) {
  y_length <- length(y)
  result <- rep(NA, y_length)
  
  for (i in 1:y_length) {
    # Determine the start and end indices for the window
    start_idx <- max(1, i - floor(n / 2))
    end_idx <- min(y_length, i + floor(n / 2))
    
    # Calculate the mean for the current window
    result[i] <- mean(y[start_idx:end_idx], na.rm = TRUE)
  }
  
  return(result)
}

downsample_fun <- function(data, ds) {
  if (is.vector(data)) {
    data[seq(1, length(data), by = ds)]
  } else {
    data[seq(1, nrow(data), by = ds), , drop = FALSE]
  }
}


# determine_tmax2
determine_tmax2 <- function(y, N=1, dt=0.1, stimulation_time=0, baseline=0, smooth=5, lwd=1.2, cex=0.6,
  tmax=NULL, y_abline=0.1, xbar=50, ybar=50, xbar_lab='ms', ybar_lab='pA') {
  if (is.null(tmax)) {
    # Calculate peak information (assumes peak.fun and abline_fun are defined)
    peak <- peak.fun(y=y, dt=dt, stimulation_time=stimulation_time, baseline=baseline, smooth=smooth)
    
    ind1 <- as.integer((stimulation_time - baseline) / dt)
    ind2 <- as.integer(stimulation_time / dt)
    y2plot <- y - mean(y[ind1:ind2])
    
    # Prepare data for plotting
    Y <- y2plot[ind1:length(y2plot)]
    X <- seq(0, dt * (length(Y) - 1), by=dt)
    
    out <- abline_fun(X, Y, N=N, y_abline=y_abline)
    A_abline <- out[1]
    avg_t.abline <- if (is.na(out[2])) max(X) else out[2]
    
    # Draw the main plot (without axes)
    plot(X, Y, col='indianred', type='l', axes=FALSE, xlab='', ylab='', lwd=lwd, main='', bty='n')
    
    # draw ablines and add text
    usr <- par('usr')
    left_axis <- usr[1]
    bottom_axis <- usr[3]
    lines(c(min(X), max(X)), c(0, 0), col='black', lwd=lwd, lty=3)
    lines(c(min(X), avg_t.abline), c(A_abline, A_abline), col='black', lwd=lwd, lty=3)
    lines(c(avg_t.abline, avg_t.abline), c(A_abline, bottom_axis), col='black', lwd=lwd, lty=3)
    ind3 <- as.integer(avg_t.abline / dt)
    text(x=max(X[ind1:ind3]) * 1.05, y=A_abline * 1.2, labels=paste0(y_abline * 100, ' %'), pos=4, cex=cex)
    text(x=max(X[ind1:ind3]) * 1.05, y=bottom_axis * 0.95, labels=paste0(avg_t.abline, ' ms'), pos=4, cex=cex)
    stim_index <- round(baseline / dt) + 1
    if (stim_index > length(X)) stim_index <- length(X)
    points(X[stim_index], Y[stim_index], pch=8, col='black', cex=1)
    x_offset <- 0.02 * diff(range(X))
    text(x=X[stim_index] + x_offset, y=Y[stim_index], labels='stim', pos=4, col='darkgray', cex=cex)
    x_limit <- avg_t.abline

  } else {
    x_limit <- tmax
  }
  
  x_limit <- x_limit + stimulation_time - baseline
    
  return(x_limit)
}

# fit_plot3
fit_plot3 <- function(traces, func=product2, lwd=1.2, cex=0.6, filter=FALSE, xbar=50, ybar=50, 
  xbar_lab='ms', ybar_lab='pA') {
  plot(traces$x, traces$y, col='gray', type='l', axes=FALSE, xlab='', ylab='',
       bty='n', lwd=lwd)
  if (filter && !is.null(traces$yfilter)) {
    lines(traces$x, traces$yfilter, col='black', type='l', lwd=lwd)
  }
  lines(traces$x, traces$yfit, col='indianred', lty=3, lwd=2 * lwd)
  if (identical(func, product2) || identical(func, product2N)) {
    lines(traces$x, traces$yfit1, col='#4C78BC', lty=3, lwd=2 * lwd)
    lines(traces$x, traces$yfit2, col='#CA92C1', lty=3, lwd=2 * lwd)
  }
  if (identical(func, product3) || identical(func, product3N)) {
    lines(traces$x, traces$yfit1, col='#F28E2B', lty=3, lwd=2 * lwd)
    lines(traces$x, traces$yfit2, col='#4C78BC', lty=3, lwd=2 * lwd)
    lines(traces$x, traces$yfit3, col='#CA92C1', lty=3, lwd=2 * lwd)
  }
  if (!is.null(traces$bl)) {
    abline(v=traces$bl, col='black', lwd=lwd, lty=3)
  }
  
  # scale bars
  usr <- par('usr')
  x_range <- usr[1:2]
  y_range <- usr[3:4]
  ybar_start <- y_range[1] + (y_range[2] - y_range[1]) / 20
  x_start <- x_range[2] - xbar - 50
  y_start <- ybar_start
  x_end <- x_start + xbar
  y_end <- y_start + ybar
  
  segments(x_start, y_start, x_end, y_start, lwd=lwd, col='black')
  segments(x_start, y_start, x_start, y_end, lwd=lwd, col='black')
  text(x=(x_start + x_end) / 2, y=y_start - ybar / 20, 
       labels=paste(xbar, xbar_lab), adj=c(0.5, 1), cex=cex)
  text(x=x_start - xbar / 4, y=(y_start + y_end) / 2, 
       labels=paste(ybar, ybar_lab), srt=90, adj=c(0.5, 0.5), cex=cex)
}

# drawPlot2
drawPlot2 <- function(traces, func=product2N, lwd=1.2, cex=1, filter=FALSE, xbar=50, ybar=50, 
                      xbar_lab='ms', ybar_lab='pA') {
  fit_plot3(traces=traces, func=func, lwd=lwd, cex=cex, filter=filter,
            xbar=xbar, ybar=ybar, xbar_lab=xbar_lab, ybar_lab=ybar_lab)
}


# MCwilcox <- function(formula, df, alternative = 'two.sided',
#                      exact = NULL, na_rm_subjects = TRUE, p_adjust = 'holm') {
#   f_str <- deparse(formula)
#   has_error <- grepl('Error', f_str)

#   if (has_error) {
#     err_part <- sub('.*Error\\((.*)\\).*', '\\1', f_str)
#     subject_var <- strsplit(err_part, '/')[[1]][1]
#     subject_var <- gsub('[[:space:]]', '', subject_var)
#     main_formula_str <- sub('\\+\\s*Error\\(.*\\)', '', f_str)
#     main_formula <- as.formula(main_formula_str)
#   } else {
#     subject_var <- NULL
#     main_formula <- formula
#   }

#   response_var <- all.vars(formula(main_formula))[1]
#   predictors <- all.vars(formula(main_formula))[-1]

#   if (na_rm_subjects && !is.null(subject_var)) {
#     df <- df[ !ave(is.na(df[[response_var]]), df[[subject_var]], FUN = any), ]
#   }

#   if (length(predictors) < 1) {
#     stop('Formula must contain at least one predictor for comparisons')
#   }

#   paired_var <- predictors[1]
#   unpaired_var <- if (length(predictors) > 1) predictors[2] else NULL
#   results <- list()

#   ### Paired ###
#   if (!is.null(subject_var) && !is.null(unpaired_var)) {
#     if (is.factor(df[[paired_var]])) {
#       levels_pair <- levels(df[[paired_var]])
#     } else {
#       levels_pair <- sort(unique(df[[paired_var]]))
#     }

#     for (lev in levels_pair) {
#       subset_df <- df[df[[paired_var]] == lev, ]
#       if (is.factor(subset_df[[unpaired_var]])) {
#         levels_unpair <- levels(subset_df[[unpaired_var]])
#       } else {
#         levels_unpair <- sort(unique(subset_df[[unpaired_var]]))
#       }
#       if (length(levels_unpair) < 2) next

#       for (i in seq_len(length(levels_unpair) - 1)) {
#         lev1 <- levels_unpair[i]
#         lev2 <- levels_unpair[i + 1]
#         d1 <- subset_df[subset_df[[unpaired_var]] == lev1, ]
#         d2 <- subset_df[subset_df[[unpaired_var]] == lev2, ]
#         common_subj <- intersect(d1[[subject_var]], d2[[subject_var]])
#         d1 <- d1[d1[[subject_var]] %in% common_subj, ]
#         d2 <- d2[d2[[subject_var]] %in% common_subj, ]
#         d1 <- d1[order(d1[[subject_var]]), ]
#         d2 <- d2[order(d2[[subject_var]]), ]
#         y1 <- d1[[response_var]]
#         y2 <- d2[[response_var]]

#         if (length(y1) > 0 && length(y1) == length(y2)) {
#           test <- wilcox.test(y1, y2, paired = TRUE, alternative = alternative, exact = exact)
#           stat_name <- if (!is.null(names(test$statistic))) names(test$statistic) else NA
#           stat_value <- as.numeric(test$statistic)
#           results[[length(results) + 1]] <- data.frame(
#             comparison   = paste('within', paired_var, lev, '(paired)'),
#             contrast     = paste(unpaired_var, ':', lev1, 'vs', lev2),
#             n            = min(sum(!is.na(y1)), sum(!is.na(y2))),
#             test         = test$method,
#             alternative  = test$alternative,
#             `test stat`  = stat_name,
#             stat         = stat_value,
#             `p value`    = test$p.value,
#             family       = 'paired',
#             stringsAsFactors = FALSE,
#             check.names = FALSE
#           )
#         }
#       }
#     }
#   }

#   ### Unpaired across levels of paired_var within unpaired_var
#   if (!is.null(unpaired_var)) {
#     if (is.factor(df[[unpaired_var]])) {
#       levels_unpair_all <- levels(df[[unpaired_var]])
#     } else {
#       levels_unpair_all <- sort(unique(df[[unpaired_var]]))
#     }

#     for (lev in levels_unpair_all) {
#       subset_df <- df[df[[unpaired_var]] == lev, ]
#       if (is.factor(subset_df[[paired_var]])) {
#         groups_pair <- levels(subset_df[[paired_var]])
#       } else {
#         groups_pair <- sort(unique(subset_df[[paired_var]]))
#       }
#       if (length(groups_pair) < 2) next
#       d1 <- subset_df[subset_df[[paired_var]] == groups_pair[1], ]
#       d2 <- subset_df[subset_df[[paired_var]] == groups_pair[2], ]

#       test <- wilcox.test(d1[[response_var]], d2[[response_var]], paired = FALSE, alternative = alternative, exact = exact)
#       stat_name <- if (!is.null(names(test$statistic))) names(test$statistic) else NA
#       stat_value <- as.numeric(test$statistic)

#       results[[length(results) + 1]] <- data.frame(
#         comparison = paste('within', unpaired_var, lev, '(unpaired)'),
#         contrast   = paste(paired_var, ':', groups_pair[1], 'vs', groups_pair[2]),
#         n          = paste(sum(!is.na(d1[[response_var]])), 'vs', sum(!is.na(d2[[response_var]]))),
#         test       = test$method,
#         alternative= test$alternative,
#         `test stat`= stat_name,
#         stat       = stat_value,
#         `p value`  = test$p.value,
#         family     = 'unpaired',
#         stringsAsFactors = FALSE,
#         check.names = FALSE
#       )
#     }
#   }

#   ### Unpaired across levels of unpaired_var within paired_var, if no subjects
#   if (is.null(subject_var) && !is.null(unpaired_var)) {
#     if (is.factor(df[[paired_var]])) {
#       levels_pair_all <- levels(df[[paired_var]])
#     } else {
#       levels_pair_all <- sort(unique(df[[paired_var]]))
#     }

#     for (lev in levels_pair_all) {
#       subset_df <- df[df[[paired_var]] == lev, ]
#       if (is.factor(subset_df[[unpaired_var]])) {
#         groups_unpair <- levels(subset_df[[unpaired_var]])
#       } else {
#         groups_unpair <- sort(unique(subset_df[[unpaired_var]]))
#       }
#       if (length(groups_unpair) < 2) next
#       d1 <- subset_df[subset_df[[unpaired_var]] == groups_unpair[1], ]
#       d2 <- subset_df[subset_df[[unpaired_var]] == groups_unpair[2], ]

#       test <- wilcox.test(d1[[response_var]], d2[[response_var]], paired = FALSE, alternative = alternative, exact = exact)
#       stat_name <- if (!is.null(names(test$statistic))) names(test$statistic) else NA
#       stat_value <- as.numeric(test$statistic)

#       results[[length(results) + 1]] <- data.frame(
#         comparison = paste('within', paired_var, lev, '(unpaired)'),
#         contrast   = paste(unpaired_var, ':', groups_unpair[1], 'vs', groups_unpair[2]),
#         n          = paste(sum(!is.na(d1[[response_var]])), 'vs', sum(!is.na(d2[[response_var]]))),
#         test       = test$method,
#         alternative= test$alternative,
#         `test stat`= stat_name,
#         stat       = stat_value,
#         `p value`  = test$p.value,
#         family     = 'unpaired',
#         stringsAsFactors = FALSE,
#         check.names = FALSE
#       )
#     }
#   }

#   ### Single-predictor unpaired comparison
#   if (is.null(subject_var) && is.null(unpaired_var)) {
#     if (is.factor(df[[paired_var]])) {
#       groups <- levels(df[[paired_var]])
#     } else {
#       groups <- sort(unique(df[[paired_var]]))
#     }
#     if (length(groups) >= 2) {
#       d1 <- df[df[[paired_var]] == groups[1], ]
#       d2 <- df[df[[paired_var]] == groups[2], ]

#       test <- wilcox.test(d1[[response_var]], d2[[response_var]], paired = FALSE, alternative = alternative, exact = exact)
#       stat_name <- if (!is.null(names(test$statistic))) names(test$statistic) else NA
#       stat_value <- as.numeric(test$statistic)

#       results[[length(results) + 1]] <- data.frame(
#         comparison = paste('between', paired_var),
#         contrast   = paste(groups[1], 'vs', groups[2]),
#         n          = paste(sum(!is.na(d1[[response_var]])), 'vs', sum(!is.na(d2[[response_var]]))),
#         test       = test$method,
#         alternative= test$alternative,
#         `test stat`= stat_name,
#         stat       = stat_value,
#         `p value`  = test$p.value,
#         family     = 'unpaired',
#         stringsAsFactors = FALSE,
#         check.names = FALSE
#       )
#     }
#   }

#   out <- if (length(results) == 1) results[[1]] else do.call(rbind, results)

#   out$`p adjusted` <- NA
#   for (fam in unique(out$family)) {
#     idx <- which(out$family == fam)
#     out$`p adjusted`[idx] <- p.adjust(out$`p value`[idx], method = p_adjust)
#   }

#   out$family <- factor(out$family, levels = c("paired", "unpaired"))
#   out <- out[order(
#     out$family,
#     sub("within (\\w+).*", "\\1", out$comparison),
#     suppressWarnings(as.numeric(sub(".*within \\w+ (\\w+) \\(.*", "\\1", out$comparison)))
#   ), ]
#   out$family <- NULL
#   return(out)
# }

MCwilcox <- function(formula, df, alternative = 'two.sided',
                     exact = NULL, na_rm_subjects = TRUE,
                     p_adjust = 'holm') {
  f_str     <- deparse(formula)
  has_error <- grepl('Error', f_str)

  if (has_error) {
    err_part         <- sub('.*Error\\((.*)\\).*', '\\1', f_str)
    subject_var      <- strsplit(err_part, '/')[[1]][1]
    subject_var      <- gsub('[[:space:]]', '', subject_var)
    main_formula_str <- sub('\\+\\s*Error\\(.*\\)', '', f_str)
    main_formula     <- as.formula(main_formula_str)
  } else {
    subject_var  <- NULL
    main_formula <- formula
  }

  response_var <- all.vars(formula(main_formula))[1]
  predictors   <- all.vars(formula(main_formula))[-1]

  # drop subjects with any missing response, if requested
  if (na_rm_subjects && !is.null(subject_var)) {
    df <- df[!ave(is.na(df[[response_var]]), df[[subject_var]], FUN = any), ]
  }

  results <- list()

  # --- single predictor + Error(subject) → paired only ---
  if (!is.null(subject_var) && length(predictors) == 1) {
    wvar <- predictors[1]
    levs <- if (is.factor(df[[wvar]])) levels(df[[wvar]]) else sort(unique(df[[wvar]]))
    if (length(levs) < 2) {
      stop("Need at least two levels of ", wvar, " for paired comparisons")
    }
    for (i in seq_len(length(levs)-1)) {
      a <- levs[i]; b <- levs[i+1]
      d1 <- df[df[[wvar]]==a, ]
      d2 <- df[df[[wvar]]==b, ]
      common <- intersect(d1[[subject_var]], d2[[subject_var]])
      d1 <- d1[d1[[subject_var]] %in% common, ]
      d2 <- d2[d2[[subject_var]] %in% common, ]
      d1 <- d1[order(d1[[subject_var]]), ]
      d2 <- d2[order(d2[[subject_var]]), ]
      y1 <- d1[[response_var]]; y2 <- d2[[response_var]]
      if (length(y1)>0 && length(y1)==length(y2)) {
        test <- wilcox.test(y1, y2, paired = TRUE,
                            alternative = alternative, exact = exact)
        snm <- if (!is.null(names(test$statistic))) names(test$statistic) else NA
        snt <- as.numeric(test$statistic)
        results[[length(results)+1]] <- data.frame(
          parameter    = response_var,
          comparison   = paste('within', wvar, '(paired)'),
          contrast     = paste(a, 'vs', b),
          n            = min(sum(!is.na(y1)), sum(!is.na(y2))),
          test         = test$method,
          alternative  = test$alternative,
          `test stat`  = snm,
          stat         = snt,
          `p value`    = test$p.value,
          family       = 'paired',
          stringsAsFactors = FALSE,
          check.names  = FALSE
        )
      }
    }
    out <- do.call(rbind, results)
    out$`p adjusted` <- p.adjust(out$`p value`, method = p_adjust)
    out$family <- NULL
    return(out)
  }

  if (is.null(subject_var) && length(predictors) == 1) {
    uvar <- predictors[1]
    levs <- if (is.factor(df[[uvar]])) levels(df[[uvar]]) else sort(unique(df[[uvar]]))
    if (length(levs) < 2) {
      stop("Need at least two levels of ", uvar, " for unpaired comparisons")
    }
    for (i in seq_len(length(levs)-1)) {
      a <- levs[i]; b <- levs[i+1]
      d1 <- df[df[[uvar]]==a, response_var]
      d2 <- df[df[[uvar]]==b, response_var]
      test <- wilcox.test(d1, d2, paired=FALSE,
                          alternative = alternative, exact = exact)
      snm <- if (!is.null(names(test$statistic))) names(test$statistic) else NA
      snt <- as.numeric(test$statistic)
      results[[length(results)+1]] <- data.frame(
        parameter    = response_var,
        comparison   = paste('within', uvar, '(unpaired)'),
        contrast     = paste(a, 'vs', b),
        n            = paste(sum(!is.na(d1)), 'vs', sum(!is.na(d2))),
        test         = test$method,
        alternative  = test$alternative,
        `test stat`  = snm,
        stat         = snt,
        `p value`    = test$p.value,
        family       = 'unpaired',
        stringsAsFactors = FALSE,
        check.names  = FALSE
      )
    }
    out <- do.call(rbind, results)
    out$`p adjusted` <- p.adjust(out$`p value`, method = p_adjust)
    out$family <- NULL
    return(out)
  }

  if (length(predictors) < 2) {
    stop('Formula must contain at least one predictor (single unpaired) or two (for paired+unpaired)')
  }
  pv <- predictors[1]; uv <- predictors[2]

  # (A) paired within pv if Error(subject) present
  if (!is.null(subject_var)) {
    lv_p <- if (is.factor(df[[pv]])) levels(df[[pv]]) else sort(unique(df[[pv]]))
    for (lev in lv_p) {
      subdf <- df[df[[pv]]==lev, ]
      lv_u  <- if (is.factor(subdf[[uv]])) levels(subdf[[uv]]) else sort(unique(subdf[[uv]]))
      if (length(lv_u)<2) next
      for (i in seq_len(length(lv_u)-1)) {
        a <- lv_u[i]; b <- lv_u[i+1]
        d1 <- subdf[subdf[[uv]]==a, ]
        d2 <- subdf[subdf[[uv]]==b, ]
        cm <- intersect(d1[[subject_var]], d2[[subject_var]])
        d1 <- d1[d1[[subject_var]] %in% cm, ]
        d2 <- d2[d2[[subject_var]] %in% cm, ]
        d1 <- d1[order(d1[[subject_var]]), ]
        d2 <- d2[order(d2[[subject_var]]), ]
        y1 <- d1[[response_var]]; y2 <- d2[[response_var]]
        if (length(y1)>0 && length(y1)==length(y2)) {
          test <- wilcox.test(y1, y2, paired = TRUE,
                              alternative = alternative, exact = exact)
          snm <- if (!is.null(names(test$statistic))) names(test$statistic) else NA
          snt <- as.numeric(test$statistic)
          results[[length(results)+1]] <- data.frame(
            parameter    = response_var,
            comparison   = paste('within', pv, lev, '(paired)'),
            contrast     = paste(a, 'vs', b),
            n            = min(sum(!is.na(y1)), sum(!is.na(y2))),
            test         = test$method,
            alternative  = test$alternative,
            `test stat`  = snm,
            stat         = snt,
            `p value`    = test$p.value,
            family       = 'paired',
            stringsAsFactors = FALSE,
            check.names  = FALSE
          )
        }
      }
    }
  }

  # (B) unpaired across uv
  lv_u_all <- if (is.factor(df[[uv]])) levels(df[[uv]]) else sort(unique(df[[uv]]))
  for (lev in lv_u_all) {
    subdf <- df[df[[uv]]==lev, ]
    gp    <- if (is.factor(subdf[[pv]])) levels(subdf[[pv]]) else sort(unique(subdf[[pv]]))
    if (length(gp)<2) next
    d1 <- subdf[subdf[[pv]]==gp[1], ]
    d2 <- subdf[subdf[[pv]]==gp[2], ]
    test <- wilcox.test(d1[[response_var]], d2[[response_var]],
                        paired=FALSE,
                        alternative = alternative, exact = exact)
    snm <- if (!is.null(names(test$statistic))) names(test$statistic) else NA
    snt <- as.numeric(test$statistic)
    results[[length(results)+1]] <- data.frame(
      parameter    = response_var,
      comparison   = paste('within', uv, lev, '(unpaired)'),
      contrast     = paste(gp[1], 'vs', gp[2]),
      n            = paste(sum(!is.na(d1[[response_var]])),
                           'vs',
                           sum(!is.na(d2[[response_var]]))),
      test         = test$method,
      alternative  = test$alternative,
      `test stat`  = snm,
      stat         = snt,
      `p value`    = test$p.value,
      family       = 'unpaired',
      stringsAsFactors = FALSE,
      check.names  = FALSE
    )
  }

  # (C) if no subject_var, also do reverse unpaired within pv
  if (is.null(subject_var)) {
    lv_p_all <- if (is.factor(df[[pv]])) levels(df[[pv]]) else sort(unique(df[[pv]]))
    for (lev in lv_p_all) {
      subdf <- df[df[[pv]]==lev, ]
      gu    <- if (is.factor(subdf[[uv]])) levels(subdf[[uv]]) else sort(unique(subdf[[uv]]))
      if (length(gu)<2) next
      d1 <- subdf[subdf[[uv]]==gu[1], ]
      d2 <- subdf[subdf[[uv]]==gu[2], ]
      test <- wilcox.test(d1[[response_var]], d2[[response_var]],
                          paired=FALSE,
                          alternative = alternative, exact = exact)
      snm <- if (!is.null(names(test$statistic))) names(test$statistic) else NA
      snt <- as.numeric(test$statistic)
      results[[length(results)+1]] <- data.frame(
        parameter    = response_var,
        comparison   = paste('within', pv, lev, '(unpaired)'),
        contrast     = paste(gu[1], 'vs', gu[2]),
        n            = paste(sum(!is.na(d1[[response_var]])),
                             'vs',
                             sum(!is.na(d2[[response_var]]))),
        test         = test$method,
        alternative  = test$alternative,
        `test stat`  = snm,
        stat         = snt,
        `p value`    = test$p.value,
        family       = 'unpaired',
        stringsAsFactors = FALSE,
        check.names  = FALSE
      )
    }
  }

  # assemble and finalize
  out <- do.call(rbind, results)
  out$`p adjusted` <- NA
  for (fam in unique(out$family)) {
    idx <- which(out$family == fam)
    out$`p adjusted`[idx] <- p.adjust(out$`p value`[idx], method = p_adjust)
  }

  # reorder exactly as before
  out$family <- factor(out$family, levels = c("paired","unpaired"))
  out <- out[ order(
    out$family,
    sub("within (\\w+).*", "\\1", out$comparison),
    suppressWarnings(as.numeric(sub(".*within \\w+ (\\w+) \\(.*", "\\1", out$comparison)))
  ), ]
  out <- out[ order(out$comparison), ]

  # **always** drop the family column
  out$family <- NULL
  rownames(out) <- seq_len(nrow(out))
  return(out)
}

save_graph <- function(svg_path, filename='graph1.svg', width=6, height=4, bg="transparent") {
  old_wd <- getwd()
  setwd(svg_path)
  dev.copy(svg, file=filename, width=width, height=height, bg=bg)
  dev.off()
  setwd(old_wd)
}


# 
analyseABFtk <- function() {

  cexVar <- if (Sys.info()["sysname"] == "Darwin") tclVar("0.6") else tclVar("1.0")

  experimentVar  <- tclVar("voltage clamp")
  unitVar        <- tclVar("")
  dataColVar     <- tclVar("")
  dtVar          <- tclVar("")
  ntracesVar     <- tclVar("")

  is.tkwin <- function(widget) {
    tryCatch({
      tclvalue(tkwinfo("exists", widget)) == "1"
    }, error = function(e) FALSE)
  }


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
    if (is.null(averaged_data) || length(averaged_data) == 0) {
      tkinsert(consoleText, 'end', 'No averaged data available.\n')
      tkyview.moveto(consoleText, 1.0)
      return()
    }

    df <- as.data.frame(do.call(cbind, averaged_data))
    colnames(df) <- as.character(seq_len(length(averaged_data)))

    download_folder <- tclvalue(folderPathVar)
    if (nchar(download_folder) == 0) {
      tkinsert(consoleText, 'end', 'No folder selected for download.\n')
      tkyview.moveto(consoleText, 1.0)
      return()
    }

    file_path <- tclvalue(tkgetSaveFile(
      initialdir = download_folder,
      defaultextension = ".csv",
      initialfile = "average.csv"
    ))

    if (nchar(file_path) == 0) {
      return()
    }

    # only reached if file_path is valid
    write.csv(df, file = file_path, row.names = FALSE)
    tkinsert(consoleText, 'end', paste0('Saved to: ', file_path, '\n'))
    tkyview.moveto(consoleText, 1.0)
  }

  abf_averages <- function(datasets, baseline = 100, stimulation_time = 350, traces2average = NULL, dataCol = 1, ylim = NULL, xlim = NULL, 
    color = 'darkgray', xbar = 100, ybar = 50, width = 5.25, height = 2.75, cex=0.6, save = FALSE, plotIt = TRUE) {
    
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
        egs_plot(x = time[[ii]], y = responses0_mean[[ii]], color = 'darkgray',
                   show_bar = FALSE, cex=cex, show_text = FALSE)
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
  # common plot settings
  tk_par_settings <- function() {
    par(mar = c(4, 5, 2, 1) + 0.1,
        mgp = c(2.5, 0.5, 0),
        tcl = -0.2)
  }

  # global plot size settings
if (Sys.info()["sysname"] == "Darwin") {
  graph_width  <<- 500
  graph_height <<- 300
  graph_hscale <<- 1.0
  graph_vscale <<- 1.0
} else {
  graph_width  <<- 750
  graph_height <<- 450
  graph_hscale <<- 1.5
  graph_vscale <<- 1.5
}

#### review_master_recordings: shows each trace with Accept/Reject/Undo buttons below ####
review_master_recordings <- function() {
  if (is.null(master_abf)) {
    tkinsert(consoleText, 'end', "No master ABF data available. Please load data first.\n")
    tkyview.moveto(consoleText, 1.0)
    return()
  }
  total_traces <<- length(master_abf$data)
  current_trace <<- 1L
  current_group_selected <<- integer(0)
  groups_list <<- list()

  children <- as.character(tkwinfo('children', plotPanel))
  if (length(children)) sapply(children, function(ch) tcl("destroy", ch))

  reviewFrame <<- tkframe(plotPanel)
  tkgrid(reviewFrame, row = 0, column = 0, sticky = 'nsew')
  tkgrid.columnconfigure(reviewFrame, 0, weight = 1)
  tkgrid.columnconfigure(reviewFrame, 1, weight = 1)
  tkgrid.columnconfigure(reviewFrame, 2, weight = 1)
  tkgrid.rowconfigure(reviewFrame,    1, weight = 1)

  infoLabel <<- tklabel(reviewFrame, text = paste('Trace', current_trace, 'of', total_traces))
  tkgrid(infoLabel, row = 0, column = 1, pady = 5)

  plotWrapper <- tkframe(reviewFrame, height = graph_height, width = graph_width)
  tkgrid(plotWrapper, row = 1, column = 1, sticky = 'nsew')
  tkgrid.rowconfigure(plotWrapper,    0, weight = 1)
  tkgrid.columnconfigure(plotWrapper, 0, weight = 1)

  reviewPlot <<- tkrplot(plotWrapper, fun = function() {
    tk_par_settings()
    cex <- as.numeric(tclvalue(cexVar))
    par(cex.lab = cex, cex.axis = cex, cex.main = cex)

    mat <- master_abf$data[[current_trace]]
    dt  <- master_abf$samplingIntervalInSec * 1000
    time <- seq(0, by = dt, length.out = nrow(mat))
    dc <- as.numeric(tclvalue(dataColVar))
    if (is.na(dc) || dc < 1 || dc > ncol(mat)) dc <- 1
    y  <- mat[, dc]
    stim_time <- as.numeric(tclvalue(stimTimeVar))
    stim_y    <- y[which.min(abs(time - stim_time))]

    plot(time, y, type = 'l', col = 'darkgray',
         xlab = 'time (ms)', ylab = tclvalue(unitVar),
         xlim = smart_axis_limits(time),
         ylim = smart_axis_limits(y),
         axes = FALSE, bty = 'l')
    axis(1); axis(2, las = 1)
    points(stim_time, stim_y, pch = 8, col = 'black')
    text(stim_time, stim_y, labels = 'stim', pos = 3, cex = cex)
  }, hscale = graph_hscale, vscale = graph_vscale)
  tkgrid(reviewPlot, row = 0, column = 0, sticky = 'nsew')

  redraw_console_master <- function() {
    tkdelete(consoleText, '1.0', 'end')
    if (length(current_group_selected) == 0) {
      tkinsert(consoleText, 'end', 'No traces selected.\n')
    } else {
      tkinsert(consoleText, 'end',
        paste0('Selected traces: ', paste(current_group_selected, collapse = ', '), '\n'))
    }
    tkyview.moveto(consoleText, 1.0)
  }

  move_next_master <- function() {
    if (current_trace < total_traces) {
      current_trace <<- current_trace + 1L
      tkconfigure(infoLabel, text = paste('Trace', current_trace, 'of', total_traces))
      tkrreplot(reviewPlot)
    } else {
      tkinsert(consoleText, 'end', 'Review complete: Approved traces stored.\n')
      tkconfigure(acceptButton, state = 'disabled')
      tkconfigure(rejectButton, state = 'disabled')
      tkconfigure(undoButton,   state = 'disabled')
      tkyview.moveto(consoleText, 1.0)
    }
  }

  navBar <- tkframe(reviewFrame)
  tkgrid(navBar, row = 2, column = 1, pady = 5)

  acceptButton <<- tkbutton(navBar, text = 'Accept', command = function() {
    current_group_selected <<- c(current_group_selected, current_trace)
    redraw_console_master()
    move_next_master()
  })
  rejectButton <<- tkbutton(navBar, text = 'Reject', command = move_next_master)
  undoButton   <<- tkbutton(navBar, text = 'Undo',   command = function() {
    if (length(current_group_selected) > 0)
      current_group_selected <<- head(current_group_selected, -1)
    if (current_trace > 1) {
      current_trace <<- current_trace - 1L
      tkconfigure(infoLabel, text = paste('Trace', current_trace, 'of', total_traces))
      tkrreplot(reviewPlot)
    }
    redraw_console_master()
  })
  averageGroupButton <<- tkbutton(navBar, text = 'Add Selected Group', command = function() {
    if (length(current_group_selected) > 0) {
      groups_list[[length(groups_list) + 1]] <<- current_group_selected
      current_group_selected <<- integer(0)
      redraw_console_master()
    }
  })
  selectionCompleteButton <<- tkbutton(navBar, text = 'Selection Complete', command = function() {
    tkinsert(consoleText, 'end', 'Review complete: Approved traces stored.\n')
    tkyview.moveto(consoleText, 1.0)
  })

  tkgrid(acceptButton, row = 0, column = 0, padx = 5)
  tkgrid(rejectButton, row = 0, column = 1, padx = 5)
  tkgrid(undoButton,   row = 0, column = 2, padx = 5)
  tkgrid(averageGroupButton, row = 1, column = 0, padx = 5, pady=5)
  tkgrid(selectionCompleteButton, row = 1, column = 1, padx = 5, pady=5)
}

review_recordings <- function() {
  children <- as.character(tkwinfo('children', plotPanel))
  if (length(children)) sapply(children, function(ch) tcl("destroy", ch))
  if (!exists('abf_analysis_result', envir = .GlobalEnv)) {
    tkinsert(consoleText, 'end', "No analysis result available for review.\n")
    tkyview.moveto(consoleText, 1.0)
    return()
  }
  result         <- get('abf_analysis_result', envir = .GlobalEnv)
  datasets       <- result$datasets
  traces2average <<- lapply(datasets, function(x) integer(0))
  current_dataset <<- 1L
  current_trace   <<- 1L

  reviewFrame <<- tkframe(plotPanel)
  tkgrid(reviewFrame, row = 0, column = 0, sticky = 'nsew')
  tkgrid.columnconfigure(reviewFrame, 0, weight = 1)
  tkgrid.columnconfigure(reviewFrame, 1, weight = 1)
  tkgrid.columnconfigure(reviewFrame, 2, weight = 1)
  tkgrid.rowconfigure(reviewFrame,    1, weight = 1)

  infoLabel <<- tklabel(reviewFrame, text = paste(names(datasets)[1], 'trace', 1))
  tkgrid(infoLabel, row = 0, column = 1, pady = 5)

  plotWrapper <- tkframe(reviewFrame, height = graph_height, width = graph_width)
  tkgrid(plotWrapper, row = 1, column = 1, sticky = 'nsew')
  tkgrid.rowconfigure(plotWrapper,    0, weight = 1)
  tkgrid.columnconfigure(plotWrapper, 0, weight = 1)

  reviewPlot <<- tkrplot(plotWrapper, fun = function() {
    tk_par_settings()
    cex <- as.numeric(tclvalue(cexVar))
    par(cex.lab = cex, cex.axis = cex, cex.main = cex)

    ds    <- datasets[[current_dataset]]
    fname <- names(datasets)[current_dataset]
    tkconfigure(infoLabel, text = paste(fname, 'trace', current_trace))

    if (current_trace > length(ds$data)) {
      plot.new()
      text(0.5, 0.5, paste('No more recordings in', fname))
    } else {
      mat <- ds$data[[current_trace]]
      dt  <- ds$samplingIntervalInSec * 1000
      time<- seq(0, by = dt, length.out = nrow(mat))
      dc  <- as.numeric(tclvalue(dataColVar))
      if (is.na(dc) || dc < 1 || dc > ncol(mat)) dc <- 1
      y   <- mat[, dc]
      stim_time <- as.numeric(tclvalue(stimTimeVar))
      stim_y    <- y[which.min(abs(time - stim_time))]

      plot(time, y, type = 'l', col = 'darkgray',
           xlab = 'time (ms)', ylab = tclvalue(unitVar),
           xlim = smart_axis_limits(time),
           ylim = smart_axis_limits(y),
           axes = FALSE, bty = 'l')
      axis(1); axis(2, las = 1)
      points(stim_time, stim_y, pch = 8, col = 'black')
      text(  stim_time, stim_y, labels = 'stim', pos = 3, cex = cex)
    }
  }, hscale = graph_hscale, vscale = graph_vscale)
  tkgrid(reviewPlot, row = 0, column = 0, sticky = 'nsew')

  redraw_console_recordings <- function() {
    tkdelete(consoleText, '1.0', 'end')
    approved <- traces2average[[current_dataset]]
    msg <- if (length(approved) == 0) {
      paste0('No approved traces for ', names(datasets)[current_dataset])
    } else {
      paste0('Approved traces for ', names(datasets)[current_dataset], ': ',
             paste(approved, collapse = ', '))
    }
    tkinsert(consoleText, 'end', paste0(msg, '\n'))
    tkyview.moveto(consoleText, 1.0)
  }

  move_next_recordings <- function() {
    ds <- datasets[[current_dataset]]
    if (current_trace < length(ds$data)) {
      current_trace <<- current_trace + 1L
      tkconfigure(infoLabel, text = paste(names(datasets)[current_dataset], 'trace', current_trace))
      tkrreplot(reviewPlot)
    } else {
      if (current_dataset < length(datasets)) {
        current_dataset <<- current_dataset + 1L
        current_trace   <<- 1L
        tkconfigure(infoLabel, text = paste(names(datasets)[current_dataset], 'trace', current_trace))
        tkrreplot(reviewPlot)
      } else {
        tkinsert(consoleText, '1.0', 'end')
        tkinsert(consoleText, 'end', 'Review complete: Approved recordings stored.\n')
        tkconfigure(acceptButton, state = 'disabled')
        tkconfigure(rejectButton, state = 'disabled')
        tkconfigure(undoButton,   state = 'disabled')
        tkyview.moveto(consoleText, 1.0)
      }
    }
  }

  navBar <- tkframe(reviewFrame)
  tkgrid(navBar, row = 2, column = 1, pady = 5)
  acceptButton <<- tkbutton(navBar, text = 'Accept', command = function() {
    traces2average[[current_dataset]] <<- c(traces2average[[current_dataset]], current_trace)
    redraw_console_recordings()
    move_next_recordings()
  })
  rejectButton <<- tkbutton(navBar, text = 'Reject', command = move_next_recordings)
  undoButton   <<- tkbutton(navBar, text = 'Undo',   command = function() {
    if (length(traces2average[[current_dataset]]) > 0)
      traces2average[[current_dataset]] <<- head(traces2average[[current_dataset]], -1)
    if (current_trace > 1) {
      current_trace <<- current_trace - 1L
      tkconfigure(infoLabel, text = paste(names(datasets)[current_dataset], 'trace', current_trace))
      try({
        if (exists('reviewPlot', inherits = TRUE)) {
          widget_id <- as.character(reviewPlot$ID)
          if (tcl('winfo', 'exists', widget_id) == '1') tkrreplot(reviewPlot)
        }
      }, silent = TRUE)

    }
    redraw_console_recordings()
  })
  tkgrid(acceptButton, row = 0, column = 0, padx = 5)
  tkgrid(rejectButton, row = 0, column = 1, padx = 5)
  tkgrid(undoButton,   row = 0, column = 2, padx = 5)
}

 



# averaging Functions
# function to average selected groups for concatenated mode.
average_selected_groups <- function() {
  if (length(groups_list) == 0) {
    tkinsert(consoleText, 'end', "No groups available for averaging. Please select groups first.\n")
    tkyview.moveto(consoleText, 1.0)
    return()
  }
  
  dt_val      <- master_abf$samplingIntervalInSec * 1000
  stim_time   <- as.numeric(tclvalue(stimTimeVar))
  base_val    <- as.numeric(tclvalue(baselineVar))
  data_column <- as.numeric(tclvalue(dataColVar))
  if (is.na(data_column) || data_column < 1) data_column <- 1

  baseline2zero <- function(y, dt, stim, baseline) {
    idx_baseline <- round(baseline / dt)
    idx_start    <- round((stim - baseline) / dt) + 1
    y0 <- y - mean(y[1:idx_baseline])
    y0[idx_start:length(y0)]
  }

  averaged_data <<- lapply(groups_list, function(indices) {
    mats <- lapply(indices, function(i)
      baseline2zero(master_abf$data[[i]][, data_column],
                    dt_val, stim_time, base_val))
    rowMeans(do.call(cbind, mats))
  })

  current_avg_index <<- 1

  children <- as.character(tkwinfo('children', plotPanel))
  if (length(children)) sapply(children, function(ch) tcl("destroy", ch))
  tkgrid.columnconfigure(plotPanel, 0, weight = 1)
  tkgrid.rowconfigure(plotPanel, 0, weight = 1)

  avgFrame <<- tkframe(plotPanel)
  tkgrid(avgFrame, row = 0, column = 0, sticky = 'nsew')
  tkgrid.columnconfigure(avgFrame, 0, weight = 1)
  tkgrid.columnconfigure(avgFrame, 1, weight = 1)
  tkgrid.columnconfigure(avgFrame, 2, weight = 1)
  tkgrid.rowconfigure(avgFrame, 1, weight = 1)

  plotWrapper <- tkframe(avgFrame, height = graph_height, width = graph_width)
  tkgrid(plotWrapper, row = 1, column = 1, sticky = 'nsew')
  tkgrid.rowconfigure(plotWrapper,    0, weight = 1)
  tkgrid.columnconfigure(plotWrapper, 0, weight = 1)

  drawAvgPlot <- function() {
    tk_par_settings()
    cex <- as.numeric(tclvalue(cexVar))
    par(cex.lab = cex, cex.axis = cex, cex.main = cex)
    
    y <- averaged_data[[current_avg_index]]
    dt_val <- master_abf$samplingIntervalInSec * 1000
    time <- seq(from = stim_time - base_val, by = dt_val, length.out = length(y))
    
    plot(time, y, type = 'l', col = 'darkgray', xlab = 'time (ms)', ylab = tclvalue(unitVar),
         xlim = smart_axis_limits(time), ylim = smart_axis_limits(y),
         axes = FALSE, bty = 'l')
    axis(1); axis(2, las = 1)
    stim_y <- y[which.min(abs(time - stim_time))]
    points(stim_time, stim_y, pch = 8, col = 'black')
    text(stim_time, stim_y, labels = 'stim', pos = 3, cex = cex)
  }

  avgPlot <<- tkrplot(plotWrapper, fun = drawAvgPlot, hscale = graph_hscale, vscale = graph_vscale)
  tkgrid(avgPlot, row = 0, column = 0, sticky = 'nsew')

  navFrame <- tkframe(avgFrame)
  tkgrid(navFrame, row = 2, column = 1, pady = 5)
  navLabel   <- tklabel(navFrame, text = paste('Average:', current_avg_index, 'of', length(averaged_data)))
  nextButton <- tkbutton(navFrame, text = 'Next', command = function() {
    current_avg_index <<- if (current_avg_index < length(averaged_data)) current_avg_index + 1 else 1
    tkconfigure(navLabel, text = paste('Average:', current_avg_index, 'of', length(averaged_data)))
    tkrreplot(avgPlot)
  })
  
  tkgrid(navLabel,   row = 0, column = 0, padx = 5)
  tkgrid(nextButton, row = 0, column = 1, padx = 5)

  tkdelete(consoleText, '1.0', 'end')
  tkinsert(consoleText, 'end', 'Averaging complete. Check the updated plot.')
  tkyview.moveto(consoleText, 1.0)
}

averageApprovedTraces_sep <- function() {
  if (length(traces2average) == 0 || all(sapply(traces2average, length) == 0)) {
    tkinsert(consoleText, 'end', "No approved traces available. Please review recordings first.\n")
    tkyview.moveto(consoleText, 1.0)
    return()
  }
  result <- abf_averages(
    datasets         = abf_analysis_result$datasets,
    traces2average   = traces2average,
    baseline         = as.numeric(tclvalue(baselineVar)),
    stimulation_time = as.numeric(tclvalue(stimTimeVar)),
    dataCol          = as.numeric(tclvalue(dataColVar)),
    color            = 'darkgray',
    xbar             = as.numeric(tclvalue(xbarVar)),
    ybar             = as.numeric(tclvalue(ybarVar)),
    cex              = as.numeric(tclvalue(cexVar)),
    plotIt           = FALSE
  )
  averaged_data <<- result$baseline_corrected_mean_data
  datasets       <- result$datasets
  current_avg_index <<- 1

  children <- as.character(tkwinfo('children', plotPanel))
  if (length(children)) sapply(children, function(ch) tcl("destroy", ch))
  tkgrid.columnconfigure(plotPanel, 0, weight = 1)
  tkgrid.rowconfigure(   plotPanel, 0, weight = 1)

  avgFrame <<- tkframe(plotPanel)
  tkgrid(avgFrame, row = 0, column = 0, sticky = 'nsew')
  tkgrid.columnconfigure(avgFrame, 0, weight = 1)
  tkgrid.columnconfigure(avgFrame, 1, weight = 1)
  tkgrid.columnconfigure(avgFrame, 2, weight = 1)
  tkgrid.rowconfigure(   avgFrame, 1, weight = 1)

  plotWrapper <- tkframe(avgFrame, height = graph_height, width = graph_width)
  tkgrid(plotWrapper, row = 1, column = 1, sticky = 'nsew')
  tkgrid.rowconfigure(plotWrapper,    0, weight = 1)
  tkgrid.columnconfigure(plotWrapper, 0, weight = 1)

  drawSingleAvg <- function() {
    tk_par_settings()
    cex <- as.numeric(tclvalue(cexVar))
    par(cex.lab = cex, cex.axis = cex, cex.main = cex)
    y    <- averaged_data[[current_avg_index]]
    dt_val <- datasets[[current_avg_index]]$samplingIntervalInSec * 1000
    time   <- seq(from = as.numeric(tclvalue(stimTimeVar)) - as.numeric(tclvalue(baselineVar)),
                  by   = dt_val,
                  length.out = length(y))

    egs_plot(x = time, y = y,
             color     = 'darkgray',
             show_bar  = TRUE,
             show_text = TRUE,
             xbar      = as.numeric(tclvalue(xbarVar)),
             ybar      = as.numeric(tclvalue(ybarVar)),
             xlim      = smart_axis_limits(time),
             ylim      = smart_axis_limits(y),
             cex       = cex)

    stim_y <- y[which.min(abs(time - as.numeric(tclvalue(stimTimeVar))))]
    points(as.numeric(tclvalue(stimTimeVar)), stim_y, pch = 8, col = 'black')
    text(  as.numeric(tclvalue(stimTimeVar)), stim_y, labels = 'stim', pos = 3, cex = cex)
  }

  avgPlot <<- tkrplot(plotWrapper, fun = drawSingleAvg,
                      hscale = graph_hscale, vscale = graph_vscale)
  tkgrid(avgPlot, row = 0, column = 0, sticky = 'nsew')

  navFrame <- tkframe(avgFrame)
  tkgrid(navFrame, row = 2, column = 1, pady = 5)
  navLabel   <- tklabel(navFrame, text = paste('Average:', current_avg_index, 'of', length(averaged_data)))
  nextButton <- tkbutton(navFrame, text = 'Next', command = function() {
    current_avg_index <<- if (current_avg_index < length(averaged_data)) current_avg_index + 1 else 1
    tkconfigure(navLabel, text = paste('Average:', current_avg_index, 'of', length(averaged_data)))
    tkrreplot(avgPlot)
  })
  tkgrid(navLabel,   row = 0, column = 0, padx = 5)
  tkgrid(nextButton, row = 0, column = 1, padx = 5)

  tkdelete(consoleText, '1.0', 'end')
  tkinsert(consoleText, 'end', 'Separate-mode averaging complete. Check the updated plot.')
  tkyview.moveto(consoleText, 1.0)
}



  # UI Setup
ABF_analysis_tk <- function() {
    tt <- tktoplevel()
    tkwm.title(tt, 'ABF Analysis')
    
    if (.Platform$OS.type == "windows") {
      hscale <- 2
      vscale <- 2
    } else {
      dpi    <- as.numeric(tclvalue(tcl('winfo','pixels', tt, '1i')))
      w_in   <- 7;  h_in  <- 7
      hscale <- (w_in * dpi) / 480
      vscale <- (h_in * dpi) / 480
    }

    sidebarFrame <- tkframe(tt)
    mainFrame   <- tkframe(tt)
    tkgrid(sidebarFrame, row = 0, column = 0, sticky = 'ns')
    tkgrid(mainFrame,   row = 0, column = 1, sticky = 'nsew')
    tkgrid.rowconfigure(tt, 0, weight = 1)
    tkgrid.columnconfigure(tt, 1, weight = 1)

    plotPanel <<- mainFrame

    ## --- folder selector ---
    folderLabel <- tklabel(sidebarFrame, text = 'Select ABF Folder:')
    tkgrid(folderLabel, row = 0, column = 0, sticky = 'w')
    folderPathVar <<- tclVar('')
    folderEntry <- tkentry(sidebarFrame, textvariable = folderPathVar, width = 30)
    tkgrid(folderEntry, row = 0, column = 1, sticky = 'w')
    browseFolderButton <- tkbutton(sidebarFrame, text = 'Browse', command = function(){
      folderPath <- tclvalue(tkchooseDirectory())
      if (nchar(folderPath) > 0) {
        tclvalue(folderPathVar) <<- folderPath
        abf_list <- list.files(path = folderPath, pattern = '\\.abf$', ignore.case = TRUE)
        if (length(abf_list) == 0) {
          tkinsert(consoleText, 'end', 'No ABF files found in the selected folder.\n')
          tkyview.moveto(consoleText, 1.0)
        } else {
          tkdelete(abfListBox, 0, 'end')
          for (f in abf_list) tkinsert(abfListBox, 'end', f)
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
    abfListBox <<- tklistbox(sidebarFrame, height = 5, selectmode = 'multiple')
    tkgrid(abfListBox, row = 2, column = 0, columnspan = 2, sticky = 'we', padx = 5, pady = 3)

    tkgrid.columnconfigure(sidebarFrame, 0, weight = 1)
    tkgrid.columnconfigure(sidebarFrame, 1, weight = 1)
    tkgrid.columnconfigure(sidebarFrame, 2, weight = 1)

    paramFrame <- tkframe(sidebarFrame)

    ## --- parameters in paramFrame ---
    tkgrid(tklabel(paramFrame, text = 'Experiment:'), row = 0, column = 0, sticky = 'w')
    experimentCombo <- ttkcombobox(paramFrame,
      textvariable = experimentVar,
      values       = c('voltage clamp','current clamp'),
      width        = 15
    )
    tkgrid(experimentCombo, row = 0, column = 1, sticky = 'w')

    tkgrid(tklabel(paramFrame, text = 'Units:'), row = 1, column = 0, sticky = 'w')
    unitEntry <- tkentry(paramFrame,
      textvariable = unitVar,
      width        = 10,
      state        = 'readonly'
    )
    tkgrid(unitEntry, row = 1, column = 1, sticky = 'w')

    tkgrid(tklabel(paramFrame, text = 'Data Column:'), row = 2, column = 0, sticky = 'w')
    dataColEntry <- tkentry(paramFrame, textvariable = dataColVar, width = 10)
    tkgrid(dataColEntry, row = 2, column = 1, sticky = 'w')
    tkbind(dataColEntry, '<FocusOut>', function(...) {
      dc <- as.numeric(tclvalue(dataColVar))
      if (!exists('abf_analysis_result', envir = .GlobalEnv)) return()
      cu <- abf_analysis_result$datasets[[1]]$channelUnits
      if (!is.na(dc) && dc>=1 && dc<=length(cu)) tclvalue(unitVar) <<- cu[dc]
    })
    tkbind(dataColEntry,'<Return>',function(...) try(tcl("focus",""),silent=TRUE))

    tkgrid(tklabel(paramFrame, text = 'dt (ms):'), row = 3, column = 0, sticky = 'w')
    dtEntry <- tkentry(paramFrame, textvariable = dtVar, width = 10)
    tkgrid(dtEntry, row = 3, column = 1, sticky = 'w')
    tkbind(dtEntry,'<Return>',function(...){})

    tkgrid(tklabel(paramFrame, text = '# traces:'), row = 4, column = 0, sticky = 'w')
    ntracesEntry <- tkentry(paramFrame, textvariable = ntracesVar, width = 10)
    tkgrid(ntracesEntry, row = 4, column = 1, sticky = 'w')
    tkbind(ntracesEntry,'<Return>',function(...){})

    tkgrid(paramFrame, row = 3, column = 0, columnspan = 2, sticky = 'we', pady = 3)
    tkgrid.columnconfigure(paramFrame, 0, weight = 1)
    tkgrid.columnconfigure(paramFrame, 1, weight = 1)

    ## --- rest of the sidebar ---
    baselineVar <<- tclVar('100')
    stimTimeVar <<- tclVar('150')
    xbarVar     <<- tclVar('100')
    ybarVar     <<- tclVar('50')
    concatMode  <<- tclVar('0')

    tkgrid(tklabel(sidebarFrame, text = 'Baseline:'), row = 4, column = 0, sticky = 'w')
    baselineEntry <- tkentry(sidebarFrame, textvariable = baselineVar, width = 10)
    tkgrid(baselineEntry, row = 4, column = 1, sticky = 'w')
    tkbind(baselineEntry,'<Return>',function(...){})

    tkgrid(tklabel(sidebarFrame, text = 'Stimulation Time:'), row = 5, column = 0, sticky = 'w')
    stimTimeEntry <- tkentry(sidebarFrame, textvariable = stimTimeVar, width = 10)
    tkgrid(stimTimeEntry, row = 5, column = 1, sticky = 'w')
    tkbind(stimTimeEntry,'<Return>',function(...){})

    tkgrid(tklabel(sidebarFrame, text = 'x-bar length:'), row = 6, column = 0, sticky = 'w')
    xbarEntry <- tkentry(sidebarFrame, textvariable = xbarVar, width = 10)
    tkgrid(xbarEntry, row = 6, column = 1, sticky = 'w')
    tkbind(xbarEntry, '<FocusOut>', function(...) {
      try(if (exists('avgPlot',inherits=TRUE)) tkrreplot(avgPlot), silent=TRUE)
    })
    tkbind(xbarEntry, '<Return>', function(...) {
      try(if (exists('avgPlot',inherits=TRUE)) tkrreplot(avgPlot), silent=TRUE)
    })

    tkgrid(tklabel(sidebarFrame, text = 'y-bar length:'), row = 7, column = 0, sticky = 'w')
    ybarEntry <- tkentry(sidebarFrame, textvariable = ybarVar, width = 10)
    tkgrid(ybarEntry, row = 7, column = 1, sticky = 'w')
    tkbind(ybarEntry, '<FocusOut>', function(...) {
      try(if (exists('avgPlot',inherits=TRUE)) tkrreplot(avgPlot), silent=TRUE)
    })
    tkbind(ybarEntry, '<Return>', function(...) {
      try(if (exists('avgPlot',inherits=TRUE)) tkrreplot(avgPlot), silent=TRUE)
    })

    tkgrid(tklabel(sidebarFrame, text = 'Text scale (cex):'), row = 8, column = 0, sticky = 'w')
    cexEntry <- tkentry(sidebarFrame, textvariable = cexVar, width = 10)
    tkgrid(cexEntry, row = 8, column = 1, sticky = 'w')
    tkbind(cexEntry, '<FocusOut>', function(...) {
      try({
        if (exists('reviewPlot', inherits = TRUE) && is.tkwin(reviewPlot$ID)) tkrreplot(reviewPlot)
        if (exists('avgPlot',    inherits = TRUE) && is.tkwin(avgPlot$ID))    tkrreplot(avgPlot)
      }, silent = TRUE)
    })

    tkbind(cexEntry, '<Return>', function(...) {
      try({
        if (exists('reviewPlot', inherits = TRUE) && is.tkwin(reviewPlot$ID)) tkrreplot(reviewPlot)
        if (exists('avgPlot',    inherits = TRUE) && is.tkwin(avgPlot$ID))    tkrreplot(avgPlot)
      }, silent = TRUE)
    })

    concatButton <- tkcheckbutton(sidebarFrame, variable = concatMode,
                                  text = 'Concatenate Imported ABFs')
    tkgrid(concatButton, row = 9, column = 0, columnspan = 3)

    tkgrid.columnconfigure(sidebarFrame, 0, weight = 1)
    tkgrid.columnconfigure(sidebarFrame, 1, weight = 1)
    tkgrid.columnconfigure(sidebarFrame, 2, weight = 1)

    consoleText <<- tktext(sidebarFrame, height = 5)
    tkgrid(consoleText, row = 10, column = 0, columnspan = 3,
           sticky = 'we', padx = 10, pady = 5)

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
    tkbind(experimentCombo, '<<ComboboxSelected>>', function(widget, ...) {
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

        folderPath <- tclvalue(folderPathVar)
        if (nchar(folderPath) == 0) {
          tkinsert(consoleText, 'end', 'Please select an ABF folder first.\n')
          tkyview.moveto(consoleText, 1.0)
          return()
        }

      folderPath <- tclvalue(folderPathVar)
      if (nchar(folderPath) == 0) {
        tkinsert(consoleText, 'end', 'Please select an ABF folder first.\n')
        tkyview.moveto(consoleText, 1.0)
        return()
      }
      selIndices <- as.integer(tkcurselection(abfListBox))
      allFiles    <- as.character(tkget(abfListBox, 0, 'end'))
      abf_files   <- if (length(selIndices) == 0) allFiles else allFiles[selIndices + 1]
      if (length(abf_files) == 0) {
        tkinsert(consoleText, 'end', 'No ABF files selected.\n')
        tkyview.moveto(consoleText, 1.0)
        return()
      }
      result <- tryCatch({
        load_abf_data(abf_files = abf_files, abf_path = folderPath)
      }, error = function(e) {
        tkinsert(consoleText, 'end', paste0('Error during data loading: ', e$message, '\n'))
        tkyview.moveto(consoleText, 1.0)
        NULL
      })
      if (!is.null(result)) {
        tkdelete(consoleText, '1.0', 'end')
        tkinsert(consoleText, 'end', paste0('Data loaded. Processed ', length(abf_files), ' file(s).\n'))
        tkyview.moveto(consoleText, 1.0)
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
          tkyview.moveto(consoleText, 1.0)
        } else {
          tkinsert(consoleText, 'end', paste0('ERROR: ', cons_msg, '\n'))
          tkyview.moveto(consoleText, 1.0)
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

    tkgrid(runAnalysisButton,       row = 11, column = 0, columnspan = 3, pady = 5)
    tkgrid(reviewButton,            row = 12, column = 0, columnspan = 3, pady = 5)
    tkgrid(avgApprovedTracesButton, row = 13, column = 0, columnspan = 3, pady = 5)
    tkgrid(tkDownloadBtn,           row = 14, column = 0, columnspan = 3, pady = 5)

    tkfocus(tt)
    tkwait.window(tt)
  }

  # launch UI
  ABF_analysis_tk()
}  



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

analysePSCtk <- function() {

  PSC_analysis_tk <- function() {
    tt <- tktoplevel()
    tkwm.title(tt, 'PSC Analysis')

    if (.Platform$OS.type == "windows") {
      hscale <- 2  
      vscale <- 2
    } else {
      # keep your old 3″×3″ DPI math on non-Windows if you like:
      dpi    <- as.numeric(tclvalue(tcl('winfo','pixels', tt, '1i')))
      w_in   <- 7; h_in <- 7
      hscale <- (w_in * dpi) / 480
      vscale <- (h_in * dpi) / 480
    }
    
    # divide window into sidebar and main panels
    sidebarFrame <- tkframe(tt)
    mainFrame <- tkframe(tt)
    plotWidget <<- NULL
    tkgrid(sidebarFrame, row=0, column=0, sticky='ns')
    tkgrid(mainFrame, row=0, column=1, sticky='nsew')
    tkgrid.rowconfigure(tt, 0, weight=0)
    tkgrid.columnconfigure(tt, 1, weight=1)
    
    # sidebar controls
    fileLabel <- tklabel(sidebarFrame, text='Upload csv or xlsx:')
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
    
    xlimVar <- tclVar('')
    tkgrid(tklabel(graphSettingsFrame, text='x limits (e.g., 25,400):'), row=7, column=0, sticky='w')
    xlimEntry <- tkentry(graphSettingsFrame, textvariable=xlimVar, width=20)
    tkgrid(xlimEntry, row=7, column=1)

    tkbind(xlimEntry, "<Return>", function() {
      if (!is.null(plotWidget)) tkrreplot(plotWidget, fun=drawPlotXlim, silent=TRUE)
    })

    cexVar <- if (Sys.info()["sysname"] == "Darwin") tclVar('0.6') else tclVar('1.4')
    tkgrid(
      tklabel(graphSettingsFrame, text='Text scale (cex):'),
      row=8, column=0, sticky='w'
    )
    cexEntry <- tkentry(graphSettingsFrame, textvariable=cexVar, width=20)
    tkgrid(cexEntry, row=8, column=1)
    tkbind(cexEntry, "<Return>", function(widget, ...) {
      if (!is.null(plotWidget)) {
        # replot initial view
        tkrreplot(plotWidget, fun=drawPlot1, silent=TRUE)
      }
    })


    # Additional sidebar controls
    userTmaxVar <- tclVar('')
    tkgrid(tklabel(sidebarFrame, text='User maximum time for fit:'), row=3, column=0, sticky='w', pady=5)
    tkgrid(tkentry(sidebarFrame, textvariable=userTmaxVar, width=10), row=3, column=1, pady=5)
    
    repeatConstraintVar <- tclVar('0')
    tkgrid(tklabel(sidebarFrame, text='Add fast constraint:'), row=4, column=0, sticky='w')
    repeatConstraintCheck <- tkcheckbutton(sidebarFrame, variable=repeatConstraintVar)
    tkgrid(repeatConstraintCheck, row=4, column=1)
    
    buttonFrame <- tkframe(sidebarFrame)
    tkgrid(buttonFrame, row=5, column=0, columnspan=3, pady=10, sticky='ew')

    tkgrid.columnconfigure(sidebarFrame, 0, weight=1)
    tkgrid.columnconfigure(sidebarFrame, 1, weight=0)
    tkgrid.columnconfigure(sidebarFrame, 2, weight=1)

    tkgrid.columnconfigure(buttonFrame, 0, weight=1)
    tkgrid.columnconfigure(buttonFrame, 1, weight=0)
    tkgrid.columnconfigure(buttonFrame, 2, weight=0)
    tkgrid.columnconfigure(buttonFrame, 3, weight=1)

    # Analysis action buttons
    runAnalysisButton <- tkbutton(buttonFrame, text='Run Initial Analysis', command=function() {
      filePath <- tclvalue(filePathVar)
      if (nchar(filePath) == 0) {
        tkinsert(consoleText, 'end', 'Please select a file first\n')
        tkyview.moveto(consoleText, 1.0)
        return()
      }
      if (nchar(tclvalue(columnVar)) == 0) {
        tkinsert(consoleText, 'end', 'Please select a column\n')
        tkyview.moveto(consoleText, 1.0)
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
      # tkrreplot(plotWidget, fun=drawPlot1)

      if (is.null(plotWidget)) {
        plotWidget <<- tkrplot(
          mainFrame,
          fun    = drawPlot1,
          hscale = hscale,
          vscale = vscale
        )
        tkgrid(plotWidget, row=0, column=0, sticky='nsew')
      } else {
        tkrreplot(plotWidget, fun=drawPlot1)
      }

    })
      
    runMainAnalysisButton <- tkbutton(buttonFrame, text='Run Main Analysis', command=function() {
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
        xlim_input <- tclvalue(xlimVar)
        if (nchar(xlim_input) > 0) {
          xlim_vals <- as.numeric(unlist(strsplit(xlim_input, ',')))
          if (length(xlim_vals) == 2) {
            out$traces <- out$traces[out$traces$x >= xlim_vals[1] & out$traces$x <= xlim_vals[2], ]
          }
        }

        tkrreplot(plotWidget, fun=function() {
          drawPlot2(traces=out$traces, func=func, lwd=lwd, cex=as.numeric(tclvalue(cexVar)), filter=filter,
                    xbar=as.numeric(tclvalue(xbarVar)), ybar=as.numeric(tclvalue(ybarVar)),
                    xbar_lab=tclvalue(xbarLabVar), ybar_lab=tclvalue(ybarLabVar))
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
      # tkdelete(consoleText, '1.0', 'end')
      # tkinsert(consoleText, 'end', 'Analysis complete.')
      tkdelete(fitOutputText, '1.0', 'end')
      tkinsert(fitOutputText, 'end', paste(capture.output(print(df_out)), collapse='\n'))
    })
    
    downloadOutputButton <- tkbutton(buttonFrame, text='Download RData', 
      command=function() {
        if (!exists('analysis_output') || is.null(analysis_output)) {
          tkinsert(consoleText, 'end', 'No analysis output available!\n')
          tkyview.moveto(consoleText, 1.0)
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
          tkinsert(consoleText, 'end', 'Output saved successfully.\n')
          tkyview.moveto(consoleText, 1.0)
        }
      }
    )
    
      downloadResultsButton <- tkbutton(buttonFrame, text='Download Output (csv/xlsx)', command=function() {
        if (!exists('analysis_output') || is.null(analysis_output)) {
          tkinsert(consoleText, 'end', 'No analysis output available!\n')
          tkyview.moveto(consoleText, 1.0)
          return()
        }
        filePath <- tclvalue(tkgetSaveFile(filetypes='{{Excel File} {.xlsx}} {{CSV File} {.csv}}'))
        if (nchar(filePath) == 0) return()
        ext <- tolower(tools::file_ext(filePath))

        data_list <- list(
          output          = analysis_output$output,
          traces          = analysis_output$traces,
          `fit criterion` = data.frame(AIC = analysis_output$AIC, BIC = analysis_output$BIC),
          `model message` = data.frame(message = analysis_output$model.message)
        )

        metadata_labels <- c(
          'Data column:','dt (ms):','Stimulation Time:','Baseline:','n:','Fit cutoff:','Function:',
          'Downsample Factor:','User maximum time for fit:','Add fast constraint:','N:','IEI:','Smooth:',
          'Method:','Weighting:','Sequential Fit:','Min interval:','Max interval:',
          'Lower bounds (comma-separated):','Upper bounds (comma-separated):','Latency limit:',
          'MLE Iterations:','Metropolis Scale:','Fit Attempts:','Random Walk Metropolis:',
          'Filter:','Filter cutoff (Hz):','Half-width fit limit:','Seed:','Decimal points:',
          'Fast constraint method:','Fast decay limit(s):','First delay constraint:'
        )

        metadata_values <- c(
          tclvalue(columnVar),tclvalue(dtVar),tclvalue(stimTimeVar),tclvalue(baselineVar),
          tclvalue(nVar),tclvalue(yAblineVar),tclvalue(funcVar),tclvalue(dsVar),
          tclvalue(userTmaxVar),as.character(as.logical(as.numeric(tclvalue(repeatConstraintVar)))),
          tclvalue(NVar),tclvalue(IEIVar),tclvalue(smoothVar),tclvalue(methodVar),
          tclvalue(weightMethodVar),as.character(as.logical(as.numeric(tclvalue(sequentialFitVar)))),
          tclvalue(intervalMinVar),tclvalue(intervalMaxVar),tclvalue(lowerVar),tclvalue(upperVar),
          tclvalue(latencyLimitVar),tclvalue(iterVar),tclvalue(metropolisScaleVar),
          tclvalue(fitAttemptsVar),as.character(as.logical(as.numeric(tclvalue(RWmVar)))),
          as.character(as.logical(as.numeric(tclvalue(filterVar)))),tclvalue(fcVar),
          tclvalue(halfWidthFitLimitVar),tclvalue(seedVar),tclvalue(dpVar),
          tclvalue(fastConstraintMethodVar),tclvalue(fastDecayLimitVar),
          as.character(as.logical(as.numeric(tclvalue(firstDelayConstraintVar))))
        )

        coerce_val <- function(v) {
          if (nzchar(v)) {
            n <- suppressWarnings(as.numeric(v))
            if (!is.na(n)) return(n)
          }
          v
        }

        if (ext == 'csv') {
          for (nm in names(data_list)) {
            write.csv(data_list[[nm]], file = sub("\\.csv$", paste0("_", nm, ".csv"), filePath), row.names = FALSE)
          }
          meta_df <- data.frame(name = metadata_labels, value = sapply(metadata_values, coerce_val), stringsAsFactors = FALSE)
          write.csv(meta_df, file = sub("\\.csv$", "_metadata.csv", filePath), row.names = FALSE)

        } else if (ext == 'xlsx') {
          wb <- createWorkbook()
          for (nm in names(data_list)) {
            addWorksheet(wb, nm)
            writeData(wb, nm, data_list[[nm]])
          }
          addWorksheet(wb, "metadata")
          writeData(wb, "metadata", c("parameter", "value"), startRow = 1, startCol = 1, colNames = FALSE)
          for (i in seq_along(metadata_labels)) {
            writeData(wb, "metadata", metadata_labels[i], startRow = i + 1, startCol = 1, colNames = FALSE)
            writeData(wb, "metadata", coerce_val(metadata_values[i]), startRow = i + 1, startCol = 2, colNames = FALSE)
          }
          saveWorkbook(wb, filePath, overwrite = TRUE)

        } else {
          tkinsert(consoleText, 'end', 'Unsupported file type. Use .csv or .xlsx\n')
          tkyview.moveto(consoleText, 1.0)
        }
      })

    exportSVGButton <- tkbutton(buttonFrame, text='Export Plot to SVG', command=function() {
      if (!exists('analysis_output') || is.null(analysis_output)) {
        tkinsert(consoleText, 'end', 'No analysis available to export!\n')
        tkyview.moveto(consoleText, 1.0)
        return()
      }
      saveFile <- tclvalue(tkgetSaveFile(filetypes='{{SVG Files} {.svg}}'))
      if (nchar(saveFile) > 0) {
        # Parse and apply xlim if provided
        traces <- analysis_output$traces
        xlim_input <- tclvalue(xlimVar)
        if (nchar(xlim_input) > 0) {
          xlim_vals <- as.numeric(unlist(strsplit(xlim_input, ',')))
          if (length(xlim_vals) == 2) {
            traces <- traces[traces$x >= xlim_vals[1] & traces$x <= xlim_vals[2], ]
          }
        }
        
        svg(filename=saveFile, width=7, height=5)
        drawPlot2(traces=traces, func=get(tclvalue(funcVar)), lwd=as.numeric(tclvalue(lwdVar)), cex=0.6,
                  filter=as.logical(as.numeric(tclvalue(filterVar))),
                  xbar=as.numeric(tclvalue(xbarVar)), ybar=as.numeric(tclvalue(ybarVar)),
                  xbar_lab=tclvalue(xbarLabVar), ybar_lab=tclvalue(ybarLabVar))
        dev.off()
        tkinsert(consoleText, 'end', 'SVG plot saved successfully.\n')
        tkyview.moveto(consoleText, 1.0)
      }
    })

    
    clearOutputButton <- tkbutton(buttonFrame, text='Clear Output', command=function() {
      analysis_output <<- NULL
      # tkdelete(consoleText, '1.0', 'end')
      tkdelete(fitOutputText, '1.0', 'end')
      # tkrreplot(plotWidget, fun=drawPlot1)
      # plotWidget <<- tkrplot(mainFrame, fun=drawPlot1)
      # tkgrid(    plotWidget,    row=0, column=0, sticky='nsew')
      if (!is.null(plotWidget)) {
        tkrreplot(plotWidget, fun=drawPlot1)
      }

    })
    
    tkgrid(runAnalysisButton,      row=0, column=0, padx=5, pady=2)
    tkgrid(runMainAnalysisButton,  row=0, column=1, padx=5, pady=2)
    tkgrid(downloadResultsButton,  row=1, column=1, padx=5, pady=2)
    tkgrid(downloadOutputButton,   row=1, column=0, padx=5, pady=2)
    tkgrid(exportSVGButton,        row=2, column=0, padx=5, pady=2)
    tkgrid(clearOutputButton,      row=2, column=1, padx=5, pady=2)

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
      cex <- as.numeric(tclvalue(cexVar))
      
      determine_tmax2(y=y_val, N=1, dt=dt, stimulation_time=stimTime, baseline=baseline, smooth=smooth, lwd=lwd,
        tmax=NULL, y_abline=y_abline, xbar=as.numeric(tclvalue(xbarVar)), ybar=as.numeric(tclvalue(ybarVar)),
        xbar_lab=tclvalue(xbarLabVar), ybar_lab=tclvalue(ybarLabVar), cex=cex)
    }
    
    drawPlotXlim <- function() {
      xlim_input <- tclvalue(xlimVar)
      cex <- as.numeric(tclvalue(cexVar))
      traces <- if (exists("analysis_output") && !is.null(analysis_output)) analysis_output$traces else NULL
      if (is.null(traces)) return()

      if (nchar(xlim_input) > 0) {
        xlim_vals <- as.numeric(unlist(strsplit(xlim_input, ",")))
        if (length(xlim_vals) == 2) {
          traces <- traces[traces$x >= xlim_vals[1] & traces$x <= xlim_vals[2], ]
        }
      }

      drawPlot2(traces=traces, func=get(tclvalue(funcVar)), lwd=as.numeric(tclvalue(lwdVar)), cex=cex,
                filter=as.logical(as.numeric(tclvalue(filterVar))),
                xbar=as.numeric(tclvalue(xbarVar)), ybar=as.numeric(tclvalue(ybarVar)),
                xbar_lab=tclvalue(xbarLabVar), ybar_lab=tclvalue(ybarLabVar))
    }
    # plotWidget <- tkrplot(tt, fun=drawPlot1)
    # tkgrid(plotWidget, row=0, column=1, sticky='nsew')
    
    # consoleText <- tktext(mainFrame, width=80, height=4)
    # tkgrid(consoleText, row=1, column=1, sticky='nsew')
    
    fitOutputLabel <- tklabel(sidebarFrame, text='Fit Output:')
    tkgrid(fitOutputLabel, row=11, column=0, columnspan=3, sticky='w', pady=c(10,2), padx=20)
    
    fitOutputText <- tktext(sidebarFrame, width=90, height=5)
    tkgrid(fitOutputText, row=12, column=0, columnspan=3, sticky='w', padx=20)
    
    tkfocus(tt)
    tkwait.window(tt)
  }

  PSC_analysis_tk()

}

analysePSCshiny <- function() {

  ui <- fluidPage(
    titlePanel('PSC Analysis'),
    sidebarLayout(
      sidebarPanel(
        fileInput('file', 'Upload csv or xlsx', accept=c('.csv', '.xlsx')),
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
                   textInput('ybar_lab', 'y-axis Units:', 'pA'),
                   textInput('xlim', 'x limits (e.g., 0,400):', '')
          )
        ),
        
        numericInput('userTmax', 'User Maximum Time for Fit:', NA),
        actionButton('run_initial', 'Run Initial Analysis'),
        actionButton('run_main', 'Run Main Analysis'),
        downloadButton('download_xlsx',  'Download Output (*.xlsx)'),
        downloadButton('download_output', 'Download RData'),
        downloadButton('download_svg', 'Download SVG Plot'),
        actionButton('clear_output', 'Clear Output')
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

        xlim_vals <- if (nchar(input$xlim) > 0) as.numeric(unlist(strsplit(input$xlim, ","))) else NULL
        traces <- state$analysis$traces
        if (!is.null(xlim_vals) && length(xlim_vals) == 2) {
          traces <- traces[traces$x >= xlim_vals[1] & traces$x <= xlim_vals[2], ]
        }

        drawPlot2(traces=traces, func=func, lwd=lwd,
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

    output$download_svg <- downloadHandler(
      filename = function() {
        paste0('PSC_plot_', Sys.Date(), '.svg')
      },
      content = function(file) {
        req(state$analysis)
        func <- switch(input$func,
                       'product1N' = product1N,
                       'product2N' = product2N,
                       'product3N' = product3N,
                       product1N)

        traces <- state$analysis$traces
        xlim_vals <- if (nchar(input$xlim) > 0) as.numeric(unlist(strsplit(input$xlim, ","))) else NULL
        if (!is.null(xlim_vals) && length(xlim_vals) == 2) {
          traces <- traces[traces$x >= xlim_vals[1] & traces$x <= xlim_vals[2], ]
        }

        svg(filename = file, width = 7, height = 5)
        drawPlot2(
          traces = traces,
          func = func,
          lwd = as.numeric(input$lwd),
          filter = input$filter,
          xbar = as.numeric(input$xbar),
          ybar = as.numeric(input$ybar),
          xbar_lab = input$xbar_lab,
          ybar_lab = input$ybar_lab
        )
        dev.off()
      }
    )

    output$download_xlsx <- downloadHandler(
      filename = function() {
        paste0(
          tools::file_path_sans_ext(basename(input$file$name)),
          "_", input$data_col,
          "_PSC_analysis.xlsx"
        )
      },
      content = function(file) {
        req(state$analysis)

        wb <- openxlsx::createWorkbook()

        # result sheets
        data_list <- list(
          output          = state$analysis$output,
          traces          = state$analysis$traces,
          `fit criterion` = data.frame(AIC = state$analysis$AIC, BIC = state$analysis$BIC),
          `model message` = data.frame(message = state$analysis$model.message)
        )
        for (nm in names(data_list)) {
          openxlsx::addWorksheet(wb, nm)
          openxlsx::writeData(wb, nm, data_list[[nm]])
        }

        # metadata sheet
        metadata_labels <- c(
          'Data column:','dt (ms):','Stimulation Time:','Baseline:','n:','Fit cutoff:','Function:',
          'Downsample Factor:','User maximum time for fit:','Add fast constraint:','N:','IEI:','Smooth:',
          'Method:','Weighting:','Sequential Fit:','Min interval:','Max interval:',
          'Lower bounds (comma-separated):','Upper bounds (comma-separated):','Latency limit:',
          'MLE Iterations:','Metropolis Scale:','Fit Attempts:','Random Walk Metropolis:',
          'Filter:','Filter cutoff (Hz):','Half-width fit limit:','Seed:','Decimal points:',
          'Fast constraint method:','Fast decay limit(s):','First delay constraint:'
        )
        metadata_values <- list(
          input$data_col,
          input$dt,
          input$stimulation_time,
          input$baseline,
          input$n,
          input$y_abline,
          input$func,
          input$ds,
          input$userTmax,
          input$fast_constraint,
          input$N,
          input$IEI,
          input$smooth,
          input$method,
          input$weight_method,
          input$sequential_fit,
          input$interval_min,
          input$interval_max,
          input$lower,
          input$upper,
          input$latency_limit,
          input$iter,
          input$metropolis_scale,
          input$fit_attempts,
          input$RWm,
          input$filter,
          input$fc,
          input$half_width_fit_limit,
          input$seed,
          input$dp,
          input$fast_constraint_method,
          input$fast_decay_limit,
          input$first_delay_constraint
        )

        numeric_labels <- c(
          'dt (ms):','Stimulation Time:','Baseline:','n:','Fit cutoff:',
          'Downsample Factor:','User maximum time for fit:','N:','IEI:','Smooth:',
          'Min interval:','Max interval:','Latency limit:','MLE Iterations:',
          'Metropolis Scale:','Fit Attempts:','Filter cutoff (Hz):',
          'Half-width fit limit:','Seed:','Decimal points:'
        )
        logical_labels <- c(
          'Add fast constraint:','Sequential Fit:','Random Walk Metropolis:','Filter:',
          'First delay constraint:'
        )

        openxlsx::addWorksheet(wb, "metadata")
        openxlsx::writeData(wb, "metadata", c("Parameter","Value"), startRow = 1, startCol = 1, colNames = FALSE)
        for (i in seq_along(metadata_labels)) {
          lbl <- metadata_labels[i]
          val <- metadata_values[[i]]
          openxlsx::writeData(wb, "metadata", lbl,       startRow = i+1, startCol = 1, colNames = FALSE)
          if (lbl %in% numeric_labels) {
            openxlsx::writeData(wb, "metadata", as.numeric(val), startRow = i+1, startCol = 2, colNames = FALSE)
          } else if (lbl %in% logical_labels) {
            openxlsx::writeData(wb, "metadata", as.logical(val), startRow = i+1, startCol = 2, colNames = FALSE)
          } else {
            openxlsx::writeData(wb, "metadata", val,      startRow = i+1, startCol = 2, colNames = FALSE)
          }
        }

        openxlsx::saveWorkbook(wb, file, overwrite = TRUE)
      }
    )
    
  }

  shinyApp(ui=ui, server=server)
}

widgetPSCtk <- function() {

  widget_PSC_tk <- function() {
    ## only one tktoplevel() here:
    tt <- tktoplevel()
    tkwm.title(tt, 'I (pA)')

    Tr            <- tclVar('5')
    Td            <- tclVar('25')
    decayUpperVar <- tclVar('0.9')
    decayLowerVar <- tclVar('0.1')

    updatePlot <- function() {
      currentTr   <- as.numeric(tclvalue(Tr))
      currentTd   <- as.numeric(tclvalue(Td))
      decay_range <- c(
        as.numeric(tclvalue(decayUpperVar)),
        as.numeric(tclvalue(decayLowerVar))
      )
      t_max  <- 250
      t_seq  <- seq(0, t_max, length.out = 2000)
      y_comb <- -(exp(-t_seq / currentTd) - exp(-t_seq / currentTr))
      y_Td   <- -exp(-t_seq / currentTd)
      y_Tr   <- -exp(-t_seq / currentTr)
      T1     <- currentTd * currentTr / (currentTd - currentTr)
      y_1mTr <- -(1 - exp(-t_seq / T1))

      peak_i <- which.min(y_comb)
      y_fall <- y_comb[peak_i:length(y_comb)]
      t_fall <- t_seq[peak_i:length(t_seq)]
      min_v  <- min(y_fall)
      up_v   <- decay_range[1] * min_v
      lo_v   <- decay_range[2] * min_v
      idx_up <- which(y_fall >= up_v)[1]
      idx_lo <- which(y_fall >= lo_v)[1]

      if (is.na(idx_lo)) {
        ext_t  <- seq(0, t_max*2, length.out = length(t_seq)*2)
        ext_y  <- -(exp(-ext_t/currentTd) - exp(-ext_t/currentTr))
        y_fall <- ext_y[peak_i:length(ext_y)]
        t_fall <- ext_t[peak_i:length(ext_t)]
        idx_lo <- which(y_fall >= lo_v)[1]
      }

      decay_time <- if (!is.na(idx_up) && !is.na(idx_lo))
        round(t_fall[idx_lo] - t_fall[idx_up], 2) else NA

      up_pct <- round(decay_range[1]*100)
      lo_pct <- round(decay_range[2]*100)
      title2 <- paste('normalised with', up_pct, '-', lo_pct, '% decay time')

      par(mfrow = c(2,1), mar = c(4,4,3,1))
      plot(t_seq, y_comb, type='l', col='grey', lwd=2, axes=FALSE,
           main=expression(e^{-t/tau[decay]} - e^{-t/tau[rise]}),
           xlab='', ylab='F(x)', xlim=c(0,t_max), ylim=c(-1,0))
      lines(t_seq, y_Td, col='indianred', lty=3, lwd=2)
      lines(t_seq, y_Tr, col='slateblue', lty=3, lwd=2)
      legend('bottomright', legend = c(
        expression(e^{-t/tau[decay]} - e^{-t/tau[rise]}),
        expression(e^{-t/tau[decay]}),
        expression(e^{-t/tau[rise]})
      ),
      col = c('grey', 'indianred', 'slateblue'), lty = c(1,3,3), lwd = 2,
      bty = 'n', inset = c(0.02, 0.15))
      axis(2, las=1, tcl=-0.2)

      plot(t_seq, y_comb/abs(min(y_comb)), type='l', col='grey', lwd=2, axes=FALSE,
           main=title2, xlab='time (ms)', ylab='normalised F(x)',
           xlim=c(0,t_max), ylim=c(-1,0))
      lines(t_seq, y_1mTr, col='slateblue', lty=3, lwd=2)
      lines(t_seq, y_Td, col='indianred', lty=3, lwd=2)
      legend('bottomright', legend = c(
        expression((e^{-t/tau[decay]} - e^{-t/tau[rise]}) / abs),
        expression(-(1 - e^{-t/tau[1]})),
        expression(-e^{-t/tau[decay]})
      ),
      col = c('gray','slateblue','indianred'), lty = c(1,3,3), lwd = 2,
      bty = 'n', inset = c(0.02, 0.15))
      if (!is.na(decay_time)) {
        abline(h=-up_v/min_v, col='darkgrey', lty=3)
        abline(h=-lo_v/min_v, col='darkgrey', lty=3)
        text(t_max*0.8, -0.3,
             paste(up_pct, '-', lo_pct, 'decay =', decay_time, 'ms'),
             col='darkgrey', cex=0.8)
      } else {
        text(t_max*0.8, -0.3, 'Decay not fully reached', col='indianred', cex=0.8)
      }
      axis(1, las=1, tcl=-0.2)
      axis(2, las=1, tcl=-0.2)
    }

    # build the UI
    img <- tkrplot(tt, fun=updatePlot, hscale=1.3, vscale=1.3)
    tkpack(img, side='top', expand=TRUE, fill='both')

    decayFrame <- tkframe(tt)
    tkpack(decayFrame, side='top', fill='x', pady=5)
    tkpack(tklabel(decayFrame, text='Decay range (high, low):'), side='left', padx=5)
    highEntry <- tkentry(decayFrame, textvariable=decayUpperVar, width=5)
    lowEntry  <- tkentry(decayFrame, textvariable=decayLowerVar, width=5)
    tkpack(highEntry, side='left', padx=2)
    tkpack(lowEntry,  side='left', padx=2)
    tkbind(highEntry, '<KeyRelease>', function() tkrplot::tkrreplot(img))
    tkbind(lowEntry,  '<KeyRelease>', function() tkrplot::tkrreplot(img))

    tkpack(tklabel(tt, text='\u03C4 rise'),  side='top')
    Tr_slider <- tkscale(tt, from=0.1, to=50, resolution=0.1,
                         showvalue=TRUE, variable=Tr, orient='horizontal',
                         command=function(...) tkrplot::tkrreplot(img))
    tkpack(Tr_slider, fill='x', padx=10, pady=5)

    tkpack(tklabel(tt, text='\u03C4 decay'), side='top')
    Td_slider <- tkscale(tt, from=0.1, to=200, resolution=1,
                         showvalue=TRUE, variable=Td, orient='horizontal',
                         command=function(...) tkrplot::tkrreplot(img))
    tkpack(Td_slider, fill='x', padx=10, pady=5)

    tkfocus(tt)
    tkwait.window(tt)
  }

  widget_PSC_tk()
}


widgetPSCshiny <- function() {

  ui <- fluidPage(
    titlePanel(
      HTML("Interactive Graphs:<br>exp(-t/&tau;<sub>decay</sub>) - exp(-t/&tau;<sub>rise</sub>)")
    ),
    sidebarLayout(
      sidebarPanel(
        numericInput(
          "decayUpper", 
          "Decay upper (fraction)", 
          value = 0.9, min = 0, max = 1, step = 0.01
        ),
        numericInput(
          "decayLower", 
          "Decay lower (fraction)", 
          value = 0.1, min = 0, max = 1, step = 0.001
        ),
        sliderInput(
          "Tr", 
          HTML("&tau;<sub>rise</sub>"), 
          min = 0.1, max = 50, value = 5, step = 0.1
        ),
        sliderInput(
          "Td", 
          HTML("&tau;<sub>decay</sub>"), 
          min = 0.101, max = 200, value = 25, step = 1
        )
      ),
      mainPanel(
        plotlyOutput("upperPlot"),
        plotlyOutput("lowerPlot")
      )
    )
  )

  server <- function(input, output, session) {
    
    # ensure Td > Tr
    observeEvent(input$Tr, {
      updateSliderInput(
        session, "Td",
        min   = input$Tr + 0.001,
        value = max(input$Td, input$Tr + 0.001)
      )
    })
    
    output$upperPlot <- renderPlotly({
      t_max  <- 250
      t_seq  <- seq(0, t_max, length.out = 2000)
      y_combined <- -(exp(-t_seq / input$Td) - exp(-t_seq / input$Tr))
      y_Td <- -exp(-t_seq / input$Td)
      y_Tr <- -exp(-t_seq / input$Tr)
      
      plot_ly(type = 'scatter', mode = 'lines') %>%
        add_trace(x = t_seq, y = y_combined, name = 'exp(-t/Td) - exp(-t/Tr)',
                  line = list(color = 'gray')) %>%
        add_trace(x = t_seq, y = y_Td, name = 'exp(-t/Td)',
                  line = list(color = 'indianred', dash = 'dot')) %>%
        add_trace(x = t_seq, y = y_Tr, name = 'exp(-t/Tr)',
                  line = list(color = 'slateblue', dash = 'dot')) %>%
        layout(
          title = 'exp(-t/τ_decay) - exp(-t/τ_rise)',
          xaxis = list(title = 't (ms)'),
          yaxis = list(title = 'F(x)')
        )
    })
    
    output$lowerPlot <- renderPlotly({
      t_max    <- 250
      t_seq    <- seq(0, t_max, length.out = 2000)
      y_comb   <- -(exp(-t_seq / input$Td) - exp(-t_seq / input$Tr))
      norm     <- y_comb / min(y_comb)
      peak_i   <- which.min(y_comb)
      y_fall   <- y_comb[peak_i:length(y_comb)]
      t_fall   <- t_seq[peak_i:length(t_seq)]
      min_val  <- min(y_fall)
      
      du <- input$decayUpper * min_val
      dl <- input$decayLower * min_val
      idx_up   <- which(y_fall >= du)[1]
      idx_lo   <- which(y_fall >= dl)[1]
      if (is.na(idx_lo)) {
        ext_t <- seq(0, t_max*2, length.out = length(t_seq)*2)
        ext_y <- -(exp(-ext_t/input$Td) - exp(-ext_t/input$Tr))
        y_fall <- ext_y[peak_i:length(ext_y)]
        t_fall <- ext_t[peak_i:length(ext_t)]
        idx_lo <- which(y_fall >= dl)[1]
      }
      decay_time <- if (!is.na(idx_up) && !is.na(idx_lo))
        round(t_fall[idx_lo] - t_fall[idx_up], 2) else NA
      
      up_pct <- round(input$decayUpper * 100)
      lo_pct <- round(input$decayLower * 100)
      title2 <- paste("normalized with", up_pct, "-", lo_pct, "% decay time")
      
      p <- plot_ly(type = 'scatter', mode = 'lines') %>%
        add_trace(x = t_seq, y = -norm, name = 'normalized F(x)',
                  line = list(color = 'gray')) %>%
        add_trace(x = t_seq,
                  y = -(1 - exp(-t_seq / (input$Td * input$Tr /
                                          (input$Td - input$Tr)))),
                  name = '-(1 - exp(-t/τ1))',
                  line = list(color = 'slateblue', dash = 'dot')) %>%
        add_trace(x = t_seq, y = -exp(-t_seq / input$Td),
                  name = '-exp(-t/τ_decay)',
                  line = list(color = 'indianred', dash = 'dot')) %>%
        layout(title = title2,
               xaxis = list(title = 'time (ms)'),
               yaxis = list(title = 'normalized F(x)'))
      
      if (!is.na(decay_time)) {
        p <- p %>%
          add_trace(x = c(0, t_max), y = rep(-du/min_val, 2),
                    showlegend = FALSE,
                    line = list(color = 'darkgrey', dash = 'dot')) %>%
          add_trace(x = c(0, t_max), y = rep(-dl/min_val, 2),
                    showlegend = FALSE,
                    line = list(color = 'darkgrey', dash = 'dot')) %>%
          add_annotations(x = t_max * 0.8, y = -0.5,
                          text = paste(up_pct, "-", lo_pct,
                                       "decay =", decay_time, "ms"),
                          showarrow = FALSE,
                          font = list(color = 'darkgrey'))
      } else {
        p <- p %>%
          add_annotations(x = t_max * 0.8, y = -0.5,
                          text = 'decay not fully reached',
                          showarrow = FALSE,
                          font = list(color = 'indianred'))
      }
      
      p
    })
  }

  shinyApp(ui, server)

}
