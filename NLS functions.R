# NLS functions
# Rewrite residFun in C++ and SS.fun in C++


# Define the C++ code as a string
cpp_code <- '
#include <Rcpp.h>
using namespace Rcpp;

// Function to calculate residuals
// [[Rcpp::export]]
NumericVector residFunCpp(NumericVector params, NumericVector y, NumericVector x, Function func) {
  int n = y.size();
  NumericVector residuals(n);
  NumericVector fitted_values = as<NumericVector>(func(params, x));
  for (int i = 0; i < n; i++) {
    residuals[i] = y[i] - fitted_values[i];
  }
  return residuals;
}

// Function to calculate sum of squares
// [[Rcpp::export]]
double SSfunCpp(NumericVector params, NumericVector x, NumericVector y, Function func) {
  NumericVector residuals = residFunCpp(params, y, x, func);
  return sum(pow(residuals, 2));
}
'

# Write the C++ code to a temporary file
cpp_file <- tempfile(fileext = ".cpp")
writeLines(cpp_code, cpp_file)

# Source the C++ file
sourceCpp(cpp_file)

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
MLE.fun <- function(logpost, par, x, y, gr=NULL, lower=NULL, upper=NULL, method='Nelder-Mead', func=product1) {
  if (method %in% c('L-BFGS-B', 'Brent')) {
    suppressWarnings(
      MLE.fit <- optim(par=par, fn=logpost, gr=gr, method=method, hessian=TRUE, 
                 x=x, y=y, lower=lower, upper=upper, func=func, 
                 control=list(fnscale=-1, maxit=2000, factr=1e7, pgtol=1e-8, parscale=rep(1, length(par))))

      # MLE.fit <- optim(par=par, fn=logpost, gr=gr, method=method, hessian=TRUE, x=x, y=y, lower=lower, upper=upper, func=func, control=list(fnscale=-1))
    )
  } else {
    suppressWarnings(
      MLE.fit <- optim(par=par, fn=logpost, gr=gr, method=method, hessian=TRUE, x=x, y=y, func=func, 
        control = list(fnscale = -1, maxit = 2000, reltol = 1e-8, parscale = rep(1, length(par))))
    )
  }

  fit <- MLE.fit$par
  h <- -solve(MLE.fit$hessian)
  p <- length(fit)
  int <- p / 2 * log(2 * pi) + 0.5 * log(det(h)) + logpost(fit, x, y, func)
  list(fit = fit, fit.se = sqrt(diag(h)), var = h, int = int, converge = MLE.fit$convergence == 0)
}

# Optimized RWmetropolis function
RWmetropolis <- function(logpost, cov.mat, scale, start, m, x, y, func) {
  pb <- length(start)
  Mpar <- matrix(0, m, pb)
  b <- matrix(t(start))
  lb <- logpost(start, x, y, func)
  a <- chol(cov.mat)
  accept <- 0
  
  random_samples <- scale * t(a) %*% matrix(rnorm(pb * m), nrow = pb)
  
  for (ii in 1:m) {
    bc <- b + random_samples[, ii]
    lbc <- logpost(t(bc), x, y, func)
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

model.selection.criteria <- function(coeffs, x, y, func) {
  res <- residFunCpp(coeffs, y, x, func)
  k <- length(coeffs)
  N <- length(res)
  # definition: log_likelihood <- -(N / 2) * log(2 * pi) - (N / 2) * log(sigma2) - (1 / (2 * sigma2)) * sum(res^2)
  # log(sigma2) == sum(res^2) / N and (1 / (2 * sigma2)) * sum(res^2) == 0.5 * N
  loglik <- -0.5 * (N * log(2 * pi) + N * log(sum(res^2) / N) + N)
  # loglik <- 0.5 * (-N * (log(2 * pi) + 1 - log(N) + log(sum(res^2))))
  df <- k + 1
  BIC <- df * log(N) - 2 * loglik
  AIC <- df * 2 - 2 * loglik
  c(AIC = AIC, BIC = BIC)
}

# Optimized fit.MLE function
fit.MLE <- function(x, y, func, N.params, sigma=5, iter=1e4, metropolis.scale=2, logpost=log.lik.post, 
  params.start=NULL, method='Nelder-Mead', lower=NULL, upper=NULL, MLE.fun.attempts=100, RWm=TRUE) {
  
  if (is.null(params.start)) {
     start <- start_optimization(x=x, y=y, func=func, N.params=N.params, SS.fun=SS.fun, lower=rep(0, N.params), upper=upper)
     params.start <- c(start, sigma)
  }
  N <- length(y)
  k <- length(params.start) - 1
  attempts <- 0
  
  l.fit <- list(int = NaN, converge = FALSE)
  while (!(!is.nan(l.fit$int) && l.fit$converge) && (attempts < MLE.fun.attempts)) {
    st1 <- if (is_product_function(func)) params.start * runif(k + 1) else params.start
    l.fit <- suppressWarnings(
      MLE.fun(logpost=logpost, par=st1, x=x, y=y, lower=lower, upper=upper, method=method, func=func)
    )
    attempts <- attempts + 1
  }
  
  if (!l.fit$converge) stop('MLE.fun does not converge')
  
  e.values <- Re(eigen.values(l.fit$var))
  if (any(e.values <= 0)) stop('MLE.fun covariance matrix not positive definite')

  if (RWm) {
    rw.fit <- RWmetropolis(logpost = logpost, cov.mat = l.fit$var, scale = metropolis.scale, start = l.fit$fit, m = iter, x = x, y = y, func = func)
    
    acceptance.rate <- rw.fit$accept
    parameters <- apply(rw.fit$par, 2, median)
    
    fits <- parameters[1:k]
    fits.se <- apply(rw.fit$par, 2, sd)[1:k]
  } else {
    fits <- l.fit$fit[1:k]
    fits.se <- sqrt(diag(l.fit$var))[1:k]
  }
  
  paramsfit <- fits
  res <- residFunCpp(paramsfit, y, x, func)
  gof.se <- (sum(res^2) / (N - k))^0.5
  msc <- model.selection.criteria(paramsfit, x, y, func)

  if (RWm) {
    list(fits = fits, fits.se = fits.se, acceptance.rate = acceptance.rate, MLE.fun.fit = l.fit$fit, MLE.fun.fit.se = l.fit$fit.se, MLE.fun.var = l.fit$var, MLE.fun.int = l.fit$int, MLE.fun.converge = l.fit$converge, MLE.fun.convergence.attempts = attempts, gof = gof.se, AIC = msc[1], BIC = msc[2])
  } else {
    list(fits = fits, fits.se = fits.se, gof = gof.se, AIC = msc[1], BIC = msc[2], model.info = l.fit$converge, model.message = l.fit$int)
  }
}

# Wrapper for fit.MLE
FIT.MLE <- function(x, y, func, N.params, sigma=5, iter=1e4, metropolis.scale=2, logpost=log.lik.post, params.start=NULL, 
  method='Nelder-Mead', lower=NULL, upper=NULL, MLE.fun.attempts=100, fit.attempts=10, RWm=TRUE) {

  
  bounds <- check_and_set_bounds(x=x, y=y, func=func, N.params=N.params, upper=upper, lower=lower)
  lower <- bounds$lower
  upper <- bounds$upper

  if (method %in% c('L-BFGS-B', 'Brent')) {
    tol_low <- 1e-08
    lower <- if (any(lower <= 0)) lower + tol_low
  }

  output.MLE <- NULL
  attempts <- 0
  while (is.null(output.MLE) && attempts < fit.attempts) {
    tryCatch({
      output.MLE <- fit.MLE(x=x, y=y, func=func, N.params=N.params, sigma=sigma, iter=iter, metropolis.scale=metropolis.scale, logpost=logpost, params.start=params.start, 
        method=method, lower=lower, upper=upper, MLE.fun.attempts=MLE.fun.attempts, RWm=RWm)
    }, error = function(e) {})
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
  pars <- ests_fun(x = x, y = y, showplot = FALSE)
  if (identical(func, product3)) {
    pars1 <- c(0.25 * pars, 15)
    if (!is.null(upper)) pars1 <- pmin(pars1, upper[1:4])
    pars2 <- c(0.5 * pars, 15)
    if (!is.null(upper)) pars2 <- pmin(pars2, upper[5:8]) 
    pars3 <- c(pars, 15)
    if (!is.null(upper)) pars3 <- pmin(pars3, upper[9:12]) 
    if (start.method=='uniform'){
      par <- (c(pars1, pars2, pars3) * runif(N.params))
    }
    else if (start.method=='lognormal'){
      par <- generate_lognormal_samples(means = c(pars1, pars2, pars3), cv = cv, n = 1)
    }

  } else if (identical(func, product2)) {
    pars1 <- c(0.25 * pars, 15)
    if (!is.null(upper)) pars1 <- pmin(pars1, upper[1:4])
    pars2 <- c(pars, 15)
    if (!is.null(upper)) pars2 <- pmin(pars2, upper[5:8]) 
    if (start.method=='uniform'){
      par <- (c(pars1, pars2) * runif(N.params))
    }
    else if (start.method=='lognormal'){
      par <- generate_lognormal_samples(means = c(pars1, pars2), cv = cv, n = 1)
    }

  } else if (identical(func, product1)) {
    pars1 <- c(pars, 15)
    if (!is.null(upper)) pars1 <- pmin(pars1, upper[1:4])
    if (start.method=='uniform'){
      par <- pars1 * runif(N.params)
    }
    else if (start.method=='lognormal'){
      par <- generate_lognormal_samples(means = pars1, cv = cv, n = 1)
    }

  }
  est <- list(par=par)
  return(est)
}

# function to check if a given function is product1 or product2
is_product_function <- function(func) {
  identical(func, product1) || identical(func, product2) || identical(func, product3)
}

check_and_set_bounds <- function(x, y, func, N.params, upper = NULL, lower = NULL) {
  if (is_product_function(func)) {
    if (is.null(upper)) {
      ests <- ests_fun(x = x, y = y)
      if (identical(func, product2)) {
        ests <- c(ests[c(1, 3, 3)], Inf, ests[c(1, 3, 3)], Inf)
        upper <- sapply(5 * ests, round_up)
        upper[c(2, 6)] <- Inf
      } else if (identical(func, product1)) {
        ests <- c(ests[c(1, 3, 3)], Inf)
        upper <- sapply(5 * ests, round_up)
        upper[2] <- Inf
      } else if (identical(func, product3)){
        ests <- c(ests[c(1, 3, 3)], Inf, ests[c(1, 3, 3)], Inf, ests[c(1, 3, 3)], Inf)
        upper <- sapply(5 * ests, round_up)
        upper[c(2, 6, 10)] <- Inf
      }  
    }
    if (is.null(lower)) lower <- rep(0, N.params)
  } else {
    upper <- rep(Inf, N.params)
    lower <- rep(-Inf, N.params)
  }
  return(list(upper = upper, lower = lower))
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

fit.LM <- function(x, y, func,  N.params, lower=NULL, upper=NULL) {
  # est <- optim(runif(N.params), SS.fun, method = 'L-BFGS-B', lower = lower, upper = upper, hessian = TRUE, y = y, x = x, func = func)
  
  start <- start_optimization(x=x, y=y, func=func, N.params=N.params, SS.fun=SS.fun, lower=lower, upper=upper)
  model <- minpack.lm::nls.lm(par=start, fn=residFunCpp, y=y, x=x, lower=lower, upper=upper, control=list(maxiter=1024, tol=1e-05, warnOnly=TRUE), func=func)
  e.values <- Re(eigen.values(model$hessian))
  if (any(e.values <= 0)) stop('Hessian not positive definite')
  fits <- summary(model)$coefficients[, 'Estimate']
  fits.se <- summary(model)$coefficients[, 'Std. Error']
  gof.se <- summary(model)$sigma
  msc <- model.selection.criteria(fits, x, y, func)
  list(fits = fits, fits.se = fits.se, gof = gof.se, AIC = msc[1], BIC = msc[2], model.info = model$info, model.message = model$message)
}

round_up <- function(x, factor=10) factor * ceiling( x / factor)
  
# Wrapper for fit.LM
FIT.LM <- function(x, y, func, N.params, lower=NULL, upper=NULL, fit.convergence.attempts = 10, fit.attempts = 10) {
  convergence.attempts <- 0
  modelinfo <- NULL
  
  bounds <- check_and_set_bounds(x=x, y=y, func=func, N.params=N.params, upper=upper, lower=lower)
  lower <- bounds$lower
  upper <- bounds$upper
  
  while ((is.null(modelinfo) || !any(modelinfo == c(1, 2, 3, 4))) && convergence.attempts < fit.convergence.attempts) {
    output.LM <- NULL
    attempts <- 0
    while (is.null(output.LM) && attempts < fit.attempts) {
      tryCatch({
        output.LM <- fit.LM(x=x, y=y, func=func, N.params=N.params, lower=lower, upper=upper)
      }, error = function(e) {})
      attempts <- attempts + 1
    }
    modelinfo <- output.LM$model.info
    convergence.attempts <- convergence.attempts + 1
  }
  c(output.LM, convergence.attempts = convergence.attempts)
}

fit.bf.LM <- function(x, y, func, N.params, lower=NULL, upper=NULL) {
  
  start <- start_optimization(x=x, y=y, func=func, N.params=N.params, SS.fun=SS.fun, lower=lower, upper=upper)

  model <- minpack.lm::nls.lm(par=start, fn=residFunCpp, y=y, x=x, lower=lower, upper=upper, 
    control = list(maxiter = 1024, tol = 1e-05, warnOnly = TRUE), func=func)
  e.values <- Re(eigen.values(model$hessian))
  if (any(e.values <= 0)) stop('Hessian not positive definite')
  fits <- summary(model)$coefficients[, 'Estimate']
  fits.se <- summary(model)$coefficients[, 'Std. Error']
  gof.se <- summary(model)$sigma
  msc <- model.selection.criteria(fits, x, y, func)
  list(fits=fits, fits.se=fits.se, gof=gof.se, AIC=msc[1], BIC=msc[2], model.info=model$info, model.message=model$message)
}

# Wrapper for fit.bf.LM
FIT.bf.LM <- function(x, y, func, N.params, lower=NULL, upper=NULL, fit.convergence.attempts = 10, fit.attempts = 10) {
  convergence.attempts <- 0
  modelinfo <- NULL

  bounds <- check_and_set_bounds(x=x, y=y, func=func, N.params=N.params, upper=upper, lower=lower)
  lower <- bounds$lower
  upper <- bounds$upper

  while ((is.null(modelinfo) || !any(modelinfo == c(1, 2, 3, 4))) && convergence.attempts < fit.convergence.attempts) {
    output.LM <- NULL
    attempts <- 0
    while (is.null(output.LM) && attempts < fit.attempts) {
      tryCatch({
        output.LM <- fit.bf.LM(x=x, y=y, func=func, N.params=N.params, lower=lower, upper=upper)
      }, error = function(e) {})
      attempts <- attempts + 1
    }
    modelinfo <- output.LM$model.info
    convergence.attempts <- convergence.attempts + 1
  }
  c(output.LM, convergence.attempts = convergence.attempts)
}

# Wrapper for FIT.bf.LM
FIT.BF.LM <- function(x, y, func, N.params, lower = NULL, upper = NULL, fit.convergence.attempts = 10, fit.attempts = 10, N.repeat = 25) {
  fits <- matrix(NA, N.repeat, N.params)
  fits.se <- matrix(NA, N.repeat, N.params)
  gof.se <- rep(NA, N.repeat)
  BIC <- rep(NA, N.repeat)
  AIC <- rep(NA, N.repeat)
  fit.convergence.attempts <- rep(NA, N.repeat)
  model.info <- rep(NA, N.repeat)
  model.message <- rep(NA, N.repeat)
  convergence.attempts <- rep(NA, N.repeat)
  
  for (iii in 1:N.repeat) {
    output1 <- FIT.bf.LM(x=x, y=y, func=func, N.params=N.params, lower=lower, upper=upper, fit.convergence.attempts=fit.convergence.attempts, fit.attempts=fit.attempts)
    if (!is.null(output1)) {
      fits[iii, ] <- output1$fits
      fits.se[iii, ] <- output1$fits.se
      gof.se[iii] <- output1$gof
      AIC[iii] <- output1$AIC
      BIC[iii] <- output1$BIC
      fit.convergence.attempts[iii] <- output1$convergence.attempts
      model.info[iii] <- output1$model.info
      model.message[iii] <- output1$model.message
      convergence.attempts[iii] <- output1$convergence.attempts
    }
  }
  
  ind <- which.min(gof.se)
  
  list(fits = fits[ind, ], fits.se = fits.se[ind, ], gof = gof.se[ind], AIC = AIC[ind], BIC = BIC[ind], model.info = model.info[ind], model.message = model.message[ind], convergence.attempts = convergence.attempts[ind])
}

# Algorithm = c('default', 'port')
fit.NLS <- function(x, y, func, form, N.params, lower=NULL, upper=NULL, algorithm='default') {
  
  param.names <- all.vars(form[[3]])

  if (is_product_function(func)){
    if (algorithm == 'port'){ 
      start <- start_optimization(x=x, y=y, func=func, N.params=N.params, SS.fun=SS.fun, lower=lower, upper=upper)
    }else{          # default is Gauss-Newton
      est <- start_optimiser(runif(N.params), SS.fun=SS.fun, lower=rep(0, N.params), upper=upper, y=y, x=x, func=func, max_attempts=10)
      start <- as.list(est$par)
    }
  }else{
    start <- start_optimization(x=x, y=y, func=func, N.params=N.params, SS.fun=SS.fun, lower=lower, upper=upper)   
  }
  names(start) <- param.names[!param.names == 'x']

  model <- suppressWarnings(nls(form, start=start, data=data.frame(x=x, y=y), algorithm=algorithm, lower=lower, upper=upper, control=list(maxiter=1024, tol=1e-05, warnOnly=TRUE)))
  
  e.values <- Re(eigen.values(summary(model)$cov.unscaled))
  if (any(e.values <= 0)) stop('Hessian not positive definite')
  
  fits <- summary(model)$coefficients[, 'Estimate']
  fits.se <- summary(model)$coefficients[, 'Std. Error']
  gof.se <- summary(model)$sigma
  msc <- model.selection.criteria(fits, x, y, func)
  list(fits=fits, fits.se=fits.se, gof=gof.se, AIC=msc[1], BIC=msc[2], model.info=model$convInfo$stopCode, model.message=model$convInfo$stopMessage)
}

# Wrapper for fit.NLS
FIT.NLS <- function(x, y, func, form, N.params, lower=NULL, upper=NULL, algorithm='default', fit.convergence.attempts=10, fit.attempts=10) {
  convergence.attempts <- 0
  modelinfo <- NULL
  
  bounds <- check_and_set_bounds(x=x, y=y, func=func, N.params=N.params, upper=upper, lower=lower)
  lower <- bounds$lower
  upper <- bounds$upper
  lower <- if (algorithm == 'port') bounds$lower else rep(-Inf, N.params)

  while ((is.null(modelinfo) || !any(modelinfo == c(0, 1, 2, 3, 4))) && convergence.attempts < fit.convergence.attempts) {
    output.NLS <- NULL
    attempts <- 0
    while (is.null(output.NLS) && attempts < fit.attempts) {
      tryCatch({
        output.NLS <- fit.NLS(x=x, y=y, func=func, form=form, N.params=N.params, lower=lower, upper=upper, algorithm=algorithm)
      }, error = function(e) {})
      attempts <- attempts + 1
    }
    modelinfo <- output.NLS$model.info
    convergence.attempts <- convergence.attempts + 1
  }
  c(output.NLS, convergence.attempts = convergence.attempts)
}

# Algorithm = c('default', 'port')
fit.NLSrobust <- function(x, y, func, form,  N.params, lower=NULL, upper=NULL, algorithm = 'default', method = 'M') {
  if (method == 'M') {
    if (is_product_function(func)){
      est <- start_optimiser(N.params=N.params, SS.fun=SS.fun, lower=lower, upper=upper, y=y, x=x, func=func, max_attempts=10)
    }else{
      est <- optim(par=runif(N.params), SS.fun, method="L-BFGS-B", lower=lower, upper=upper, hessian=TRUE,  y=y, x=x, func=func) 
    }
    start <- as.list(est$par)
    param.names <- all.vars(form[[3]])
    names(start) <- param.names[!param.names == 'x']
    mod <- suppressWarnings(robustbase::nlrob(form=form, start=start, data=data.frame(x=x, y=y), algorithm=algorithm, method=method, control=robustbase:::nlrob.control(method), lower=lower, upper=upper))
  } else {
    mod <- suppressWarnings(robustbase::nlrob(form=form, data = data.frame(x=x, y=y), method=method, control=robustbase:::nlrob.control(method), lower=lower, upper=upper))
  }

  e.values <- Re(eigen.values(vcov(mod)))
  if (any(e.values <= 0)) stop('Hessian not positive definite')
  
  fits <- summary(mod)$coefficients[, 'Estimate']
  fits.se <- summary(mod)$coefficients[, 'Std. Error']
  gof.se <- summary(mod)$Scale
  msc <- model.selection.criteria(fits, x, y, func)
  list(fits = fits, fits.se = fits.se, gof = gof.se, AIC = msc[1], BIC = msc[2], model.message = summary(mod)$status)
}

# Wrapper for fit.NLSrobust
FIT.NLSrobust <- function(x, y, func, form, N.params, lower=NULL, upper=NULL, algorithm = 'default', method = 'M', fit.convergence.attempts = 10, fit.attempts = 10) {
  
  convergence.attempts <- 0
  modelinfo <- NULL
  
  bounds <- check_and_set_bounds(x=x, y=y, func=func, N.params=N.params, upper=upper, lower=lower)
  lower <- bounds$lower
  upper <- bounds$upper

  while ((is.null(modelinfo) || modelinfo != 'converged') && convergence.attempts < fit.convergence.attempts) {
    output.NLS <- NULL
    attempts <- 0
    while (is.null(output.NLS) && attempts < fit.attempts) {
      tryCatch({
        output.NLS <- fit.NLSrobust(x=x, y=y, func=func, form=form, N.params=N.params, lower=lower, upper=upper, algorithm=algorithm, method=method)
      }, error = function(e) {})
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
product_rise_and_decay_percent <- function(tau1, tau2, interval=c(0.1, 0.9), showplot=FALSE) {
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
      cat("Failed to find roots for p =", p, "\n")
      cat("Interval1:", interval1, "\n")
      cat("Interval2:", interval2, "\n")
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

product3 <- function(params, x) {
  product1(params[1:4], x) + product1(params[5:8], x) + product1(params[9:12], x)
}


sign_fun <- function(y, direction_method = c('smooth', 'regression', 'cumsum'), k=5) {
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
    # Find the peak value
    peak_value <- max(abs(smoothed_signal), na.rm = TRUE)
    peak_value <- smoothed_signal[which.max(abs(smoothed_signal))]
    
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
    peak_value <- max(abs(fitted_values), na.rm = TRUE)
    peak_value <- fitted_values[which.max(abs(fitted_values))]
    
  } else if (method == 'cumsum') {
    # Calculate the cumulative sum of the signal
    cumsum_signal <- cumsum(y)
    # Find the peak value of the cumulative sum
    peak_value <- max(abs(cumsum_signal), na.rm = TRUE)
    peak_value <- cumsum_signal[which.max(abs(cumsum_signal))]
  }
  
  # Determine if the peak is positive or negative
  peak_direction <- ifelse(peak_value > 0, 1, -1)
  
  # Return the sign of the peak direction
  return(peak_direction[[1]])
}

out.fun <- function(params, interval = c(0.1, 0.9), dp = 3, sign=1) {
  
  # Calculate derived metrics
  A1 <- sign*params[[1]]
  tau.rise <- tau_rise(params[[2]], params[[3]])
  tau.decay <- params[[3]]
  tpeak <- product_tpeak(params[[2]], params[[3]])
  
  diffs <- product_rise_and_decay_percent(tau1 = params[[2]], tau2 = params[[3]], interval = interval, showplot = FALSE)
  trise.percent <- diffs[1]
  tdecay.percent <- diffs[2]
  
  area1 <- product_area(params[[1]], params[[2]], params[[3]])
  delay1 <- params[[4]]
  
  # Calculate rise and decay interval labels
  r_percent_start <- interval[1] * 100
  r_percent_end <- interval[2] * 100
  r_label <- paste0('r', r_percent_start, '_', r_percent_end)
  
  d_label <- paste0('d', interval[2] * 100, '_', interval[1] * 100)
  
  # Create data frame with calculated values
  df_output <- data.frame(
    amp = round(A1, dp),
    τrise = round(tau.rise, dp),
    τdecay = round(tau.decay, dp),
    tpeak = round(tpeak, dp),
    delay = round(delay1, dp),
    area = round(area1, dp)
  )
  
  # Add dynamic columns for rise and decay percentages
  df_output[[r_label]] <- round(trise.percent, dp)
  df_output[[d_label]] <- round(tdecay.percent, dp)
    # Reorder columns to match the desired output format
  df_output <- df_output[c('amp', 'τrise', 'τdecay', 'tpeak', r_label, d_label, 'delay', 'area')]
  
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
      bounds[2] <- adjust_tau(bounds[2], bounds[3], upper=upper)
    } else if (identical(func, product2)) {
      bounds[2] <- adjust_tau(bounds[2], bounds[3], upper=upper)
      bounds[6] <- adjust_tau(bounds[6], bounds[7], upper=upper)
    } else if (identical(func, product3)) {
      bounds[2] <- adjust_tau(bounds[2], bounds[3], upper=upper)
      bounds[6] <- adjust_tau(bounds[6], bounds[7], upper=upper)
      bounds[10] <- adjust_tau(bounds[10], bounds[11], upper=upper)
    }
  }
  return(bounds)
}

FIT <- function(response, dt=0.1, func=product2, method= c('LM', 'BF.LM', 'GN', 'port', 'robust', 'MLE'), 
  stimulation_time=0, baseline=0, fast.decay.limit=NULL, latency.limit=NULL, lower=NULL, upper=NULL, filter=FALSE, fc=1000, interval=c(0.1, 0.9), 
  MLEsettings=list(iter=1e4, metropolis.scale=1.5, fit.attempts=100, RWm=FALSE), 
  MLE.method=c('L-BFGS-B', 'Nelder-Mead', 'BFGS','CG', 'SANN', 'Brent'),
  response_sign_method = c('smooth', 'regression', 'cumsum'), 
  dp=3, lwd=1.2, xlab='time (ms)', ylab='PSC (pA)', width=5, height=5, 
  return.output=FALSE, show.output=TRUE, show.plot=TRUE){

  dx <- dt 
  method <- match.arg(method)
  y <- response
  if (all(is.na(y[(which(!is.na(y))[length(which(!is.na(y)))] + 1):length(y)]))) {
    y <- y[!is.na(y)]
  }

  sign <- sign_fun(y, direction_method=response_sign_method) 
  y <- sign * y

  x <- seq(0, (length(y) - 1) * dx, by = dx)

  if (filter){
    ind = 20
    fc = fc; fs = 1/dx*1000; bf <- butter(2, fc/(fs/2), type='low')
    yfilter <- signal::filter(bf, y)
  } else {
    ind=1
    yfilter=y
  }
  x.orig <- x

  ind1 <- (stimulation_time - baseline)/dx
  ind2 <- baseline/dx
  
  yorig <- y[ind1:length(y)]
  yfilter <- yfilter[ind1:length(yfilter)]
  xorig <- seq(0, dx * (length(yorig) - 1), by = dx)

  if (is_product_function(func)) {
    yorig <- yorig - mean(yorig[1:ind2])
    yfilter <- yfilter - mean(yfilter[1:ind2])

    y2fit <- yfilter[ind2:length(yfilter)]
    x2fit <- seq(0, dx * (length(y2fit) - 1), by = dx)

    # check if both upper and fast.decay.limit are specified
    if ( !is.null(upper) && (!is.null(fast.decay.limit) || !is.null(latency.limit)) ) {
      warning("both 'upper' boundaries and at either 'fast decay limit' and/or 'latency limit' are specified:\ninputs 'fast.decay.limit' and/or 'latency.limit' will be ignored")
      fast.decay.limit <- NULL
      latency.limit <- NULL
    }

    # adjusts bounds for tau_rise to correct form involving tau1
    upper <- adjust_product_bounds(bounds=upper, func=func, upper=TRUE)
    lower <- adjust_product_bounds(bounds=lower, func=func)
      
    if (is.null(upper) && (!is.null(fast.decay.limit) || !is.null(latency.limit))){

      if (identical(func, product1)){
        upper <- c(sign * Inf, Inf, Inf, Inf) 
      } else if (identical(func, product2)) {
        upper <- c(sign * Inf, Inf, Inf, Inf, sign * Inf, Inf, Inf, Inf)
      } else if (identical(func, product3)) {
        upper <- c(sign * Inf, Inf, Inf, Inf, sign * Inf, Inf, Inf, Inf, sign * Inf, Inf, Inf, Inf)
      }

      if (!is.null(fast.decay.limit)){
        upper[3] <- fast.decay.limit[1]
        if (identical(func, product3)){
          upper[7] <-  if (length(fast.decay.limit)==1) fast.decay.limit[1] else fast.decay.limit[2]
        }
      }

      if (!is.null(latency.limit)){
        if (identical(func, product1)){
          upper[4] <- latency.limit
        } else if (identical(func, product2)){
          upper[c(4,8)] <- latency.limit
        } else if (identical(func, product3)){
          upper[c(4,8,12)] <- latency.limit
        }
      }

    }
  
    if (identical(func, product1)){
      param_names <- c('a', 'b', 'c', 'd')
      N.params <- length(param_names)
      form.2.fit <- y ~ (a / (((b / (b + c)) ^ (b / c)) * c / (b + c))) * (1 - exp(-(x - d) / b)) * exp(-(x - d) / c) * (x >= d)

      if (!is.null(lower) && length(lower) != N.params) stop("number of lower boundaries should be equal to the number of parameters")
      if (!is.null(upper) && length(upper) != N.params) stop("number of upper boundaries should be equal to the number of parameters")
      if (!is.null(lower)) lower[1] <- sign * lower[1] 
      if (!is.null(upper)) upper[1] <- sign * upper[1] 

    } else if (identical(func, product2)){
      param_names <- c('a', 'b', 'c', 'd', 'e', 'f', 'g', 'h')
      N.params <- length(param_names)
      form.2.fit <- y ~ (a / (((b / (b + c)) ^ (b / c)) * c / (b + c))) * (1 - exp(-(x - d) / b)) * exp(-(x - d) / c) * (x >= d) + 
                        (e / (((f / (f + g)) ^ (f / g)) * g / (f + g))) * (1 - exp(-(x - h) / f)) * exp(-(x - h) / g) * (x >= h)
      if (!is.null(lower) && length(lower) != N.params) stop("number of lower boundaries should be equal to the number of parameters")
      if (!is.null(upper) && length(upper) != N.params) stop("number of upper boundaries should be equal to the number of parameters")
      if (!is.null(lower)) { lower[1] <- sign * lower[1]; lower[5] <- sign * lower[5] }
      if (!is.null(upper)) { upper[1] <- sign * upper[1]; upper[5] <- sign * upper[5] }
    
    } else if (identical(func, product3)){
      param_names <- c('a', 'b', 'c', 'd', 'e', 'f', 'g', 'h', 'i', 'j', 'k', 'l')
      N.params <- length(param_names)
      form.2.fit <- y ~ (a / (((b / (b + c)) ^ (b / c)) * c / (b + c))) * (1 - exp(-(x - d) / b)) * exp(-(x - d) / c) * (x >= d) + 
                        (e / (((f / (f + g)) ^ (f / g)) * g / (f + g))) * (1 - exp(-(x - h) / f)) * exp(-(x - h) / g) * (x >= h) +
                        (i / (((j / (j + k)) ^ (j / k)) * g / (j + k))) * (1 - exp(-(x - l) / j)) * exp(-(x - l) / k) * (x >= l)
      if (!is.null(lower) && length(lower) != N.params) stop("number of lower boundaries should be equal to the number of parameters")
      if (!is.null(upper) && length(upper) != N.params) stop("number of upper boundaries should be equal to the number of parameters")
      if (!is.null(lower)) { lower[1] <- sign * lower[1]; lower[5] <- sign * lower[5]; lower[9] <- sign * lower[9] }
      if (!is.null(upper)) { upper[1] <- sign * upper[1]; upper[5] <- sign * upper[5]; upper[9] <- sign * upper[9] }
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

    output <- FIT.MLE(x=x2fit, y=y2fit, func=func, N.params=N.params, sigma=sd.est, iter=iter, metropolis.scale=metropolis.scale, 
      logpost=log.lik.post, params.start=NULL, method=MLE.method, lower=lower, upper=upper, MLE.fun.attempts=100, fit.attempts=fit.attempts, RWm=RWm) 

  } else if (method == 'LM'){
     output <- FIT.LM(x=x2fit, y=y2fit, func=func, N.params=N.params, lower=lower, upper=upper, fit.convergence.attempts=10, fit.attempts=10)

  } else if (method == 'BF.LM'){
    output <- FIT.bf.LM(x=x2fit, y=y2fit, func=func, N.params=N.params, lower=lower, upper=upper, fit.convergence.attempts=10, fit.attempts=10)

  } else if (method == 'GN'){
    output <- FIT.NLS(x=x2fit, y=y2fit, func=func, form=form.2.fit, N.params=N.params, lower=lower, upper=upper)

  } else if (method == 'port'){
    output <- FIT.NLS(x=x2fit, y=y2fit, func=func, form=form.2.fit, N.params=N.params, lower=lower, upper=upper, algorithm='port')

  } else if (method == 'robust'){
    algorithm <- if (!is.null(lower) || !is.null(upper)) 'port' else 'default'
    # currently only 'M' working
    algorithm <-  'port'
    output <- FIT.NLSrobust(x=x2fit, y=y2fit, func=func, form=form.2.fit, N.params=N.params, lower=lower, upper=upper, algorithm=algorithm, method='M')
  }
    

  if (is_product_function(func)) {
    if (identical(func, product1)){
      df_output <- out.fun(params=output$fits[1:4], interval=interval, dp=dp, sign=sign)
      output$fits[1] <- sign * output$fits[1]
    } else if (identical(func, product2)){  
      df_output1 <- out.fun(params=output$fits[1:4], interval=interval, dp=dp, sign=sign)
      df_output2 <- out.fun(params=output$fits[5:8], interval=interval, dp=dp, sign=sign)
      df_output <- if (df_output1[3] < df_output2[3]) rbind('fast' = df_output1, 'slow' = df_output2) else rbind('fast' = df_output2, 'slow' = df_output1)
      output$fits[1] <- sign * output$fits[1]; output$fits[5] <- sign * output$fits[5]
      if (output$fits[3] > output$fits[7]){
        output$fits <- c(output$fits[5:8], output$fits[1:4])
        output$fits.se <- c(output$fits.se[5:8], output$fits.se[1:4])
      }
    } else if (identical(func, product3)){  
      # Compute the output data frames
      df_output1 <- out.fun(params=output$fits[1:4],  interval=interval, dp=dp, sign=sign)
      df_output2 <- out.fun(params=output$fits[5:8],  interval=interval, dp=dp, sign=sign)
      df_output3 <- out.fun(params=output$fits[9:12], interval=interval, dp=dp, sign=sign)

      # Create a list of the outputs
      output_list <- list(df_output1, df_output2, df_output3)

      # Extract the third elements from each output and determine the order
      order_indices <- order(sapply(output_list, function(x) x[[3]]))

      # Reorder the output list based on the third element
      output_ordered <- output_list[order_indices]

      # Combine the outputs in the correct order
      df_output <- rbind('fast' = output_ordered[[1]], 'medium' = output_ordered[[2]], 'slow' = output_ordered[[3]])

      # Adjust the order of the fits and fits.se based on the sorted order
      output$fits <- c(output$fits[(4*order_indices[1]-3):(4*order_indices[1])],
                       output$fits[(4*order_indices[2]-3):(4*order_indices[2])],
                       output$fits[(4*order_indices[3]-3):(4*order_indices[3])])

      output$fits.se <- c(output$fits.se[(4*order_indices[1]-3):(4*order_indices[1])],
                          output$fits.se[(4*order_indices[2]-3):(4*order_indices[2])],
                          output$fits.se[(4*order_indices[3]-3):(4*order_indices[3])])

      # adjust signs
      output$fits[1] <- sign * output$fits[1]
      output$fits[5] <- sign * output$fits[5]
      output$fits[9] <- sign * output$fits[9]
    }

    traces=data.frame(x=xorig, y= sign * yorig, yfilter= sign * yfilter)
    fits <- output$fits
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


fit_plot <- function(traces, func=product2, xlab='time (ms)', ylab='PSC amplitude (pA)', xlim=NULL, ylim=NULL, bl=NULL, lwd=1.2, filter=FALSE, width=4, height=4, filename='trace.svg', save=FALSE) {
  
  if (save) {
    # Open SVG device
    svg(file=filename, width=width, height=height, bg='transparent')
  } else {
    dev.new(width=width, height=height, noRStudioGD=TRUE)
  }

  plot(traces$x, traces$y, col='gray', xlab=xlab, ylab=ylab, xlim=xlim, ylim=ylim, type='l', bty='l', las=1, lwd=lwd, main='')
  
  if (filter) {
    lines(traces$x, traces$yfilter, col='black', type='l', lwd=lwd)
  }
  
  lines(traces$x, traces$yfit, col='indianred', lty=3, lwd=2 * lwd)
  
  if (identical(func, product2)) {
    lines(traces$x, traces$yfit1, col='#4C78BC', lty=3, lwd=2 * lwd)
    lines(traces$x, traces$yfit2, col='#CA92C1', lty=3, lwd=2 * lwd)
  }

  if (identical(func, product3)) {
    lines(traces$x, traces$yfit1, col='#F28E2B', lty=3, lwd=2 * lwd)
    lines(traces$x, traces$yfit2, col='#4C78BC', lty=3, lwd=2 * lwd)
    lines(traces$x, traces$yfit3, col='#CA92C1', lty=3, lwd=2 * lwd)
  }

  if (!is.null(bl)) abline(v=bl, col='black', lwd=lwd, lty=3)

  if (save) {
    dev.off()
  }
}

nFIT <- function(response, n=30, dt=0.1, func=product2, method= c('LM', 'BF.LM', 'GN', 'port', 'robust', 'MLE'), 
  stimulation_time=0, baseline=0, fast.decay.limit=NULL, fast.constraint=FALSE, fast.constraint.method=c('rise', 'peak'), first.delay.constraint=FALSE,
  latency.limit=NULL, lower=NULL, upper=NULL, filter=FALSE, fc=1000, interval=c(0.1, 0.9), 
  MLEsettings=list(iter=1e4, metropolis.scale=1.5, fit.attempts=100, RWm=FALSE),  MLE.method=c('L-BFGS-B', 'Nelder-Mead', 'BFGS','CG', 'SANN', 'Brent'), 
  response_sign_method = c('smooth', 'regression', 'cumsum'), dp=3, lwd=1.2, xlab='time (ms)', ylab='PSC (pA)', width=5, height=5, return.output=FALSE, show.output=TRUE, show.plot=TRUE, seed=42) {
  
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
      
      FIT(response=response, dt=dt, func=func, method=method, stimulation_time=stimulation_time, baseline=baseline, fast.decay.limit=fast.decay.limit, 
                      latency.limit=latency.limit, lower=lower, upper=upper, filter=filter, fc=fc, interval=interval, MLEsettings=MLEsettings, 
                      MLE.method=MLE.method, response_sign_method=response_sign_method, dp=dp, lwd=lwd, xlab=xlab, ylab=ylab, width=width, height=height, 
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
      
      if (identical(func, product2)){

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

      }else if (identical(func, product3)){
      
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
  
  if (show.plot) fit_plot(traces=traces, func=func, xlab=xlab, ylab=ylab, lwd=lwd, filter=filter, width=width, height=height)
  
  if (show.output) print(df_output)
  
  if (return.output) {
    return(list(output = df_output, fits = output$fits, fits.se = output$fits.se, gof = output$gof, AIC = output$AIC, BIC = output$BIC, model.message = output$model.message, sign=sign, traces=traces, fit_results=fit_results))
  }
}

# use start fun to generate initial estimated values
# ests_fun <- function(x, y, k = 10, showplot=FALSE) {
#   n <- length(y)
  
#   ysmooth <- rep(NA, n)
#   for (i in 1:(n - k + 1)) {
#     ysmooth[i + floor(k/2)] <- mean(y[i:(i + k - 1)])
#   }

#   # fc = 300, dx = 0.1
#   # Smooth the signal using a Butterworth filter
#   # fs <- 1 / dx * 1000
#   # bf <- butter(5, fc / (fs / 2), type = 'low')
#   # ysmooth <- filter(bf, y)
  
#   # Find the peak value
#   idx.peak <- which.max(abs(ysmooth))
#   t.peak <- x[idx.peak]
#   A <- ysmooth[idx.peak]
  
#   # Calculate 20% and 80% of the peak amplitude
#   A.20 <- 0.20 * A
#   A.80 <- 0.80 * A
  
#   # Find all indices where y crosses 20% and 80% amplitude levels during rising phase
#   rise_idx.20 <- which(y[1:idx.peak] >= A.20)
#   rise_idx.80 <- which(y[1:idx.peak] >= A.80)
  
#   # Find all indices where y crosses 20% and 80% amplitude levels during decaying phase
#   decay_idx.80 <- which(y[idx.peak:length(y)] <= A.80) + idx.peak - 1
#   decay_idx.20 <- which(y[idx.peak:length(y)] <= A.20) + idx.peak - 1
  
#   # Calculate times for these points
#   t.rise.20 <- x[rise_idx.20]
#   t.rise.80 <- x[rise_idx.80]
#   t.decay.80 <- x[decay_idx.80]
#   t.decay.20 <- x[decay_idx.20]
  
#   # Average the times if there are multiple points
#   avg_t.rise.20 <- median(t.rise.20)
#   avg_t.rise.80 <- median(t.rise.80)
#   avg_t.decay.80 <- median(t.decay.80)
#   avg_t.decay.20 <- median(t.decay.20)
  
#   # Calculate 20-80% rise time and 80-20% decay time
#   rise.time <- avg_t.rise.80 - avg_t.rise.20
#   decay.time <- avg_t.decay.20 - avg_t.decay.80

#   # Plot the original and smoothed signals
#     # Find roots for both p values and plot if required
#   if (showplot) {
#     plot(x, y, col = 'indianred', xlab = 'x', type = 'l', bty = 'l', las = 1, main = '')
#     lines(x, ysmooth, col = 'lightgray')
#     abline(h = A.20, col = 'black', lty = 3)
#     abline(h = A.80, col = 'black', lty = 3)
#   }

#   return(c(A, rise.time, decay.time))
# }

# use start fun to generate initial estimated values
ests_fun <- function(x, y, showplot = FALSE) {
  n <- length(y)
  
  dx = x[2]-x[1]
  k <- round(ceiling(1/dx))

  ysmooth <- rep(NA, n)
  for (i in 1:(n - k + 1)) {
    ysmooth[i + floor(k/2)] <- mean(y[i:(i + k - 1)])
  }

  # Find the peak value
  idx.peak <- which.max(abs(ysmooth))
  t.peak <- x[idx.peak]
  A <- ysmooth[idx.peak]
  
  # Calculate 20% and 80% of the peak amplitude
  A.20 <- 0.20 * A
  A.80 <- 0.80 * A
  
  # Find all indices where y crosses 20% and 80% amplitude levels during rising phase
  rise_idx.20 <- which(diff(sign(y[1:idx.peak] - A.20)) != 0) # which(y[1:idx.peak] >= A.20)
  rise_idx.80 <- which(diff(sign(y[1:idx.peak] - A.80)) != 0) # which(y[1:idx.peak] >= A.80)
  
  # Find all indices where y crosses 20% and 80% amplitude levels during decaying phase
  decay_idx.80 <- which(diff(sign(y[idx.peak:length(y)] - A.80)) != 0) + idx.peak # which(y[idx.peak:length(y)] <= A.80) + idx.peak - 1
  decay_idx.20 <- which(diff(sign(y[idx.peak:length(y)] - A.20)) != 0) + idx.peak # which(y[idx.peak:length(y)] <= A.20) + idx.peak - 1
  
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
      t.decay.20 <- x2 + (A.20 - y2) / slope
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
    plot(x, y, col = 'indianred', xlab = 'x', type = 'l', bty = 'l', las = 1, main = '')
    lines(x, ysmooth, col = 'lightgray')
    abline(h = A.20, col = 'black', lty = 3)
    abline(h = A.80, col = 'black', lty = 3)
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

load_data2 <- function(wd, name) {
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
    data_list[[sheet]] <- openxlsx::read.xlsx(file_path, sheet = sheet)
  }
  
  return(data_list)
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

determine_tmax <- function(y, dt=0.1, stimulation_time=0, baseline=0, smooth=5, tmax=NULL, y_abline=0.1, height=5, width=5){ 
  if (is.null(tmax)){
    # Plot the data
    # plot(x, y, col='indianred', xlab='x', type='l', bty='l', las=1, main=paste('trace', ii))
    peak <- peak.fun(y=y, dt=dt, stimulation_time=stimulation_time, baseline=baseline, smooth=smooth)

    # Get one fit and assess
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
  }else{
    x_limit <- tmax
  }

  x_limit <- x_limit + stimulation_time - baseline
  return(x_limit)
}

sequential_fit <- function(response, n=30, dt=0.1, func=product2, method= c('LM', 'BF.LM', 'GN', 'port', 'robust', 'MLE'), 
    stimulation_time=0, baseline=0, tmax=NULL, y_abline=0.1, fast.decay.limit=NULL, fast.constraint=FALSE, fast.constraint.method=c('rise', 'peak'), first.delay.constraint=FALSE,
    latency.limit=NULL, lower=NULL, upper=NULL,  filter=FALSE, fc=1000, interval=c(0.1, 0.9), MLEsettings=list(iter=1e4, metropolis.scale=1.5, fit.attempts=100, RWm=FALSE), 
    MLE.method=c('L-BFGS-B', 'Nelder-Mead', 'BFGS','CG', 'SANN', 'Brent'), smooth=5, response_sign_method = c('smooth', 'regression', 'cumsum'), dp=3, 
    lwd=1.2, xlab='time (ms)', ylab='PSC (pA)', width=5, height=5, return.output=FALSE, show.output=TRUE, show.plot=TRUE, seed=42){
  
  y <- response
  x_limit <- determine_tmax(y=y, dt=dt, stimulation_time=stimulation_time, baseline=baseline, smooth=smooth, tmax=tmax, y_abline=y_abline, width=width, height=height)

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
  out <- nFIT(response=y2fit, n=n, dt=dt, func=func, method=method, stimulation_time=stimulation_time, baseline=baseline, 
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
      out <- nFIT(response=y2fit, n=n, dt=dt, func=func, method=method, stimulation_time=stimulation_time, baseline=baseline, 
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

nFIT_sequential <- function(response, n=30, dt=0.1, func=product2, method= c('LM', 'BF.LM', 'GN', 'port', 'robust', 'MLE'), 
    stimulation_time=0, baseline=0, fit.limits=NULL, fast.decay.limit=NULL, fast.constraint=FALSE, fast.constraint.method=c('rise', 'peak'), 
    first.delay.constraint=FALSE, latency.limit=NULL, lower=NULL, upper=NULL, filter=FALSE, fc=1000, interval=c(0.1, 0.9), 
    MLEsettings=list(iter=1e4, metropolis.scale=1.5, fit.attempts=100, RWm=FALSE), MLE.method=c('L-BFGS-B', 'Nelder-Mead', 'BFGS','CG', 'SANN', 'Brent'), 
    response_sign_method = c('smooth', 'regression', 'cumsum'), dp=3, lwd=1.2, xlab='time (ms)', ylab='PSC (pA)', width=5, height=5, return.output=FALSE, 
    show.output=TRUE, show.plot=TRUE, seed=42){

  fast.constraint.method <- match.arg(fast.constraint.method)

  y <- response
  if (identical(func, product1)){
    functions <- list(product1)
  }else if (identical(func, product2)){
    functions <- list(product2, product1)
  }else if (identical(func, product3)){
    functions <- list(product3, product2, product1)
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
    outputs[[ii]] <- sequential_fit(response=y2fit, n=n, dt=dt, func=functions[[ii]], method=method, stimulation_time=stimulation_time, baseline=baseline, tmax=tmax, 
                    fast.decay.limit=fast.decay.limit, fast.constraint=fast.constraint, fast.constraint.method=fast.constraint.method, first.delay.constraint=first.delay.constraint,
                    latency.limit=latency.limit, lower=lower, upper=upper, filter=filter, fc=fc, interval=interval, MLEsettings=MLEsettings, MLE.method=MLE.method, 
                    response_sign_method=response_sign_method, dp=dp, lwd=lwd, xlab=xlab, ylab=ylab, width=width, height=height, return.output=TRUE, show.output=output_logic, show.plot=output_logic, seed=seed)

    if (is.null(tmax)) { # Prompt user if they want to repeat
      cat('Do you want to repeat fit with new time base? (y/n): ')
      repeat_fit <- tolower(readLines(n = 1))
      while (repeat_fit == 'y') {
          dev.off()
          outputs[[ii]] <- sequential_fit(response=y2fit, n=n, dt=dt, func=functions[[ii]], method=method, stimulation_time=stimulation_time, baseline=baseline, tmax=tmax, 
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

  if (identical(func, product1)){
    df_output <- out.fun(params=fits[1:4], interval=interval, dp=dp, sign=1)
  } else if (identical(func, product2)){  
    df_output1 <- out.fun(params=fits[1:4], interval=interval, dp=dp, sign=1)
    df_output2 <- out.fun(params=fits[5:8], interval=interval, dp=dp, sign=1)
    df_output <- rbind('fast' = df_output1, 'slow' = df_output2)
  } else if (identical(func, product3)){  
    # Compute the output data frames
    df_output1 <- out.fun(params=fits[1:4],  interval=interval, dp=dp, sign=1)
    df_output2 <- out.fun(params=fits[5:8],  interval=interval, dp=dp, sign=1)
    df_output3 <- out.fun(params=fits[9:12], interval=interval, dp=dp, sign=1)
    df_output <- rbind('fast' = df_output1, 'medium' =  df_output2, 'slow' =  df_output3)
  }

  traces <- traces_fun(y=y, fits=fits, dt=dt, stimulation_time=stimulation_time, baseline=baseline, func=func, filter=filter, fc=fc, height=height, width=width)

  # lwd=1.2
  # xlab='time (ms)'
  # ylab='PSC (pA)'

  if (show.plot) fit_plot(traces=traces, func=func, xlab=xlab, ylab=ylab, lwd=lwd, filter=filter, width=width, height=height)
  if (show.output) print(df_output)
    
  idx1 <- baseline/dt
  idx2 <- max(t_limits)/dt

  k = length(fits)
  gof.se <- sqrt(sum((traces$y[idx1:idx2] - traces$yfit[idx1:idx2])^2) / (length(traces$y[idx1:idx2])-k))
  msc <- model.selection.criteria(coeffs=fits, x=traces$x[idx1:idx2]-traces$x[idx1], y=traces$y[idx1:idx2], func=func)

  if (return.output) {
    out <- list(output=df_output, fits=fits, gof=gof.se, AIC=msc[1], BIC=msc[2], traces=traces, fit.limits=t_limits - stimulation_time + baseline)
    return(out)
  }
}

analyse_PSC <- function(response, dt=0.1, n=30, stimulation_time=150, baseline=50, smooth=5, func=product2, method='LM', sequential.fit=FALSE, fit.limits=NULL, 
  MLEsettings=list(iter=1e3, metropolis.scale=1.5, fit.attempts=10, RWm=FALSE), filter=FALSE, fc=1000, interval=c(0.1, 0.9), lower=NULL, upper=NULL,  fast.decay.limit=NULL, 
  fast.constraint=FALSE, fast.constraint.method=c('rise', 'peak'), first.delay.constraint=FALSE, latency.limit=NULL, y_abline=0.1, dp=3, lwd=1.2, xlab='time (ms)', ylab='PSC (pA)', 
  return.output=TRUE, height=5, width=5, seed=42) {
  
  y <- response
  if (all(is.na(y[(which(!is.na(y))[length(which(!is.na(y)))] + 1):length(y)]))) {
    y <- y[!is.na(y)]
  }

  x <- seq(0, (length(y) - 1) * dt, by = dt)

  if (!sequential.fit){
    
    tmax <- fit.limits
    x_limit <- determine_tmax(y=y, dt=dt, stimulation_time=stimulation_time, baseline=baseline, smooth=smooth, tmax=tmax, y_abline=y_abline, width=width, height=height)   
    adjusted_response <- y[x < x_limit]
    
    # Execute nFIT
    out <- nFIT(response=adjusted_response, n=n, dt=dt, func=func, method=method, MLEsettings=MLEsettings, stimulation_time=stimulation_time, baseline=baseline, 
      filter=filter, fc=fc, interval=interval, fast.decay.limit=fast.decay.limit, fast.constraint=fast.constraint, fast.constraint.method=fast.constraint.method, 
      first.delay.constraint=first.delay.constraint, lower=lower, upper=upper, latency.limit=latency.limit, return.output=TRUE, show.plot=FALSE, dp=dp, height=height, width=width, seed=seed)

    out$traces <- traces_fun(y=y, fits=out$fits, dt=dt,  stimulation_time=stimulation_time, baseline=baseline, func=func, filter=filter, fc=fc, height=height, width=width)
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
        out <- nFIT(response=adjusted_response, n=n, dt=dt, func=func, method=method, MLEsettings=MLEsettings, stimulation_time=stimulation_time, baseline=baseline, 
          filter=filter, fc=fc, interval=interval, fast.decay.limit=fast.decay.limit, fast.constraint=TRUE, fast.constraint.method=fast.constraint.method, 
          first.delay.constraint=first.delay.constraint, lower=lower, upper=upper, latency.limit=latency.limit, return.output=TRUE, show.plot=FALSE, dp=dp, height=height, width=width, seed=seed)

        out$traces <- traces_fun(y=y, fits=out$fits, dt=dt,  stimulation_time=stimulation_time, baseline=baseline, func=func, filter=filter, fc=fc, height=height, width=width)
        fit_plot(traces=out$traces, func=func, xlab=xlab, ylab=ylab, lwd=lwd, filter=filter, width=width, height=height)
    
      }
    }
    
  }else{
    
    out <- nFIT_sequential(response=y, n=n, dt=dt, func=func, method=method, stimulation_time=stimulation_time, baseline=baseline, fit.limits=fit.limits, 
      fast.decay.limit=fast.decay.limit, fast.constraint=fast.constraint, fast.constraint.method=fast.constraint.method, first.delay.constraint=first.delay.constraint,
      latency.limit=latency.limit, lower=lower, upper=upper, filter=filter, fc=fc, interval=interval, MLEsettings=MLEsettings, MLE.method=MLE.method, dp=dp, 
      lwd=lwd, xlab=xlab, ylab=ylab, width=width, height=height, return.output=TRUE, show.output=TRUE, show.plot=TRUE, seed=seed)
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


boxplot_calculator <- function(data, type = 6) {
  unique_x <- unique(data$x)
  result <- data.frame(x = numeric(), Q1 = numeric(), Q3 = numeric(), Median = numeric(), Min = numeric(), Max = numeric(), MAD = numeric())
  
  for (i in 1:length(unique_x)) {
    current_x <- unique_x[i]
    d <- data$y[data$x == current_x]
    
    q1 <- quantile(d, probs = 0.25, type = type)
    q3 <- quantile(d, probs = 0.75, type = type)
    iqr <- q3 - q1  # Calculate IQR
    
    lower_bound <- q1 - 1.5 * iqr  # Lower bound for outliers
    upper_bound <- q3 + 1.5 * iqr  # Upper bound for outliers
    
    # Exclude outliers
    d_filtered <- d[d >= lower_bound & d <= upper_bound]
    
    median_val <- median(d)
    min_val <- min(d_filtered)
    max_val <- max(d_filtered)
    
    # Calculate MAD
    mad <- median(abs(d - median_val))
    
    result <- rbind(result, data.frame(x = current_x, Q1 = q1, Q3 = q3, Median = median_val, Min = min_val, Max = max_val, MAD = mad))
  }
  
  rownames(result) <- NULL  # Remove row names
  return(result)
}


WBplot <- function(data, wid = 0.2, cap = 0.05, xlab = '', ylab = 'PSP amplitude (mV)', 
                   xrange = c(0.75, 2.25), yrange = c(0, 400), main = '', tick_length = 0.02, 
                   x_tick_interval = NULL, y_tick_interval = 100, lwd = 0.8, type = 6) {
  
  boxplot_values <- boxplot_calculator(data, type)
  
  if (is.null(x_tick_interval)){
    x_ticks <- unique(data$x)
  }else{
    x_ticks <- seq(xrange[1], xrange[2], by = x_tick_interval)
  }
  xrange <- xrange + c(-wid, wid)
  
  # Ensure background is off and plot area is clear
  par(bg = NA)
  plot(1, type = 'n', ylim = yrange, xlim = xrange, xlab = xlab, ylab = ylab, 
       main = main, xaxt = 'n', yaxt = 'n', bty = 'n', lwd = lwd)
  
  for (i in 1:nrow(boxplot_values)) {
    current_x <- boxplot_values$x[i]
    
    rect(current_x - wid, boxplot_values$Q1[i], current_x + wid, boxplot_values$Q3[i], col = 'white', lwd = lwd)
    segments(current_x, boxplot_values$Q1[i], current_x, boxplot_values$Min[i], lwd = lwd)
    segments(current_x, boxplot_values$Q3[i], current_x, boxplot_values$Max[i], lwd = lwd)
    segments(current_x - cap, boxplot_values$Min[i], current_x + cap, boxplot_values$Min[i], lwd = lwd)
    segments(current_x - cap, boxplot_values$Max[i], current_x + cap, boxplot_values$Max[i], lwd = lwd)
    segments(current_x - wid * 1.1, boxplot_values$Median[i], current_x + wid * 1.1, boxplot_values$Median[i], col = 'black', lwd = 3 * lwd)
  }
  
  # Set the x-axis ticks
  axis(1, at = x_ticks, labels = x_ticks, tcl = -tick_length)
  
  # Set the y-axis ticks
  y_ticks <- seq(yrange[1], yrange[2], by = y_tick_interval)
  axis(2, at = y_ticks, tcl = -tick_length, las = 1)
}


BoxPlot <- function(data, wid=0.2, cap=0.05, xlab='', ylab='PSP amplitude (mV)', main='', xrange=c(0.75,2.25), yrange=c(0, -400), 
  tick_length=0.02, x_tick_interval = NULL, y_tick_interval=100, lwd=0.8, type=6, amount=0.05, p.cex=0.5, filename='boxplot.svg', 
  height=2.5, width=4, save=False){
  if (save) {
    # Open SVG device
    svg(file=filename, width=width, height=height, bg='transparent')
  } else {
    dev.new(width=width, height=height, noRStudioGD=TRUE)
  }

  WBplot(data=data, wid=wid, cap=cap, xlab=xlab, ylab=ylab, main=main, xrange=xrange, yrange=yrange, 
    tick_length=tick_length, x_tick_interval=x_tick_interval, y_tick_interval=y_tick_interval, lwd=lwd, type=type)
  set.seed(42)
  data$x_jitter <- jitter(data$x, amount=amount)
  points(data$x_jitter, data$y, pch=19, bg='transparent', col='darkgray', lwd=lwd/2, cex=p.cex)

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


scatter_plot <- function(scatter, xlim=c(0, 400), ylim=c(0, 400), x_tick_interval=100, y_tick_interval=100, height=4, width=4, main='',
                         colors=c("black", "indianred"), open_symbols=FALSE, p.cex=0.5, filename='scatter.svg', save=FALSE) {
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
  axis(1, at=x_ticks, labels=x_ticks, tcl=tick_length)

  # Customize y-axis with horizontal labels
  axis(2, at=y_ticks, labels=y_ticks, tcl=tick_length, las=1)

  if (save) {
    dev.off()
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

single_fit_egs <- function(single_example, xlim=NULL, ylim=NULL, lwd=1, show_text=FALSE, normalise=FALSE, height=4, width=2.5, xbar=100, ybar=50, save=FALSE, log_y=FALSE, colors=c('#4C77BB', '#CA92C1'), filename='plot.svg') {
  
  if (save) {
    # Open SVG device
    svg(file=filename, width=width, height=height, bg='transparent')
  } else {
    dev.new(width=width, height=height, noRStudioGD=TRUE)
  }
  
  x <- single_example$x
  y <- single_example$y
  
  fit1 <- single_example$y_fit1
  fit2 <- single_example$y_fit2
  
  if (is.null(xlim)) xlim <- c(min(x), max(x))
  if (is.null(ylim)) {
    if (log_y) {
      y <- -y
      fit1 <- -fit1
      fit2 <- -fit2
      ylim <- c(log(1), log(max(y[y > 0], na.rm=TRUE)))
    } else {
      ylim <- c(-max(y, na.rm=TRUE), 0)
    }
  }

  if (log_y) {
    y <- ifelse(y > 0, log(pmax(y, .Machine$double.eps)), NA) # ifelse(y > 0, log(y), NA)
    fit1 <- ifelse(fit1 > 0, log(fit1), NA)
    fit2 <- ifelse(fit2 > 0, log(fit2), NA)
  } 

  idx1 <- which.min(abs(x - xlim[1]))
  idx2 <- which.min(abs(x - xlim[2]))

  plot(x[idx1:idx2], y[idx1:idx2], type='l', col='#808285', xlim=xlim, ylim=ylim, bty='n', lwd=lwd, lty=1, axes=FALSE, frame=FALSE, xlab='', ylab='')

  fits <- cbind(fit1, fit2)
  
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
    mtext('PSC amplitude (pA)', side=2, line=2.5)
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



traces_smoothfits <- function(y, fits, dt=0.1,  stimulation_time=150, baseline=50, func=product1, filter=FALSE, fc=1000, upsample.fit = c(upsample=TRUE, factor=100)){
  
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


  if (identical(func, product1)){
    fits[4] <- fits[4] + baseline
  } else if (identical(func, product2)){
    fits[4] <- fits[4] + baseline; fits[8] <- fits[8] + baseline
  } else if (identical(func, product3)){
    fits[4] <- fits[4] + baseline; fits[8] <- fits[8] + baseline; fits[12] <- fits[12] + baseline
  }    
  traces2$yfit <- func(fits,traces2$x+dx)  
  if (identical(func, product2)){
    traces2$yfit1 <- product1(fits[1:4],traces2$x+dx) 
    traces2$yfit2 <- product1(fits[5:8],traces2$x+dx) 
  } 
  if (identical(func, product3)){
    traces2$yfit1 <- product1(fits[1:4],traces2$x+dx) 
    traces2$yfit2 <- product1(fits[5:8],traces2$x+dx) 
    traces2$yfit3 <- product1(fits[9:12],traces2$x+dx) 
  } 

  return(list(original=traces1, fits=traces2))
}



fit_plot2 <- function(traces, func=product1, xlab='time (ms)', ylab='PSC amplitude (pA)', xlim=NULL, ylim=NULL, bl=NULL, lwd=1.2, filter=FALSE, width=4, height=4, filename='trace.svg', save=FALSE) {
  
  if (save) {
    # Open SVG device
    svg(file=filename, width=width, height=height, bg='transparent')
  } else {
    dev.new(width=width, height=height, noRStudioGD=TRUE)
  }

  plot(traces$original$x, traces$original$y, col='gray', xlab=xlab, ylab=ylab, xlim=xlim, ylim=ylim, type='l', bty='l', las=1, lwd=lwd, main='')
  
  if (filter) {
    lines(traces$original$x, traces$original$yfilter, col='black', type='l', lwd=lwd)
  }
  
  lines(traces$fit$xfit, traces$fit$yfit, col='#CD5C5C', lty=3, lwd=2 * lwd)
  
  if (identical(func, product2)) {
    lines(traces$fit$x, traces$fit$yfit1, col='#4C78BC', lty=3, lwd=2 * lwd)
    lines(traces$fit$x, traces$fit$yfit2, col='#CA92C1', lty=3, lwd=2 * lwd)
  }

  if (identical(func, product3)) {
    lines(traces$fit$x, traces$fit$yfit1, col='#F28E2B', lty=3, lwd=2 * lwd)
    lines(traces$fit$x, traces$fit$yfit2, col='#4C78BC', lty=3, lwd=2 * lwd)
    lines(traces$fit$x, traces$fit$yfit3, col='#CA92C1', lty=3, lwd=2 * lwd)
  }

  if (!is.null(bl)) abline(v=bl, col='black', lwd=lwd, lty=3)

  if (save) {
    dev.off()
  }
}



smooth.plots <- function(y, fits, dt=0.1,  stimulation_time=150, baseline=50, func=product1, filter=FALSE, fc=1000, upsample.fit = c(upsample=FALSE, factor=100),
  xlab='time (ms)', ylab='', xlim=NULL, ylim=NULL, lwd=1.2, width=5, height=5, filename='trace.svg', save=FALSE){

  traces <- traces_smoothfits(y=y, fits=fits, dt=dt,  stimulation_time=stimulation_time, baseline=baseline, func=func, filter=filter, fc=fc, upsample.fit = upsample.fit)

  fit_plot2(traces=traces, func=func, xlab=xlab, ylab=ylab, xlim=xlim, ylim=ylim, lwd=lwd, filter=filter, width=width, height=height, filename=filename, save=save) 
  
}


