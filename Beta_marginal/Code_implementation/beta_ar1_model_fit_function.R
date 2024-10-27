##########################################################
##
## The driver function that fits our method
##   assuming a Beta marginal distribution with an 
##   underlying AR(1) correlation structure. 
## Three versions of the model exist
##   - one with a mean shift, constant precision
##   - one with both shape parameters potentially changing
##   - one with a constrant trend but changing intercept
##     (based on Beta-regression)

source("beta_ar1_model_helper_functions.R")
source("beta_ar1_model_loglik_functions.R")

##########################################################
## This is the call function to compute the maximum 
##   likelihood estimator for the mean model and
##   AR(1) coefficient. 
## Least Squares and Yule-Walker estimates are found as initial
##   estimates before numerically finding the likelihood
##
## In the case of approx=TRUE, the likelihood function is 
##   calculated uses the LSE and YW estimates
##   This speeds up the GA search (substantially) to find
##   the appropriate segmentation

fit_beta_ar1_chpt <- function(
    series, chpts=NULL, type=c("mean", "trend", "shape"), approx=FALSE ) {
  
  series_name <- deparse(substitute(series))
  type <- match.arg(type)
  
  m <- length(chpts)
  
  ## Segment the series up based on the change point configuration
  INT_SHIFTS <- as.factor(findInterval(1:length(series), chpts))
  GROUPINGS <- split(series, INT_SHIFTS)
  Regime_lengths <- sapply(GROUPINGS, length)
  names(Regime_lengths) <- paste("Regime", 0:m)
  
  ###############################################
  ## Initial alpha & beta values - this depends 
  ##   on the model structure
  if(type=="shape") {
    log_lik_fun <- log_likelihood_beta_ar1_shape
    
    ## Get MMEs for alpha & beta within each segment
    params <- sapply(GROUPINGS, get_mme_shape_parameter)
    alpha_vec <- rep(params[1,], Regime_lengths)
    beta_vec <- rep(params[2,], Regime_lengths)
    ## On a log-scale for the optimization routine
    beta_parms <- log(c(params[1,], params[2,] ) )
    opt_parscale <- c(rep(0.01, (length(beta_parms))), 0.001)
  } else{
    ## For both the change in mean and trend models
    ##   we use generalized linear model ideas
    logit_link <- make.link("logit")
    series_logit <- logit_link$linkfun(series)
    
    if(type=="trend") {
      TIME_TREND <- 1:length(series)
      log_lik_fun <- log_likelihood_beta_ar1_trend
      form_so_far <- " ~ TIME_TREND"
      num_terms <- 2 + m
    } else if(type=="mean") {
      log_lik_fun <- log_likelihood_beta_ar1_mean
      form_so_far <- " ~ 1"
      num_terms <- 1 + m
    } 
    
    ## If change points exist, add it as a factor variable
    INT_NAMES <- NULL
    if(m > 0){
      form_so_far <- paste(form_so_far, "+ INT_SHIFTS")
      INT_NAMES <- paste("Regime", 1:m)
    }
    
    lm_fit <- lm(as.formula(paste("series_logit", form_so_far)))
    X_mat <- model.matrix(as.formula(form_so_far), data=as.data.frame(series_logit))
    
    yhat <- logit_link$linkinv(lm_fit$fitted.values)
    dlink <- 1/logit_link$mu.eta(lm_fit$fitted.values)
    res <- lm_fit$residuals
    N <- length(series)
    sigma2 <- sum(res^2)/((N - (num_terms) ) * (dlink)^2)
    phi_y <- yhat * (1 - yhat)/(N * sigma2) - 1/length(series)
    #precision <- sum(phi_y)     #suppressWarnings(log(sum(phi_y)))
    precision <- suppressWarnings(log(sum(phi_y)))
    ## Convert to standard alpha & beta shape parameters (vectors)
    alpha_vec <- yhat*exp(precision)
    beta_vec <- (1-yhat)*exp(precision)
    beta_parms <- c(lm_fit$coefficients, precision)
    opt_parscale <- c(rep(0.1, (length(beta_parms)-1)), 0.01, 0.001)
  }
  
  
  ## Approximate the underlying Std. Normal series
  ##   and get an estimate for the underlying correlation structure
  Z_series <- qnorm(pbeta(series, shape1=alpha_vec, shape2=beta_vec ) )
  phi <- ar(Z_series, order.max=1, aic=FALSE, demean=FALSE)$ar
  
  ## Now put everything into a vector and feed into the log-liklihood
  init_pars <- c(beta_parms, asin(phi) )
  
  if(approx) {
    opt_out <- list(par = init_pars,
                    value = as.numeric(log_lik_fun(init_pars, series=series, chpts=chpts)) )
  } else {
  opt_out <- optim(init_pars, fn=log_lik_fun, 
                   series=series, chpts=chpts,
                   method="BFGS",
                   control=list(parscale=opt_parscale,
                                reltol=1.0e-07))
  }
  
  ## Now convert the mean-parameters (on a logit scale) to proportions)
  
  if(type=="trend") {
    Coeffs <- opt_out$par[1:(m+2)]
    names(Coeffs) <- c("Intercept", "Trend", INT_NAMES)
    Parameters <- list(Coefficients = Coeffs,
                       Precision = exp(opt_out$par[m+3]),
                       AR_coef = sin(opt_out$par[m+4]) )
    Fitted <- logit_link$linkinv(X_mat%*%opt_out$par[1:(m+2)])
  } else if(type=="mean") {
    Means <- as.vector(logit_link$linkinv(unique(X_mat)%*%opt_out$par[1:(m+1)]))
    names(Means) <- paste("Regime", 0:m)
    Parameters <- list(Means = Means,
                       Precision = exp(opt_out$par[m+2]),
                       AR_coef = sin(opt_out$par[m+3]) ) 
    Fitted <- logit_link$linkinv(X_mat%*%opt_out$par[1:(m+1)])
  } else {
    Alpha <- exp(opt_out$par[1:(m+1)])
    Beta <- exp(opt_out$par[(m+2):(2*m+2)])
    names(Alpha) <- names(Beta) <- paste("Regime", 0:m)
    alpha_vec <- rep(Alpha, Regime_lengths)
    beta_vec <- rep(Beta, Regime_lengths)
    Fitted <- as.vector(alpha_vec/(alpha_vec + beta_vec))
    Parameters <- list(Alpha=Alpha, Beta=Beta,
                       AR_coef = sin(opt_out$par[2*(m+1)+1]) )
  }
  
  out <- list(series = series,
              Series_name = series_name,
              Model_type = type,
              chpts = chpts,
              num_chpts = m,
              Regime_lengths = Regime_lengths,
              fitted.values = Fitted,
              Regimes = INT_SHIFTS,
              Parameters = Parameters,
              loglik = -opt_out$value
  )
  
  penalties <- calc_penalty_segment_beta_ar1(out)
  out$MDL <- -2*out$loglik + penalties$MDL_penalty
  out$AIC <- -2*out$loglik + penalties$AIC_penalty
  out$BIC <- -2*out$loglik + penalties$BIC_penalty

  out$Full_optim_output <- opt_out
  class(out) <- "BetaAR.Chpt"
  
  return(out)
}


