##############################################################
## 
## This code provides variations of the change point
##   methods (MDL or BIC based) with the assumption
##   of iid Beta distributed data.
##
## We still include the three variations
##   - mean and precision
##   - shape parameters
##   - Beta regression for trend.

library(fitdistrplus)
library(betareg)

beta_fit_mean_chpt <- function(series, chpts=NULL, quick_eval=FALSE) {
  
  INT_SHIFTS <- as.factor(findInterval(1:length(series), chpts))
  GROUPINGS <- split(series, INT_SHIFTS)
  
  if(quick_eval) reltol <- 1.0e-4
  else reltol <- 1.0e-10
  
  mod_form <- "series ~ 1"
  
  if(length(chpts)> 0) {
    mod_form <- paste(mod_form, "+ INT_SHIFTS")
  }
  
  fit <- betareg(as.formula(mod_form), control=betareg.control(reltol=reltol))
  loglik <- fit$loglik

  mdl_changepoint_penalty <- ifelse(length(chpts)>0, 
                                    sum(log(sapply(GROUPINGS, length)))/2 +  ##intercept parameters for each segment
                                      log(length(series))/2 +                ## precision parameter
                                      log(length(chpts)) +                   ## Number of change points
                                      sum(log(chpts)),                       ## chpt locations
                                    log(length(series)) )                    ## No chpt penalty, 2 parameters
                                    
  bic_changepoint_penalty <- ifelse(length(chpts)>0,
                                    2*length(chpts)+2,                       ## (m+1) means + precision + m-chpts & determining m
                                    2)                                       ## Two parameters
  list(loglik=loglik, 
       MDL=-2*loglik + mdl_changepoint_penalty,
       BIC=-2*loglik + log(length(series))*bic_changepoint_penalty,
       Fit = fit)
}

beta_fit_shape_chpt <- function(series, chpts=NULL, method="mme") {
  
  INT_SHIFTS <- as.factor(findInterval(1:length(series), chpts))
  GROUPINGS <- split(series, INT_SHIFTS)
  
  log_lik_values <- sapply(GROUPINGS, function(x) {fitdist(x,"beta", method=method)$loglik})

  loglik <- sum(log_lik_values)
  
  mdl_changepoint_penalty <- ifelse(length(chpts)>0, 
                                    2*sum(log(sapply(GROUPINGS, length)))/2 +  ## two shape parameters for each segment
                                      log(length(chpts)) +                     ## Number of change points
                                      sum(log(chpts)),                         ## chpt locations
                                    sum(log(sapply(GROUPINGS, length))))       ## two shape parameters
  
  list(loglik=loglik, 
       MDL=-2*loglik + mdl_changepoint_penalty,
       BIC=-2*loglik + log(length(series))*(2*(length(chpts)+1)+length(chpts)+1))
}


ga_search_function_beta_chpt_mdl <- function(x, series, model_type="mean", min.regime=8, max.chpts=10, quick_eval=TRUE) {
  change.pts <- which(x==1) + 1
  if(sum(x) > max.chpts) {  # no more than max.chpts
    return(-1111111)
  } else if(min(diff(c(1,change.pts, (length(series)+1)))) < min.regime) {  # each regime must be at least min.regime obs
    return(-9999999)
  } else {
    if(model_type=="mean") fit_series <- beta_fit_mean_chpt(series=series, chpts=change.pts, quick_eval=quick_eval)
    else fit_series <- beta_fit_shape_chpt(series=series, chpts=change.pts, method="mme")
    -(fit_series$MDL )
  }
}

ga_search_function_beta_chpt_bic <- function(x, series, model_type="mean", min.regime=8, max.chpts=10, quick_eval=TRUE) {
  change.pts <- which(x==1) + 1
  if(sum(x) > max.chpts) {  # no more than max.chpts
    return(-1111111)
  } else if(min(diff(c(1,change.pts, (length(series)+1)))) < min.regime) {  # each regime must be at least min.regime obs
    return(-9999999)
  } else {
    if(model_type=="mean") fit_series <- beta_fit_mean_chpt(series=series, chpts=change.pts, quick_eval=quick_eval)
    else fit_series <- beta_fit_shape_chpt(series=series, chpts=change.pts, method="mme")
    -(fit_series$BIC )
  }
}

