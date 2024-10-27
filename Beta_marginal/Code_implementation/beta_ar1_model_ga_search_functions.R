#####################################################
##
## Functions used by the Genetic Algorithm
##  to search on.  It is a modification of our
##  fit_beta_ar1_chpt() function where we 
##  optimize over MDL, AIC or BIC
##

source("beta_ar1_model_fit_function.R")


ga_search_function_beta_ar1_chpt_mdl <- function(x, series, model_type, approx=FALSE, min.regime=8, max.chpts=10) {
  change.pts <- which(x==1) + 1
  if(sum(x) > max.chpts) {  # no more than max.chpts
    return(-1111111)
  } else if(min(diff(c(1,change.pts, (length(series)+1)))) < min.regime) {  # each regime must be at least min.regime obs
    return(-9999999)
  } else {
    fit_series <- fit_beta_ar1_chpt(series=series, chpts=change.pts, type=model_type, approx=approx)
    -(fit_series$MDL)
  }
}

ga_search_function_beta_ar1_chpt_aic <- function(x, series, model_type, approx=FALSE, min.regime=8, max.chpts=10) {
  change.pts <- which(x==1) + 1
  if(sum(x) > max.chpts) {  # no more than max.chpts
    return(-1111111)
  } else if(min(diff(c(1,change.pts, (length(series)+1)))) < min.regime) {  # each regime must be at least min.regime obs
    return(-9999999)
  } else {
    fit_series <- fit_beta_ar1_chpt(series=series, chpts=change.pts, type=model_type, approx=approx)
    -(fit_series$AIC )
  }
}


ga_search_function_beta_ar1_chpt_bic <- function(x, series, model_type, approx=FALSE, min.regime=8, max.chpts=10) {
  change.pts <- which(x==1) + 1
  if(sum(x) > max.chpts) {  # no more than max.chpts
    return(-1111111)
  } else if(min(diff(c(1,change.pts, (length(series)+1)))) < min.regime) {  # each regime must be at least min.regime obs
    return(-9999999)
  } else {
    fit_series <- fit_beta_ar1_chpt(series=series, chpts=change.pts, type=model_type, approx=approx)
    -(fit_series$BIC )
  }
}