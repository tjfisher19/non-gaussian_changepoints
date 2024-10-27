#####################################################
##
## Functions used by the Genetic Algorithm
##  to search on.  It is a modification of our
##  fit_beta_ar1_chpt() function where we 
##  optimize over MDL, AIC or BIC
##

source("mv_pois_iid_fit_function.R")


ga_search_function_mvpois_iid_chpt_mdl <- function(x, series, approx=FALSE, min.regime=5, max.chpts=5) {
  change.pts <- which(x==1) + 1
  if(sum(x) > max.chpts) {  # no more than max.chpts
    return(-1111111)
  } else if(min(diff(c(1,change.pts, (nrow(series)+1)))) < min.regime) {  # each regime must be at least min.regime obs
    return(-9999999)
  } else {
    fit_series <- fit_mvpois_iid(series=series, chpts=change.pts)
    -(fit_series$MDL)
  }
}


ga_search_function_mvpois_iid_chpt_bic <- function(x, series, approx=FALSE, min.regime=5, max.chpts=5) {
  change.pts <- which(x==1) + 1
  if(sum(x) > max.chpts) {  # no more than max.chpts
    return(-1111111)
  } else if(min(diff(c(1,change.pts, (nrow(series)+1)))) < min.regime) {  # each regime must be at least min.regime obs
    return(-9999999)
  } else {
    fit_series <- fit_mvpois_iid(series=series, chpts=change.pts)
    -(fit_series$BIC )
  }
}

