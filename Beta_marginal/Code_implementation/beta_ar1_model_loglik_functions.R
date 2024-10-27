#############################################################################
##
## These functions calculate the underlying log-likelihood function
##
##  There is the main function - loglik_beta_ar1
##    which calculcates the log likelihood based on our derivation
##
## and then three wrapper functions corresponding to the three implementations
##    one with a mean and precision
##    one with shape parameters
##    one with a beta-regression (trend)
###



loglik_beta_ar1 <- function(x, alpha_vec, beta_vec, phi) {
  N <- length(x)
  
  Z_series <- qnorm(pbeta(x, shape1=alpha_vec, shape2=beta_vec) )
  pred_errors <- Z_series[2:N] - phi*Z_series[1:(N-1)]
  
  -1/2*log(2*pi) - Z_series[1]^2/2 -                                     ## Condition on obs 1
    (N-1)/2*log(2*pi*(1-phi^2) ) - (1/2)*sum(pred_errors^2)/(1-phi^2) +  ## Gaussian AR(1) part
    N/2*log(2*pi) + (1/2)*sum(Z_series^2) + sum(dbeta(x, shape1=alpha_vec, shape2=beta_vec, log=TRUE))    ## Jacobian part
  
}


######################################################
## For a given time series (series) and changepoint
##   configuration (chpts) with parameter values (par)
##   compute the log likelihood
## Note: par is a vector of length 2(m+1)+1 where
##   there are m changepoints.  Thus m+1 regimes
##   and two Beta shape parameters in each; thus 2(m+1)
##   the additional term is the AR(1) coefficient.

log_likelihood_beta_ar1_shape <- function(par, series, chpts) {
  
  m <- length(chpts)
  N <- length(series)
  
  groupings <- split(series, as.factor(findInterval(1:length(series), chpts)))
  alpha_vec <- rep(exp(par[1:(m+1)]), sapply(groupings, length))
  beta_vec <- rep(exp(par[(m+2):(length(par)-1)]), sapply(groupings, length))
  phi <- sin(par[length(par)])
  -loglik_beta_ar1(x=series, alpha_vec=alpha_vec, beta_vec=beta_vec, phi=phi)
}

######################################################
## For a given time series (series) and changepoint
##   configuration (chpts) with parameter values (par)
##   compute the log likelihood for the model with a mean shift
## Note: par is a vector of length (m+1)+2 where
##   there are m changepoints (shifts in intercept).
##   Thus m+1 regimes (intercept & m shifts). 
##   We also have a precision parameter and 
##   the additional term is the AR(1) coefficient.

log_likelihood_beta_ar1_mean <- function(par, series, chpts) {
  
  m <- length(chpts)
  mod_formula <- "~ 1"
  INT_SHIFTS <- as.factor(findInterval(1:length(series), chpts))
  if(m > 0 ) {
    mod_formula <- paste(mod_formula, "+ INT_SHIFTS")
  }
  X_mat <- model.matrix(as.formula(mod_formula), data=as.data.frame(series))
  
  Beta_terms <- par[1:(m+1)]
  logit_link <- make.link("logit")
  mu_terms <- logit_link$linkinv(X_mat%*%(Beta_terms) )
  #mu_terms <- 1/(1 + exp(-X_Mat%*%log(Beta_terms) ) )
  alpha_vec <- mu_terms*exp(par[m+2])        ## par[m+2] is the precision term
  beta_vec <- (1 - mu_terms)*exp(par[m+2])
  phi <- sin(par[length(par)])
  
  -loglik_beta_ar1(x=series, alpha_vec=alpha_vec, beta_vec=beta_vec, phi=phi)
}


######################################################
## For a given time series (series) and changepoint
##   configuration (chpts) with parameter values (par)
##   compute the log likelihood for the model with a trend
## Note: par is a vector of length (m+1)+3 where
##   there are m changepoints (shifts in intercept).
##   Thus m+1 regimes (intercept & m shifts). 
##   We also have a trend component, precision parameter
##   the additional term is the AR(1) coefficient.

log_likelihood_beta_ar1_trend <- function(par, series, chpts) {
  
  m <- length(chpts)
  TIME_TREND <- 1:length(series)
  INT_SHIFTS <- as.factor(findInterval(1:length(series), chpts))
  mod_formula <- "~ TIME_TREND"
  
  if(m > 0 ) {
    mod_formula <- paste(mod_formula, "+ INT_SHIFTS")
  }
  X_mat <- model.matrix(as.formula(mod_formula))
  Beta_terms <- par[1:(m+2)]
  logit_link <- make.link("logit")
  mu_terms <- logit_link$linkinv(X_mat%*%(Beta_terms) )
  #mu_terms <- 1/(1 + exp(-X_mat%*%Beta_terms))
  
  alpha_vec <- mu_terms*exp(par[m+3])        ## par[m+3] is the precision term
  beta_vec <- (1 - mu_terms)*exp(par[m+3])
  phi <- sin(par[length(par)])
  
  -loglik_beta_ar1(x=series, alpha_vec=alpha_vec, beta_vec=beta_vec, phi=phi)
}
