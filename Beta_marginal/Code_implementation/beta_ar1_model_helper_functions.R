
######################################
## Get Method of Moment Estimators
##    for a sequence (x) of Beta
##    random variables. Returns the
##    two shape parameters

get_mme_shape_parameter <- function(x) {
  x_bar <- mean(x)
  v_bar <- var(x)
  c(x_bar*(x_bar*(1-x_bar)/v_bar - 1), (1-x_bar)*(x_bar*(1-x_bar)/v_bar - 1) )
}


#############################################################
## For a given fitted Beta-AR(1) Change point model
##    we calculate the penalty terms for AIC, BIC and MDL

calc_penalty_segment_beta_ar1 <- function(fitted_model) {
  
  if(fitted_model$Model_type=="mean") {
    
    mdl_changepoint_penalty <- ifelse(fitted_model$num_chpts>0, 
                                      sum(log(fitted_model$Regime_lengths))/2 + ## mean for each segment estimation
                                        log(length(fitted_model$series))/2 +    ## Precision parameter estimation
                                        log(fitted_model$num_chpts) +           ## Number of change points
                                        sum(log(fitted_model$chpts)),           ## chpt locations
                                      log(length(fitted_model$series) ) )       ## No penalty on chpts, two parameters
    
    aic_changepoint_terms <- ifelse(fitted_model$num_chpts>0,
                                    (fitted_model$num_chpts+1) +                ## (m+1) means + precision
                                      (fitted_model$num_chpts+1),               ## (m change points + determining m)
                                    2)                                          ## no penalty on chpts, two parameters
  } else if(fitted_model$Model_type=="trend") {
    mdl_changepoint_penalty <- ifelse(fitted_model$num_chpts>0, 
                                      sum(log(fitted_model$Regime_lengths))/2 + ## intercept term in each segment
                                        2*log(length(fitted_model$series))/2 +  ## trend and precision estimate
                                        log(fitted_model$num_chpts) +           ## Number of change points
                                        sum(log(fitted_model$chpts)),           ## chpt locations
                                      log(length(fitted_model$series) )*3/2 )   ## intercept + slope + precision
    
    aic_changepoint_terms <- ifelse(fitted_model$num_chpts>0,
                                    (fitted_model$num_chpts+1) + 1              ## (m+1) intercepts + trend + precision
                                    (fitted_model$num_chpts+1),                 ## (m change points + determining m)
                                    3)                                          ## intercept + slope + precision
    
  } else if(fitted_model$Model_type=="shape") {
    mdl_changepoint_penalty <- ifelse(fitted_model$num_chpts>0, 
                                      2*sum(log(fitted_model$Regime_lengths))/2 +  ## 2 shape parameters for each segment
                                        log(fitted_model$num_chpts) +              ## Number of change points
                                        sum(log(fitted_model$chpts)),              ## chpt locations
                                      log(length(fitted_model$series) ) )          ## No chpt penalty, 2 shape terms
    
    aic_changepoint_terms <- ifelse(fitted_model$num_chpts>0,
                                    2*(fitted_model$num_chpts+1) +              ## 2*(m+1) shape parameters
                                      (fitted_model$num_chpts+1),               ## (m change points + determining m)
                                    2)                                          ## 2 shape parameters
  }
  
  ### This part does not really matter for our purposes since we are not varying the AR order
  ##   but included for completeness
  
  mdl_AR_penalty <- 1*log(length(fitted_model$series))/2 +    ## Cost of phi_i for i=1
    log(1) +                                                    ## Cost of estimating p, AR(p)
    log(length(fitted_model$series))/2                          ## Cost of estimating ACVF_0
  
  
  aic_AR_terms <- 1   # Just one parameter is estimate for AIC & BIC
  
  list(MDL_penalty = mdl_changepoint_penalty + mdl_AR_penalty,
       AIC_penalty = 2*(aic_changepoint_terms + aic_AR_terms),      # 2*number of terms
       BIC_penalty = log(length(fitted_model$series))*(aic_changepoint_terms + aic_AR_terms) )
}

#############################################################
## For a given fitted Beta-AR(1) Change point model
##    a convenience function print() that outputs 
##  Most of the relavent information, including:
##   the series name, the model type (trend, mean or shape),
##   the change point locations, the lengths of the regimes,
##   the fitted model parameters and measure of fit which
##   include the log-likelihood, AIC, BIC and MDL penalized
##   measures of fit

print.BetaAR.Chpt <- function(out) {
  
  if(out$Model_type=="shape") {
    cat("Change Points in the Shape parameters\n\n")
  } else if(out$Model_type=="trend") {
    cat("Change Points in the intercept for a Beta-regression with trend\n\n")
  } else if(out$Model_type=="mean") {
    cat("Change Points in the Mean with a constant Precision parameter\n\n")
  } else {
    stop("Error: Invalid model type")
  }
  cat(paste0("Series: ", out$Series_name, "\n\n") )
  cat("Change Point Locations\n")
  print(out$chpts)
  cat("\nRegime Lengths\n")
  print(out$Regime_lengths)
  cat("\nFitted Parameters\n")
  print(out$Parameters)
  tmp <- c(out$loglik, out$AIC, out$BIC, out$MDL)
  names(tmp)<- c("Log Likelihood", "AIC", "BIC", "MDL")
  cat("\n")
  cat("Measures of fit\n")
  print(tmp)
}


