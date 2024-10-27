
##################################################
##
##  Functions that run our Beta simulations
##    when there is not change point present


library(GA)

source("beta_ar1_model_ga_search_functions.R")
source("beta_iid_shape_functions.R")

##########################################
## CONSTANTS
##########################################

## Genetic algorithm
cores <- 11
min.regime <- 8
PopSize <- 300
maxit <- 500
max.run <- 50

rng_seed <- 1    # For data generation

## Series - sample size, AR-coef, shape1 & shape2
## These values are based on the MLE for the Baseball data
##   as if no change points are present
N <- 100
rep_series <- 100
phi <- 0.30       ## Reset for each simulation
Alpha <- 8.68
Beta <- 365
## Note:
Alpha/(Alpha+Beta)   ## Mean
Alpha+Beta           ## Precision

### The GA search function is here, it needs N to be declared
source("simulate_search_functions.R")

###############################################
## No correlation
phi <- 0
set.seed(rng_seed)
innov <- replicate(rep_series, rnorm(N), simplify=FALSE)

fake_series <- lapply(innov, function(z) { qbeta(pnorm(z), shape1=Alpha, shape2=Beta)})


## This counter is included so we can see progress
cnt <- 0
full_fit_fake_series <- lapply(fake_series, fit_fake_series)

parameters <- list(Alpha=Alpha, Beta=Beta, phi=phi)
save(full_fit_fake_series, fake_series, parameters, 
     file="results/beta_no_change_sims_noCorrelation.RData")

###############################################
## Moderate correlation (close to baseball)

phi <- 0.30
set.seed(rng_seed)
innov <- replicate(rep_series, arima.sim(model=list(ar=phi), n=N, sd=sqrt(1-phi^2)), simplify = FALSE )

fake_series <- lapply(innov, function(z) { qbeta(pnorm(z), shape1=Alpha, shape2=Beta)})

cnt <- 0
full_fit_fake_series <- lapply(fake_series, fit_fake_series)

parameters <- list(Alpha=Alpha, Beta=Beta, phi=phi)
save(full_fit_fake_series, fake_series, parameters, 
     file="results/beta_no_change_sims_baseballParams.RData")



###############################################
## Strong Positive correlation

phi <- 0.60
set.seed(rng_seed)
innov <- replicate(rep_series, arima.sim(model=list(ar=phi), n=N, sd=sqrt(1-phi^2)), simplify = FALSE )

fake_series <- lapply(innov, function(z) { qbeta(pnorm(z), shape1=Alpha, shape2=Beta)})

cnt <- 0
full_fit_fake_series <- lapply(fake_series, fit_fake_series)

parameters <- list(Alpha=Alpha, Beta=Beta, phi=phi)
save(full_fit_fake_series, fake_series, parameters, 
     file="results/beta_no_change_sims_strongCorrelation.RData")


