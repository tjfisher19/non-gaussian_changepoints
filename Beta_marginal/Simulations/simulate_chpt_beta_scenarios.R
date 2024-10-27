


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
## Single Change Point - Mean
phi <- 0.30
set.seed(rng_seed)
innov <- replicate(rep_series, arima.sim(model=list(ar=phi), n=N, sd=sqrt(1-phi^2)), simplify = FALSE )

mu <- c(rep(0.023, 40), rep(0.030, 60))
precision <- 500
Alpha <- mu*precision
Beta <- (1-mu)*precision

fake_series <- lapply(innov, function(z) { qbeta(pnorm(z), shape1=Alpha, shape2=Beta)})
# 
# 
## This counter is included so we can see progress
cnt <- 0
full_fit_fake_series <- lapply(fake_series, fit_fake_series)

parameters <- list(Alpha=Alpha, Beta=Beta, phi=phi)
save(full_fit_fake_series, fake_series, parameters,
     file="results/beta_change_sims_singleChangePoint.RData")

###############################################
## Change in Precision

# phi <- 0.30
# set.seed(rng_seed)
# innov <- replicate(rep_series, arima.sim(model=list(ar=phi), n=N, sd=sqrt(1-phi^2)), simplify = FALSE )
# 
# mu <- 0.023
# precision <- c(rep(400, 40), rep(700, 60) )
# Alpha <- mu*precision
# Beta <- (1-mu)*precision
# 
# fake_series <- lapply(innov, function(z) { qbeta(pnorm(z), shape1=Alpha, shape2=Beta)})
# 
# cnt <- 0
# full_fit_fake_series <- lapply(fake_series, fit_fake_series)
# 
# parameters <- list(Alpha=Alpha, Beta=Beta, phi=phi)
# save(full_fit_fake_series, fake_series, parameters, 
#      file="beta_change_sims_precisionChangePoint.RData")

###############################################
## Change in mean and precision

phi <- 0.30
set.seed(rng_seed)
innov <- replicate(rep_series, arima.sim(model=list(ar=phi), n=N, sd=sqrt(1-phi^2)), simplify = FALSE )

mu <- c(rep(0.023, 40), rep(0.030, 60))
precision <- c(rep(300, 40), rep(700, 60) )
Alpha <- mu*precision
Beta <- (1-mu)*precision

fake_series <- lapply(innov, function(z) { qbeta(pnorm(z), shape1=Alpha, shape2=Beta)})

cnt <- 0
full_fit_fake_series <- lapply(fake_series, fit_fake_series)

parameters <- list(Alpha=Alpha, Beta=Beta, phi=phi)
save(full_fit_fake_series, fake_series, parameters,
     file="results/beta_change_sims_meanPrecisionChangePoint.RData")

###############################################
## Two Changes Mean

phi <- 0.30
set.seed(rng_seed)
innov <- replicate(rep_series, arima.sim(model=list(ar=phi), n=N, sd=sqrt(1-phi^2)), simplify = FALSE )

mu <- c(rep(0.023, 40), rep(0.030, 30), rep(0.022, 30))
precision <- 600
Alpha <- mu*precision
Beta <- (1-mu)*precision

fake_series <- lapply(innov, function(z) { qbeta(pnorm(z), shape1=Alpha, shape2=Beta)})

cnt <- 0
full_fit_fake_series <- lapply(fake_series, fit_fake_series)

parameters <- list(Alpha=Alpha, Beta=Beta, phi=phi)
save(full_fit_fake_series, fake_series, parameters,
     file="results/beta_change_sims_twoChangePoints.RData")

###############################################
## Two Changes Mean and Precision

phi <- 0.30
set.seed(rng_seed)
innov <- replicate(rep_series, arima.sim(model=list(ar=phi), n=N, sd=sqrt(1-phi^2)), simplify = FALSE )

mu <- c(rep(0.023, 40), rep(0.030, 30), rep(0.022, 30))
precision <- c(rep(600, 40), rep(300, 30), rep(600, 30))
Alpha <- mu*precision
Beta <- (1-mu)*precision

fake_series <- lapply(innov, function(z) { qbeta(pnorm(z), shape1=Alpha, shape2=Beta)})

cnt <- 0
full_fit_fake_series <- lapply(fake_series, fit_fake_series)

parameters <- list(Alpha=Alpha, Beta=Beta, phi=phi)
save(full_fit_fake_series, fake_series, parameters,
    file="results/beta_change_sims_twoMeanPrecisionChangePoints.RData")

###############################################
## Four changes Mean Only

phi <- 0.30
set.seed(rng_seed)
innov <- replicate(rep_series, arima.sim(model=list(ar=phi), n=N, sd=sqrt(1-phi^2)), simplify = FALSE )
 
mu <- c(rep(0.021, 25), rep(0.027, 25), rep(0.021, 15), rep(0.029, 25), rep(0.034, 10) )
precision <- 1500
Alpha <- mu*precision
Beta <- (1-mu)*precision
 
fake_series <- lapply(innov, function(z) { qbeta(pnorm(z), shape1=Alpha, shape2=Beta)})
# 
# 
cnt <- 0
full_fit_fake_series <- lapply(fake_series, fit_fake_series)
# 
parameters <- list(Alpha=Alpha, Beta=Beta, phi=phi)
save(full_fit_fake_series, fake_series, parameters, 
     file="results/beta_change_sims_fourChangePoints.RData")

###############################################
## Four changes Mean and precision

phi <- 0.30
set.seed(rng_seed)
innov <- replicate(rep_series, arima.sim(model=list(ar=phi), n=N, sd=sqrt(1-phi^2)), simplify = FALSE )

mu <- c(rep(0.021, 25), rep(0.027, 25), rep(0.021, 15), rep(0.029, 25), rep(0.034, 10) )
precision <- c(rep(500, 25), rep(600, 25), rep(700, 15), rep(1000, 25), rep(1250, 10) )
Alpha <- mu*precision
Beta <- (1-mu)*precision

fake_series <- lapply(innov, function(z) { qbeta(pnorm(z), shape1=Alpha, shape2=Beta)})


cnt <- 0
full_fit_fake_series <- lapply(fake_series, fit_fake_series)

parameters <- list(Alpha=Alpha, Beta=Beta, phi=phi)
save(full_fit_fake_series, fake_series, parameters, 
     file="results/beta_change_sims_fourMeanPrecisionChangePoints.RData")



