
##################################################
##
##  Functions that run our Beta simulations
##    when there is not change point present


library(GA)

source("mv_pois_ar1_fit_functions.R")
source("mv_pois_ar1_generation_testing.R")

##########################################
## CONSTANTS
##########################################

## Genetic algorithm
cores <- 27
min.regime <- 5
PopSize <- 100
maxit <- 250
max.run <- 25

rng_seed <- 1    # For data generation


### The GA search function is here, it needs N to be declared
source("ga_search_functions.R")

###############################################

load("data/nochange_pos_pos.RData")

cor_config <- "pos_pos"
## This counter is included so we can see progress
cnt <- 0
full_fit_fake_series <- lapply(X_series, fit_fake_series)


save(full_fit_fake_series, X_series, cor_config,
     file="results/mvpois_no_change_sims_pos_pos.RData")



