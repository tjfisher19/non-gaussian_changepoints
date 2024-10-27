

library(tidyverse)
load("./data/tropicalStormCounts.RData")

source("joint_distro_code.R")

my.fun.mv.pois.mdl <- function(x, series, min.regime=0, max.chpts=10) {
  change.pts <- which(x==1) + 1
  if(sum(x) > max.chpts) {  # no more than max.chpts
    return(-1111111)
  } else if(min(diff(c(1,change.pts, (nrow(series)+1)))) < min.regime) {  # each regime must be at least min.regime obs
    return(-9999999)
  } else {
    fit_series <- mle_mv_pois_chpt(X_mat=series, chpts=change.pts, approx=TRUE)
    -(fit_series$MDL)
  }
}

my.fun.mv.pois.bic <- function(x, series, min.regime=0, max.chpts=10) {
  change.pts <- which(x==1) + 1
  if(sum(x) > max.chpts) {  # no more than max.chpts
    return(-1111111)
  } else if(min(diff(c(1,change.pts, (nrow(series)+1)))) < min.regime) {  # each regime must be at least min.regime obs
    return(-9999999)
  } else {
    fit_series <- mle_mv_pois_chpt(X_mat=series, chpts=change.pts, approx=TRUE)
    -(fit_series$BIC)
  }
}


library(parallel)
library(GA)

N <- dim(storm_counts_wide)[1]
PopSize <- 100
min.regime <- 5
maxit <- 500
max.run <- 50
set.seed(123)   ## For consistent result

suggest_matrix <- matrix(0,
                         ncol= N-1,
                         nrow= PopSize)
for(i in 2:(NCOL(suggest_matrix)-2*min.regime+3) ) {
  suggest_matrix[i, min.regime + i-2] <- 1
}
for(i in (NCOL(suggest_matrix)-min.regime+2):PopSize) {
  suggest_matrix[i, sample(min.regime:((NCOL(suggest_matrix)-min.regime+1)/2), size=1)] <- 1
  suggest_matrix[i, sample(((NCOL(suggest_matrix)-min.regime+1)/2+1):(NCOL(suggest_matrix)-min.regime+1), size=1)] <- 1
}


####################################
## Any changes to the joint series
####################################
joint.pois.ga.mdl <- ga(type="binary", 
                        fitness=my.fun.mv.pois.mdl,
                        nBits=N-1,
                        series = storm_counts_wide[,-1],
                        min.regime=min.regime,
                        max.chpts=floor(N/min.regime),
                        popSize=PopSize,
                        suggestions = suggest_matrix,
                        parallel = 8,
                        maxiter=5000,
                        run=500)
joint.pois.ga.bic <- ga(type="binary", 
                        fitness=my.fun.mv.pois.bic,
                        nBits=N-1,
                        series = storm_counts_wide[,-1],
                        min.regime=min.regime,
                        max.chpts=floor(N/min.regime),
                        popSize=PopSize,
                        suggestions = suggest_matrix,
                        parallel = 8,
                        maxiter=5000,
                        run=500)


joint.pois.ga.mdl@solution
joint.pois.ga.bic@solution

which(joint.pois.ga.mdl@solution==1)+1980
which(joint.pois.ga.bic@solution==1)+1980

save(joint.pois.ga.bic, joint.pois.ga.mdl,
     file="tropicalStorm_gaFindings.RData")
