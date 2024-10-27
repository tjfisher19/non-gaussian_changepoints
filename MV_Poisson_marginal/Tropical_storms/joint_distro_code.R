
library(mvtnorm)
library(mnormt)

###############################
## Our multivariate Normal
##   with the Poisson stuff embedded
##
## b is the X terms from the inverse transformation
## a is the X-1 terms (lower terms)
##
## R is the underlying correlation
##
## This function returns the pdf for that vector and those values
##
distro_function <- function(a, b, R) {
  tryCatch({sadmvn(lower=a, upper=b, mean=rep(0, length(a)),
         varcov=R) },
         error=function(cond) {
           0
         } )
}


###########################################
## This internal function is what will be optimized...
##   it accepts a single vector with all the parameters
##   and calculates the log likelihood
## Parameters are transformed to help in optimization
internal_fun <- function(pars, d, X_mat, chpts) {
  
  m <- length(chpts)
  regimes <- as.factor(findInterval(1:nrow(X_mat), chpts) )
  
  
  A <- as.array(ltMatrices(pars[((m+1)*d+1):length(pars)], diag=FALSE ))[,,1]
  A <- A - diag(diag(A))
  A <- A + diag(apply(A, 1, function(x) sqrt(1 - sum(x*x))))
  R <- tcrossprod(A,A)
  
  lambda_mat <- NULL
  for(i in 0:m) {
    lambda_mat <- rbind(
      lambda_mat,
      matrix(rep(exp(pars[(i*d+1):(i*d+d)]),
                 each=sum(regimes==i)), ncol=d)
    )
  }
  #print(pars)
  b_series <- qnorm(ppois(X_mat, lambda=lambda_mat))
  a_series <- qnorm(ppois(X_mat-1, lambda=lambda_mat))

  sum(log(sapply(1:nrow(X_mat), function(i) { distro_function(a=a_series[i,], b=b_series[i,], R=R)}) ))
  
}


##############################################
## Here is the main function for MV Poisson
##   with no serial (temporal) correlation
##

mle_mv_pois_chpt <- function(X_mat, chpts=NULL, factr=1e12, approx=TRUE) {

  d <- dim(X_mat)[2]
  m <- length(chpts)
  X_mat <- as.matrix(X_mat)
  
  lambda_init <- c()
  regimes <- as.factor(findInterval(1:nrow(X_mat), chpts) )
  GROUPINGS <- split(1:nrow(X_mat), regimes)
  Regime_lengths <- sapply(GROUPINGS, length)
  if(m > 0) {
    for(i in 0:m) {
      lambda_init <- c(lambda_init, colMeans(X_mat[regimes==i, ]) )
    }
  } else {
    lambda_init <- colMeans(X_mat)
  }
  lambda_mat <- NULL
  for(i in 0:m) {
    lambda_mat <- rbind(
      lambda_mat,
      matrix(rep(lambda_init[(i*d+1):(i*d+d)],
                 each=sum(regimes==i)), ncol=d)
    )
  }
  
  b_series <- qnorm(ppois(as.matrix(X_mat), lambda=lambda_mat))
  R <- cor(b_series)
  a_series <- qnorm(ppois(as.matrix(X_mat)-1, lambda=lambda_mat))
  
  ## Ir approx==TRUE - Just compute the likelihood with the MME-type estimates
  if(approx==TRUE) {
    loglik <- sum(log(sapply(1:nrow(X_mat), function(i) { distro_function(a=a_series[i,], b=b_series[i,], R=R)}) ))
    MDL <- ifelse(m > 0,
                  -2*loglik + d*sum(log(Regime_lengths))/2 + log(m) + sum(log(chpts)) + d*(d-1)/2*log(nrow(X_mat))/2,
                  -2*loglik + d*sum(log(nrow(X_mat)))/2 + d*(d-1)/2*log(nrow(X_mat))/2 )
    BIC <- -2*loglik + log(nrow(X_mat))*(d*(m+1) + m + d*(d-1)/2)
    out <- list(loglik = loglik,
                MDL = MDL,
                BIC = BIC,
                chpts = chpts,
                Lambda = lambda_init,
                R=R)
  } else {  ## Here optimize of the initial values
    init_pars <- c()
    for(i in 0:m) {
      init_pars <- c(init_pars,
                     log(colMeans(X_mat[findInterval(1:nrow(X_mat), chpts)==i, ])) )
    }
    init_pars <- c(init_pars,
                   t(chol(R))[lower.tri(t(chol(R)), diag=FALSE)] )
    
    opt_out <- optim(par=init_pars, fn=internal_fun, d=d, X_mat=X_mat, chpts=chpts,
                     method="BFGS",
                     control=list(trace=0,
                                  parscale=c(rep(1.0,d*(m+1)), rep(0.1,length(init_pars)-(m+1)*d) ),
                                  factr=factr,
                                  fnscale=-1))
    A <- as.array(ltMatrices(opt_out$par[((m+1)*d+1):length(init_pars)], diag=FALSE ))[,,1]
    A <- A - diag(diag(A))
    A <- A + diag(apply(A, 1, function(x) sqrt(1 - sum(x*x))))
    R <- tcrossprod(A,A)
    
    loglik <- opt_out$value
    MDL <- ifelse(m > 0,
                  -2*loglik + d*sum(log(Regime_lengths))/2 + log(m) + sum(log(chpts)) + d*(d-1)/2*log(nrow(X_mat))/2,
                  -2*loglik + d*sum(log(nrow(X_mat)))/2 + d*(d-1)/2*log(nrow(X_mat))/2 )
    BIC <- -2*loglik + log(nrow(X_mat))*(d*(m+1) + m + d*(d-1)/2)
    out <- list(opt_out = opt_out, 
                loglik = opt_out$value,
                MDL = MDL,
                BIC = BIC,
                chpts = chpts,
                Lambda = exp(opt_out$par[1:(d*(m+1))]),
                R = R)
  }
  
  out
}


