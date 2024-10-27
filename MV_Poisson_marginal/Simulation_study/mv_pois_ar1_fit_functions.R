library(mnormt)
library(mvtnorm)
library(expm)
library(matrixStats)


likelihood_mvpois_ar1 <- function(pars, series, chpts=NULL, N) { 
  
  ## First we un-transform all the variables and get all the bookkeeping
  d <- dim(series)[2]
  
  m <- length(chpts)
  lambda_vec <- exp(pars[1:((m+1)*d)])
  regimes <- findInterval(1:nrow(series), chpts)
  
  ##  Get the lambda matrix
  lambda_mat <- NULL
  for(i in 0:m) {
    lambda_mat <- rbind(
      lambda_mat,
      matrix(rep(lambda_vec[(i*d+1):(i*d+d)],
                 each=sum(regimes==i)), ncol=d)
    )
  }
  
  ## Get our Gamma1 matrix
  #Gamma1 <- matrix((atan(pars[((m+1)*d+1):((m+1)*d+d*d)])) / pi * 2, nrow=d, ncol=d)
  Gamma1 <- matrix(pars[((m+1)*d+1):((m+1)*d+d*d)], nrow=d, ncol=d)
  
  ## Get the Gamma0 matrix (Cholesky lower values passed in)
  A <- as.array(ltMatrices(pars[(d*d+(m+1)*d+1):length(pars)], diag=TRUE ))[,,1]
  Gamma0 <- tcrossprod(A,A)
  
  ## From Gamma0 and Gamma1 we can calculate everything
  
  ## ACF matrix at lag 0 - the Z-distribution matrix
  Sigma_Z <- diag((1/sqrt(diag(Gamma0))))%*%Gamma0%*%diag((1/sqrt(diag(Gamma0))))
  
  ## Now we can get the Phi & variance-covariance for the epsilon
  ##   as they are both functions of Gamma0 & Gamma1
  Phi <- Gamma1%*%solve(Gamma0)

  ## From Phi and Sigma_Z we can get the Sigma for Epsilon
  A_mat <- diag(rep(1,d*d) ) - kronecker(Phi, Phi)
  Sigma_eps <- matrix(as.vector(Sigma_Z)%*%A_mat, nrow=d )
  Sigma_eps <- Sigma_Z - Phi%*%Sigma_Z%*%t(Phi)

  ## Initialization of some matrices
  prev_Z <- matrix(rep(0, d*N), nrow=d )
  epsilon <- array(rep(0, d*N), dim=c(d, N ) )

  ## Get the bounds based on the observed poisson
  a <- qnorm(ppois(series-1, lambda=lambda_mat) )
  b <- qnorm(ppois(series, lambda=lambda_mat) )
 
  ## Time point t=1, conditioning step
  
  ## Initial Weight
  w_Z <- mnormt::sadmvn(lower=a[1,], upper=b[1,], 
                        mean=rep(0,2), varcov = Sigma_Z)

  ## Accumulator for Weight
  log_w <- w_Z <- rep(log(w_Z), N)
  
  ## Generate a set of initial Z terms from
  ##   truncated Normal, here we use the Sigma_Z
  prev_Z <- t(mnormt::rmtruncnorm(N, mean=rep(0,d), varcov=Sigma_Z, 
                                  lower=a[1,], upper=b[1,]) )

  ## Now everything is built off of time point t=1
  for (t in 2:dim(series)[1]) {

    ##step 1, find Z_hat (based on previous Z)
    Z_hat <- Phi %*% prev_Z
    
    ## step 2, find epsilon
    ##   Both these terms are effectively residuals
    ##   and have mean zero, with Sigma_eps covariance
    lowbound <- (a[t,] - Z_hat)
    upbound  <- (b[t,] - Z_hat)

    for(i in 1:N) {

      ## Generate MV Truncated Normal with epsilon covariance
      epsilon[,i] <- mnormt::rmtruncnorm(1, mean=rep(0,d), varcov=Sigma_eps,
                                         lower=lowbound[,i], upper=upbound[,i])
      w_Z[i] <- mnormt::sadmvn(lower=lowbound[,i], upper=upbound[,i], 
                            mean=rep(0,2), varcov = Sigma_Z)
    }

    ## step 3, update Z for next iteration
    ##   The epsilon terms have variance Sigma_eps
    ##   so the Z terms have our desired covariance structure
    prev_Z <- Z_hat + epsilon

    ## step 4, update w
    ## w_Z is the current weight and is calculated above
    ##    for efficiency (one nested loop)
    ## then add to the accumulator

    log_w <- log_w + log(w_Z)
    
  }
  
  -log(N) + logSumExp(log_w)
  
}

fit_mvpois_ar1 <- function(series, chpts=NULL,
                           N=500, approx=TRUE) {
  approx=TRUE ## Override user for now, optim not working
  
  ### Initial values
  m <- length(chpts)
  d <- dim(series)[2]
  
  ##  Get lambda values based on sample means in each regime
  lambda_init <- c()
  regimes <- as.factor(findInterval(1:nrow(series), chpts) )
  GROUPINGS <- split(1:nrow(series), regimes)
  Regime_lengths <- sapply(GROUPINGS, length)
  if(m > 0) {
    for(i in 0:m) {
      lambda_init <- c(lambda_init, colMeans(series[regimes==i, ]) )
    }
  } else {
    lambda_init <- colMeans(series)
  }
  lambda_mat <- NULL
  for(i in 0:m) {
    lambda_mat <- rbind(
      lambda_mat,
      matrix(rep(lambda_init[(i*d+1):(i*d+d)],
                 each=sum(regimes==i)), ncol=d)
    )
  }
  
  ## From the determined lambda values, figure out the Z_series
  Z_series <- qnorm(ppois(series, lambda=lambda_mat) )
  
  ## Now determine an estimate for the variance-covariance matrix
  ##    for the underlying white noise process (truncated MV Normal)
  ## EVERYTHING is a function of Gamma0 and Gamma1
  ##   Calculation is relegated to the likelihood function
  ##   so we can eventually get optim working
  
  Gamma <- acf(Z_series, lag.max=1, type="covariance", plot=FALSE, demean=FALSE)$acf
  Gamma0 <- Gamma[1,,]
  Gamma1 <- Gamma[2,,]
  
  ## The follow are the Variance of the Z series
  ##   the AR coefficient matrix and the variance for the epsilons
  ##   they are computed here because we output the parameters
  Sigma_Z <- diag((1/sqrt(diag(Gamma0))))%*%Gamma0%*%diag((1/sqrt(diag(Gamma0))))
  Phi <- Gamma1%*%solve(Gamma0)

  A_mat <- diag(rep(1,d*d) ) - kronecker(Phi, Phi)
  Sigma_epsilon <- matrix(as.vector(Sigma_Z)%*%A_mat, nrow=d )
  Sigma_epsilon <- Sigma_Z - Phi%*%Sigma_Z%*%t(Phi)

  ## Construct the vector of parameters 
  ##   (transformed to eventually use in optim() )
  
  init_pars <- c(log(lambda_init), 
                 #tan(as.vector(Gamma1*pi/2)), 
                 as.vector(Gamma1),
                 t(chol(Gamma0))[lower.tri(t(chol(Gamma0)), diag=TRUE)] )
  
  ## Matrix of probability values -- currently not used
  #Rp = matrix(runif(dim(series)[1] * N), nrow = N)
  
  if(approx==TRUE) {
    ## Get Likelihood value from Particle Filtering
    loglik <- likelihood_mvpois_ar1(init_pars, series=series, chpts=chpts,
                                    N=N)
    ## Extra parameters for output
    Lambda <- matrix(lambda_init, nrow=m+1, byrow=TRUE)
    colnames(Lambda) <- colnames(series)
    rownames(Lambda) <- paste0("Regime ", 1:(m+1))
    MDL <- ifelse(m > 0,
                  -2*loglik + d*sum(log(Regime_lengths ) )/2 + log(m) + sum(log(chpts)),
                  -2*loglik + d*sum(log(Regime_lengths))/2)
    BIC <- -2*loglik + log(nrow(series))*(d*(m+1) + m)
  } else { ## This part is currently disabled since approx=TRUE
    opt_out <- optim(init_pars,
                     likelihood_mvpois_ar1,
                     series=series, chpts=chpts, Rp=Rp, 
                     method="BFGS", 
                     control = list(factr = 1e12, fnscale=-1,
                                    parscale=c(rep(1,d*(m+1)), rep(0.01,d*d), rep(0.01, length(init_pars)-(m+1)*d-d*d) ),
                                    trace=3) )
    loglik <- opt_out$value
  }
   list(loglik = loglik,
        Lambda = Lambda,
        Phi = Phi,
        Sigma_Z = Sigma_Z,
        Sigma_epsilon = Sigma_epsilon,
        MDL = MDL,
        BIC = BIC)
}


