




fit_mvpois_iid <- function(series, chpts=NULL) {
  
  ### Initial values
  m <- length(chpts)
  d <- dim(series)[2]
  
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

  loglik <- sum(dpois(series, lambda=lambda_mat, log=TRUE))
  
  Lambda <- matrix(lambda_init, nrow=m+1, byrow=TRUE)
  rownames(Lambda) <- paste("Regime", 1:(m+1) )
  colnames(Lambda) <- colnames(series)
  
  MDL <- ifelse(m > 0,
                -2*loglik + sum(log(Regime_lengths ) )/2*d + log(m) + sum(log(chpts)),
                -2*loglik + sum(log(Regime_lengths))/2*d)
  BIC <- -2*loglik + log(nrow(series))*(d*(m+1) + m)
  
  list(loglik = loglik,
       Lambda = Lambda,
       MDL = MDL,
       BIC = BIC)
}
