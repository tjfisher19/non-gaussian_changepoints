library(tsDyn)

gen_mv_pois_series_no_change <- function(
    n=50, n_series=1,
    ar_cor = "pos", sigma_cor="neg",
    lambda = c(13.5, 17.8) ) {
  
  if(ar_cor=="pos") {
    Phi <- matrix(c(0.40, 0.10, 0.10, 0.40), 2)
  } else {
    Phi <- matrix(c(0.40, -0.10, -0.10, 0.40), 2)
  }
  if(sigma_cor=="neg") {
    Gamma0 <- cbind(c(1,-0.30), c(-0.30,1) )
  } else {
    Gamma0 <- cbind(c(1, 0.30), c(0.30,1) )
  }
  
  A_mat <- diag(rep(1,NROW(Phi)*NROW(Phi)) ) - kronecker(Phi, Phi)
  WN_Sigma <- matrix(as.vector(Gamma0)%*%A_mat, nrow=2 )
  
  ## Generate n_series + 5 worth of noise (burn in)
  Z_series <- VAR.sim(B=Phi, n=n*(n_series+5), 
                      include="none", varcov=WN_Sigma)
  Z_series <- Z_series[-(1:(n*5)),]   ## Remove burn in rows
  
  X_series <- list()

  for(i in 1:n_series) {
    X_series[[i]] <- matrix(nrow=n, ncol=length(lambda))
    for(j in 1:length(lambda)) {
      X_series[[i]][,j] <- qpois(pnorm(Z_series[(n*(i-1)+1):(n*i),j]), lambda=lambda[j])
    }
  }
  X_series
}




gen_mv_pois_series_changes <- function(
    n=50, n_series=1,
    ar_cor = "pos", sigma_cor="neg",
    lambda = rbind(c(10, 18), c(16.5, 17) ),
    chpts=c(20) ) {
  
  if(ar_cor=="pos") {
    Phi <- matrix(c(0.40, 0.10, 0.10, 0.40), 2)
  } else {
    Phi <- matrix(c(0.40, -0.10, -0.10, 0.40), 2)
  }
  if(sigma_cor=="neg") {
    Gamma0 <- cbind(c(1,-0.30), c(-0.30,1) )
  } else {
    Gamma0 <- cbind(c(1, 0.30), c(0.30,1) )
  }
  
  A_mat <- diag(rep(1,NROW(Phi)*NROW(Phi)) ) - kronecker(Phi, Phi)
  WN_Sigma <- matrix(as.vector(Gamma0)%*%A_mat, nrow=2 )
  
  ## Generate n_series + 5 worth of noise (burn in)
  Z_series <- VAR.sim(B=Phi, n=n*(n_series+5), 
                      include="none", varcov=WN_Sigma)
  Z_series <- Z_series[-(1:(n*5)),]   ## Remove burn in rows
  
  
  ### Get my Lambda values
  
  regimes <- as.factor(findInterval(1:n, chpts) )
  lambda_mat <- NULL
  for(i in 0:length(chpts)) {
    lambda_mat <- rbind(
      lambda_mat,
      matrix(rep(lambda[(i+1),],
                 each=sum(regimes==i)), ncol=ncol(lambda))
    )
  }

  X_series <- list()  
  for(i in 1:n_series) {
    X_series[[i]] <- matrix(nrow=n, ncol=ncol(lambda))
    for(j in 1:ncol(lambda)) {
      X_series[[i]][,j] <- qpois(pnorm(Z_series[(n*(i-1)+1):(n*i),j]), lambda=lambda_mat[,j])
    }
  }
  X_series
}

