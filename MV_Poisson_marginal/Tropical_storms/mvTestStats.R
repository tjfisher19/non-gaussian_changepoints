

RobbinsBlockSingle <- function(x, y, lag.spot )
{
  x <- as.matrix(x);
  y <- as.matrix(y);
  d.x <- NCOL(x);
  d.y <- NCOL(y);
  k <- d.x + d.y;
  n <- NROW(x);  # We can error check later

  out <- matrix(0, nrow = k, ncol = k)
  out <- diag(k);
  X <- cbind(x,y);
  Accmat <- acf(X, lag.max = abs(lag.spot), plot = FALSE, type = "covariance", na.action=na.pass)$acf
  C0 <- Accmat[1, ,];
  C01 <- C0[1:d.x, 1:d.x];
  C02 <- C0[(d.x+1):k, (d.x+1):k];
  C0[(d.x+1):k, 1:(d.x)] <- 0;
  C0[1:(d.x), (d.x+1):k] <- 0;
  L01 <- t(chol(solve(C01)));
  L02 <- t(chol(solve(C02)));
  C0.inv <- solve(C0);
  
  L <- t(chol(C0.inv));
  if(lag.spot >= 0)
    Ch <- Accmat[abs(lag.spot)+1, 1:d.x, (d.x+1):k]
  else
    Ch <- t(Accmat[abs(lag.spot)+1, (d.x+1):k, 1:d.x] );
#  print(Ch)
#  print(dim(L02))
#  print(dim(Ch))
#  print(dim(L01))
  #tmp <- sqrt(n/(n-abs(lag.spot)))*crossprod(t(crossprod(L01,Ch)), L02);
  tmp <- crossprod(t(crossprod(L01,Ch)), L02)
  #tmp <- sqrt(n/(n-abs(lag.spot)))*crossprod(t(crossprod(L, Accmat[abs(lag.spot)+1,,])), L);
#  if(lag.spot >= 0)
#    tmp <- tmp[(1:d.x), (d.x+1):k]
#  else
#    tmp <- tmp[(d.x+1):k, (1:d.x)];
#  out[(1:d.x),(d.x+1):k] <- tmp;
#  out[(d.x+1):k, (1:d.x)] <- t(tmp);
  out <- diag(k);
  out[1:d.x,(d.x+1):k]=tmp;
  out[(d.x+1):k,1:(d.x)]=t(tmp);
  return(out);
}



sumLogMatrixTest <- function(x, lag.max=5, fitdf, ljung=FALSE, weighted=FALSE) {
  x <- as.matrix(x);
  n <- NROW(x);
  d.x <- NCOL(x);
  
  log.terms <- rep(0, lag.max)
  for(i in 1:lag.max) {
    mat <- RobbinsBlockSingle(x,x, lag.spot=i);
    log.terms[i] <- -n*log(abs(det(mat)));
  }
  if(ljung) weights <- n/(n-1:lag.max)
  else weights <- rep(1, lag.max);
  if(weighted) more.weights <- (lag.max:1)/lag.max
  else more.weights <- rep(1, lag.max);
  stat <- sum(more.weights*weights*log.terms)
  if(weighted) {
    K1 <- d.x*d.x*sum(more.weights) - d.x*d.x*fitdf
    K2 <- 2*d.x*d.x*sum(more.weights^2) - 2*d.x*d.x*fitdf
    shape <- K1^2/K2
    scale <- K2/K1
    pval <- pgamma(stat, shape=shape, scale=scale, lower.tail=FALSE);
  } else {
    pval <- pchisq(stat, df=d.x*d.x*(lag.max-fitdf), lower.tail=FALSE);
  }
  c(stat, pval);
}



geoSumLogMatrixTest <- function(x, lag.max=5, fitdf, ljung=FALSE, weight=0.9) {
  x <- as.matrix(x);
  n <- NROW(x);
  d.x <- NCOL(x);
  
  log.terms <- rep(0, lag.max)
  for(i in 1:lag.max) {
    mat <- RobbinsBlockSingle(x,x, lag.spot=i);
    log.terms[i] <- -n*log(abs(det(mat)));
  }
  if(ljung) v <- n/(n-1:lag.max)
  else v <- rep(1, lag.max);
  
  if(fitdf==0)
    weights <- weight^(0:(lag.max-1))
  else
    weights <- fitdf*weight^(0:(lag.max-1));
  
  stat <- sum(weights*v*log.terms)
  
  cumulant1 <- d.x*d.x*sum(weights);
  cumulant2 <- 2*d.x*d.x*sum(weights^2) - 2*d.x*d.x*fitdf;
  shape <- cumulant1^2/cumulant2;
  scale <- cumulant2/cumulant1;
  pval <- pgamma(stat, shape=shape, scale=scale, lower.tail=FALSE);
  
  c(stat, pval);
}

HoskingTest <- function(x, lag, fitdf, weighted=FALSE) {
  acfMat <- acf(x, lag.max=lag, plot=FALSE, type="correlation", na.action=na.pass)$acf;
  invR0 <- solve(acfMat[1,,]);
  n <- dim(x)[1];
  k <- dim(x)[2];
  
  foo <- function(j) {
    tvecR <- t(as.vector(acfMat[j+1, , ]));
    crossprod(t(tvecR), crossprod(t(kronecker(invR0, invR0)), t(tvecR)));
  }
  ind <- 1:lag
  
  prodVec <- sapply(ind, foo);
  if(weighted) weights <- (lag:1)/lag
  else weights <- rep(1,lag)
  Hosking <- n*sum(weights*prodVec);
  if(weighted) {
    K1 <- k*k*sum(weights )-k*k*fitdf
    K2 <- 2*k*k*(sum(weights^2)-fitdf)
    shape <- K1^2/K2
    scale <- K2/K1
    pval <- pgamma(Hosking, shape=shape, scale=scale, lower.tail=FALSE);
  } else {
    pval <- pchisq(Hosking, df=k*k*(lag-fitdf), lower.tail=FALSE);
  }
  c(Hosking, pval)
}

LBHoskingTest <- function(x, lag, fitdf, weighted=FALSE) {
  acfMat <- acf(x, lag.max=lag, plot=FALSE, type="correlation", na.action=na.pass)$acf;
  invR0 <- solve(acfMat[1,,]);
  n <- dim(x)[1];
  k <- dim(x)[2];
  
  foo <- function(j) {
    tvecR <- t(as.vector(acfMat[j+1, , ]));
    crossprod(t(tvecR), crossprod(t(kronecker(invR0, invR0)), t(tvecR)));
  }
  ind <- 1:lag
  
  prodVec <- sapply(ind, foo);
  if(weighted) weights <- (lag:1)/lag
  else weights <- rep(1, lag)
  Hosking <- n*sum(n/(n-ind)*weights*prodVec);
  if(weighted) {
    K1 <- k*k*sum(weights)-k*k*fitdf
    K2 <- 2*k*k*(sum(weights^2)-fitdf)
    shape <- K1^2/K2
    scale <- K2/K1
    pval <- pgamma(Hosking, shape=shape, scale=scale, lower.tail=FALSE);
  } else {
    pval <- pchisq(Hosking, df=k*k*(lag-fitdf), lower.tail=FALSE);
  }
  c(Hosking, pval)
}


