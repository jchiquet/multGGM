library(pbmcapply)

rgraph <- function(p,nedges) {

  A <- matrix(0,p,p)
  A[sample(which(upper.tri(A)), nedges)] <- 1

  return(A + t(A))
}

oneRun <- function(G, n, K, hardeness, pb) {
  increase(pb)
  X <- rmultivar(G, n, K, hardeness)
  lambda.u <- getlambda(X, "univar")
  res.AUC <-  c()
  for (k in 1:K) {
    res.AUC <- c(res.AUC,AUC(univariate(X[[k]], lambda.u, verbose=FALSE), G))
  }
  res.AUC <- c(res.AUC, AUC(multivariate(X, verbose=FALSE), G))
  names(res.AUC) <- c(paste("univar", 1:K), "multivar")
  return(res.AUC)
}

learn.networks <- function(nsim, p, n, K, hardeness=1, nedges=p, verbose=TRUE) {

  ## cat("\nSimulations will end once the Progression Bar will have printed 100 '.' or 10 '|'\n\n")
  pb <- ProgressBar(100, stepLength=100/nsim)
  return(replicate(nsim, oneRun(rgraph(p,nedges), n, K, hardeness, pb)))  
}

rmultivar <- function(G, n, K = 2, hardness = 1, eps=.Machine$double.eps, scale=TRUE) {

  p <- ncol(G) # number of nodes per graph
  ## split the adjacency matrix to multivariate space with the kronecker product
  A <- (G + diag(rep(1,p)) ) %x% matrix(1,K,K)

  ## Remove axes from with null eigen values to make it invertible
  eigenA <- eigen(A)
  D <- eigenA$values
  U <- eigenA$vectors
  D[D < eps] <- eps
  B <- U %*% diag(D) %*% t(U)

  ## The concentration matrix is built from B
  ## with hardness manage thourgh diagonal smoothing
  
  ## Then, generate the associate multivariate Gaussian sample
  X <- rmvnorm(n,sigma=solve(B + hardness * diag(diag(B))), method="chol")
  
  Xs <- sapply(1:K,function(k) {X[,seq(k,(K*p),by=K)]},simplify=FALSE)
  
  return(lapply(Xs, scale,scale,scale))
}

perf.roc <- function(theta.hat, theta.star) {
  
  roc <- function(theta) {
    nzero <- which(theta != 0)
    zero  <- which(theta == 0)

    true.nzero <- which(theta.star != 0)
    true.zero  <- which(theta.star == 0)

    TP <- sum(nzero %in% true.nzero)
    TN <- sum(zero %in%  true.zero)
    FP <- sum(nzero %in% true.zero)
    FN <- sum(zero %in%  true.nzero)

    recall    <- TP/(TP+FN) ## also recall and sensitivity
    fallout   <- FP/(FP+TN) ## also 1 - specificit

    res <-  round(c(fallout,recall),3)
    res[is.nan(res)] <- 0
    names(res) <- c("fallout","recall")
    return(res)
  }

  if (is.list(theta.hat)) {
    return(as.data.frame(do.call(rbind, lapply(theta.hat, roc))))
  } else {
    return(roc(theta.hat))
  }
}

perf.auc <- function(roc) {
  fallout <- c(0,roc$fallout,1)
  recall  <- c(0,roc$recall, 1)
  dx <- diff(fallout)
  return(sum(c(recall[-1]*dx, recall[-length(recall)]*dx))/2)
}
