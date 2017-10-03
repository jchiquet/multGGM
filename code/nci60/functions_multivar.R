
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

group.norm <- function(beta,group,norm=2) {
  return(switch(norm,
                "1"   = c(tapply(abs(beta), group,sum)),
                "2"   = c(sqrt(tapply(beta^2, group,sum))),
                "inf" = c(tapply(abs(beta), group, max))))
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

getlambda <- function(X, type="univar", lmin=1e-3, lsize=50, norm=2) {
  ## get the smallest lambda to produce null vector for the
  ## lasso/elastic net fit among the K samples

  getlmax <- switch(type,
                   "univar"   = getlmax.univar,
                   "multivar" = getlmax.multivar)
  
  lmax <- getlmax(X, norm)
  lsize <- 50 # number of penalties
  return(10^seq(from=log10(max(lmax)),log10(lmin*lmax),len=lsize))  
}

getlmax.univar <- function(Xs, norm) {
    
  K <- length(Xs)
  n <- nrow(Xs[[1]])
  p <- ncol(Xs[[1]])
  
  return(max(sapply(1:p, function(i) {
    y <- c(sapply(Xs, function(x) x[, i]))
    x <- matrix(0,K*n,K*(p-1))
    for (k in 1:K) {
      x[((k-1)*n+1):(k*n),((k-1)*(p-1)+1):(k*(p-1))] <- Xs[[k]][, -i]
    }
    return(max(abs(crossprod(y,x))))})))
}

getlmax.multivar <- function(Xs, norm) {

  K <- length(Xs)
  n <- nrow(Xs[[1]])
  p <- ncol(Xs[[1]])
  
  return(max(sapply(1:p, function(i) {  
    y <- c(sapply(Xs, function(x) x[, i]))
    x <- matrix(0,K*n,K*(p-1))
    B <- c()
    for (k in 1:K) {
      B <- cbind(B, Xs[[k]][, -i])
    }
    
    ## Blockwise contruciton of the huge X matrix K (block diagonal with
    ## K block equal to
    ## X[[1]] ... X[[K]
    x <- matrix(0,K*n,K*K*(p-1))
    for (k in 1:K) {
      x[((k-1)*n+1):(k*n),((k-1)*K*(p-1)+1):(k*K*(p-1))] <- B
    }
    
    group <- rep(rep(1:(p-1),K),K)
    pk <- tabulate(group)
    penscale <- switch(norm, "2" = sqrt(pk), "inf" = rep(1,length(pk)))
    
    return(max(group.norm(crossprod(y,x),group,switch(norm, "inf" = 1, 2))/penscale))})))

}

multivariate <- function(X, lambda=getlambda(X, "multivar"), norm=2, gamma=0.01, verbose=TRUE) {

  require(scoop)
  
  stopifnot(is.list(X))
  K <- length(X)
  p <- ncol(X[[1]])
  n <- nrow(X[[1]])
  if (verbose) {cat("\ncurrent block is")}
  return(sapply(1:p, function(i) {
    if(verbose){cat("",i)}
    y <- c(sapply(X, function(x) x[, i]))
    x <- matrix(0,K*n,K*(p-1))
    B <- c()
    for (k in 1:K) {
      B <- cbind(B, X[[k]][, -i])
    }
    
    ## Blockwise contruction of the huge X matrix K (block diagonal with
    ## K block equal to
    ## X[[1]] ... X[[K]
    x <- matrix(0,K*n,K*K*(p-1))
    for (k in 1:K) {
      x[((k-1)*n+1):(k*n),((k-1)*K*(p-1)+1):(k*K*(p-1))] <- B
    }
    group.rem <- rep(rep(1:(p-1),K),K)
    o <- order(group.rem,decreasing=FALSE)
    x       <- x[,o]
    group.rem   <- group.rem[o]
    
    out <- group.lasso(x,y,group.rem,lambda=lambda,intercept=FALSE,normalize=FALSE)
    
    out@coefficients  <- out@coefficients[, order(o)]
    
    betai <- t(apply(out@coefficients,1,group.norm,group.rem[order(o)],norm))
    
    beta  <- matrix(0, length(lambda), p)
    beta[, -i] <- betai !=0
    return(t(beta))
  }))
}

## The resulting network is stock is a matrix with nlambda x p and
## column, so as the first 30 rows a the network for 

univariate <- function(X, lambda=getlambda(list(X), "univar"), gamma=0.01, verbose=TRUE) {
  require(quadrupen)
  
  p <- ncol(X)
  if (verbose) {cat("\ncurrent variable is")}
  return(sapply(1:p, function(i) {
    if (verbose) {cat("",i)}
    y <- X[, i]
    x <- X[,-i]
    out <- elastic.net(x,y,lambda1=lambda,lambda2=gamma,intercept=FALSE,normalize=FALSE)
    betai <- as.matrix(out@coefficients) != 0
    beta  <- matrix(0, length(lambda), p)
    beta[, -i] <- betai
    return(t(beta))
  }))
}

extract.net <- function(nets, sym.rule="AND") {

  p <- ncol(nets)
  lambda.size <- nrow(nets) %/% p

  getOneNet <- function(l,nets) {
    net <- nets[(p*(l-1)+1):(l*p), ]
    net <- switch(sym.rule,
                  "AND" = pmin(net,t(net)),
                  "OR"  = pmax(net,t(net)))
  }
  
  return(sapply(1:lambda.size, getOneNet, nets, simplify=FALSE))
  
}

AUC <- function(nets, ref, sym.rule="AND", plot=FALSE, main="ROC curve") {

  p <- ncol(nets)
  lambda.size <- nrow(nets) %/% p
  one.ROC.point <- function(l, nets, ref) {
    net <- nets[(p*(l-1)+1):(l*p), ]
    net <- switch(sym.rule,
                  "AND" = pmin(net,t(net)),
                  "OR"  = pmax(net,t(net)))
    
    ## Removing diagonal terms...
    diag(ref) <- 0
    diag(net) <- 0
    
    ## List of infered edges
    inferred <- which(abs(net[upper.tri(net)]) > 0)
  
    ## List of true edges
    true <- which(abs(ref[upper.tri(ref)]) > 0)

    ## Computing TP, FP, ...
    TP <- sum(inferred %in% true)
    FP <- length(inferred) - TP
    TN <- length(ref[upper.tri(ref)]) - length(true) - FP
    FN <- length(true) - TP
    
    recall  <- TP/length(true)
    fallout <- FP/(FP+TN)
        
    return(c(fallout,recall))
  }
  ROC <- sapply(1:lambda.size, one.ROC.point, nets, ref)
  specificity <- c(0, ROC[1, ], 1)
  sensitivity <- c(0, ROC[2, ], 1)
  dx <- diff(specificity)
  AUC <- sum(c(sensitivity[-1]*dx, sensitivity[-length(sensitivity)]*dx))/2

  if (plot) {
    plot(specificity, sensitivity, xlab="specificity", ylab="sensitivity", type="l", main=main)
    abline(0,1, lty=3)
    text(0.7,0.05, paste("AUC=",round(AUC,3)))
  }

  return(AUC)
}
