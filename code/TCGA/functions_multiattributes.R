library(Matrix)
library(glmnet)
library(gglasso)
library(pbmcapply)

lasso <- function(x, y, group, lambda, cv.choice) {
  if (cv.choice == "none") {
    out <- glmnet(x, y, lambda=lambda, intercept=FALSE, standardize=FALSE)
    beta <- coef.glmnet(out)[-1,]
  } else {
    out <- cv.glmnet(x, y, lambda=lambda, intercept=FALSE, standardize=FALSE)
    beta <- matrix(coef.cv.glmnet(out, s=paste("lambda",cv.choice, sep="."))[-1], ncol=1)
  }
  abs(beta)
}

grplasso <- function(x, y, group, lambda, cv.choice) {
  if (cv.choice == "none") {
    out <- gglasso(x, y, group, lambda=lambda, intercept=FALSE)
    beta <- coef(out)[-1,]
  } else {
    out <- cv.gglasso(x, y, group, lambda=lambda, intercept=FALSE)
    beta <- matrix(coef(out, s=paste("lambda",cv.choice, sep="."))[-1], ncol=1)
  }
  beta <- apply(beta, 2, function(b) sqrt(rowsum(b^2,group)))
}

MultiVariateNS <- function(X, sym.rule="AND", cv.choice=c("1se", "min", "none"),
                           nlambda=50, min.ratio=1e-3, mc.cores=1) {

  cv.choice <- match.arg(cv.choice)
  LAPPLY <- ifelse(mc.cores > 1, pbmclapply, mclapply)
  
  stopifnot(is.list(X) | is.matrix(X)) 
  if (is.matrix(X)) {
    warning("even if a single omic is provided, should be in a list")
    X <- list(X)
  }
  
  ## scaling to make everything comparable
  X <- lapply(X, scale)
  
  ## problem dimension
  K <- length(X)
  p <- ncol(X[[1]])
  n <- nrow(X[[1]])
  
  ## define group indexes
  group <- rep(rep(1:(p-1),K), K)

  ## create the list of data sets corresponding to each variable
  training_sets <- vector("list", p)
  for (i in 1:p) {
    training_sets[[j]] <- list(
      y   = c(sapply(X, function(x) x[, j])),
      x   = as.matrix(bdiag(rep(list(do.call(cbind, lapply(X, function(x) x[, -j]))), K))),
      var = j
    )
  }

  ## Get the largest lambda value across all variable that produce a null vector as estimate
  lambda.max <- max(sapply(training_sets, function(data) {
     return(max(sqrt(rowsum(crossprod(data$x, data$y)^2, group)/(2*n))))
  }))

  lambda <- 10^seq(from=log10(lambda.max),log10(min.ratio*lambda.max),len=nlambda)

  ## compute all coefficients
  f_infer <- ifelse(K==1, lasso, grplasso)
  
  coefficients <- do.call(rbind,LAPPLY(training_sets, function(d) {
      coef <- f_infer(d$x, d$y, group, lambda, cv.choice)
      beta  <- matrix(0, p, ncol(coef))
      beta[-d$j, ] <- as.matrix(coef)
      beta
  }, mc.cores=mc.cores))
  
  networks <- apply(coefficients, 2, function(beta) {
    net <- matrix(beta != 0,p,p)
    net <- Matrix(switch(sym.rule,
                   "AND" = pmin(net,t(net)),
                   "OR"  = pmax(net,t(net))))
    net
  })

  return(list(coefficients = coefficients, networks = networks, penalties = lambda))
}
