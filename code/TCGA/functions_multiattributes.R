library(parallel)
library(Matrix)
library(glmnet)
library(grpreg)

inferMultinet <- function(X, sym.rule="AND", cv.choice="1se", verbose=TRUE, mc.cores=10) {
    
    stopifnot(is.list(X))
    K <- length(X)
    p <- ncol(X[[1]])
    n <- nrow(X[[1]])
    if (verbose) {cat("\ncurrent block is")}
    res <- do.call(rbind, mclapply(1:p,
       function(i) {
           if(verbose){cat("",i)}
           y <- c(sapply(X, function(x) x[, i]))
           x <- bdiag(rep(list(do.call(cbind, lapply(X, function(x) x[, -i]))), K))
           group.rem <- rep(rep(1:(p-1),K),K)
           
           out <- cv.grpreg(as.matrix(x),y,group.rem)
           beta <- rep(0,p)
           lambda.cv <- switch(cv.choice,
                               "min" = out$lambda.min,
                               "1se" = max(out$lambda[out$cve < min(out$cve)+out$cvse]))
           beta[-i] <- predict(out, type="norm", lambda=lambda.cv)
           
           return(beta)
       }, mc.cores=mc.cores))
    
    return(Matrix(switch(sym.rule, "AND" = pmin(res,t(res)), "OR"  = pmax(res,t(res))),
                  dimnames=list(colnames(X[[1]]),colnames(X[[1]]))))
    
}

inferNet <- function(X, sym.rule="AND", cv.choice="1se", verbose=TRUE, mc.cores=10) {

    stopifnot(is.matrix(X))
    p <- ncol(X)
    n <- nrow(X)
    if (verbose) {cat("\ncurrent block is")}
    res <- do.call(rbind, mclapply(1:p,
       function(i) {
           if(verbose){cat("",i)}
           out <- cv.glmnet(X[, -i], X[, i], intercept=FALSE, standardize=FALSE)
           lambda.cv <- switch(cv.choice, "min" = out$lambda.min, "1se" = out$lambda.1se)
           beta <- rep(0,p); beta[-i] <- coef(out, ls=lambda.cv)[-1] # remove intercept, still added even when intercept=FALSE
           beta
       }, mc.cores=mc.cores))
    
    return(Matrix(switch(sym.rule, "AND" = pmin(res,t(res)), "OR"  = pmax(res,t(res))),
                  dimnames=list(colnames(X[[1]]),colnames(X[[1]]))))
    
}


