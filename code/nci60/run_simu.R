rm(list=ls())
library(R.utils)
library(mvtnorm)
source("functions_multivar.R")

## fixed coefficients
nsim <- 50
p <- 20
hardeness <- 1
out <- list()

k <- 0
for (K in c(2,4,8)) {
  cat("\n\n===========================================================")
  cat("\nK =",K)
  
  out[[k]] <- sapply(c(p/2,p,2*p), function(n) {
    cat("\n - n =", n," ")
    return(learn.networks(nsim, p, n, K, hardeness))
  })
}

noms <- c("univar.1","univar.1","multivariate")
  
par(mfrow=c(1,3))
boxplot(matrix(out[[1]][, 1], ncol=3, byrow=TRUE), las=3, names=noms, main="n=p/2", ylim=c(0.35, 1))
boxplot(matrix(out[[1]][, 2], ncol=3, byrow=TRUE), las=3, names=noms, main="n=p",  ylim=c(0.35, 1))
boxplot(matrix(out[[1]][, 3], ncol=3, byrow=TRUE), las=3, names=noms, main="n=2p", ylim=c(0.35, 1))

##boxplot(t(out), main="AUC")


