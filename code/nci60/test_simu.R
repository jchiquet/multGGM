rm(list=ls())
library(simone)
library(quadrupen)
library(scoop)
library(mvtnorm)

source("functions_multivar.R")

p <- 20    # size of the network
n <- 100   # number of experiments
K <- 2     # number of task

## a random network
G <- graph.affiliation(p,c(1/3,1/3,1/3),0.1,0.1)$x
## a list of K Gaussian samples
X <- rmultivar(G, n, K)

cat("\nUNIVARIATE INFERENCE - FIRST ATTRIBUTE")
outunivar1 <- univariate(X[[1]], getlambda(X, "univar"))

cat("\nUNIVARIATE INFERENCE - SECOND ATTRIBUTE")
outunivar2 <- univariate(X[[2]], getlambda(X, "univar"))

cat("\nMULTIVARIATE INFERENCE")
outmultivar <- multivariate(X)

par(mfrow=c(1,3))
AUC.univar1  <- AUC(outunivar1 , G, plot=TRUE, main="Univariate Inference onn Task 1")
AUC.univar2  <- AUC(outunivar2 , G, plot=TRUE, main="Univariate Inference onn Task 2")
AUC.multivar <- AUC(outmultivar, G, plot=TRUE, main="Multivariate Inference")
