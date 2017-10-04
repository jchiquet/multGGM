rm(list=ls())
library(mvtnorm)
source("functions_utils.R")
source("functions_multiattributes.R")

p <- 20
n <- 100
hardeness <- 1
nedges <- p

G <- rgraph(p, nedges)

X <- rmultivar(G, n, K=4, hardeness)

univar.1se <- MultiVariateNS(X[[1]], cv.choice="1se")
univar.min <- MultiVariateNS(X[[1]], cv.choice="min")
univar.all <- MultiVariateNS(X[[1]], cv.choice="none")

bivar.1se <- MultiVariateNS(X, cv.choice="1se")
bivar.min <- MultiVariateNS(X, cv.choice="min")
bivar.all <- MultiVariateNS(X, cv.choice="none")

perf.roc(univar.1se$networks, G)
perf.roc(univar.min$networks, G)
perf.roc(bivar.1se$networks, G)
perf.roc(bivar.min$networks, G)

perf.auc(perf.roc(univar.all$networks, G))
perf.auc(perf.roc(bivar.all$networks, G))

