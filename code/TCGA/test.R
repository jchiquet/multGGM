rm(list=ls())
library(mvtnorm)
source("functions_utils.R")
source("functions_multiattributes.R")

p <- 20
n <- 50
hardeness <- 1
nedges <- p

G <- rgraph(p, nedges)

X <- rmultivar(G, n, K=4, hardeness)

Xuniv <- Reduce("rbind", X)
univar.1se <- MultiVariateNS(Xuniv, cv.choice="1se")
univar.min <- MultiVariateNS(Xuniv, cv.choice="min")
univar.all <- MultiVariateNS(Xuniv, cv.choice="none")

bivar.1se <- MultiVariateNS(X, cv.choice="1se")
bivar.min <- MultiVariateNS(X, cv.choice="min")
bivar.all <- MultiVariateNS(X, cv.choice="none")

perf.roc(univar.1se$networks, G)
perf.roc(univar.min$networks, G)
perf.roc(bivar.1se$networks, G)
perf.roc(bivar.min$networks, G)

perf.auc(perf.roc(univar.all$networks, G))
perf.auc(perf.roc(bivar.all$networks, G))

