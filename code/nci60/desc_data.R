rm(list=ls())
library(R.utils)
library(mvtnorm)
library(quadrupen)
library(scoop)
library(gRc)

source("functions_multivar.R")
load("gpconsensus.RData")

X1 <-  scale(t(dat$gene.select), TRUE,TRUE)
X2 <-  scale(t(dat$protein.select), TRUE,TRUE)

cols.type <- as.numeric(as.factor(gsub(":.*", "",dat$sites)))

expr.PCA <- PCA(X1, graph = FALSE)
plot(expr.PCA, choix="ind", col.ind=cols.type, axes=c(1,2))

prot.PCA <- PCA(X2, graph = FALSE)
plot(prot.PCA, choix="ind", col.ind=cols.type, axes=c(1,2))

plot(expr.PCA, choix="var")
plot(prot.PCA, choix="var")

p <- ncol(X1)
n <- nrow(X1)

luniv1  <- getlambda(list(X1), "univar", lmin=0.1, lsize=1000)
luniv2  <- getlambda(list(X2), "univar", lmin=0.1, lsize=1000)
lbivar  <- getlambda(list(X1,X2), "multivar", lmin=0.1, lsize=1000)

cat("\nMULTIVARIATE INFERENCE")
cat("\ncurrent variable is")
outbivar <- multivariate(list(X1,X2), lambda=lbivar)

cat("\nUNIVARIATE INFERENCE - FIRST ATTRIBUTE")
cat("\ncurrent variable is")
outunivar1 <- univariate(X1, lambda=luniv1)

cat("\nUNIVARIATE INFERENCE - SECOND ATTRIBUTE")
cat("\ncurrent variable is")
outunivar2 <- univariate(X2, lambda=luniv2)

nets.univar1 <- extract.net(outunivar1)
nets.univar2 <- extract.net(outunivar2)
nets.bivar   <- extract.net(outbivar)

ind <- 10
expr.var <- glasso(cor(X1), 0, zero = which(nets.univar1[[ind]] == 0, arr.ind = TRUE), penalize.diagonal = FALSE)$w
prot.var <- glasso(cor(X2), 0, zero = which(nets.univar2[[ind]] == 0, arr.ind = TRUE), penalize.diagonal = FALSE)$w

T.expr <- X1 %*% svd(cov2cor(expr.var))$v
T.prot <- X2 %*% svd(cov2cor(prot.var))$v

par(mfrow=c(1,2))
plot(T.expr[,1], T.expr[,2], col=cols.type, pch=20)
plot(T.prot[,2], T.prot[,3], col=cols.type, pch=20)

expr.bivar <- glasso(cor(X1), 0, zero = which(nets.bivar[[ind]] == 0, arr.ind = TRUE), penalize.diagonal = FALSE)$w
prot.bivar <- glasso(cor(X2), 0, zero = which(nets.bivar[[ind]] == 0, arr.ind = TRUE), penalize.diagonal = FALSE)$w

T.expr.bivar <- X1 %*% svd(cov2cor(expr.bivar))$v
T.prot.bivar <- X2 %*% svd(cov2cor(prot.bivar))$v

par(mfrow=c(1,2))
plot(T.expr.bivar[,1], T.expr.bivar[,2], col=cols.type, pch=20)
plot(T.prot.bivar[,2], T.prot.bivar[,3], col=cols.type, pch=20)

save(nets.univar1, nets.univar2, nets.bivar, file="res_desc.RData")

