rm(list=ls())
library(tidyverse)
library(viridis)
library(Matrix)

load(file="RATHER_prot+expr-2017-10-13.RData")

image(protNet.1se$networks[[1]])
image(exprNet.1se$networks[[1]])
image(bivarNet.1se$networks[[1]])

nedges.prot <- sapply(protNet.all$networks, sum)
nedges.expr <- sapply(exprNet.all$networks, sum)
nedges.mult <- sapply(bivarNet.all$networks, sum)

max.edges <- 1000
J12 <-vector("numeric", max.edges)
J13 <-vector("numeric", max.edges)
J23 <-vector("numeric", max.edges)
for (nedges in 1:max.edges) {
  prot <- which(protNet.all$networks[[match(TRUE, nedges.prot >= 2*nedges)]] !=0)
  expr <- which(exprNet.all$networks[[match(TRUE, nedges.expr >= 2*nedges)]] !=0)
  mult <- which(bivarNet.all$networks[[match(TRUE, nedges.mult >= 2*nedges)]] !=0)

  J12[nedges] <- length(intersect(prot,expr)) / length(union(prot,expr))

  J13[nedges] <- length(intersect(prot,mult)) / length(union(prot,mult))

  J23[nedges] <- length(intersect(expr,mult)) / length(union(expr,mult))

}

dp <- data.frame(jaccard = c(J12, J13, J23),
           couple = rep(c("protein/expr.","protein/expr+protein","expr/expr+protein"), each=max.edges),
           edges = rep(1:max.edges, 3))
p <- ggplot(dp, aes(x=edges, y=jaccard, group=couple, colour=couple)) + geom_line() +
  theme_bw(base_size = 20) + scale_color_viridis(discrete=TRUE) + ylim(0.,0.4)


