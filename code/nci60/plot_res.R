rm(list=ls())
library(ggplot2)
library(reshape2)

load(file="res_desc.RData")

nedges.univar1 <- sapply(nets.univar1, function(x) sum(x!=0))
nedges.univar2 <- sapply(nets.univar2, function(x) sum(x!=0))
nedges.bivar   <- sapply(nets.bivar, function(x) sum(x!=0))

nedges <- 50

net1 <- nets.univar1[[match(TRUE, nedges.univar1 >= 2*nedges)]]
net2 <- nets.univar2[[match(TRUE, nedges.univar2 >= 2*nedges)]]
mnet <- nets.bivar[[match(TRUE, nedges.bivar >= 2*nedges)]]


data <- rbind(melt(net1), melt(net2), melt(net1 | net2), melt(net1 & net2),
              melt(mnet))
data$net <- rep(c("net1","net2","net1 | net2","net2 & net2","mnet"), each=91*91)

d <- ggplot(data, aes(x=Var1, y=Var2, fill=value)) + geom_tile()
d+facet_grid(.~net)

J12 <- c()
J13 <- c()
J23 <- c()

for (nedges in 1:500) {
  net1 <- which(nets.univar1[[match(TRUE, nedges.univar1 >= 2*nedges)]] !=0)
  net2 <- which(nets.univar2[[match(TRUE, nedges.univar2 >= 2*nedges)]] !=0)
  mnet <- which(nets.bivar[[match(TRUE, nedges.bivar >= 2*nedges)]] !=0)

  J12[nedges] <- length(intersect(net1,net2)) / length(union(net1,net2))
  
  J13[nedges] <- length(intersect(net1,mnet)) / length(union(net1,mnet))

  J23[nedges] <- length(intersect(net2,mnet)) / length(union(net2,mnet))
  
}

plot(J12, type="l", ylim=c(0,0.5))
lines(J13)
lines(J23)

