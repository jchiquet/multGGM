rm(list=ls())
if (Sys.info()["user"] == 'jchiquet') {
  setwd("~/git/multiggm/code/")
} else {
  setwd("~/")
}
data.dir <- "../data/brca/tcga/"

source("functions_multiattributes.R")
library(Matrix)

## ===============================================================
## PART 0. FORMATING THE DATA
readANDformat <- FALSE

if (readANDformat) {
## get transcriptomic data (take a while..)
expr <- read.table(paste0(data.dir,"data_RNA_Seq_v2_expression_median.txt"), header=TRUE, row.names = 1)

## get proteomic data
rppa <- read.table(paste0(data.dir,"data_rppa.txt"), header=TRUE, sep="\t", row.names = 1)
rppa <- rppa[-which(rowSums(is.na(rppa)) > 0), ] # remove NA - of course better could be done

## removing doublon in terms of genes
vrppa <- sapply(strsplit(as.character(Reduce(rbind, strsplit(rownames(rppa), "|", fixed=TRUE))[, 1]), " ", fixed=TRUE), function(v) v[1])
rppa <- rppa[!duplicated(vrppa), ]; rownames(rppa) <- vrppa[!duplicated(vrppa)]

## stick to those genes in the expression set
genes <- intersect(rownames(rppa),  rownames(expr))
expr <- expr[match(genes, rownames(expr)), ]
rppa <- rppa[match(genes, rownames(rppa)), ]

## stick top common samples
samples <- intersect(colnames(expr), colnames(rppa))
expr <- t(expr[, match(samples, colnames(expr))])
rppa <- t(rppa[, match(samples, colnames(rppa))])

save(expr, rppa, file="BRCA_rppa+genes.RData")

} else {
  load(file="BRCA_rppa+genes.RData")
}

# heatmap(cor(expr), symm=TRUE)
# heatmap(cor(rppa), symm=TRUE)

## ===============================================================
## PART 1. MULTIATTRIBUTE NETWORK INFERENCE

## infer the multiattributes networks for CIT and NHU data sets
mnet <- inferMultinet(list(scale(expr), scale(rppa))); adj <- 1*(mnet != 0)

## ===============================================================
## PART 3. NETWORK VIZUALIZATION

library(blockmodels)
## Adjust a SBM to the inferred CIT network
my_model <- BM_bernoulli("SBM_sym",as.matrix(adj))
my_model$estimate(); SBM <- list()
SBM$cl <- apply(my_model$memberships[[which.max(my_model$ICL)]]$Z, 1, which.max)
SBMT$pi <- my_model$model_parameters[[which.max(my_model$ICL)]]$pi

## image(mnet_CIT[order(SBM_CIT$cl),order(SBM_CIT$cl)])
## image(mnet_NHU[order(SBM_NHU$cl),order(SBM_NHU$cl)])

library(igraph)
library(RColorBrewer)
pal   <- brewer.pal(10, "Set3")

g <- graph_from_adjacency_matrix(mnet, mode = "undirected", weighted = TRUE, diag = FALSE)
V(g)$class <- SBM$cl
V(g)$size <- 5
V(g)$frame.color <- "white"
V(g)$color <- pal[V(g)$class]
V(g)$label <- ""
E(g)$arrow.mode <- 0

pdf(file="CIT_multiattributenetwork.pdf")
plot(g, edge.width=E(g)$weight)
plot(g, edge.width=E(g)$weight, layout=layout_in_circle(g,order(V(g)$class)))
dev.off()
