rm(list=ls())
working_directory <-
  switch(Sys.info()["user"],
    "jchiquet" = '~/git/multGGM/code/RATHER/'
## put your own directory here
  )
setwd(working_directory)

library(multivarNetwork)

rppa <- read.table('RPPA_RATHER_ILC.csv', sep=",")
expr <- read.table('expr_data_processed.csv', sep="\t")
mc.cores <- 10

cat("\nMULTIVARIATE INFERENCE")
cat("\ncurrent variable is")
bivarNet.min <- multivarNetwork(list(expr,prot), cv.choice="min", mc.core=mc.cores)
bivarNet.1se <- multivarNetwork(list(expr,prot), cv.choice="1se", mc.core=mc.cores)
bivarNet.all <- multivarNetwork(list(expr,prot),mc.core=mc.cores, nlambda=1000)

cat("\nUNIVARIATE INFERENCE - FIRST ATTRIBUTE")
cat("\ncurrent variable is")
exprNet.min <- multivarNetwork(expr, cv.choice="min", mc.core=mc.cores)
exprNet.1se <- multivarNetwork(expr, cv.choice="1se", mc.core=mc.cores)
exprNet.all <- multivarNetwork(expr, mc.core=mc.cores, nlambda=1000)

cat("\nUNIVARIATE INFERENCE - SECOND ATTRIBUTE")
cat("\ncurrent variable is")
protNet.min <- multivarNetwork(prot, cv.choice="min", mc.core=mc.cores)
protNet.1se <- multivarNetwork(prot, cv.choice="1se", mc.core=mc.cores)
protNet.all <- multivarNetwork(prot, mc.core=mc.cores, nlambda=1000)

save(protNet.min, protNet.1se, protNet.all,
     exprNet.min, exprNet.1se, exprNet.all,
     bivarNet.min, bivarNet.1se, bivarNet.all,  file=paste0("nci60_prot+expr-",Sys.Date(),".RData") )

