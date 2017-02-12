rm(list=ls())
gc()

## some fast well understood examples:

library(pls)	## to get their example data.
library(ExPosition)
source('Utilities_PLS_CompTests.R')
data(oliveoil)

X <- expo.scale(as.matrix(oliveoil[,"chemical"]))
Y <- expo.scale(as.matrix(oliveoil[,"sensory"]))


res <- simple.pls(X,Y)	
prettyScree(res$eigs)
bs.res <- broken.stick(res$eigs)
k.res <- kaiser(res$eigs)
boot.res <- boot.comps(X,Y,F,F,eigs=res$eigs)
perm.between.res <- perm.between.comps(X,Y,F,F,eigs=res$eigs)
perm.within.res <- perm.within.comps(X,Y,F,F,eigs=res$eigs)


