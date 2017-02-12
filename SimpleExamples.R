rm(list=ls())
gc()

## some fast well understood examples:


library(ExPosition)
source('Utilities_PLS_CompTests.R')
data(beer.tasting.notes)

## This is a good example of exactly 1 component.
dat.X <- expo.scale(beer.tasting.notes$data)
# dat.X_center <- matrix(attributes(dat.X)$`scaled:center`,nrow(dat.X),ncol(dat.X),byrow=T)
# dat.X_scale <- matrix(attributes(dat.X)$`scaled:scale`,nrow(dat.X),ncol(dat.X),byrow=T)
# dat.X_recon <- dat.X * dat.X_scale + dat.X_center


res <- epPCA(dat.X,F,F,graphs=F)
## Now reconstruct dat.X with exactly 1 underlying source component, but put it in different places with the same (fake) singular values
dat.X1 <- cbind(	
	res$ExPosition.Data$pdq$p[,1],
	0,
	0
) %*% diag(c(3,2,1)) %*% t(res$ExPosition.Data$pdq$q[,1:3])
# dat.X1_recon <- dat.X1 * dat.X_scale + dat.X_center
# expo.scale(dat.X1_recon) / dat.X1

dat.X2 <- cbind(	
	0,
	res$ExPosition.Data$pdq$p[,1],
	0
) %*% diag(c(3,2,1)) %*% t(res$ExPosition.Data$pdq$q[,1:3])
#dat.X2_recon <- dat.X2 * dat.X_scale + dat.X_center
#expo.scale(dat.X2_recon) / dat.X2

dat.X3 <- cbind(	
	0,
	0,
	res$ExPosition.Data$pdq$p[,1]
) %*% diag(c(3,2,1)) %*% t(res$ExPosition.Data$pdq$q[,1:3])
#dat.X3_recon <- dat.X3 * dat.X_scale + dat.X_center

dat.X1_3 <- cbind(	
	res$ExPosition.Data$pdq$p[,1],
	0,
	res$ExPosition.Data$pdq$p[,3]
) %*% diag(c(3,2,1)) %*% t(res$ExPosition.Data$pdq$q[,1:3])




	# ## between data and itself
# XX.res <- simple.pls(dat.X,dat.X)	
# #prettyScree(XX.res$eigs)
# bs.res <- broken.stick(XX.res$eigs)
# k.res <- kaiser(XX.res$eigs)
# boot.res <- boot.comps(dat.X,dat.X,F,F,eigs=XX.res$eigs)
# perm.between.res <- perm.between.comps(dat.X,dat.X,F,F,eigs=XX.res$eigs)
# perm.within.res <- perm.within.comps(dat.X,dat.X,F,F,eigs=XX.res$eigs)



	## between data and a single low variance source of itself
# XX3.res <- simple.pls(dat.X,dat.X3)	
# #prettyScree(XX3.res$eigs)
# #bs.res <- broken.stick(XX3.res$eigs)
# #k.res <- kaiser(XX3.res$eigs)
# boot.res <- boot.comps(dat.X,dat.X3,F,F,eigs=XX3.res$eigs)
# perm.between.res <- perm.between.comps(dat.X,dat.X3,F,F,eigs=XX3.res$eigs)
# perm.within.res <- perm.within.comps(dat.X,dat.X3,F,F,eigs=XX3.res$eigs)

XX1_3.res <- simple.pls(dat.X,dat.X1_3)	
#prettyScree(XX1_3.res$eigs)
bs.res <- broken.stick(XX1_3.res$eigs)
k.res <- kaiser(XX1_3.res$eigs)
boot.res <- boot.comps(dat.X,dat.X1_3,F,F,eigs=XX1_3.res$eigs)
perm.between.res <- perm.between.comps(dat.X,dat.X1_3,F,F,eigs=XX1_3.res$eigs)
perm.within.res <- perm.within.comps(dat.X,dat.X1_3,F,F,eigs=XX1_3.res$eigs)






