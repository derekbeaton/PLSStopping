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

broken.boot <- apply(boot.res$boot.eigs,1,broken.stick)
colSums(t(broken.boot))

boot_lower.ci_above.boot.mean <- apply(boot.res$boot.eigs,2,function(x){ sort(x)[nrow(boot.res$boot.eigs)*.025] }) > mean(res$eigs)
boot_lower.ci_above.fixed.mean <- apply(boot.res$boot.eigs,2,function(x){ sort(x)[nrow(boot.res$boot.eigs)*.025] }) > mean(boot.res$boot.eigs)



### NOTE: The deflation version is trickier than I thought.
	## I had originally assumed I could just regress $lx or $ly from each column in its respective matrix
	## I was wrong, and that does not actually account for what I want to take out
	## So I have to use the deflation step from PLSR but for each X & Y
		## This is going to be insanely expensive.
	
## now check the regressed version then perform tests on it
X.new <- X
Y.new <- Y
perm.within_deflate <- perm.between_deflate <- boot.eigs_deflate <- matrix(NA,100,length(res$d))
boot.eigs_deflate[,1] <- boot.res$boot.eigs[,1]
perm.between_deflate[,1] <- perm.between.res$perm.eigs[,1]
perm.within_deflate[,1] <- perm.within.res$perm.eigs[,1]

for(j in 1:(length(res$d)-1)){
		
	X.new <- X.new - (res$lx[,j] %o% res$u[,j])
	Y.new <- Y.new - (res$ly[,j] %o% res$v[,j])		
	# reg.res <- simple.pls(X.new,Y.new)
	# print(reg.res$d)
	# pause()
	
	boot.eigs_deflate[,(j+1)] <- boot.comps(X,Y,F,F,eigs=res$eigs)$boot.eigs[,1]
	perm.between_deflate[,(j+1)] <- perm.between.comps(X,Y,F,F,eigs=res$eigs)$perm.eigs[,1]
	perm.within_deflate[,(j+1)] <- perm.within.comps(X,Y,F,F,eigs=res$eigs)$perm.eigs[,1]		

	print(j)
	
}


	## bootstrap is behaving much more differently -- converging to the middle.
colSums(boot.eigs_deflate < matrix(rowMeans(boot.eigs_deflate,na.rm=T),100,ncol(boot.eigs_deflate),byrow=F)) / 100
(res$eigs - colMeans(boot.eigs_deflate)) / apply(boot.eigs_deflate,2,sd)
(colMeans(boot.eigs_deflate) / apply(boot.eigs_deflate,2,sd))



	## these are not like the permutation ones above -- the above ones really more so reflect the idea of an "omnibus"
		## probably because the first component is so big.
colSums(perm.between_deflate > matrix(res$eigs,100,ncol(perm.between_deflate),byrow=T)) / 100
colSums(perm.within_deflate > matrix(res$eigs,100,ncol(perm.within_deflate),byrow=T)) / 100

