simple.pls <- function(X,Y,k=0,tol=.Machine$double.eps*1000){

	## X and Y better be preprocessed!
	R <- t(X) %*% Y

	if(k==0){
		k <- min(dim(R))
	}

	svd.res <- svd(R,nu=k,nv=k)
	comps.to.keep <- which(!(svd.res$d^2 < tol))
		
	svd.res$u <- as.matrix(svd.res$u[,comps.to.keep])
	svd.res$v <- as.matrix(svd.res$v[,comps.to.keep])
	svd.res$d <- svd.res$d[comps.to.keep]		
	
	svd.res$eigs <- svd.res$d^2
	svd.res$tau <- svd.res$eigs / sum(svd.res$eigs)
	if(length(svd.res$d)>1){
		svd.res$fi <- svd.res$u %*% diag(svd.res$d)
		svd.res$fj <- svd.res$v %*% diag(svd.res$d)	
	}
	else{
		svd.res$fi <- svd.res$u * svd.res$d
		svd.res$fj <- svd.res$v * svd.res$d
	}
	svd.res$lx <- X %*% svd.res$u
	svd.res$ly <- Y %*% svd.res$v

	return(svd.res)
	
}

broken.stick <- function(eigs){
	
	exp.var <- eigs/sum(eigs) * 100
	broken.stick.distribution <- unlist(lapply(X = 1:length(eigs),FUN = function(x, n) { return(tail((cumsum(1/x:n))/n, n = 1))}, n = length(eigs))) * 100
	broken.stick.comps <- exp.var > broken.stick.distribution
	broken.stick.comps[head(which(!broken.stick.comps), n = 1):length(eigs)] <- rep(FALSE,length(head(which(!broken.stick.comps), n = 1):length(eigs)))
	
	return(broken.stick.comps)
}

kaiser <- function(eigs){
	return(eigs > mean(eigs))
}

boot.comps <- function(X,Y,center=T,scale=T,eigs=NA,k=0,iters=100){
	
	if(is.na(eigs)){
		boot.eigs <-	 matrix(NA,iters,max(c(ncol(X),ncol(Y))))
	}else{
		boot.eigs <-	 matrix(NA,iters,length(eigs))
	}
	
	pb <- txtProgressBar(min = 0, max = iters, style = 3)
	for(i in 1:iters){
		boot.samps <- sample(nrow(X),nrow(X),replace=T)
		boot.X <- expo.scale(X[boot.samps,],center=center,scale=scale)
		boot.Y <- expo.scale(Y[boot.samps,],center=center,scale=scale)	
		boot.res_eigs <- simple.pls(boot.X,boot.Y,k=k)$eigs
		boot.eigs[i,1:min(ncol(boot.eigs),length(boot.res_eigs))] <- boot.res_eigs
	 setTxtProgressBar(pb, i)	
	}
		## rudimentary p values to estimate how often each eigenvalue is above its respective mean
	boot.ps <- colSums(boot.eigs < matrix(rowMeans(boot.eigs,na.rm=T),iters,ncol(boot.eigs),byrow=F)) / iters
		## we can do eigenvalue bootstrap ratio tests, too...
	if(!is.na(eigs)){
		boot.bsrs <- (eigs - colMeans(boot.eigs)) / apply(boot.eigs,2,sd)
		return( list(boot.eigs= boot.eigs,boot.ps=boot.ps, boot.bsrs_from.eigs= boot.bsrs,boot.bsrs=colMeans(boot.eigs) / apply(boot.eigs,2,sd)) ) 		
	}
	return( list(boot.eigs= boot.eigs,boot.ps=boot.ps,boot.bsrs=colMeans(boot.eigs) / apply(boot.eigs,2,sd)) )
}

perm.between.comps <- function(X,Y,center=T,scale=T,eigs=NA,k=0,iters=100){
	
	if(is.na(eigs)){
		stop("We need eigenvalues to test against!")
	}
	
	perm.eigs <- matrix(NA,iters,length(eigs))
	pb <- txtProgressBar(min = 0, max = iters, style = 3)	
	for(i in 1:iters){
		perm.samps <- sample(nrow(X),nrow(X),replace=F)
		perm.X <- expo.scale(X[perm.samps,],center=center,scale=scale)
		Y <- expo.scale(Y,center=center,scale=scale)	
		perm.res_eigs <- simple.pls(perm.X,Y,k=k)$eigs
		perm.eigs[i,1:min(ncol(perm.eigs),length(perm.res_eigs))] <- perm.res_eigs
		setTxtProgressBar(pb, i)
	}
    perm.ps <- colSums(perm.eigs > matrix(eigs,iters,ncol(perm.eigs),byrow=T)) / iters	
	return(list(perm.eigs = perm.eigs,perm.ps=perm.ps))
}

perm.within.comps <- function(X,Y,center=center,scale=scale,eigs=NA,k=0,iters=100){
	
	if(is.na(eigs)){
		stop("We need eigenvalues to test against!")
	}
		
	perm.eigs <- matrix(NA,iters,length(eigs))		
	pb <- txtProgressBar(min = 0, max = iters, style = 3)			
	for(i in 1:iters){
		perm.X <- expo.scale(apply(X,2,sample),center=center,scale=scale)
		Y <- expo.scale(Y,center=center,scale=scale)
		perm.res_eigs <- simple.pls(perm.X,Y,k=k)$eigs
		perm.eigs[i,1:min(ncol(perm.eigs),length(perm.res_eigs))] <- perm.res_eigs		
		setTxtProgressBar(pb, i)		
	}
    perm.ps <- colSums(perm.eigs > matrix(eigs,iters,ncol(perm.eigs),byrow=T)) / iters	
	return(list(perm.eigs = perm.eigs,perm.ps=perm.ps))
}


## Split-half becomes a bit tricky. What exactly are we trying to get from it?
	## maximum replication of e.g., cor(u1, u2) or cor(lx1, lx2)?
	## minimum value between cor(lx1,ly2) - observed SVs (i.e., each split is highly similar)?
	## also, should data be normed before or within?
	## 		and, should we create null distributions within the splits?
## Not easy for now mostly because there are many, many possibilities.
	
# split.half_simple <- function(X,Y,center=center,scale=scale,k=0){
	
	# sh1 <- sort(sample(nrow(X), floor(nrow(X)/2), replace=F ))
	# sh1.X <- expo.scale(X[sh1,],center=center,scale=scale)
	# sh1.Y <- expo.scale(Y[sh1,],center=center,scale=scale)	
	# sh1_res <- simple.pls(sh1.X,sh1.Y,k=k)
	
	
	# sh2 <- setdiff(1:nrow(X),sh1)
	# sh2.X <- expo.scale(X[sh2,],center=center,scale=scale)
	# sh2.Y <- expo.scale(Y[sh2,],center=center,scale=scale)	
	# sh2_res <- simple.pls(sh2.X,sh2.Y,k=k)
		
	
	
# }

# split.half_withNULL <- function(X,Y,center=center,scale=scale,k=0){
	
	# sh1 <- sort(sample(nrow(X), floor(nrow(X)/2), replace=F ))
	# sh1.X <- expo.scale(X[sh1,],center=center,scale=scale)
	# sh1.Y <- expo.scale(Y[sh1,],center=center,scale=scale)	
	# sh1_res <- simple.pls(sh1.X,sh1.Y,k=k)
	
	
	# sh2 <- setdiff(1:nrow(X),sh1)
	# sh2.X <- expo.scale(X[sh2,],center=center,scale=scale)
	# sh2.Y <- expo.scale(Y[sh2,],center=center,scale=scale)	
	# sh2_res <- simple.pls(sh2.X,sh2.Y,k=k)
		
	
	
# }
