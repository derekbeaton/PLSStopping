rm(list=ls())
gc()

## load from SlimPosition.
  ## should eventually load from URL
source("/Volumes/JOHNNYFIVE/Professional/Software/ExPosition-Family/Code/R/Development/SlimPosition/R/gsvd.R")
source("/Volumes/JOHNNYFIVE/Professional/Software/ExPosition-Family/Code/R/Development/SlimPosition/R/tolerance.svd.R")
source("/Volumes/JOHNNYFIVE/Professional/Software/ExPosition-Family/Code/R/Development/SlimPosition/R/invert.rebuild_matrix.R")
source("/Volumes/JOHNNYFIVE/Professional/Software/ExPosition-Family/Code/R/Development/SlimPosition/R/power.rebuild_matrix.R")
source("/Volumes/JOHNNYFIVE/Professional/Software/ExPosition-Family/Code/R/Development/SlimPosition/R/isDiagonal.matrix.R")
source("/Volumes/JOHNNYFIVE/Professional/Software/ExPosition-Family/Code/R/Development/SlimPosition/R/sp.cca.R")
source("/Volumes/JOHNNYFIVE/Professional/Software/ExPosition-Family/Code/R/Development/SlimPosition/R/sp.pls.R")
source("/Volumes/JOHNNYFIVE/Professional/Software/ExPosition-Family/Code/R/Development/SlimPosition/R/sp.rrr.R")
source("/Volumes/JOHNNYFIVE/Professional/Software/ExPosition-Family/Code/R/Development/SlimPosition/R/expo.scale.R")
source("/Volumes/JOHNNYFIVE/Professional/Software/ExPosition-Family/Code/R/Development/SlimPosition/R/sp.component_plot.R")
source("/Volumes/JOHNNYFIVE/Professional/Software/ExPosition-Family/Code/R/Development/SlimPosition/R/sp.latentvar_plot.R")
source("/Volumes/JOHNNYFIVE/Professional/Software/ExPosition-Family/Code/R/Development/SlimPosition/R/sp.scree.R")
source("/Volumes/JOHNNYFIVE/Professional/Software/ExPosition-Family/Code/R/Development/SlimPosition/R/sp.pca.R")

load('/Volumes/JOHNNYFIVE/Professional/Software/ExPosition-Family/Code/R/Development/SlimPosition/data/two.table.wine.rda')


broken.stick.test <- function(eigs){
  broken.stick.comps <- (eigs/sum(eigs) * 100) > (unlist(lapply(X = 1:length(eigs), FUN = function(x, n) {return(tail((cumsum(1/x:n))/n, n = 1))}, n = length(eigs))) * 100)
  broken.stick.comps[head(which(!broken.stick.comps), n = 1):length(eigs)] <- rep(FALSE, length(head(which(!broken.stick.comps), n = 1):length(eigs)))
  return( max(which(broken.stick.comps)) )
}

kaiser <- function(eigs){
  return(eigs > mean(eigs))
}

## I still don't even know what this is exactly but I found it here:
  ## from https://stats.stackexchange.com/questions/33917/how-to-determine-significant-principal-components-using-bootstrapping-or-monte-c
  ## I think it's supposed to use the last N components to build a regression line, and anything above the line is "a keeper".
      ## it's weird and I don't entirely get it.
multi.part.reg <- function(eigs){
  
  comp.seq <- 0:(round(length(eigs)/2))
  reg.comps <- matrix(0,length(comp.seq),length(eigs))
  
  for(i in comp.seq){
    ntrail <- length(eigs)-i
    x <- seq(length(eigs)-ntrail+1, length(eigs))
    y <- log(tail(eigs, ntrail))
    #fit <- lm(y ~ x)
    #fit.pred <- predict(lm(y ~ x))
    new <- data.frame(x = 1:length(eigs))
    pred.fit <- predict(lm(y ~ x), new)
    
    reg.comps[i+1,which((log(eigs) - pred.fit) > 0)] <- 1
  }
  return(reg.comps)
}


X <- wine$objective
Y <- wine$subjective

  ## it is important to have these as references.
X.pca <- sp.pca(X, scale = "SS1", graphs=F, compact=F)
  X.broken.comps <- broken.stick.test(X.pca$d.orig^2)
  X.kaiser.comps <- kaiser(X.pca$d.orig^2)
  X.reg.comps <- multi.part.reg(X.pca$d.orig^2)
  
Y.pca <- sp.pca(Y,scale = "SS1", graphs=F, compact=F)
  Y.broken.comps <- broken.stick.test(Y.pca$d.orig^2)
  Y.kaiser.comps <- kaiser(Y.pca$d.orig^2)
  Y.reg.comps <- multi.part.reg(Y.pca$d.orig^2)

  
## the tests to do: 
  # as-is
    # broken stick
    # kaiser
    # multi-part regression
  # resampling
    # permutation between
    # permutation within
    # bootstrap
      # lower CI > mean
      # lower CI > boot mean
      # broken stick
      # kaiser?
    # split-half; each below have many variations.
      # cor fss
        # cor(fi)
        # cor(fj)
      # cor lvs
        # cor(lx)
        # cor(ly)
      # perm
      # Rv/R2 psuedo-step down
  
  
pls.res <- sp.pls(X, Y, graphs = F, compact = F)
  pls.broken.comps <- broken.stick.test(pls.res$d.orig^2)
  pls.kaiser.comps <- kaiser(pls.res$d.orig^2)
  
  
cca.res <- sp.cca(X, Y, graphs = F, compact = F)
  pls.broken.comps <- broken.stick.test(cca.res$d.orig^2)
  pls.kaiser.comps <- kaiser(cca.res$d.orig^2)
  
  
rrr.res <- sp.rrr(X, Y, graphs = F, compact = F)
  rrr.broken.comps <- broken.stick.test(rrr.res$d.orig^2)
  rrr.kaiser.comps <- kaiser(rrr.res$d.orig^2)

  
iters <- 100
rrr.sing.permb <- cca.sing.permb <- pls.sing.permb <- rrr.sing.permw <- cca.sing.permw <- pls.sing.permw <- rrr.sing.boots <- cca.sing.boots <- pls.sing.boots <- matrix(NA,iters,max( length(pls.res$d.orig),length(cca.res$d.orig),length(rrr.res$d.orig) )) 
for(i in 1:iters){
  
  boot.order <- sample(nrow(X),nrow(X),replace = T)
  
  pls.boot <- sp.pls(X[boot.order,], Y[boot.order,], graphs = F, compact = F)
    pls.sing.boots[i, 1:min( length(pls.boot$d.orig), ncol(pls.sing.boots)) ] <- pls.boot$d.orig[1:min( length(pls.boot$d.orig), ncol(pls.sing.boots)) ]
  cca.boot <- sp.cca(X[boot.order,], Y[boot.order,], graphs = F, compact = F)
    cca.sing.boots[i, 1:min( length(cca.boot$d.orig), ncol(cca.sing.boots)) ] <- cca.boot$d.orig[1:min( length(cca.boot$d.orig), ncol(cca.sing.boots)) ]  
  rrr.boot <- sp.rrr(X[boot.order,], Y[boot.order,], graphs = F, compact = F)
    rrr.sing.boots[i, 1:min( length(rrr.boot$d.orig), ncol(rrr.sing.boots)) ] <- rrr.boot$d.orig[1:min( length(rrr.boot$d.orig), ncol(rrr.sing.boots)) ]        
    
  
  X.win.perm <- apply(X,2,sample)
  pls.win <- sp.pls(X.win.perm, Y, graphs = F, compact = F)
    pls.sing.permw[i, 1:min( length(pls.win$d.orig), ncol(pls.sing.permw)) ] <- pls.win$d.orig[1:min( length(pls.win$d.orig), ncol(pls.sing.permw)) ]  
  cca.win <- sp.cca(X.win.perm, Y, graphs = F, compact = F)
    cca.sing.permw[i, 1:min( length(cca.win$d.orig), ncol(cca.sing.permw)) ] <- cca.win$d.orig[1:min( length(cca.win$d.orig), ncol(cca.sing.permw)) ]    
  rrr.win <- sp.rrr(X.win.perm, Y, graphs = F, compact = F)    
    rrr.sing.permw[i, 1:min( length(rrr.win$d.orig), ncol(rrr.sing.permw)) ] <- rrr.win$d.orig[1:min( length(rrr.win$d.orig), ncol(rrr.sing.permw)) ]    
  
  
  X.bw.perm <- X[sample(nrow(X),nrow(X),replace=F),]
  pls.bw <- sp.pls(X.bw.perm, Y, graphs = F, compact = F)
    pls.sing.permb[i, 1:min( length(pls.bw$d.orig), ncol(pls.sing.permb)) ] <- pls.bw$d.orig[1:min( length(pls.bw$d.orig), ncol(pls.sing.permb)) ]    
  cca.bw <- sp.cca(X.bw.perm, Y, graphs = F, compact = F)
    cca.sing.permb[i, 1:min( length(cca.bw$d.orig), ncol(cca.sing.permb)) ] <- cca.bw$d.orig[1:min( length(cca.bw$d.orig), ncol(cca.sing.permb)) ]      
  rrr.bw <- sp.rrr(X.bw.perm, Y, graphs = F, compact = F)      
    rrr.sing.permb[i, 1:min( length(rrr.bw$d.orig), ncol(rrr.sing.permb)) ] <- rrr.bw$d.orig[1:min( length(rrr.bw$d.orig), ncol(rrr.sing.permb)) ]      
  
}  



pls.broken.boot <- apply(pls.sing.boots,1,broken.stick.test)
cca.broken.boot <- apply(cca.sing.boots,1,broken.stick.test)
rrr.broken.boot <- apply(rrr.sing.boots,1,broken.stick.test)


pls.boot_lower.ci_above.boot.mean <- apply(pls.sing.boots,2,function(x){ sort(x)[iters*.025] }) > mean(pls.res$d.orig)
cca.boot_lower.ci_above.boot.mean <- apply(cca.sing.boots,2,function(x){ sort(x)[iters*.025] }) > mean(cca.res$d.orig)
rrr.boot_lower.ci_above.boot.mean <- apply(rrr.sing.boots,2,function(x){ sort(x)[iters*.025] }) > mean(rrr.res$d.orig)


pls.boot_lower.ci_above.fixed.mean <- apply(pls.sing.boots,2,function(x){ sort(x)[iters*.025] }) > mean(pls.sing.boots)
cca.boot_lower.ci_above.fixed.mean <- apply(cca.sing.boots,2,function(x){ sort(x)[iters*.025] }) > mean(cca.sing.boots)
rrr.boot_lower.ci_above.fixed.mean <- apply(rrr.sing.boots,2,function(x){ sort(x)[iters*.025] }) > mean(rrr.sing.boots)



