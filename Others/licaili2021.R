######### Fixed Effects Estimation #########
library(Matrix)
library(MASS)
library(glmnet)
library(scalreg)
library(tictoc)

CV_fixed_effects <- function(xlist, ylist, zlist, grp, a.seq=seq(0.5, 10, 0.5), inf.coord,
                             standardize, intercept,
                             std_response,
                             train_grps, mi){
    n <- length(unique(grp))
    N <- length(unlist(ylist))
    n.tr <- round(0.75*n)#four fold to train and one fold to test
    best.err <- sum(unlist(ylist)^2)
    best.a <- a.seq[1]
    
    for(ia in 1:length(a.seq)){
        est.re <- Fixed_effects_estimation(xlist[1:n.tr], ylist[1:n.tr],
                                           zlist[1:n.tr], grp = grp[1:sum(mi[train_grps][1:n.tr])],
                                           a = a.seq[ia], inf.coord = inf.coord,
                                           standardize = standardize, intercept = intercept,
                                           std_response = std_response,
                                           train_grps = train_groups[1:n.tr],
                                           mi = mi[train_grps][1:n.tr])
        pred.err <- sum((unlist(ylist[-(1:n.tr)]) - drop(do.call(rbind, xlist[-(1:n.tr)]) %*% est.re$beta.hat.lasso))^2)
        if(pred.err < best.err){
            best.err <- pred.err
            best.a <- a.seq[ia]
        }
        
    }
    list(best.a = best.a)
}

Fixed_effects_estimation <- function(xlist, ylist, zlist, grp, a, inf.coord = NULL,
                                     standardize, intercept, std_response,
                                     train_grps, mi){
    N <- length(unlist(ylist))
    n <- length(unique(grp))
    p <- ncol(xlist[[1]])
    
    
    x <- do.call(rbind, xlist)
    y <- unlist(ylist)
    
    #mixed effect model fitting
    x.a <- x
    y.a <- y
    tr.a <- 0
    
    for (i in 1:n){
        cur.mem <- which(grp == train_grps[i])
        sigmai <- a*tcrossprod(zlist[[i]]) + diag(rep(1, nrow(zlist[[i]])))
        sigmai.svd <- svd(sigmai)
        Sig.a.inv.half <- sigmai.svd$u %*% diag(1 / sqrt(sigmai.svd$d)) %*%
            t(sigmai.svd$u)
        x.a[cur.mem,] <-  Sig.a.inv.half %*% x[cur.mem, ] 
        y.a[cur.mem] <-  Sig.a.inv.half %*% as.matrix(y[cur.mem])
        tr.a <- tr.a + sum(1 / sigmai.svd$d) # compute the effective sample size
    }
    sig.init <- scalreg(x.a, y.a)$hsigma # scaled-lasso to get tuning parameter
    lasso.fit <- glmnet(x.a, y.a,
                        lambda = sig.init * sqrt(2 * log(p)/N),
                        standardize = standardize, intercept = intercept, standardize.response = std_response)
    beta.hat.lasso <- as.vector(coef.glmnet(lasso.fit))[-1]
    beta.hat.all <- as.vector(coef(lasso.fit))
    ### debiased Lasso
    beta.hat.db.sd <- rep(NA,length=length(inf.coord))
    beta.hat.db <- rep(NA,length=length(inf.coord))
    
    if(is.null(inf.coord)) {
        return(list(beta.hat.lasso = beta.hat.lasso,
                    beta.hat.db = beta.hat.db, tr.a = tr.a))
    }
    
    res <- y.a - x.a %*% beta.hat.lasso
    
    for(j in 1:length(inf.coord)){
        col.j <- j
        sig.x <- scalreg(x.a[, -col.j], x.a[, col.j])$hsigma
        kappa.hat.mlm <- glmnet(x.a[, -col.j], x.a[, col.j],
                                lambda = sig.x * sqrt(2*log(p)/N),
                                standardize = standardize, intercept = intercept, standardize.response = std_response)
        gam.j <- as.vector(coef.glmnet(kappa.hat.mlm))[-1]
        wj.mlm <- x.a[, j] - x.a[, -j] %*% gam.j
        beta.hat.db[j] <- beta.hat.lasso[j] + sum(wj.mlm * res)/sum(wj.mlm * x.a[,j])
        
        num <- 0
        for(i in 1:n){
            cur.mem <- which(grp == train_grps[i])
            num <- num + sum(wj.mlm[cur.mem] * res[cur.mem])^2
        }
        
        beta.hat.db.sd[j] <- sqrt(num)/sum(wj.mlm * as.vector(x.a[,j]))
    }
    
    list(beta.hat.lasso = beta.hat.lasso, beta.hat.db = beta.hat.db,
         beta.hat.db.sd = beta.hat.db.sd, tr.a = tr.a,
         beta.hat.all = beta.hat.all,
         xtrain_mean = colMeans(do.call("rbind", xlist)),
         xtrain_sd = apply(do.call("rbind", xlist), 2, sd))
}

