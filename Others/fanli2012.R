############### Fan and Li (2012) ###############
######## For Simulation Study ########
######## Function to calculate the number of times the correct fixed effects are selected ########


calc_CF <- function(ylist, xlist, zlist,
                    capM, beta.true, xtrain_sds,
                    standardize, intercept, std_response,
                    newylist, newxlist, newzlist) {
    library(Matrix)
    library(glmnet)
    zmat <- as.matrix(bdiag(zlist)) # N x nq
    xmat <- do.call(rbind, xlist) # N x p
    y <- matrix(unlist(ylist), ncol = 1) # N x 1
    
    Pz <- chol2inv(chol(diag(length(y)) + tcrossprod(zmat %*% capM, zmat)))
    Pz_svd <- svd(Pz)
    Pz_half <- Pz_svd$u %*% tcrossprod(diag(sqrt(Pz_svd$d)), Pz_svd$u)
    y_scaled <- Pz_half %*% y
    xmat_scaled <- Pz_half %*% xmat
    
    cv.out <- cv.glmnet(xmat_scaled, y_scaled, alpha = 1,
                        intercept = intercept, standardize = standardize, standardize.response = std_response)
    
    # lambda value chosen by 10-fold CV
    beta.hat.lmin <- coef(cv.out, s = cv.out$lambda.min)[-1]
    
    rmse.lmin <- sqrt(mean((beta.hat.lmin - beta.true) ^ 2))
    beta.hat.lmin.orig.scale <- beta.hat.lmin / xtrain_sds
    rmse.lmin.scale <- sqrt(mean((beta.hat.lmin.orig.scale - beta.true) ^ 2))
    
    true_indx.lmin <- which(beta.true != 0)
    est_indx.lmin <- which(beta.hat.lmin != 0)
    cf.lmin <- length(intersect(est_indx.lmin, true_indx.lmin)) 
    tpr.lmin <- cf.lmin / length(true_indx.lmin)
    fpr.lmin <- length(intersect(est_indx.lmin, which(beta.true == 0))) / length(which(beta.true == 0))

    ## Prediction
    
    y.test <- as.matrix(unlist(newylist), ncol = 1)
    x.test <- do.call(rbind, newxlist)
    y.hat.lmin <- x.test %*% beta.hat.lmin
    
    mspe.lmin <- mean((y.hat.lmin - y.test) ^ 2)
    
    list(beta.hat.lmin = beta.hat.lmin, cf.lmin = cf.lmin,
         tpr.lmin = tpr.lmin, fpr.lmin = fpr.lmin,
         rmse.lmin = rmse.lmin,
         y.hat.lmin = y.hat.lmin, ytrue.test = y.test, mspe.lmin = mspe.lmin,
         beta.hat.lmin.orig.scale = beta.hat.lmin.orig.scale,
         rmse.lmin.scale = rmse.lmin.scale,
         xmat_scaled = xmat_scaled, y_scaled = y_scaled,
         Pz = Pz, Pz_half = Pz_half)
}


######## Function to calculate coverage of confidence intervals for beta ########
########(following Debiased estimator by LiCaiLi2021) ########

covg_fanli <- function(lasso_est, xtrain_sds, beta.true, n, m, intercept, standardize, std_response) {
    library(Matrix)
    library(glmnet)
    p <- length(beta.true)
    y_scaled <- lasso_est$y_scaled
    xmat_scaled <- lasso_est$xmat_scaled
    
    beta.hat.db.sd.lmin <- rep(NA, p)
    beta.hat.db.lmin <- rep(NA, p)
    
    res.lmin <- y_scaled - xmat_scaled %*% lasso_est$beta.hat.lmin
    
    for(j in 1:p){
        kappa.hat.mlm <- cv.glmnet(as.matrix(xmat_scaled[, -j], ncol = p-1), as.vector(xmat_scaled[, j]),
                                   alpha = 1, nfolds = 10,
                                   intercept = intercept, standardize = standardize, standardize.response = std_response)
        gam.j.lmin <- coef(kappa.hat.mlm, s = kappa.hat.mlm$lambda.min)[-1]
        #gam.j.l1se <- coef(kappa.hat.mlm, s = kappa.hat.mlm$lambda.1se)[-1]
        wj.mlm.lmin <- xmat_scaled[, j] - xmat_scaled[, -j] %*% gam.j.lmin
        #wj.mlm.l1se <- xmat_scaled[, j] - xmat_scaled[, -j] %*% gam.j.l1se
        beta.hat.db.lmin[j] <- lasso_est$beta.hat.lmin[j] + sum(wj.mlm.lmin * res.lmin)/sum(wj.mlm.lmin * xmat_scaled[,j])
        #beta.hat.db.l1se[j] <- lasso_est$beta.hat.l1se[j] + sum(wj.mlm.l1se * res.l1se)/sum(wj.mlm.l1se * xmat_scaled[,j])
        
        num.lmin <- 0
        for(i in 1:n){
            cur.mem <- ((i-1)*m+1):(i*m)
            num.lmin <- num.lmin + sum(wj.mlm.lmin[cur.mem] * res.lmin[cur.mem])^2
        }
        
        beta.hat.db.sd.lmin[j] <- sqrt(num.lmin)/sum(wj.mlm.lmin * as.vector(xmat_scaled[,j]))
    }
    
    beta.hat.db.lmin.scaled <- beta.hat.db.lmin / xtrain_sds
    beta.hat.db.sd.lmin.scaled <- beta.hat.db.sd.lmin / xtrain_sds
    z_95 <- qnorm(0.975); z_90 <- qnorm(0.95); z_80 <- qnorm(0.90)
    
    covgs_95.lmin <- rep(NA, p)
    covgs_90.lmin <- rep(NA, p)
    covgs_80.lmin <- rep(NA, p)
    
    for(pp in 1:p) {
        covgs_95.lmin[pp] <- as.numeric(beta.hat.db.lmin[pp] - z_95 * beta.hat.db.sd.lmin[pp] <= beta.true[pp] & beta.true[pp] <= beta.hat.db.lmin[pp] + z_95 * beta.hat.db.sd.lmin[pp])
        covgs_90.lmin[pp] <- as.numeric(beta.hat.db.lmin[pp] - z_90 * beta.hat.db.sd.lmin[pp] <= beta.true[pp] & beta.true[pp] <= beta.hat.db.lmin[pp] + z_90 * beta.hat.db.sd.lmin[pp])
        covgs_80.lmin[pp] <- as.numeric(beta.hat.db.lmin[pp] - z_80 * beta.hat.db.sd.lmin[pp] <= beta.true[pp] & beta.true[pp] <= beta.hat.db.lmin[pp] + z_80 * beta.hat.db.sd.lmin[pp])
    }
    
    covgs_95.lmin.scaled <- rep(NA, p)
    covgs_90.lmin.scaled <- rep(NA, p)
    covgs_80.lmin.scaled <- rep(NA, p)
    
    for(pp in 1:p) {
        covgs_95.lmin.scaled[pp] <- as.numeric(beta.hat.db.lmin.scaled[pp] - z_95 * beta.hat.db.sd.lmin.scaled[pp] <= beta.true[pp] & beta.true[pp] <= beta.hat.db.lmin.scaled[pp] + z_95 * beta.hat.db.sd.lmin.scaled[pp])
        covgs_90.lmin.scaled[pp] <- as.numeric(beta.hat.db.lmin.scaled[pp] - z_90 * beta.hat.db.sd.lmin.scaled[pp] <= beta.true[pp] & beta.true[pp] <= beta.hat.db.lmin.scaled[pp] + z_90 * beta.hat.db.sd.lmin.scaled[pp])
        covgs_80.lmin.scaled[pp] <- as.numeric(beta.hat.db.lmin.scaled[pp] - z_80 * beta.hat.db.sd.lmin.scaled[pp] <= beta.true[pp] & beta.true[pp] <= beta.hat.db.lmin.scaled[pp] + z_80 * beta.hat.db.sd.lmin.scaled[pp])
    }
    
    list(beta.hat.db.lmin = beta.hat.db.lmin,
         beta.hat.db.sd.lmin = beta.hat.db.sd.lmin,
         covgs_95.lmin = covgs_95.lmin, covgs_90.lmin = covgs_90.lmin, covgs_80.lmin = covgs_80.lmin,
         beta.hat.db.lmin.scaled = beta.hat.db.lmin.scaled,
         beta.hat.db.sd.lmin.scaled = beta.hat.db.sd.lmin.scaled,
         covgs_95.lmin.scaled = covgs_95.lmin.scaled, covgs_90.lmin.scaled = covgs_90.lmin.scaled, covgs_80.lmin.scaled = covgs_80.lmin.scaled)
}


########### For Data Analysis ############
fanli_fes <- function(ylist, xlist, zlist, capM) {
    library(Matrix)
    library(glmnet)
    zmat <- as.matrix(bdiag(zlist)) # N x nq
    xmat <- do.call(rbind, xlist) # N x p
    y_train <- matrix(unlist(ylist), ncol = 1) # N x 1
    
    Pz <- chol2inv(chol(diag(length(y_train)) + tcrossprod(zmat %*% capM, zmat)))
    Pz_svd <- svd(Pz)
    Pz_half <- Pz_svd$u %*% tcrossprod(diag(sqrt(Pz_svd$d)), Pz_svd$u)
    y_scaled <- Pz_half %*% y_train
    xmat_scaled <- Pz_half %*% xmat
    
    cv.result <- cv.glmnet(xmat_scaled, y_scaled, alpha = 1)
    
    # lambda value chosen by 10-fold CV
    beta.hat.lmin <- as.vector(coef(cv.result, s = cv.result$lambda.min))
    
    ## Predicted values of y_scaled (training data)
    pred_y_scaled_lmin <- cbind(rep(1, nrow(xmat_scaled)), xmat_scaled) %*% beta.hat.lmin
    
    list(beta.hat.lmin = beta.hat.lmin, 
         xmat_scaled = xmat_scaled, y_scaled = y_scaled,
         Pz = Pz, Pz_half = Pz_half,
         cv.result = cv.result,
         pred_y_scaled_lmin = pred_y_scaled_lmin)
}



fanli_pred <- function(ylist, xlist, zlist, capM, standardize, intercept,
                         newylist, newxlist, newzlist) {
    library(Matrix)
    library(glmnet)
    library(MASS)
    zmat <- as.matrix(bdiag(zlist)) # N x nq
    xmat <- do.call(rbind, xlist) # N x p
    y_train <- matrix(unlist(ylist), ncol = 1) # N x 1
    
    Pz <- chol2inv(chol(diag(length(y_train)) + tcrossprod(zmat %*% capM, zmat)))
    Pz_svd <- svd(Pz)
    Pz_half <- Pz_svd$u %*% tcrossprod(diag(sqrt(Pz_svd$d)), Pz_svd$u)
    y_scaled <- Pz_half %*% y_train
    xmat_scaled <- Pz_half %*% xmat
    
    cv.result <- cv.glmnet(xmat_scaled, y_scaled, alpha = 1, standardize = standardize, intercept = intercept)
    
    # lambda value chosen by 10-fold CV
    beta.hat.lmin <- as.vector(coef(cv.result, s = cv.result$lambda.min))
    
    ## Predicted values of y_scaled (training data)
    pred_y_scaled_lmin <- cbind(rep(1, nrow(xmat_scaled)), xmat_scaled) %*% beta.hat.lmin
    
    ## Out-of-sample Prediction
    y.test <- as.matrix(unlist(newylist), ncol = 1)
    x.test <- cbind(rep(1, nrow(y.test)), do.call(rbind, newxlist))
    y.hat.lmin <- drop(x.test %*% beta.hat.lmin)
    mspe.lmin <- mean((y.hat.lmin - drop(y.test)) ^ 2)
    
    list(beta.hat.lmin = beta.hat.lmin, 
         xmat_scaled = xmat_scaled, y_scaled = y_scaled,
         Pz = Pz, Pz_half = Pz_half,
         cv.result = cv.result,
         pred_y_scaled_lmin = pred_y_scaled_lmin,
         y.hat.lmin = y.hat.lmin,
         y_true = drop(y.test), x.test = x.test,
         mspe.lmin = mspe.lmin)
}


