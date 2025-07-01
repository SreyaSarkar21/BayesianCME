# ylist: List of response vectors where the i-th component is a m_i observations for the i-th subject.
# xlist: List of fixed effect covariates where the i-th component contains m_i observations for the i-th subject.
# zlist: List of random effect covariates where the i-th component contains m_i observations for the i-th subject.
# a0: shape hyperparameter for inverse-gamma prior for error variance.
# b0: scale hyperparameter for inverse-gamma prior for error variance.
# gammaVar0: variance hyperparameter for Gaussian prior for the vectorized compressed covariance parameter.
# S: k_1 x q random projection matrix with entries iid from Gaussian distribution with mean 0 and variance 1/k_1.
# R: k_2 x q random projection matrix with entries iid from Gaussian distribution with mean 0 and variance 1/k_2.
# niter: number of MCMC iterations
# nburn: number of burn-in samples
# nthin: thinning size specifying every nthin-th sample is retained.

sampleCME <- function (ylist, xlist, zlist, a0, b0, gammaVar0,
                       R, S, niter, nburn, nthin) {
    
    p <- ncol(xlist[[1]]) ## number of fixed effects
    q <- ncol(zlist[[1]]) ## number of random effects
    n <- length(xlist) ## number of subjects
    N <- sum(unlist(lapply(ylist, length))) ## total number of observations
    k1 <- nrow(S); k2 <- nrow(R) ## compression dimensions
    
    ################## For storing final samples #####################
    errVarSamp <- rep(0.0, (niter - nburn) / nthin)
    betaSamp <- matrix(0.0, (niter - nburn) / nthin, p)
    gammaSamp <- matrix(NA, (niter - nburn) / nthin, k1 * k2)
    lambda2Samp <- matrix(0.0, (niter - nburn) / nthin, p) ## local shrinkage
    nuSamp <- matrix(0.0, (niter - nburn) / nthin, p)
    delta2Samp <- rep(0.0, (niter - nburn) / nthin) ## global shrinkage
    xiSamp <- rep(0.0, (niter - nburn) / nthin)
    
    Sigma_gamma <- gammaVar0 * diag(k1 * k2) ## prior covariance matrix of compressed covariance parameter
    
    ####################### Initialization ######################
    beta <- numeric(length = p)
    errVar <- 1 / rgamma(1, shape = a0, rate = b0)
    Gamma <- matrix(rnorm(k1 * k2, 0, sd = sqrt(gammaVar0)), k1, k2)
    nu <- 1 / rgamma(p, shape = 0.5, rate = 1)
    xi <- 1 / rgamma(1, shape = 0.5, rate = 1)
    delta2 <- 1 / rgamma(1, shape = 0.5, rate = 1/xi)
    lambda2 <- 1 / rgamma(p, shape = 0.5, rate = 1/nu)
    
    #### Compress Zi's ####
    zlist <- lapply(zlist, function(z) tcrossprod(z, S))
    
    cts <- 0
    startTime <- proc.time()
    for(its in 1:niter) {
        if(its %% 1000 == 0) cat("iteration: ", its, "\n")
        cycle1Samp <- sampleCycle_1(ylist, xlist, zlist,
                                    beta, errVar, Gamma, Sigma_gamma, R, S)
        Gamma <- matrix(cycle1Samp$gammaSamp, k1, k2)
        
        cycle2Samp <- sampleCycle_2(ylist, xlist, zlist,
                                    beta, errVar, a0, b0, lambda2, delta2, nu, xi,
                                    cycle1Samp$gammaSamp, R, S)
        beta <- cycle2Samp$betaSamp
        errVar <- cycle2Samp$errVarSamp
        lambda2 <- cycle2Samp$lambda2Samp
        delta2 <- cycle2Samp$delta2Samp
        nu <- cycle2Samp$nuSamp
        xi <- cycle2Samp$xiSamp
        
        if (its > nburn & its %% nthin == 0) {
            cts <- cts + 1
            errVarSamp[cts] <- cycle2Samp$errVarSamp
            betaSamp[cts, ] <- cycle2Samp$betaSamp
            gammaSamp[cts, ] <- cycle1Samp$gammaSamp
            lambda2Samp[cts, ] <- cycle2Samp$lambda2Samp
            nuSamp[cts, ] <- cycle2Samp$nuSamp
            delta2Samp[cts] <- cycle2Samp$delta2Samp
            xiSamp[cts] <- cycle2Samp$xiSamp
        }
    }
    endTime <- proc.time()
    
    list(
        errVarSamp = errVarSamp, betaSamp = betaSamp, gammaSamp = gammaSamp,
        lambda2Samp = lambda2Samp, delta2Samp = delta2Samp,
        sampler_time = endTime - startTime
    )
}


sampleCycle_1 <- function(ylist, xlist, zlist, beta, errVar, Gamma, Sigma_gamma, R, S) {
    n <- length(ylist)
    p <- ncol(xlist[[1]])
    k1 <- nrow(S)
    k2 <- nrow(R)
    
    ## Imputation Step: Drawing the unobserved random effects d_i's from their full-conditional distribution
    ranSamplist <- vector("list", n)
    for (gg in 1:n) {
        zstGR <- zlist[[gg]] %*% Gamma %*% R
        Vyd <- errVar * tcrossprod(zstGR, R)
        VyyInv <- chol2inv(chol(errVar * tcrossprod(zstGR) + errVar * diag(1, nrow(zlist[[gg]]))))
        VydTvyyInv <- crossprod(Vyd, VyyInv)
        postRanVar <- errVar * tcrossprod(R) - drop(VydTvyyInv %*% Vyd)
        postRanMean <- drop(VydTvyyInv %*% (ylist[[gg]] - xlist[[gg]] %*% beta))
        ranSamplist[[gg]] <- postRanMean + drop(crossprod(chol(postRanVar), rnorm(length(postRanMean)))) #+ diag(1e-6, k2)
    }
    
    ## Sampling of gamma, vec of k_1 x k_2 Gamma matrix, the compressed covariance parameter ##
    residi <- list(); zchecki <- list()
    for(gg in 1:n) {
        residi[[gg]] <- ylist[[gg]] - xlist[[gg]] %*% beta
        zchecki[[gg]] <- kronecker(matrix(ranSamplist[[gg]], nrow = 1), zlist[[gg]]) ## di^{\T} \otimes Zi
    }
    
    resid <- unlist(residi)
    zcheck <- do.call("rbind", zchecki)
    postVar_gamma <- chol2inv(chol(crossprod(zcheck) / errVar +
                                       chol2inv(chol(Sigma_gamma))))
    postMean_gamma <- drop(postVar_gamma %*% crossprod(zcheck, resid) / errVar)
    
    gammaSamp <- postMean_gamma + drop(crossprod(chol(postVar_gamma), rnorm(length(postMean_gamma))))
    
    list(gammaSamp = gammaSamp, ranSamplist = ranSamplist)
    
}

sampleCycle_2 <- function(ylist, xlist, zlist,
                          beta, errVar, a0, b0, lambda2, delta2, nu, xi,
                          gammaSamp, R, S) {
    
    n <- length(ylist)
    N <- sum(unlist(lapply(ylist, length)))
    p <- ncol(xlist[[1]])
    k1 <- nrow(S)
    k2 <- nrow(R)
    
    gammaMat <- matrix(gammaSamp, k1, k2)
    
    ## Marginalization ##
    VarErri <- lapply(zlist,
                      function(zi) {tcrossprod(zi %*%
                                                   gammaMat %*% R) + diag(1, nrow(zi))})
    VarErri_svd <- lapply(VarErri, svd)
    VarErri_inv_half <- vector("list", n)
    for(i in 1:n) {
        if(nrow(zlist[[i]]) == 1) {
            VarErri_inv_half[[i]] <- (VarErri_svd[[i]]$u)*(VarErri_svd[[i]]$u)/sqrt(VarErri_svd[[i]]$d)
        } else {
            VarErri_inv_half[[i]] <- VarErri_svd[[i]]$u %*% tcrossprod(diag(1 / sqrt(VarErri_svd[[i]]$d)), VarErri_svd[[i]]$u)
        }
    }
    scaled_xi <- lapply(1:n, function(ii) {VarErri_inv_half[[ii]] %*% xlist[[ii]]})
    scaled_yi <- lapply(1:n, function(ii) {VarErri_inv_half[[ii]] %*% ylist[[ii]]})
    xmat <- do.call(rbind, scaled_xi)
    yvec <- unlist(scaled_yi)
    
    vinv_beta_post <- crossprod(xmat) + diag(1 / (delta2 * lambda2))
    covBeta <- errVar * chol2inv(chol(vinv_beta_post))
    xty <- crossprod(xmat, yvec)
    muBeta <- drop(covBeta %*% xty) / errVar
    beta <- muBeta + drop(crossprod(chol(covBeta), rnorm(p))) 
    
    resids <- yvec - drop(xmat %*% beta)
    shp_err <- a0 + 0.5 * (N + p) 
    scl_err <- b0 + 0.5 * (sum(resids^2) + sum((beta * beta) / (delta2 * lambda2))) 
    errVar <- 1 / rgamma(1, shape = shp_err, rate = scl_err)
    
    for (pp in 1:p) {
        lambda2[pp] <- 1 / rgamma(1, shape = 1,
                                  rate = 1 / nu[pp] + (beta[pp]^2) / (2 * delta2 * errVar))
        nu[pp] <- 1 / rgamma(1, shape = 1, rate = 1 + 1 / (lambda2[pp]))
    }
    delta2 <- 1 / rgamma(1, shape = 0.5 * (p + 1),
                         rate = 1 / xi + sum((beta * beta) / (lambda2)) / (2 * errVar))
    xi <- 1 / rgamma(1, shape = 1, rate = 1 + 1 / delta2)
    
    list(betaSamp = beta, errVarSamp = errVar, delta2Samp = delta2,
         lambda2Samp = lambda2, nuSamp = nu, xiSamp = xi)
    
}

