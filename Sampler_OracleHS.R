# ylist: List of response vectors where the i-th component is a m_i observations for the i-th subject.
# xlist: List of fixed effect covariates where the i-th component contains m_i observations for the i-th subject.
# zlist: List of random effect covariates where the i-th component contains m_i observations for the i-th subject.
# Sigma: true value of Sigma used to generate the data.
# a0: shape hyperparameter for inverse-gamma prior for error variance.
# b0: scale hyperparameter for inverse-gamma prior for error variance.
# niter: number of MCMC iterations
# nburn: number of burn-in samples
# nthin: thinning size specifying every nthin-th sample is retained.

oracleBVSHS <- function (ylist, xlist, zlist, Sigma, a0, b0,
                         niter, nburn, nthin) {
    
    p <- ncol(xlist[[1]])
    q <- ncol(zlist[[1]])
    n <- length(xlist)
    N <- sum(unlist(lapply(ylist, length)))
    
    ################## For storing final samples #####################
    errVarSamp <- rep(0.0, (niter - nburn) / nthin)
    betaSamp <- matrix(0.0, (niter - nburn) / nthin, p)
    lambda2Samp <- matrix(0.0, (niter - nburn) / nthin, p)
    nuSamp <- matrix(0.0, (niter - nburn) / nthin, p)
    delta2Samp <- rep(0.0, (niter - nburn) / nthin)
    xiSamp <- rep(0.0, (niter - nburn) / nthin)
    
    ####################### Initialization ######################
    beta <- numeric(length = p)
    errVar <- 1 / rgamma(1, shape = a0, rate = b0) 
    nu <- 1 / rgamma(p, shape = 0.5, rate = 1) 
    xi <- 1 / rgamma(1, shape = 0.5, rate = 1) 
    delta2 <- 1 / rgamma(1, shape = 0.5, rate = 1/xi) 
    lambda2 <- 1 / rgamma(p, shape = 0.5, rate = 1/nu) 
    
    cts <- 0
    startTime <- proc.time()
    for(its in 1:niter) {
        if (its %% 1000 == 0) cat("iteration: ", its, "\n")
        samp <- sampleHS(ylist, xlist, zlist, Sigma,
                         beta, errVar, delta2, lambda2, xi, nu,
                         a0, b0)
        errVar <- samp$errVarSamp
        beta <- samp$betaSamp
        delta2 <- samp$delta2Samp
        lambda2 <- samp$lambda2Samp
        xi <- samp$xiSamp
        nu <- samp$nuSamp
        if (its > nburn & its %% nthin == 0) {
            cts <- cts + 1
            errVarSamp[cts] <- samp$errVarSamp
            betaSamp[cts, ] <- samp$betaSamp
            lambda2Samp[cts, ] <- samp$lambda2Samp
            nuSamp[cts, ] <- samp$nuSamp
            delta2Samp[cts] <- samp$delta2Samp
            xiSamp[cts] <- samp$xiSamp
        }
    }
    endTime <- proc.time()
    
    list(errVarSamp = errVarSamp, betaSamp = betaSamp,
         lambda2Samp = lambda2Samp, delta2Samp = delta2Samp,
         sampler_time = endTime - startTime)
}

sampleHS <- function(ylist, xlist, zlist,
                     Sigma, beta, errVar, delta2, lambda2, xi, nu,
                     a0, b0) {
    
    p <- ncol(xlist[[1]])
    q <- ncol(zlist[[1]])
    n <- length(xlist)
    N <- sum(unlist(lapply(ylist, length)))
    
    VarErri <- lapply(zlist,
                      function(zi) {zi %*% tcrossprod(Sigma, zi) + diag(1, nrow(zi))})
    VarErri_svd <- lapply(VarErri, svd)
    VarErri_inv_half <- lapply(VarErri_svd, function(svdi) {svdi$u %*%
            tcrossprod(diag(1 / sqrt(svdi$d)), svdi$u)})
    scaled_xi <- lapply(1:n, function(ii) {VarErri_inv_half[[ii]] %*% xlist[[ii]]})
    scaled_yi <- lapply(1:n, function(ii) {VarErri_inv_half[[ii]] %*% ylist[[ii]]})
    xmat <- do.call(rbind, scaled_xi)
    yvec <- as.matrix(unlist(scaled_yi), ncol = 1)
    
    vinv_beta_post <- crossprod(xmat) + diag(1 / (delta2 * lambda2))
    covBeta <- errVar * chol2inv(chol(vinv_beta_post))
    xty <- crossprod(xmat, yvec)
    muBeta <- drop(covBeta %*% xty) / errVar
    beta <- muBeta + drop(crossprod(chol(covBeta), rnorm(p))) 
    
    resids <- yvec - drop(xmat %*% beta)
    shp_err <- a0 + 0.5 * (N + p) ## a0 = prior shape of errVar
    scl_err <- b0 + 0.5 * (sum(resids^2) + sum((beta * beta) / (delta2 * lambda2))) ## b0 = prior scale of errVar
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

