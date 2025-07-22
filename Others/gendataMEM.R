n <- 48 ## First 36 subjects used for training, next 12 subjects used for testing
m <- c(4, 8, 12)
p <- 300
q <- 14
A <- as.matrix(expand.grid(Rep = 1:50, m = c(4, 8, 12)), ncol = 2)

genData <- function(N, n, p, q, tau2, ranef_covtype = "diag_psd", fixef_corr = "Indep", seedval) {
    library(matrixStats)
    library(Matrix)
    set.seed(seedval)
    
    m <- N / n ## cluster size
    
    ## true beta coefficients 
    nonzero_fixef <- c(1, 0.5, 0.2, 0.1, 0.05)
    fixef <- c(nonzero_fixef, rep(0, p - length(nonzero_fixef)))
    
    ## generate group
    group <- rep(1:n, each = m)
    
    ## generate fixed effect covariates 
    xlist <- list()
    if(fixef_corr == "Indep") {
        for(gg in 1:n) {
            xlist[[gg]] <- matrix(rnorm(m * p, 0, 1), m, p)
        }
    } else if(fixef_corr == "Toep") {
        for(gg in 1:n) {
            X_ind <- matrix(rnorm(m * p, 0, 1), ncol = p)
            xlist[[gg]] <- t(apply(X_ind, 1, function(u) crossprod(chol(toeplitz(0.5 ^ (0:(p-1)))), u)))
        }
    }
    
    ## generate random effect covariates
    zlist <- list()
    for(gg in 1:n) {
        zlist[[gg]] <- matrix(rnorm(m * q, 0, 1), m, q)
    }
    
    ## random effects covariance matrix
    if(ranef_covtype == "diag_psd") {
        ranefCov <- diag(c(rep(0.5, floor(q/2)), rep(0, q - floor(q/2)))) # diagonal psd Sigma based on Li et al. 2021
        svdS <- svd(ranefCov)
        ranefCovSqrt <- svdS$u %*% diag(sqrt(svdS$d)) # L
    } else if(ranef_covtype == "Toep_pd") {
        ranefCov <- toeplitz(0.5 ^ (0:(q-1))) # Toeplitz pd Sigma
        ranefCovSqrt <- chol(ranefCov) # L
    } else if(ranef_covtype == "bdiag_psd") {
        k0 <- 3 # round(log(q)); q <- 14
        nr <- c(5, 5, 4)
        ranefCovSqrt <- matrix(0, q, k0)
        ranefCovSqrt[1:5, 1] <- runif(nr[1], 0, 3)
        ranefCovSqrt[6:10, 2] <- runif(nr[2], 0, 3)
        ranefCovSqrt[11:14, 3] <- runif(nr[3], 0, 3)
        for(i in 1:k0) {
            ranefCovSqrt[, i] <- ranefCovSqrt[, i] / sqrt(sum(ranefCovSqrt[, i] * ranefCovSqrt[, i]))
        }
        ranefCov <- tcrossprod(ranefCovSqrt)
    }
    
    ## compute linear predictors and generate observations
    ylist <- list()
    for(gg in 1:n) {
        mu <- drop(xlist[[gg]] %*% fixef)
        ci <- rnorm(q, 0, sd = sqrt(tau2)) # var(ci) = tau2 * I_q
        ylist[[gg]] <- mu + drop(zlist[[gg]] %*% ranefCovSqrt %*% ci) +
            rnorm(N / n, 0, sd = sqrt(tau2))
    }
    
    list(fixef = fixef,
         ranefCov = ranefCov,
         L = ranefCovSqrt,
         ranef_covtype = ranef_covtype,
         errVar = tau2,
         group = group,
         xlist = xlist,
         zlist = zlist,
         ylist = ylist
    )
}

