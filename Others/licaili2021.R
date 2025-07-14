######### Fixed Effects Estimation #########
library(MASS)
library(glmnet)
library(scalreg)
library(tictoc)
Psi.def<- function(eta,q, cov='diag'){ ##generating different Psi
    if (cov=='diag') {psi = diag(eta,q)}
    if (cov=='sym') {psi=toeplitz(eta^(0:(q-1))); d=1}
    if (cov=='semi') {psi = diag(c(rep(eta[1],floor(q/2)),rep(0,q - floor(q/2))));d=q}
    psi
}

CV_fixed_effects <- function(x, y, z, grp, m, a.seq = seq(0.5, 10, 0.5), inf.coord,
                             standardize, intercept, std_response) {
    subject_ids <- unique(grp)
    n <- length(subject_ids)
    n.tr <- round(0.75 * n)
    best.err <- sum(y^2)
    best.a <- a.seq[1]
    
    for (ia in seq_along(a.seq)) {
        # Subset by subject
        train_ids <- subject_ids[1:n.tr]
        test_ids <- subject_ids[(n.tr + 1):n]
        
        train_idx <- which(grp %in% train_ids)
        test_idx <- which(grp %in% test_ids)
        
        est.re <- Fixed_effects_estimation(
            x = x[train_idx, ],
            y = y[train_idx],
            z = z[train_idx, ],
            m = m,
            grp = grp[train_idx],
            a = a.seq[ia],
            inf.coord = inf.coord,
            standardize, intercept, std_response
        )
        
        pred.err <- sum((y[test_idx] - x[test_idx, ] %*% est.re$beta.hat.lasso)^2)
        
        if (pred.err < best.err) {
            best.err <- pred.err
            best.a <- a.seq[ia]
        }
    }
    
    list(best.a = best.a)
}

Fixed_effects_estimation <- function(x, y, z, m, grp, a, inf.coord,
                                     standardize, intercept, std_response) {
    N <- length(y)
    p <- ncol(x)
    group_levels <- sort(unique(grp))
    n <- length(group_levels)
    
    ###preprocessing
    x <- scale(x)
    y <- y - mean(y)
    
    x.a <- x
    y.a <- y
    tr.a <- 0
    
    for (i in seq_along(group_levels)) {
        cur.mem <- which(grp == group_levels[i])
        mi <- length(cur.mem)
        if (mi == 0) next
        
        zi <- as.matrix(z[cur.mem, ])
        sigmai <- a * tcrossprod(zi) + diag(1, mi)
        sigmai.svd <- svd(sigmai)
        Sig.a.inv.half <- sigmai.svd$u %*% diag(1 / sqrt(sigmai.svd$d)) %*% t(sigmai.svd$u)
        x.a[cur.mem, ] <- Sig.a.inv.half %*% as.matrix(x[cur.mem, ])
        y.a[cur.mem] <- Sig.a.inv.half %*% as.matrix(y[cur.mem])
        tr.a <- tr.a + sum(1 / sigmai.svd$d)
    }
    
    sig.init <- scalreg(x.a, y.a)$hsigma
    lasso.fit <- glmnet(x.a, y.a, lambda = sig.init * sqrt(2 * log(p) / N),
                        standardize = standardize, intercept = intercept, standardize.response = std_response)
    beta.hat.lasso <- as.vector(coef(lasso.fit)[-1])
    beta.hat.all <- as.vector(coef(lasso.fit))
    
    beta.hat.db <- rep(NA, length(inf.coord))
    beta.hat.db.sd <- rep(NA, length(inf.coord))
    
    if (is.null(inf.coord)) {
        return(list(beta.hat.lasso = beta.hat.lasso, beta.hat.db = beta.hat.db, tr.a = tr.a))
    }
    
    res <- y.a - x.a %*% beta.hat.lasso
    
    for (j in seq_along(inf.coord)) {
        col.j <- inf.coord[j]
        sig.x <- scalreg(x.a[, -col.j], x.a[, col.j])$hsigma
        kappa.hat.mlm <- glmnet(x.a[, -col.j], x.a[, col.j], lambda = sig.x * sqrt(2 * log(p) / N),
                                standardize = standardize, intercept = intercept, standardize.response = std_response)
        gam.j <- as.vector(coef(kappa.hat.mlm)[-1])
        wj.mlm <- x.a[, col.j] - x.a[, -col.j] %*% gam.j
        beta.hat.db[j] <- beta.hat.lasso[col.j] + sum(wj.mlm * res) / sum(wj.mlm * x.a[, col.j])
        
        num <- 0
        for (i in seq_along(group_levels)) {
            cur.mem <- which(grp == group_levels[i])
            num <- num + sum(wj.mlm[cur.mem] * res[cur.mem])^2
        }
        
        beta.hat.db.sd[j] <- sqrt(num) / sum(wj.mlm * x.a[, col.j])
    }
    
    list(
        beta.hat.lasso = beta.hat.lasso,
        beta.hat.db = beta.hat.db,
        beta.hat.db.sd = beta.hat.db.sd,
        tr.a = tr.a,
        beta.hat.all = beta.hat.all
    )
}

