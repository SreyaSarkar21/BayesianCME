source("Others/gendataMEM.R")

n <- 48; m <- 4; p <- 300; q <- 14
dat <- genData(N = n * m, n = n, p = p, q = q, tau2 = 0.01, ranef_covtype = "diag_psd", fixef_corr = "Indep", seedval = 1)
k1 <- round(log(q)); k2 <- round(log(q))
S <- matrix(rnorm(k1 * q, 0, 1 / sqrt(k1)), k1, q)
R <- matrix(rnorm(k2 * q, 0, 1 / sqrt(k2)), k2, q)

n_train <- 36
train_grps <- 1:n_train
test_grps <- (n_train+1):length(dat$ylist)
library(cme)
res <- sampleCME(ylist = dat$ylist[train_grps], xlist = dat$xlist[train_grps],
                 zlist = dat$zlist[train_grps], a0 = 0.01, b0 = 0.01, gammaVar0 = 1,
                 R = R, S = S, niter = 15000, nburn = 5000, nthin = 1)

betaPostMed <- apply(res$betaSamp, 2, median) ## posterior median
betaTrue <- dat$fixef
q025 <- numeric(length = p)
q975 <- numeric(length = p)
covgs_95 <- numeric(length = p)
for(pp in 1:p) {
    q025[pp] <- quantile(res$betaSamps[, pp], 0.025)
    q975[pp] <- quantile(res$betaSamps[, pp], 0.975)
    covgs_95[pp] <- as.numeric(q025[pp] <= betaTrue[pp] & betaTrue[pp] <= q975[pp])
}
ciwidth_95 <- q975 - q025

modes <- function(x) {
    ux <- unique(x)
    tab <- tabulate(match(x, ux))
    ux[tab == max(tab)]
}

source("Others/S2M.R")
s2m_res <- S2M(res$betaSamp,
               lower = abs(min(betaPostMed[betaPostMed != 0]))/100,
               upper = abs(min(betaPostMed[betaPostMed != 0])),
               l = 25)
s2m_vs <- S2M.vs(s2m_res, H = min(modes(s2m_res$H.b.i)))
tpr <- length(intersect(1:5, s2m_vs)) / 5 ## First 5 indices contain the true non-zero betas
fpr <- length(intersect(6:300, s2m_vs)) / 295

pred_res <- predictCME(sampler_res = res, ylist_test = dat$ylist[test_grps],
                       xlist_test = dat$xlist[test_grps], zlist_test = dat$zlist[test_grps],
                       S = S, R = R, nom.level = 0.95)

ytrue <- unlist(dat$ylist[test_grps])
covgs_pi_95 <- numeric(length = length(ytrue))
for(ij in 1:length(ytrue)) {
    covgs_pi_95[ij] <- as.numeric(pred_res$lower_pi[ij] <= ytrue[ij] & ytrue[ij] <= pred_res$upper_pi[ij])
}

mean(covgs_pi_95)
mean(pred_res$upper_pi - pred_res$lower_pi) ## mean width of prediction interval

