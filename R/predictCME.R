#' @title predictCME
#'
#' @description This function performs posterior predictive sampling using the compressed mixed-effects (CME) model.
#'
#' @param sampler_res a list of posterior samples of \eqn{(\beta, \tau^2, \gamma)} returned by the function \code{\link{sampleCME}}.
#' @param ylist_test a list of response vectors where the \eqn{i^{th}} component contains \eqn{m_i} observations for the \eqn{i^{th}} subject.
#' @param xlist_test a list of fixed effect covariates where the \eqn{i^{th}} component contains \eqn{m_i} observations for the \eqn{i^{th}} subject.
#' @param zlist_test a list of random effect covariates where the \eqn{i^{th}} component contains \eqn{m_i} observations for the \eqn{i^{th}} subject.
#' @param S the \eqn{k_1 \times q} random projection matrix used to implement the collapsed Gibbs sampler, with entries iid from a Gaussian distribution with mean 0 and variance \eqn{1/k_1}.
#' @param R the \eqn{k_2 \times q} random projection matrix used to implement the collapsed Gibbs sampler, with entries iid from a Gaussian distribution with mean 0 and variance \eqn{1/k_2}.
#' @param nom.level nominal level for constructing prediction intervals.
#' @return a list containing the following components:
#' \describe{
#' \item{ypred}{predicted values of response.}
#' \item{mspe}{mean square prediction error.}
#' \item{lower_pi}{lower limits of the prediction intervals.}
#' \item{upper_pi}{upper limits of the prediction intervals.}
#' \item{yhat_samples}{posterior samples of the response.}
#' }
#' @importFrom stats rnorm rgamma quantile
#' @export

predictCME <- function(sampler_res, ylist_test, xlist_test, zlist_test,
                       R, S, nom.level) {
    k1 <- nrow(S); k2 <- nrow(R)
    n_test <- length(ylist_test)
    resids <- vector("list", n_test)
    preds <- vector("list", n_test)
    yhat_samples <- vector("list", n_test)
    lower_pi <- vector("list", n_test)
    upper_pi <- vector("list", n_test)
    
    comp_zlist <- lapply(zlist_test, function(foo) {tcrossprod(foo, S)})

    for (gg in 1:n_test) {
        yihat_iters <- list()
        for (tt in 1:length(sampler_res$errVarSamp)) {
            predCov <- sampler_res$errVarSamp[tt] * (tcrossprod(comp_zlist[[gg]] %*% matrix(sampler_res$gammaSamp[tt, ], k1, k2) %*% R) +
                                                 diag(1, length(ylist_test[[gg]])))
            yihat_iters[[tt]] <- xlist_test[[gg]] %*% sampler_res$betaSamp[tt, ] + drop(crossprod(chol(predCov), rnorm(length(ylist_test[[gg]]))))
        }

        yhat_samples[[gg]] <- do.call("cbind", yihat_iters)
        preds[[gg]] <- rowMeans(yhat_samples[[gg]])
        resids[[gg]] <- ylist_test[[gg]] - preds[[gg]]

        grpvec <- 1:length(ylist_test[[gg]])

        lower_pi[[gg]] <- apply(yhat_samples[[gg]], 1, function(u) quantile(u, (1 - nom.level)/2))
        upper_pi[[gg]] <- apply(yhat_samples[[gg]], 1, function(u) quantile(u, 1 - (1 - nom.level)/2))
    }

    mspe <- mean(unlist(lapply(resids, function(foo) {mean(foo ^ 2)})))

    result_Pred <- list(ypred = unlist(preds), mspe = mspe,
                        lower_pi = unlist(lower_pi), upper_pi = unlist(upper_pi),
                        yhat_samples = yhat_samples)
    result_Pred
}
