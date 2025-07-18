% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/sampleCME.R
\name{sampleCME}
\alias{sampleCME}
\title{sampleCME}
\usage{
sampleCME(ylist, xlist, zlist, a0, b0, gammaVar0, R, S, niter, nburn, nthin)
}
\arguments{
\item{ylist}{a list of response vectors where the \eqn{i^{th}} component contains \eqn{m_i} observations for the \eqn{i^{th}} subject.}

\item{xlist}{a list of fixed effect covariates where the \eqn{i^{th}} component contains \eqn{m_i \times p} design matrix for the \eqn{i^{th}} subject.}

\item{zlist}{a list of random effect covariates where the \eqn{i^{th}} component contains \eqn{m_i \times q} design matrix for the \eqn{i^{th}} subject.}

\item{a0}{shape hyperparameter for inverse-gamma prior for error variance.}

\item{b0}{scale hyperparameter for inverse-gamma prior for error variance.}

\item{gammaVar0}{variance hyperparameter \eqn{\sigma_{\gamma}^2} for Gaussian prior on the columnwise vectorization of the compressed covariance parameter \eqn{\Gamma}.}

\item{R}{a \eqn{k_2 \times q} random projection matrix with entries iid from a Gaussian distribution with mean 0 and variance \eqn{1/k_2}.}

\item{S}{a \eqn{k_1 \times q} random projection matrix with entries iid from a Gaussian distribution with mean 0 and variance \eqn{1/k_1}.}

\item{niter}{number of MCMC iterations.}

\item{nburn}{number of burn-in samples.}

\item{nthin}{thinning size specifying every nthin-th sample is retained.}
}
\value{
a list containing the following components:
\describe{
\item{betaSamp}{a matrix with \eqn{\bigl(\text{niter}-\text{nburn}\bigr)/\text{nthin}} rows and \eqn{p} columns containing posterior samples of fixed effect parameter \eqn{\beta}.}
\item{errVarSamp}{a numeric vector of length \eqn{\bigl(\text{niter}-\text{nburn}\bigr)/\text{nthin}} containing posterior samples of error variance \eqn{\tau^2}.}
\item{gammaSamp}{a matrix with \eqn{\bigl(\text{niter}-\text{nburn}\bigr)/\text{nthin}} rows and \eqn{k_1 k_2} columns containing posterior samples of columnwise vectorization of the compressed covariance parameter \eqn{\Gamma}.}
\item{sampler_time}{time taken by the sampler to run.}
}
}
\description{
This function implements the collapsed Gibbs sampler for fitting the compressed mixed-effects (CME) model.
}
