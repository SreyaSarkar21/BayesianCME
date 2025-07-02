# Bayesian Compressed Mixed-Effects Models (CME)

This repositary contains the R code for the simulation studies and riboflavinV100 data analysis.

## Files

* `gendataMEM.R` contains the code to simulate clustered data with n subjects with m observations each, for different choices of random effects covariance matrix.
* `Sampler_CME.R` contains the code for the collapsed Gibbs sampler devised for posterior sampling of the parameters using CME model.
* `Sampler_OracleHS.R` contains the code for OracleHS, the Bayesian oracle competitor of CME, where the random effects covariance matrix is set to its true value and a Horseshoe prior is assigned on the fixed effects coefficient.
* `FanLi2012.R` contains the code for implementing the penalized quasi-likelihood method for fixed effects selection by Fan and Li (2012) [1].
* `LiCaiLi2021.R` contains the code implementing the penalized quasi-likelihood estimation and inference procedures for fixed effects selection by Li et al. (2021) [2], using their published supplementary code as a reference.

## References

[1] Yingying Fan and Runze Li (2012). "Variable Selection in Linear Mixed Effects Models." *The Annals of Statistics*, Vol. 40, No. 4, pp. 2043–2068. DOI: [10.1214/12-AOS1028](https://doi.org/10.1214/12-AOS1028).
[2] Sai Li, T. Tony Cai, and Hongzhe Li (2021). "Inference for High-Dimensional Linear Mixed-Effects Models: A Quasi-Likelihood Approach." *Journal of the American Statistical Association*, 117(540), 1835–1846. DOI: [10.1080/01621459.2021.1888740](https://doi.org/10.1080/01621459.2021.1888740).

