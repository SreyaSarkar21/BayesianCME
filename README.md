# Bayesian Compressed Mixed-Effects Models (CME)

This repositary contains:

1. **R package cme**
    - Implements the Bayesian Compressed Mixed-Effects (CME) Model.
    - Install this library with:
    ```r
    devtools::install_github("SreyaSarkar21/BayesianCME")
    library(cme)
    ```


2. **Others folder**
    * contains codes to implement the competing methods.
    * `gendataMEM.R` contains the code used to simulate clustered data with n subjects with m observations each, for different choices of random effects covariance matrix.
    * `sampler_OracleHS.R` contains the code for OracleHS, the Bayesian oracle competitor of CME, where the random effects covariance matrix is set to its true value and a Horseshoe prior is assigned on the fixed effects coefficient.
    * `fanli2012.R` contains the code for implementing the penalized quasi-likelihood method for fixed effects selection by [1].
    * `licaili2021.R` contains the code implementing the penalized quasi-likelihood estimation and inference procedures for fixed effects selection by [2], using their published supplementary code as a reference.
    * `S2M.R` contains the [code](https://github.com/hlstat/Sequential-2-Means/blob/master/S2M.R) implementing the sequential-2-means algorithm by [3] for fixed effects selection.
    * `ExampleCME.R` contains example on how to implement the CME model.
    * The riboflavinV100 data[4] used for benchmarking is available to download [here](https://www.annualreviews.org/content/journals/10.1146/annurev-statistics-022513-115545#supplementary_data).

## References

[1] Yingying Fan and Runze Li (2012). "Variable Selection in Linear Mixed Effects Models." *The Annals of Statistics*, Vol. 40, No. 4, pp. 2043–2068. DOI: [10.1214/12-AOS1028](https://doi.org/10.1214/12-AOS1028).
[2] Sai Li, T. Tony Cai, and Hongzhe Li (2021). "Inference for High-Dimensional Linear Mixed-Effects Models: A Quasi-Likelihood Approach." *Journal of the American Statistical Association*, 117(540), 1835–1846. DOI: [10.1080/01621459.2021.1888740](https://doi.org/10.1080/01621459.2021.1888740).
[3] Hanning Li and Debdeep Pati (2017). "Variable Selection Using Shrinkage Priors." *Computational Statistics & Data Analysis*, 107, 107–119. DOI: [10.1016/j.csda.2016.10.008](https://www.sciencedirect.com/science/article/abs/pii/S0167947316302353?via%3Dihub)
[4] Peter Bühlmann, Markus Kalisch, and Lukas Meier (2014). "High-Dimensional Statistics with a View Toward Applications in Biology." *Annual Review of Statistics and Its Application*, 1, 255–278. DOI: [10.1146/annurev-statistics-022513-115545](https://www.annualreviews.org/content/journals/10.1146/annurev-statistics-022513-115545).
