# simtrait

The `simtrait` package enables simulation of complex traits with user-set number of causal loci and the desired heritability of the trait (the proportion of variance due to genetic effects).

In all cases a genotype matrix is required, which may be simulated with other packages (such as `bnpsd`) or may consist of real data.
Also, in order to scale the locus effect sizes to yield the desired the heritability, you will need to provide the true ancestral allele frequencies (most accurate, though only available for simulations) or the mean kinship (can be estimated for real data, though this is less accurate).

## Installation

NOT TRUE YET:
You can install the released version of simtrait from [CRAN](https://CRAN.R-project.org) with:

``` r
install.packages("simtrait")
```

Install from private github repo.
You'll need your GitHub API key setup!
``` r
library(devtools) # install if needed
devtools::install_github("OchoaLab/simtrait")
```


## Example

This is a basic example which shows you how to solve a common problem:

Should show how to simulate from BNPSD.
``` r
# construct a dummy genotype matrix
X <- matrix(
           data=c(0,1,2,1,2,1,0,0,1),
           nrow=3,
           byrow=TRUE
           )
# made up mean kinship for example
mean_kinship <- 0.1
# create simulated trait and associated data
obj <- sim_trait(X=X, m_causal=1, herit=0.8, mean_kinship=mean_kinship)
# trait vector
obj$y
# randomly-picked causal locus index
obj$i
# locus effect size vector
obj$beta

# covariance of simulated traits
V <- cov_trait(kinship=Phi, herit=0.8)
```
