# simtrait

The `simtrait` package enables simulation of complex traits with user-set number of causal loci and the desired heritability of the trait (the proportion of variance due to genetic effects).

In all cases a genotype matrix is required, which may be simulated with other packages (such as `bnpsd`) or may consist of real data.
Additionally, in order to scale the locus effect sizes to yield the desired the heritability, you will need to provide the true ancestral allele frequencies (most accurate, though only available for simulations) or the kinship matrix (can be estimated for real data, though this is less accurate).

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

### Simulate an admixed population

The first step is to simulate genotypes from an admixed population, to have an example where there is population structure (a non-trivial kinship matrix) and known ancestral allele frequencies.

``` r
devtools::install_github("StoreyLab/bnpsd") # install development version
library(bnpsd) # to simulate an admixed population

# dimensions of data/model
m <- 10000 # number of loci
n_ind <- 30 # number of individuals, smaller than usual for easier visualizations
k <- 3 # number of intermediate subpops

# define population structure
F <- 1:k # FST values for k=3 subpopulations
s <- 0.5 # bias coeff of standard Fst estimator
Fst <- 0.3 # desired final Fst
obj <- q1d(n_ind, k, s=s, F=F, Fst=Fst) # data
# in this case return value is a named list with three items:
Q <- obj$Q # admixture proportions
F <- obj$F # rescaled Fst vector for intermediate subpops

# get pop structure parameters of the admixed individuals
Theta <- coanc(Q,F) # the coancestry matrix
Phi <- coanc_to_kinship(Theta) # kinship matrix

# draw allele freqs and genotypes
out <- rbnpsd(Q, F, m, wantP=FALSE, wantB=FALSE, noFixed=TRUE) # exclude variables not of interest
X <- out$X # genotypes
p_anc <- out$Pa # ancestral AFs
```

### Simulate a random trait


``` r
library(simtrait) # load this package

# parameters of simulation
m_causal <- 100
herit <- 0.8

# create simulated trait and associated data

# version 1: known p_anc (prefered, only applicable to simulated data)
obj <- sim_trait(X=X, m_causal=m_causal, herit=herit, p_anc=p_anc)
# version 2: known kinship (can be estimated with popkin, so more broadly applicable but less precise control of heritability)
obj <- sim_trait(X=X, m_causal=m_causal, herit=herit, kinship=Phi)

# outputs in both versions:
# trait vector
obj$y
# randomly-picked causal locus index
obj$i
# locus effect size vector
obj$beta

# theoretical covariance of the simulated traits
V <- cov_trait(kinship=Phi, herit=herit)
```
