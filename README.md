# simtrait

The `simtrait` package enables simulation of complex traits with user-set number of causal loci and the desired heritability of the trait (the proportion of variance due to genetic effects).

The main function requires a simulated genotype matrix, including the true ancestral allele frequencies.
These parameters are necessary to correctly specify the desired correlation structure.
See the package `bnpsd` for simulating genotypes for admixed individuals (example below).

Simulating a trait from real genotypes is possible with a good kinship matrix estimate.
See the package `popkin` for accurate kinship estimation.

## Installation

<!-- 
You can install the released version of simtrait from [CRAN](https://CRAN.R-project.org) with:
``` r
install.packages("simtrait")
``` 
-->

Install the latest development version from GitHub:
``` r
install.packages("devtools") # if needed
library(devtools)
devtools::install_github("OchoaLab/simtrait", build_opts=c())
```

You can see the package vignette, which has additional documentation, by typing this into your R session:
``` r
vignette('simtrait')
```


## Example

The code below has two parts: (1) simulate genotypes, and (2) simulate the trait.

### Simulate an admixed population

The first step is to simulate genotypes from an admixed population, to have an example where there is population structure and known ancestral allele frequencies.
We use the external package `bnpsd` to achieve this.

``` r
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

Here we apply our package to this simulated genotype data.

``` r
library(simtrait) # load this package

# parameters of simulation
m_causal <- 100
herit <- 0.8

# create simulated trait and associated data

# version 1: known p_anc (prefered, only applicable to simulated data)
obj <- sim_trait(X=X, m_causal=m_causal, herit=herit, p_anc=p_anc)
# version 2: known kinship (more broadly applicable but fewer guarantees)
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
