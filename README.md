# simtrait

The `simtrait` R package enables simulation of complex traits with user-set number of causal loci and the desired heritability of the trait (the proportion of variance due to genetic effects).

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
```R
install.packages("devtools") # if needed
library(devtools)
install_github("OchoaLab/simtrait", build_opts = c())
```

You can see the package vignette, which has more detailed documentation, by typing this into your R session:
```R
vignette('simtrait')
```


## Example

The code below has two parts: (1) simulate genotypes, and (2) simulate the trait.

### Simulate an admixed population

The first step is to simulate genotypes from an admixed population, to have an example where there is population structure and known ancestral allele frequencies.
We use the external package `bnpsd` to achieve this.

```R
library(bnpsd) # to simulate an admixed population

# dimensions of data/model
# number of loci
m_loci <- 10000
# number of individuals, smaller than usual for easier visualizations
n_ind <- 30
# number of intermediate subpops
k_subpops <- 3

# define population structure
# FST values for k = 3 subpopulations
inbr_subpops <- 1 : k_subpops
# bias coeff of standard Fst estimator
bias_coeff <- 0.5
# desired final Fst of admixed individuals
Fst <- 0.3
obj <- admix_prop_1d_linear(
    n_ind,
    k_subpops,
    bias_coeff = bias_coeff,
    coanc_subpops = inbr_subpops,
    fst = Fst
)
admix_proportions <- obj$admix_proportions
# rescaled Fst vector for intermediate subpops
inbr_subpops <- obj$coanc_subpops

# get pop structure parameters of the admixed individuals
concestry <- coanc_admix(admix_proportions, inbr_subpops)
kinship <- coanc_to_kinship(concestry)

# draw allele freqs and genotypes
out <- draw_all_admix(admix_proportions, inbr_subpops, m_loci)
# genotypes
X <- out$X
# ancestral allele frequencies
p_anc <- out$p_anc
```

### Simulate a random trait

Here we apply our package to this simulated genotype data.

```R
library(simtrait) # load this package

# parameters of simulation
m_causal <- 100
herit <- 0.8

# create simulated trait and associated data

# version 1: known p_anc (prefered, only applicable to simulated data)
obj <- sim_trait(X = X, m_causal = m_causal, herit = herit, p_anc = p_anc)
# version 2: known kinship (more broadly applicable but fewer guarantees)
obj <- sim_trait(X = X, m_causal = m_causal, herit = herit, kinship = kinship)

# outputs in both versions:
# trait vector
obj$trait
# randomly-picked causal locus index
obj$causal_indexes
# locus effect size vector
obj$causal_coeffs

# theoretical covariance of the simulated traits
V <- cov_trait(kinship = kinship, herit = herit)
```
