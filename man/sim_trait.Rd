% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/sim_trait.R
\name{sim_trait}
\alias{sim_trait}
\title{Simulate a complex trait from genotypes}
\usage{
sim_trait(
  X,
  m_causal,
  herit,
  p_anc = NULL,
  kinship = NULL,
  mu = 0,
  sigma_sq = 1,
  labs = NULL,
  labs_sigma_sq = NULL,
  maf_cut = NA,
  loci_on_cols = FALSE,
  m_chunk_max = 1000,
  fes = FALSE
)
}
\arguments{
\item{X}{The \code{m}-by-\code{n} genotype matrix (if \code{loci_on_cols = FALSE}, transposed otherwise), or a BEDMatrix object.
This is a numeric matrix consisting of reference allele counts (in \code{c(0, 1, 2, NA)} for a diploid organism).}

\item{m_causal}{The desired number of causal loci.}

\item{herit}{The desired heritability (proportion of trait variance due to genetics).}

\item{p_anc}{The length-\code{m} vector of true ancestral allele frequencies.
Optional but recommended for simulations.
Either this or \code{kinship} must be specified.}

\item{kinship}{The mean kinship value of the individuals in the data.
The \code{n}-by-\code{n} kinship matrix of the individuals in the data is also accepted.
Optional but recommended for real data.
Either this or \code{p_anc} must be specified.}

\item{mu}{The desired parametric mean value of the trait (scalar, default 0).}

\item{sigma_sq}{The desired parametric variance factor of the trait (scalar, default 1).
Corresponds to the variance of an outbred individual.}

\item{labs}{Optional labels assigning individuals to groups, to simulate environment group effects.
Values can be numeric or strings, simply assigning the same values to individuals in the same group.
If vector (single environment), length must be number of individuals.
If matrix (multiple environments), individuals must be along rows, and environments along columns.
The environments are not required to be nested.
If this is non-\code{NULL}, then \code{labs_sigma_sq} must also be given!}

\item{labs_sigma_sq}{Optional vector of environment effect variance proportions, one value for each environment given in \code{labs} (a scalar if \code{labs} is a vector, otherwise its length should be the number of columns of \code{labs}).
Ignored unless \code{labs} is also given.
As these are variance proportions, each value must be non-negative and \code{sum(labs_sigma_sq) + herit <= 1} is required so residual variance is non-negative.}

\item{maf_cut}{The optional minimum allele frequency threshold (default \code{NA}, no threshold).
This prevents rare alleles from being causal in the simulation.
Threshold is applied to the \emph{sample} allele frequencies and not their true parametric values (\code{p_anc}), even if these are available.}

\item{loci_on_cols}{If \code{TRUE}, \code{X} has loci on columns and individuals on rows; if \code{FALSE} (the default), loci are on rows and individuals on columns.
If \code{X} is a BEDMatrix object, loci are always on the columns (\code{loci_on_cols} is ignored).}

\item{m_chunk_max}{BEDMatrix-specific, sets the maximum number of loci to process at the time.
If memory usage is excessive, set to a lower value than default (expected only for extremely large numbers of individuals).}

\item{fes}{If \code{TRUE}, causal coefficients are inversely proportional to the square root of \code{p_anc * ( 1 - p_anc )} (estimated when \code{p_anc} is unavailable), which ensures \emph{fixed effect sizes} (FES) per causal locus.
Signs (+/-) are drawn randomly with equal probability.
If \code{FALSE} (the default), \emph{random coefficients} (RC) are drawn from a standard Normal distribution.
In both cases coefficients are rescaled to result in the desired heritability.}
}
\value{
A named list containing:
\itemize{
\item \code{trait}: length-\code{n} vector of the simulated trait
\item \code{causal_indexes}: length-\code{m_causal} vector of causal locus indexes
\item \code{causal_coeffs}: length-\code{m_causal} vector of coefficients at the causal loci
\item \code{group_effects}: length-\code{n} vector of simulated environment group effects, or 0 (scalar) if not simulated
}

However, if \code{herit = 0} then \code{causal_indexes} and \code{causal_coeffs} will have zero length regardless of \code{m_causal}.
}
\description{
Simulate a complex trait given a SNP genotype matrix and model parameters, which are minimally: the number of causal loci, the heritability, and either the true ancestral allele frequencies used to generate the genotypes or the mean kinship of all individuals.
An optional minimum marginal allele frequency for the causal loci can be set.
The output traits have by default a zero mean and unit variance (for outbred individuals), but those parameters can be modified.
The code selects random loci to be causal, constructs coefficients for these loci (scaled appropriately) and random Normal independent non-genetic effects and random environment group effects if specified.
There are two models for constructing causal coefficients: random coefficients (RC; default) and fixed effect sizes (FES; i.e., coefficients roughly inversely proportional to allele frequency; use \code{fes = TRUE}).
Suppose there are \code{m} loci and \code{n} individuals.
}
\details{
To center and scale the trait and locus coefficients vector correctly to the desired parameters (mean, variance, heritability), the parametric ancestral allele frequencies (\code{p_anc}) must be known.
This is necessary since in the heritability model the genotypes are random variables (with means given by \code{p_anc} and a covariance structure given by \code{p_anc} and the kinship matrix), so these genotype distribution parameters are required.
If \code{p_anc} are known (true for simulated genotypes), then the trait will have the specified mean and covariance matrix in agreement with \code{\link[=cov_trait]{cov_trait()}}.
To simulate traits using real genotypes, where \code{p_anc} is unknown, a compromise that works well in practice is possible if the mean \code{kinship} is known (see package vignette).
We recommend estimating the mean kinship using the \code{popkin} package!
}
\examples{
# construct a dummy genotype matrix
X <- matrix(
    data = c(
        0, 1, 2,
        1, 2, 1,
        0, 0, 1
    ),
    nrow = 3,
    byrow = TRUE
)
# made up ancestral allele frequency vector for example
p_anc <- c(0.5, 0.6, 0.2)
# made up mean kinship
kinship <- 0.2
# desired heritability
herit <- 0.8

# create simulated trait and associated data
# default is *random coefficients* (RC) model
obj <- sim_trait(X = X, m_causal = 2, herit = herit, p_anc = p_anc)

# trait vector
obj$trait
# randomly-picked causal locus indexes
obj$causal_indexes
# regression coefficients vector
obj$causal_coeffs

# *fixed effect sizes* (FES) model
obj <- sim_trait(X = X, m_causal = 2, herit = herit, p_anc = p_anc, fes = TRUE)

# either model, can apply to real data by replacing `p_anc` with `kinship`
obj <- sim_trait(X = X, m_causal = 2, herit = herit, kinship = kinship)

}
\seealso{
\code{\link[=cov_trait]{cov_trait()}}, \code{\link[=sim_trait_mvn]{sim_trait_mvn()}}
}
