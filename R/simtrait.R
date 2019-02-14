#' simtrait: a package for simulating complex traits from genotypes
#'
#' This package enables simulation of complex (polygenic and continuous) traits from a simulated or real genotype matrix.
#' The focus is on controling the mean and covariance structure of the data to yield the desired heritability under arbitrary population structures (any underlying kinship matrix).
#' The main function is \code{\link{sim_trait}}, which returns the simulated trait and the vector of causal loci (randomly selected) and their effect sizes (randomly drawn and scaled appropriately).
#' The function \code{\link{cov_trait}} computes the expected covariance matrix of the trait given the model parameters (namely the desired heritability and the true kinship matrix).
#'
#' The recommended inputs are simulated genotypes with known ancestral allele frequencies.
#' The `bnpsd` package simulates genotypes for admixed individuals, resulting in a complex population structure.
#'
#' For real data it is necessary to estimate the kinship matrix.
#' The `popkin` package provides high-accuracy kinship estimates.
#' 
#' @examples
#' # construct a dummy genotype matrix
#' X <- matrix(
#'            data=c(0,1,2,1,2,1,0,0,1),
#'            nrow=3,
#'            byrow=TRUE
#'            )
#' # made up ancestral allele frequency vector for example
#' p_anc <- c(0.5, 0.6, 0.2)
#' # desired heritability
#' herit <- 0.8
#' # create simulated trait and associated data
#' obj <- sim_trait(X=X, m_causal=2, herit=herit, p_anc=p_anc)
#' # trait vector
#' obj$y
#' # randomly-picked causal locus indeces
#' obj$i
#' # locus effect size vector
#' obj$beta
#'
#' # create a dummy kinship matrix for example
#' Phi <- matrix(
#'              data=c(0.6,0.1,0, 0.1,0.6,0.1, 0,0.1,0.6),
#'              nrow=3,
#'              byrow=TRUE
#'              )
#' # covariance of simulated traits
#' V <- cov_trait(kinship=Phi, herit=herit)
#' 
#' @docType package
#' @name simtrait
"_PACKAGE"

