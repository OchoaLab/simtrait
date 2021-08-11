#' simtrait: simulate complex traits from genotypes
#'
#' This package enables simulation of complex (polygenic and continuous) traits from a simulated or real genotype matrix.
#' The focus is on constructing the mean and covariance structure of the data to yield the desired heritability.
#' The main function is [sim_trait()], which returns the simulated trait and the vector of causal loci (randomly selected) and their coefficients.
#' The causal coefficients are constructed under two models: *random coefficients* (RC) and *fixed effect sizes* (FES).
#' The function [cov_trait()] computes the expected covariance matrix of the trait given the model parameters (namely the desired heritability and the true kinship matrix).
#' Infinitesimal traits (without causal loci) can also be simulated using [sim_trait_mvn()].
#'
#' Package also provides some functions for evaluating genetic association approaches.
#' [pval_srmsd()] and [pval_infl()] quantify null p-value accuracy, while [pval_aucpr()] quantifies predictive power.
#'
#' The recommended inputs are simulated genotypes with known ancestral allele frequencies.
#' The `bnpsd` package simulates genotypes for admixed individuals, resulting in a complex population structure.
#'
#' For real data it is necessary to estimate the kinship matrix.
#' [popkin::popkin()]` provides high-accuracy kinship estimates.
#' 
#' @examples
#' # construct a dummy genotype matrix
#' X <- matrix(
#'     data = c(
#'         0, 1, 2,
#'         1, 2, 1,
#'         0, 0, 1
#'     ),
#'     nrow = 3,
#'     byrow = TRUE
#' )
#' # made up ancestral allele frequency vector for example
#' p_anc <- c(0.5, 0.6, 0.2)
#' # desired heritability
#' herit <- 0.8
#' # create a dummy kinship matrix for example
#' # make sure it is positive definite!
#' kinship <- matrix(
#'     data = c(
#'         0.6, 0.1, 0.0,
#'         0.1, 0.5, 0.0,
#'         0.0, 0.0, 0.5
#'     ),
#'     nrow = 3
#' )
#'
#' # create simulated trait and associated data
#' # default is *random coefficients* (RC) model
#' obj <- sim_trait(X = X, m_causal = 2, herit = herit, p_anc = p_anc)
#' # trait vector
#' obj$trait
#' # randomly-picked causal locus indeces
#' obj$causal_indexes
#' # regression coefficient vector
#' obj$causal_coeffs
#'
#' # *fixed effect sizes* (FES) model
#' obj <- sim_trait(X = X, m_causal = 2, herit = herit, p_anc = p_anc, fes = TRUE)
#'
#' # either model, can apply to real data by replacing `p_anc` with `kinship`
#' obj <- sim_trait(X = X, m_causal = 2, herit = herit, kinship = kinship)
#'
#' # covariance of simulated traits
#' V <- cov_trait(kinship = kinship, herit = herit)
#' 
#' # draw simulated traits (matrix of replicates) from infinitesimal model
#' traits <- sim_trait_mvn( rep = 10, kinship = kinship, herit = herit )
#' traits
#'
#' # Metrics for genetic association approaches
#' 
#' # simulate truly null p-values, which should be uniform
#' pvals <- runif(10)
#' # for toy example, take these p-value to be truly causal
#' causal_indexes <- c(1, 5, 7)
#'
#' # calculate desired measures
#' # this one quantifies p-value uniformity
#' pval_srmsd( pvals, causal_indexes )
#' # related, calculates inflation factors
#' pval_infl( pvals )
#' # this one quantifies predictive power
#' pval_aucpr( pvals, causal_indexes )
#'
#' 
#' @docType package
#' @name simtrait
"_PACKAGE"

