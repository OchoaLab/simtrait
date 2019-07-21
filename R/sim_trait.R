#' Simulate a complex trait from genotypes
#'
#' Simulate a complex trait \eqn{y} given a SNP genotype matrix and model parameters (the desired heritability and the true ancestral allele frequencies used to generate the genotypes, or alternatively the kinship matrix of the individuals).
#' Users can choose the number of causal loci and minimum marginal allele frequency requirements for the causal loci.
#' The code selects random loci to be causal, draws random Normal effect sizes for these loci (scaled appropriately) and random independent non-genetic effects.
#' Below let there be \eqn{m} loci and \eqn{n} individuals.
#'
#' In order to center and scale the trait and locus effect size vector correctly to the desired parameters (mean, variance factor, and heritability), the parametric ancestral allele frequencies (\code{p_anc}) must be known.
#' This is necessary since in the context of Heritability the genotypes are themselves random variables (with means given by \code{p_anc} and a covariance structure given by \code{p_anc} and the kinship matrix), so the parameters of the genotypes must be taken into account.
#' If \code{p_anc} are indeed known (true for simulated genotypes), then the trait will have the specified mean and covariance matrix in agreement with \code{\link{cov_trait}}.
#'
#' If the desire is to simulate a trait using real genotypes, where \code{p_anc} is unknown, a compromise that works well in practice is possible if the kinship matrix (\code{kinship}) is known (see package vignette).
#' The kinship matrix can be estimated accurately using the \code{popkin} package!
#' 
#' @param X The \eqn{m \times n}{m-by-n} genotype matrix (if `loci_on_cols = FALSE`, transposed otherwise), or a BEDMatrix object`.
#' This is a numeric matrix consisting of reference allele counts (in \code{c(0,1,2,NA)} for a diploid organism).
#' @param m_causal The number of causal loci desired.
#' @param herit The desired heritability (proportion of trait variance due to genetics).
#' @param p_anc The length-\eqn{m} vector of true ancestral allele frequencies.
#' Recommended way to adjust the simulated trait to achieve the desired heritability and covariance structure.
#' Either this or \code{kinship} must be specified.
#' @param kinship The \eqn{n \times n}{n-by-n} kinship matrix of the individuals in the data.
#' This offers an alternative way to adjust the simulated parameters parameters to achieve the desired covariance structure for real genotypes, since \code{p_anc} is only known for simulated data.
#' Either this or \code{p_anc} must be specified.
#' @param mu The desired parametric mean value of the trait (default zero).
#' The sample mean of the trait will not be exactly zero, but instead have an expectation of \code{mu} (with potentially large variance depending on the kinship matrix and the heritability).
#' @param sigma_sq The desired parametric variance factor of the trait (default 1).
#' This factor corresponds to the variance of an outbred individual (see \code{\link{cov_trait}}).
#' @param maf_cut The optional minimum allele frequency threshold (default 5\%).
#' This prevents rare alleles from being causal in the simulation.
#' Note that this threshold is applied to the sample allele frequencies and not their true parametric values (\code{p_anc}), even if these are available.
#' @param loci_on_cols If \code{TRUE}, \eqn{X} has loci on columns and individuals on rows; if false (the default), loci are on rows and individuals on columns.
#' If \eqn{X} is a BEDMatrix object, loci are taken to be on the columns (regardless of the value of \code{loci_on_cols}).
#'
#' @return A list containing the simulated \code{trait} (length \eqn{n}), the vector of causal locus indexes \code{causal_indexes} (length \eqn{m_causal}), and the locus effect size vector \code{causal_coeffs} (length \eqn{m_causal}) at the causal loci.
#' However, if \code{herit = 0} then \code{causal_indexes} and \code{causal_coeffs} will have zero length regardless of \code{m_causal}.
#'
#' @examples
#' # construct a dummy genotype matrix
#' X <- matrix(
#'            data = c(0,1,2,1,2,1,0,0,1),
#'            nrow = 3,
#'            byrow = TRUE
#'            )
#' # made up ancestral allele frequency vector for example
#' p_anc <- c(0.5, 0.6, 0.2)
#' 
#' # create simulated trait and associated data
#' obj <- sim_trait(X = X, m_causal = 2, herit = 0.8, p_anc = p_anc)
#' 
#' # trait vector
#' obj$trait
#' # randomly-picked causal locus indexes
#' obj$causal_indexes
#' # locus effect size vector
#' obj$causal_coeffs
#' 
#' @export
sim_trait <- function(X, m_causal, herit, p_anc, kinship, mu = 0, sigma_sq = 1, maf_cut = 0.05, loci_on_cols = FALSE) {
    # check for missing parameters
    if (missing(X))
        stop('genotype matrix `X` is required!')
    if (missing(m_causal))
        stop('the number of causal loci `m_causal` is required!')
    if (missing(herit))
        stop('the heritability `herit` is required!')
    if (missing(p_anc) && missing(kinship))
        stop('either the true ancestral allele frequency vector `p_anc` or `kinship` are required!')
    
    # other checks
    if (length(mu) != 1)
        stop('`mu` must be a scalar! (input has length ', length(mu), ')')
    if (length(sigma_sq) != 1)
        stop('`sigma_sq` must be a scalar! (input has length ', length(sigma_sq), ')')
    if (length(herit) != 1)
        stop('`herit` must be a scalar! (input has length ', length(herit), ')')
    if (herit < 0)
        stop('`herit` must be non-negative!')
    if (herit > 1)
        stop('`herit` cannot be greater than 1!')
    if (sigma_sq <= 0)
        stop('`sigma_sq` must be positive!')

    # simplifies subsetting downstream
    if (class(X) == 'BEDMatrix')
        loci_on_cols <- TRUE

    # get m_loci to check m_causal before more heavy computations
    m_loci <- if (loci_on_cols) ncol(X) else nrow(X)
    if (m_causal > m_loci)
        stop('m_causal (', m_causal, ') exceeds the number of loci (', m_loci, ')!')
    
    if (herit == 0) {
        # lots of work can be avoided in this edge case
        # the index and coefficients vectors are empty
        causal_indexes <- c()
        causal_coeffs <- c()
        G <- 0 # construct a trivial genotype effect of zero (becomes vector automatically later)
    } else {
        
        ###################
        ### CAUSAL LOCI ###
        ###################
        
        # compute marginal allele frequencies
        p_anc_hat <- allele_freqs(X, loci_on_cols = loci_on_cols)
        
        # select random SNPs! this performs the magic...
        # also runs additional checks
        causal_indexes <- select_loci(p_anc_hat, m_causal, maf_cut)
        
        # draw random SNP coefficients for selected loci
        causal_coeffs <- stats::rnorm(m_causal, 0, 1)

        # subset data to consider causal loci only
        p_anc_hat <- p_anc_hat[causal_indexes]
        if (!missing(p_anc))
            p_anc <- p_anc[causal_indexes] # subset if available
        # the subset of causal data
        # (drop = FALSE for keeping as a matrix even if m_causal == 1)
        if (loci_on_cols) {
            # also transpose for consistent behavior downstream
            X <- t( X[, causal_indexes, drop = FALSE] )
        } else{
            X <- X[causal_indexes, , drop = FALSE]
        }
        
        ###############
        ### KINSHIP ###
        ###############
        
        if (!missing(kinship)) {
            # precompute some things when this is present
            mean_kinship <- mean(kinship)
        }
        
        #############
        ### SCALE ###
        #############
        
        # to scale causal_coeffs to give correct heritability, we need to estimate the pq = p(1-p) vector
        # calculate pq = p_anc * (1 - p_anc) in one of two ways
        if ( !missing(p_anc) ) { # this takes precedence, it should be most accurate
            # direct calculation
            pq <- p_anc * (1 - p_anc)
        } else if ( !missing(kinship) ) {
            # indirect, infer from genotypes and mean kinship
            # recall E[ p_anc_hat * (1 - p_anc_hat) ] = pq * (1 - mean_kinship), so we solve for pq:
            pq <- p_anc_hat * (1 - p_anc_hat) / (1 - mean_kinship)
        } else {
            # a redundant check (if this were so, it should have died earlier)
            stop('either `p_anc` or `kinship` must be specified!')
        }
        
        # the initial genetic variance is
        sigma_0 <- sqrt( 2 * sum( pq * causal_coeffs^2 ) )
        # adjust causal_coeffss so final variance is sigma_sq*herit as desired!
        causal_coeffs <- causal_coeffs * sqrt( sigma_sq * herit ) / sigma_0 # scale by standard deviations
        
        # construct genotype signal
        if (any(is.na(X))) {
            # if any of the causal loci are missing, let's treat them as zeroes
            # this isn't perfect but we must do something to apply this to real data
            X[is.na(X)] <- 0
        }
        G <- drop( causal_coeffs %*% X ) # this is a vector
        # NOTE by construction:
        # Cov(G) = 2 * herit * kinship
        
        ##############
        ### CENTER ###
        ##############

        # calculate the mean of the genetic effect
        if ( !missing(p_anc) ) {
            # parametric solution
            muXB <- 2 * drop( causal_coeffs %*% p_anc )
        } else {
            # works very well assuming causal_coeffs and p_anc are uncorrelated!
            muXB <- 2 * sum( causal_coeffs ) * mean( p_anc_hat )
        }
        # in all cases:
        # - remove the mean from the genotypes (muXB)
        # - add the desired mean
        G <- G - muXB + mu
    }

    if (herit == 1) {
        E <- 0 # in this edge case there is no "noise", just genotype effects
    } else {
        # length of E
        n_ind <- ncol(X)
        # draw noise
        E <- stats::rnorm(n_ind, 0, (1 - herit) * sigma_sq ) # noise has mean zero but variance ((1-herit) * sigma_sq)
        # NOTE by construction:
        # Cov(E) = (1-herit) * sigma_sq * I
    }

    # lastly, here's the trait:
    trait <- G + E

    # return all these things
    list(trait = trait, causal_indexes = causal_indexes, causal_coeffs = causal_coeffs)
}
