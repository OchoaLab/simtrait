#' Per-locus heritability contribution from allele frequency and causal coefficient
#'
#' Calculates the vector of per-locus heritability values, with each causal locus `i` calculated as
#' `h_i^2 = 2 * p_i * ( 1 - p_i ) * beta_i^2 / sigma_sq`,
#' where `p_i` is the ancestral allele frequency,
#' `beta_i` is the causal coefficient,
#' and `sigma_sq` is the trait variance scale.
#' These are all assumed to be true parameters (not estimated).
#'
#' @param p_anc The ancestral allele frequency vector.
#' @param causal_coeffs The vector of causal coefficients.
#' @param causal_indexes The optional vector of causal indexes.
#' If `NULL` (default), `p_anc` and `causal_coeffs` are assumed to be for causal loci only (must be the same length).
#' If non-`NULL`, `causal_loci` is used to subset both `p_anc` and `causal_coeffs` as needed: if each of these vectors is longer than `causal_loci`, then it is subset; otherwise they must have equal lengths as `causal_loci` or an error is thrown.
#' @param sigma_sq The parametric variance factor of the trait (default 1).
#' This factor corresponds to the variance of an outbred individual.
#'
#' @return The vector of per-locus heritability contributions.
#' The sum of these values gives the overall heritability.
#' This value can be greater than one (or wrong, more generally) if `sigma_sq` is misspecified.
#'
#' @examples
#' # create toy random data
#' m_loci <- 10
#' # ancestral allele frequencies
#' p_anc <- runif( m_loci )
#' # causal loci
#' causal_coeffs <- rnorm( m_loci ) / m_loci
#' # resulting heritability contributions vector
#' herit_loci( p_anc, causal_coeffs )
#'
#' @seealso
#' [sim_trait()] generates random traits by drawing causal loci and their coefficients to fit a desired heritability.
#' [cov_trait()] calculates the covariance structure of the random traits.
#' 
#' @export
herit_loci <- function(
                       p_anc,
                       causal_coeffs,
                       causal_indexes = NULL,
                       sigma_sq = 1
                       ){
    if ( missing( p_anc ) )
        stop( '`p_anc` (vector of ancestral allele frequencies) is required!' )
    if ( missing( causal_coeffs ) )
        stop( '`causal_coeffs` vector is required!' )

    # make sure sigma_sq has valid values
    if ( length( sigma_sq ) != 1 )
        stop( '`sigma_sq` must have length 1 (instead has length ', length( sigma_sq ), ')!' )
    if ( is.na( sigma_sq ) )
        stop( '`sigma_sq` cannot be `NA`!' )
    if ( sigma_sq <= 0 )
        stop( '`sigma_sq` must be positive (passed ', sigma_sq ,')!' )

    # subset data if desired
    if ( !is.null( causal_indexes ) ) {
        # first determine how many loci are causal
        m_loci <- length( causal_indexes )
        # if p_anc has more elements, subset
        if ( length( p_anc ) > m_loci ) 
            p_anc <- p_anc[ causal_indexes ]
        # if causal_coeffs has more elements, subset
        if ( length( causal_coeffs ) > m_loci ) 
            causal_coeffs <- causal_coeffs[ causal_indexes ]
    }

    # at this point, both main vectors should be same length, check!
    if ( length( p_anc ) != length( causal_coeffs ) )
        stop( 'lengths of `p_anc` (', length( p_anc ), ') and `causal_coeffs` (', length( causal_coeffs ), ') disagree!' )

    # calculate heritability vector, one value per causal locus
    herit_vec <- 2 * p_anc * ( 1 - p_anc ) * causal_coeffs^2 / sigma_sq

    return( herit_vec )
}
