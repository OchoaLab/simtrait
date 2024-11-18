#' Calculate maximum likelihood estimates of allele frequencies from Beta model
#'
#' For each sample allele frequency `x` in the input, this function calculates the maximum likelihood estimate of the true allele frequency p` assuming that the sample was drawn from a Beta distribution with mean `p` and variance p*(1-p)*kinship`, where `kinship` is the mean kinship of the sample and is known.
#' Note that, oddly, this function estimates the unknown `p` from a single data point `x`.
#' Nevertheless, this procedure results in favorable shrinkage of estimate towards 0.5.
#' 
#' @param p_anc_est A vector of sample allele frequencies (each refered to as `x` above)
#' @param kinship The mean kinship coefficient of the data
#'
#' @return A vector of maximum likelihood estimates of allele frequencies.
#' The length of the vector is the same as the input `p_anc_est`, and each element in the output was estimated from each corresponding element in the input only.
#' 
#' @examples
#' # select a relatively high value for example
#' kinship <- 0.1
#' 
#' # try a grid of values, including edge cases (0,1)
#' m_loci <- 1000
#' p_anc_est <- 0 : (m_loci - 1) / (m_loci - 1)
#'
#' # calculate the desired estimates!
#' p_anc_mle <- p_anc_est_beta_mle( p_anc_est, kinship )
#' 
#' # notice values tend to be shrunk towards 0.5,
#' # except values near the edges of the range are shrunk less
#' plot( p_anc_est, p_anc_mle )
#' abline(0,1)
#'
#' @seealso
#' [inv_var_est_bayesian()] for another way to get unbiased estimates of a specific class of inverse variance terms of interest in this project.
#'
#' @export
p_anc_est_beta_mle <- function( p_anc_est, kinship ){
    # validate inputs
    if ( missing( p_anc_est ) )
        stop( '`p_anc_est` is required!' )
    if ( missing( kinship ) )
        stop( '`kinship` is required!' )
    # detailed validations of p_anc_est
    if ( !is.numeric( p_anc_est ) )
        stop( '`p_anc_est` must be numeric!' )
    if ( any( p_anc_est < 0, na.rm = TRUE ) )
        stop( '`p_anc_est` cannot have negative values!' )
    if ( any( p_anc_est > 1, na.rm = TRUE ) )
        stop( '`p_anc_est` cannot exceed 1!' )
    # detailed validations of kinship
    if ( !is.numeric( kinship ) )
        stop( '`kinship` must be numeric!' )
    if ( length( kinship ) > 1 )
        stop( '`kinship` must be a scalar!' )
    if ( kinship <= 0 )
        stop( '`kinship` must be positive!' )
    if ( kinship >= 1 )
        stop( '`kinship` must be less than one!' )
    
    nu <- 1 / kinship - 1
    # optimize works on scalar functions, and we want to do it separately for each value in the vector `p_anc_est`
    p_mle <- vapply(
        p_anc_est,
        FUN = function( x ) p_anc_est_beta_mle_single( x, nu ),
        FUN.VALUE = 0.5,
        USE.NAMES = FALSE
    )
    return( p_mle )
}

# wrapper that handles edge cases elegantly
# no validations here, they are already done outside
p_anc_est_beta_mle_single <- function( x, nu ) {
    # preserve NAs without errors (test before all other cases)
    if ( is.na( x ) ) return( NA )
    # these are the edge cases that return infinities for dbeta, by continuity we want them to be fixed:
    if ( x == 0 || x == 1 ) return( x )
    # code for all other cases
    stats::optimize(
               function( p ) stats::dbeta( x, p * nu, ( 1 - p ) * nu ),
               lower = 0,
               upper = 1,
               maximum = TRUE
           )$maximum
}
