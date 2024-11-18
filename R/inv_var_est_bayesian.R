## p_anc_est is the empirical estimate of p, (1/2n)sum_{j=1}^n x_{ij}
## g is the value of the gamma exponent, g=0.5 for regression coefficients
## kinship is the mean kinship value
# the function transforms kinship into nu=((1/kinship)-1)
# nested function pr.ph.p() gives the Pr(p-hat|p)
# the numerator of the expectation is calculated, the integral of w*pr.ph.p() from 0 to 0.5 
# the denominator of the expectation is calculated, the integral of pr.ph.p() from 0 to 0.5
# the expectation is returned

#' Calculate expectations of inverse variance terms over the posterior distribution of the ancestral allele frequency
#'
#' This function assumes that each sample allele frequency `x` in the input has a Beta distribution with mean `p`  and variance p*(1-p)*kinship`, where `p` is the true allele frequency and is unknown, while `kinship` is the mean kinship of the sample and is known.
#' Using this assumption, an expectation over the posterior distribution of `p` is calculated for what we call an inverse variance term, namely `1 / ( p * ( 1 - p ) )^g`, where `p * ( 1 - p )` corresponds to the Bernoulli variance, and `g` is some desired exponent.
#' Note that this function estimates these posteriors each from a single data point `x`.
#'
#' @inheritParams p_anc_est_beta_mle
#' @param g The exponent of the variance term
#'
#' @return A vector of expectations of `1 / ( p * ( 1 - p ) )^g` over the marginalized, unknown true ancestral allelel frequencies `p`, weighted by how likely they are in the posterior given each input sample allele frequency.
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
#' # here g=0.5 estimates `1 / sqrt( p * ( 1 - p ) )`
#' inv_var_est <- inv_var_est_bayesian( p_anc_est, kinship, 0.5 )
#' 
#' @seealso
#' [p_anc_est_beta_mle()] for another way to get less biased estimates of the ancestral allele frequencies and related terms.
#' 
#' @export
inv_var_est_bayesian <- function( p_anc_est, kinship, g ){
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
    inv_var_est <- vapply(
        p_anc_est,
        FUN = function( x ) inv_var_est_bayesian_single( x, nu, g ),
        FUN.VALUE = 0.5,
        USE.NAMES = FALSE
    )
    return( inv_var_est )
}

inv_var_est_bayesian_single <- function( x, nu, g ){
    # preserve NAs without errors (test before all other cases)
    if ( is.na( x ) ) return( NA )
    # these are the edge cases that return infinities for dbeta, here just return infinity without errors
    if ( x == 0 || x == 1 ) return( Inf )
    # code for all other cases
    numer <- stats::integrate(
        function( p ) ( 1 / ( p * ( 1 - p ) )^g ) * stats::dbeta( x, p * nu, ( 1 - p ) * nu ),
        lower = 0,
        upper = 1,
        subdivisions = 2000
    )
    denom <- stats::integrate(
        function( p ) stats::dbeta( x, p * nu, ( 1 - p ) * nu ),
        lower = 0,
        upper = 1
    )
    expec <- numer$value / denom$value
    return(expec)
}
