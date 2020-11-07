#' Calculate inflation factor from p-values
#'
#' The inflation factor is defined as the median association test statistic divided by the expected median under the null hypothesis, which is typically assumed to have a Chi square distribution.
#' This function takes a p-value distribution and maps its median back to the chi-square value (using the quantile function) in order to compute the inflation factor in the Chi square scale.
#' The full p-value distribution (a mix of null and alternative cases) is used to calculate the desired median value (the true `causal_loci` is not needed, unlike `\link[pval_srmsd]`).
#'
#' @param pvals The vector of association p-values to analyze.
#' This function assumes all p-values are provided (a mix of null and alternative tests).
#' `NA` values are allowed in input, are removed in calculating the median.
#' Non-`NA` values outside of \[0,1\] will trigger an error.
#' @param df The degrees of freedom of the assumed Chi square distribution (default 1).
#'
#' @return The inflation factor
#'
#' @examples
#' # simulate truly null p-values, which should be uniform
#' pvals <- runif(10)
#' # calculate desired measure
#' pval_infl( pvals )
#'
#' @seealso
#' `\link[pval_srmsd]`, a more robust measure of null p-value accuracy, but which requires knowing the true causal loci.
#' 
#' @export
pval_infl <- function( pvals, df = 1 ) {
    if ( missing( pvals ) )
        stop( '`pvals` is required!' )
    # check range of data here, to complain if it was bad
    if ( any( pvals < 0, na.rm = TRUE ) )
        stop( 'Input p-values included negative values!' )
    if ( any( pvals > 1, na.rm = TRUE ) )
        stop( 'Input p-values included values exceeding 1!' )
    
    # observed median p-value
    pval_median_obs <- stats::median( pvals, na.rm = TRUE )
    # expected median p-value
    pval_median_exp <- 0.5
    
    # convert observed and expected median p-values to Chi-square statistics, using the quantile function
    chisq_median_obs <- stats::qchisq( 1 - pval_median_obs, df = df )
    chisq_median_exp <- stats::qchisq( 1 - pval_median_exp, df = df )
    
    # the inflation factor is this ratio
    lambda <- chisq_median_obs / chisq_median_exp
    return( lambda )
}
