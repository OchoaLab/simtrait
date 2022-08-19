#' Estimate calibrated power
#'
#' Given a significance level `alpha` and p-values with known causal status, this function estimates the calibrated power.
#' First it estimates the p-value threshold at which the desired type I error of `alpha` is achieved, then it uses this p-value threshold (not `alpha`) to estimate statistical power.
#' Note that these simple empirical estimates are likely to be inaccurate unless the number of p-values is much larger than `1/alpha`.
#'
#' @inheritParams pval_type_1_err
#' @param causal_indexes The vector of causal indexes, defining the true classes used for calibrated power estimation.
#' Values of `causal_indexes` as returned by `sim_trait` work.
#' There must be at least one causal index and at least one non-causal case.
#'
#' @return The calibrated power estimates at each `alpha`
#'
#' @examples
#' # simulate truly null p-values, which should be uniform
#' pvals <- runif(10)
#' # for toy example, take the first two p-values to be truly causal
#' causal_indexes <- 1:2
#' # estimate desired measure
#' pval_power_calib( pvals, causal_indexes )
#'
#' @seealso
#' [pval_aucpr()], a robust proxy for calibrated power that integrates across significance thresholds.
#' 
#' @export
pval_power_calib <- function(pvals, causal_indexes, alpha = 0.05) {
    if ( missing( pvals ) )
        stop( '`pvals` is required!' )
    if ( missing( causal_indexes ) )
        stop( '`causal_indexes` is required!' )
    if ( length( causal_indexes ) == 0 )
        stop( '`causal_indexes` must have at least one index!' )
    
    # in some cases there is nothing to do (LMM has singular information matrix)
    # NA is best value to return in that case (scalar)
    if ( is.null( pvals ) )
        return( NA )
    
    # check range of data here, to complain if it was bad
    if ( any( pvals < 0, na.rm = TRUE ) )
        stop( 'Input p-values included negative values!' )
    if ( any( pvals > 1, na.rm = TRUE ) )
        stop( 'Input p-values included values exceeding 1!' )

    # separate p-values for both classes
    pvals_alt <- pvals[ causal_indexes ]
    pvals_nul <- pvals[ -causal_indexes ]

    # make sure both lists are non-empty
    # pvals_alt was tested through testing causal_indexes
    # just test the other one
    if ( length( pvals_nul ) == 0 )
        stop( 'There were no null (non-causal) cases!' )

    # first thing is to determine the p-value threshold at which the type I error rate is as desired
    # this handy function returns that value, the p-value threshold below which there is an `alpha` proportion of the null data
    # NOTE: this gives an answer no matter how short `pvals_nul` is, even when it's a scalar!  (If quantile not observed, the minimum appears to be used.)  So we don't have to worry about tricky cases there resulting in missing values or errors.
    p_cut <- stats::quantile( pvals_nul, probs = alpha, na.rm = TRUE, names = FALSE )
    
    # now calculate power at that threshold
    power <- sapply( p_cut, function( x ) mean( pvals_alt <= x, na.rm = TRUE ) )
    
    return( power )

}
