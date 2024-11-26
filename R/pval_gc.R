#' Apply genomic control transformation to a list of p-values
#'
#' Genomic control (GC) is a procedure designed to correct the calibration of p-values calculated from test statistics that are supposed to have a Chi-squared distribution but are misspecified and have a different median than desired, in which case the test statistics are divided by the inflation factor (the ratio between the observed and desired medians).
#' This function applies the correction using p-values only, calculating the underlying test statistics assuming they are Chi-squared ditributed with known degrees of freedom.
#' The function guarantees that the corrected p-values have a median of 0.5, as desired if they were calibrated, but otherwise there is no guarantee that these corrected p-values will have a reasonable distribution.
#' 
#' @inheritParams pval_infl
#'
#' @return A list with the following named elements:
#' - `pvals`: the vector of GC-corrected p-values
#' - `lambda`: the inflation factor of the input p-values, which was used to calculate the corrected p-values
#'
#' @examples
#' # a simulated set of highly inflated p-values, skewed toward zero
#' pvals <- rbeta( 100, 1, 10 )
#' hist( pvals )
#'
#' # calculate the GC-corrected p-values
#' obj <- pval_gc( pvals )
#' pvals2 <- obj$pvals
#' # and get the inflation factor of the original data:
#' obj$lambda
#'
#' # note GC-corrected p-values are often not uniform:
#' hist( pvals2 )
#' # but they have a median of 0.5
#' median( pvals2 )
#' 
#' @seealso
#' [pval_infl()], which is used internally to calculate the inflation factor.
#'
#' @export
pval_gc <- function( pvals, df = 1 ) {
    if ( missing( pvals ) )
        stop( '`pvals` is required!' )

    # in some cases there is nothing to do (LMM has singular information matrix)
    if ( is.null( pvals ) )
        return( NULL )
    
    # let `pval_infl` test p-values otherwise (range tests)
    
    # calculate inflation factor
    lambda <- pval_infl( pvals, df )
    # calculate chi squared statistics from observed p-values assuming that is their true distribution
    chisqs <- stats::qchisq( pvals, df = df, lower.tail = FALSE )
    # apply genomic control correction
    chisqs <- chisqs / lambda
    # convert back to p-values
    pvals <- stats::pchisq( chisqs, df = df, lower.tail = FALSE )
    
    # done, return those!
    return( list( pvals = pvals, lambda = lambda ) )
}
