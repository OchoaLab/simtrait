#' Calculate type I error rate
#'
#' Given a significance level and p-values with known causal status, this function calculates the type I error rate, defined as the proportion of null p-values that are below or equal to the threshold.
#'
#' @inheritParams pval_srmsd
#' @param alpha The desired significance level (default 0.05).
#'
#' @return The type I error rate
#'
#' @examples
#' # simulate truly null p-values, which should be uniform
#' pvals <- runif(10)
#' # for toy example, take the first p-value to be truly causal (will be ignored below)
#' causal_indexes <- 1
#' # calculate desired measure
#' pval_type_1_err( pvals, causal_indexes )
#'
#' @seealso
#' [pval_srmsd()] to directly quantify null p-value uniformity, a more robust alternative to type I error rate. 
#'
#' [pval_infl()] for the more traditional inflation factor, which focuses on the median of the full distribution (combination of causal and null cases).
#'
#' @export
pval_type_1_err <- function( pvals, causal_indexes, alpha = 0.05 ) {
    if ( missing( pvals ) )
        stop( '`pvals` is required!' )
    if ( missing( causal_indexes ) )
        stop( '`causal_indexes` is required!' )
    
    # in some cases there is nothing to do (LMM has singular information matrix)
    if (is.null(pvals))
        return(NA)

    # check range of data here, to complain if it was bad
    if ( any( pvals < 0, na.rm = TRUE ) )
        stop( 'Input p-values included negative values!' )
    if ( any( pvals > 1, na.rm = TRUE ) )
        stop( 'Input p-values included values exceeding 1!' )
    
    if ( is.null( causal_indexes ) ) {
        # only case where null p-values are all p-values
        pvals_null <- pvals
    } else {
        # though we could handle this as the NULL case, I think it's safest to do this instead.
        if ( length( causal_indexes ) == 0 )
            stop( 'non-NULL `causal_indexes` must have at least one index!' )
        
        # remove causal loci from vector
        pvals_null <- pvals[ -causal_indexes ]
        
        if ( length( pvals_null ) == 0 )
            stop( 'No loci were null (non-causal)!  (`pvals[ -causal_indexes ]` had length 0)' )
    }

    # this calculates the desired type I error rate
    tie <- mean( pvals_null <= alpha, na.rm = TRUE )
    
    return( tie )
}

