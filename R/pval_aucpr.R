#' Area under the precision-recall curve
#'
#' Calculates the Precision-Recall (PR) Area Under the Curve (AUC) given a vector of p-values and the true classes (causal (alternative) vs non-causal (null)).
#' This is a wrapper around the function `pr.curve` from the `PRROC` package, which actually calculates the AUC (see that for details).
#'
#' @param pvals The vector of association p-values to analyze.
#' `NA` values are allowed in input, are internally set to 1 (worst score) prior to AUC calculation (to prevent methods to get good AUCs by setting more cases to NA).
#' Non-`NA` values outside of \[0,1\] will trigger an error.
#' @param causal_indexes The vector of causal indexes, defining the true classes used for AUC calculation.
#' Values of `causal_indexes` as returned by `sim_trait` work.
#' There must be at least one causal index and at least one non-causal case.
#'
#' @return The PR AUC scalar value.
#'
#' However, if the input `pvals` is `NULL` (taken for case of singular association test, which is rare but may happen), then the returned value is `NA`.
#'
#' @examples
#' # simulate truly null p-values, which should be uniform
#' pvals <- runif(10)
#' # for toy example, take the first two p-values to be truly causal
#' causal_indexes <- 1:2
#' # calculate desired measure
#' pval_aucpr( pvals, causal_indexes )
#'
#' @seealso
#' The function `pr.curve` from the `PRROC` package, which is used internally by `pval_aucpr`.
#' 
#' @export
pval_aucpr <- function(pvals, causal_indexes) {
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
    
    # sometimes p-values are missing (GCATest!), treat as worst possible value (p=1)
    if ( anyNA( pvals ) )
        pvals[ is.na( pvals ) ] <- 1
    
    # check range of data here, to complain if it was bad
    if ( any( pvals < 0 ) )
        stop( 'Input p-values included negative values!' )
    if ( any( pvals > 1 ) )
        stop( 'Input p-values included values exceeding 1!' )
    
    # turn p-values into "scores" (for some reason this is needed to have the correct PR curve; is it just the order flip?)
    scores <- - log( pvals )

    # separate scores for both classes
    scores_alt <- scores[ causal_indexes ]
    scores_nul <- scores[ -causal_indexes ]

    # make sure both lists are non-empty
    # scores_alt was tested through testing causal_indexes
    # just test the other one
    if ( length( scores_nul ) == 0 )
        stop( 'There were no null (non-causal) cases!' )
    
    # generate data, skip curve (default)
    pr <- PRROC::pr.curve( scores_alt, scores_nul )
    
    # return the AUC only!
    return( pr$auc.integral )
}
