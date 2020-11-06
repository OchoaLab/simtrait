#' Signed RMSD measure of null p-value uniformity
#'
#' Quantifies null p-value uniformity by computing the RMSD (root mean square deviation) between the sorted observed null (truly non-causal) p-values and their expected quantiles under a uniform distribution.
#' Meant as a more robust alternative to the "inflation factor" common in the GWAS literature, which compares median values only and uses all p-values (not just null p-values).
#' Our signed RMSD, to correspond with the inflaction factor, includes a sign that depends on the median null p-value:
#' Positive if this median is <= 0.5 (corresponds with test statistic inflation), negative otherwise (test statistic deflation).
#' Zero corresponds to uniform null p-values, which arises in expectation only if test statistics have their assumed null distribution (there is no misspecification, including inflation).
#'
#' @param pvals The vector of association p-values to analyze.
#' This function assumes all p-values are provided (a mix of null and alternative tests), which is subset using the following parameter.
#' `NA` values are allowed in input, are removed in calculating the signed RMSD.
#' Non-`NA` values outside of \[0,1\] will trigger an error.
#' @param causal_indexes The vector of causal indexes, whose p-values will be omitted.
#' Values of `causal_indexes` as returned by `sim_trait` work.
#' There must be at least one causal index.
#' This parameter is required to prevent use of this function except when the true status of every test (null vs alternative) is known.
#' @param detailed If `FALSE` (default) only SRMSD is returned.
#' If `TRUE`, sorted null p-values without NAs and their expectations are returned (useful for plots).
#'
#' @return If `detailed` is `FALSE`, returns the signed RMSD between the observed p-value order statistics and their expectation under true uniformity.
#'
#' If `detailed` is `TRUE`, returns a named list containing:
#' - `srmsd`: The signed RMSD between the observed p-value order statistics and their expectation under true uniformity.
#' - `pvals_null`: Sorted null p-values (observed order statistics).  If any input null p-values were `NA`, these have been removed here (removed by `sort`).
#' - `pvals_unif`: Expected order statistics assuming uniform distribution, same length as `pvals_null`.
#'
#' The detailed data is returned as it is useful for plots.
#' 
#' However, if the input `pvals` is `NULL` (taken for case of singular association test, which is rare but may happen), then the returned value is `NA` if `detailed` was `FALSE`, or otherwise the list contains `NA`, `NULL` and `NULL` for the above three items.
#'
#' @examples
#' # simulate truly null p-values, which should be uniform
#' pvals <- runif(10)
#' # for toy example, take the first p-value to be truly causal (will be ignored below)
#' causal_indexes <- 1
#' # calculate desired measure
#' pval_srmsd( pvals, causal_indexes )
#'
#' @seealso
#' `\link[rmsd]` for the generic root-mean-square deviation function.
#'
#' @export
pval_srmsd <- function(pvals, causal_indexes, detailed = FALSE) {
    if ( missing( pvals ) )
        stop( '`pvals` is required!' )
    if ( missing( causal_indexes ) )
        stop( '`causal_indexes` is required!' )
    if ( length( causal_indexes ) == 0 )
        stop( '`causal_indexes` must have at least one index!' )
    
    # in some cases there is nothing to do (LMM has singular information matrix)
    if (is.null(pvals))
        if (detailed) {
            return(
                list(
                    srmsd = NA, # NA is best value to return in that case (scalar)
                    pvals_null = NULL, # NULL for vectors
                    pvals_unif = NULL
                )
            )
        } else
            return(NA)
    
    # remove causal loci from vector
    pvals_null <- pvals[ -causal_indexes ]

    if ( length( pvals_null ) == 0 )
        stop( 'No loci were null (non-causal)!  (`pvals[ -causal_indexes ]` had length 0)' )
    
    # sort p-values, which has the added benefit of removing NAs
    pvals_null <- sort( pvals_null )

    # check range of data here, to complain if it was bad
    if ( any( pvals_null < 0 ) )
        stop( 'Input p-values included negative values!' )
    if ( any( pvals_null > 1 ) )
        stop( 'Input p-values included values exceeding 1!' )
    
    # expected quantiles
    pvals_unif <- stats::qunif( stats::ppoints( length( pvals_null ) ) )
    
    # compute error metric (RMSD), will add sign next
    srmsd <- rmsd( pvals_unif, pvals_null )
    
    # add sign depending on median null p-value
    # negative means conservative (good)
    if ( stats::median( pvals_null ) > 0.5 )
        srmsd <- - srmsd

    if ( detailed ) {
        # return items of interest in a list
        return(
            list(
                srmsd = srmsd,
                pvals_null = pvals_null, # sorted!
                pvals_unif = pvals_unif
            )
        )
    } else
        return( srmsd )
}

