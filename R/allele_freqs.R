#' Compute locus allele frequencies
#'
#' On a regular matrix, this is essentially a wrapper for [colMeans()] or [rowMeans()] depending on `loci_on_cols`.
#' On a BEDMatrix object, the locus allele frequencies are computed keeping memory usage low.
#'
#' @param X The genotype matrix (regular R matrix or BEDMatrix object).
#' Missing values are ignored in averages.
#' @param loci_on_cols If `TRUE`, `X` has loci on columns and individuals on rows; if false (the default), loci are on rows and individuals on columns.
#' If `X` is a BEDMatrix object, code assumes loci on columns (`loci_on_cols` is ignored).
#' @param fold If `TRUE`, allele frequencies are converted to minor allele frequencies.
#' Default is to return frequencies for the given allele counts in `X` (regardless of whether it is the minor or major allele).
#' @param m_chunk_max BEDMatrix-specific, sets the maximum number of loci to process at the time.
#' If memory usage is excessive, set to a lower value than default (expected only for extremely large numbers of individuals).
#' @param subset_ind Optionally subset individuals by providing their indexes (negative indexes to exclude) or a boolean vector (in other words, the usual ways to subset matrices).
#' Most useful for BEDMatrix inputs, to subset chunks and retain low memory usage.
#'
#' @return The vector of allele frequencies, one per locus.
#' Names are set to the locus names, if present.
#'
#' @examples
#' # Construct toy data
#' X <- matrix(
#'     c(0, 1, 2,
#'       1, 0, 1,
#'       1, NA, 2),
#'     nrow = 3,
#'     byrow = TRUE
#' )
#' 
#' # row means
#' allele_freqs(X)
#' c(1/2, 1/3, 3/4)
#'
#' # row means, in minor allele frequencies
#' allele_freqs(X, fold = TRUE)
#' c(1/2, 1/3, 1/4)
#'
#' # col means
#' allele_freqs(X, loci_on_cols = TRUE)
#' c(1/3, 1/4, 5/6)
#'
#' @export
allele_freqs <- function(
                         X,
                         loci_on_cols = FALSE,
                         fold = FALSE,
                         m_chunk_max = 1000, # optimal value not tested directly
                         subset_ind = NULL
                         ) {
    # behavior depends on class
    if ('BEDMatrix' %in% class(X)) {
        # extract data dimensions
        m_loci <- ncol(X)
        n_ind <- nrow(X)
        # allocate desired vector of allele frequencies
        p_anc_est <- vector('numeric', m_loci)
        
        # navigate chunks
        i_chunk <- 1 # start of first chunk (needed for matrix inputs only; as opposed to function inputs)
        while (TRUE) { # start an infinite loop, break inside as needed
            # this means all SNPs have been covered!
            if (i_chunk > m_loci)
                break
            
            # indexes to extract loci, and also so save to FstTs and FstBs vectors
            indexes_loci_chunk <- i_chunk : min(i_chunk + m_chunk_max - 1, m_loci) # range of SNPs to extract in this chunk

            # Only loci_on_cols==TRUE cases is supported here, this is only for BEDMatrix
            # DO NOT transpose for our usual setup (this is faster)
            Xi <- X[, indexes_loci_chunk, drop = FALSE]
            
            # apply subsetting to this chunk of matrix
            if ( !is.null( subset_ind ) )
                Xi <- Xi[ subset_ind, , drop = FALSE ]
            
            # compute and store the values we want!
            # because we didn't transpose, use colMeans here istead of rowMeans! 
            p_anc_est[ indexes_loci_chunk ] <- colMeans(Xi, na.rm = TRUE)/2

            # fold allele frequencies if requested
            # NOTE: perform on each chunk, further preserving memory
            if ( fold ) 
                p_anc_est[ indexes_loci_chunk ] <- fold_allele_freqs( p_anc_est[ indexes_loci_chunk ] )
            
            # update starting point for next chunk! (overshoots at the end, that's ok)
            i_chunk <- i_chunk + m_chunk_max
        }

        # set names to be SNP names
        # this is done by col/rowMeans when X is not BEDMatrix below, so let's match that here
        names(p_anc_est) <- colnames(X)
    } else if (is.matrix(X)) {
        # apply subsetting to whole matrix in this case
        if ( !is.null( subset_ind ) )
            X <- if ( loci_on_cols ) X[ subset_ind, , drop = FALSE ] else X[ , subset_ind, drop = FALSE ]
        
        # compute allele frequencies directly, all at once
        if (loci_on_cols) {
            p_anc_est <- colMeans(X, na.rm = TRUE)/2
        } else{
            p_anc_est <- rowMeans(X, na.rm = TRUE)/2
        }
        
        # either way, fold allele frequencies if requested
        if ( fold ) 
            p_anc_est <- fold_allele_freqs( p_anc_est )
    } else
        stop('Only matrix and BEDMatrix supported!  Instead, class was: ', toString( class(X) ) )
    
    # return the vector whichever way it was computed
    return( p_anc_est )
}

# internal function that folds allele frequencies that have already been calculated from genotypes or otherwise
fold_allele_freqs <- function( ps ) {
    if ( missing( ps ) )
        stop( 'Allele frequency vector `ps` is required!' )
    # identify loci to flip
    indexes <- ps > 0.5
    # actually flip
    ps[ indexes ] <- 1 - ps[ indexes ]
    # return modified vector
    return( ps )
}
