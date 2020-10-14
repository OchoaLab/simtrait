#' Compute locus allele frequencies
#'
#' On a regular matrix, this is essentially a wrapper for colMeans or rowMeans depending on \code{loci_on_cols}.
#' On a BEDMatrix object, the locus allele frequencies are computed keeping memory usage low.
#'
#' @param X The genotype matrix (regular R matrix or BEDMatrix object).
#' Missing values are ignored in averages.
#' @param loci_on_cols If \code{TRUE}, \eqn{X} has loci on columns and individuals on rows; if false (the default), loci are on rows and individuals on columns.
#' If \eqn{X} is a BEDMatrix object, columns are averaged to yield locus allele frequencies (regardless of the value of \code{loci_on_cols}).
#' @param m_chunk_max BEDMatrix-specific, sets the maximum number of loci to process at the time.
#' If memory usage is excessive, set to a lower value than default (expected only for extremely large numbers of individuals).
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
#' # col means
#' allele_freqs(X, loci_on_cols = TRUE)
#' c(1/3, 1/4, 5/6)
#'
#' @export
allele_freqs <- function(
                         X,
                         loci_on_cols = FALSE,
                         m_chunk_max = 1000 # optimal value not tested directly
                         ) {
    # behavior depends on class
    if ('BEDMatrix' %in% class(X)) {
        # extract data dimensions
        m_loci <- ncol(X)
        n_ind <- nrow(X)
        # allocate desired vector of allele frequencies
        p_anc_hat <- vector('numeric', m_loci)
        
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
            
            # compute and store the values we want!
            # because we didn't transpose, use colMeans here istead of rowMeans! 
            p_anc_hat[ indexes_loci_chunk ] <- colMeans(Xi, na.rm = TRUE)/2
            
            # update starting point for next chunk! (overshoots at the end, that's ok)
            i_chunk <- i_chunk + m_chunk_max
        }

        # set names to be SNP names
        # this is done by col/rowMeans when X is not BEDMatrix below, so let's match that here
        names(p_anc_hat) <- colnames(X)
    } else if (is.matrix(X)) {
        if (loci_on_cols) {
            p_anc_hat <- colMeans(X, na.rm = TRUE)/2
        } else{
            p_anc_hat <- rowMeans(X, na.rm = TRUE)/2
        }
    } else
        stop('Only matrix and BEDMatrix supported!  Instead, class was: ', toString( class(X) ) )
    
    # return the vector whichever way it was computed
    return( p_anc_hat )
}
