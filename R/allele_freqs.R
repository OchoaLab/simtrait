#' Compute locus allele frequencies
#'
#' On a regular matrix, this is essentially a wrapper for colMeans or rowMeans depending on \code{loci_on_cols}.
#' On a BEDMatrix object, the locus allele frequencies are computed keeping memory usage low.
#'
#' @param X The genotype matrix (regular R matrix or BEDMatrix object).
#' Missing values are ignored in averages.
#' @param loci_on_cols If \code{TRUE}, \eqn{X} has loci on columns and individuals on rows; if false (the default), loci are on rows and individuals on columns.
#' If \eqn{X} is a BEDMatrix object, columns are averaged to yield locus allele frequencies (regardless of the value of \code{loci_on_cols}).
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
allele_freqs <- function(X, loci_on_cols = FALSE) {
    # behavior depends on class
    if (class(X) == 'BEDMatrix') {
        # extract the number of loci
        m_loci <- ncol(X)
        # allocate desired vector of allele frequencies
        p_anc_hat <- vector('numeric', m_loci)
        # navigate loci
        for (i in 1 : m_loci) {
            # compute mean one locus at the time
            p_anc_hat[i] <- mean(X[, i], na.rm = TRUE)/2
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
        stop('Only matrix and BEDMatrix supported!  class: ', class(X))
    
    # return the vector whichever way it was computed
    return( p_anc_hat )
}
