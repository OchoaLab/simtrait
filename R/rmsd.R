#' Root mean square deviation
#'
#' Calculates the euclidean distance between two vectors `x` and `y` divided by the square root of the lengths of the vectors.
#' `NA` values are ignored by default when calculating the mean squares (so the denominator is the number of non-`NA` differences).
#'
#' @param x The first vector to compare (required).
#' @param y The second vector to compare (required).
#' Lengths of `x` and `y` must be equal.
#' @param na.rm If `TRUE` (default), `NA` values are removed before calculating the mean square difference.
#' If `FALSE`, any missing values in either `x` or `y` result in `NA` returned.
#' Passed to [mean()], see that for more info.
#'
#' @return the square root of the mean square difference between `x` and `y`, after removing `NA` comparisons (cases where either is `NA`).
#'
#' @examples
#' x <- rnorm(10)
#' y <- rnorm(10)
#' rmsd( x, y )
#' 
#' @export
rmsd <- function( x, y, na.rm = TRUE ) {
    if ( missing( x ) || is.null( x ) )
        stop( '`x` is required!' )
    if ( missing( y ) || is.null( y ) )
        stop( '`y` is required!' )
    if ( length( x ) != length( y ) )
        stop( 'Lengths of `x` (', length( x ), ') and `y` (', length( y ), ') do not agree!' )
    return( sqrt( mean( ( x - y )^2, na.rm = na.rm ) ) )
}
