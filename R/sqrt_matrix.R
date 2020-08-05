# internal function
sqrt_matrix <- function( V, tol = 1e-06 ){
    if ( missing( V ) )
        stop('`V` is missing!')
    # perform eigendecomposition (slowest step for large matrices)
    obj <- eigen( V )
    eigenvals <- obj$values
    eigenvecs <- obj$vectors
    # make sure matrix square root is defined properly
    if ( any( eigenvals < -tol * abs( eigenvals[ 1L ] ) ) )
        stop("'V' is not positive definite")
    # replace any negative values with zero
    eigenvals[ eigenvals < 0 ] <- 0
    # this is the desired matrix square root
    V_sqrt <- eigenvecs %*% diag( sqrt( eigenvals ) )
    return( V_sqrt )
}
