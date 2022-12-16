# used internally by: sim_trait, cov_trait
# always run because it checks heritability, and calculates residual variance correctly and consistently!
check_herit_labs <- function(
                             herit,
                             labs,
                             labs_sigma_sq,
                             n_ind
                             ) {
    # always check herit
    if ( length( herit ) != 1)
        stop( '`herit` must be a scalar! (input has length ', length( herit ), ')' )
    if ( herit < 0 )
        stop( '`herit` must be non-negative!' )
    if ( herit > 1 )
        stop( '`herit` cannot be greater than 1!' )
    
    # residual variance when there are no group effects
    # the previous checks guarantee this is non-negative and smaller or equal than 1
    sigma_sq_residual <- 1 - herit
    
    # if there are group effects, check those
    if ( !is.null( labs ) ) {
        if ( is.null( labs_sigma_sq ) )
            stop( '`labs_sigma_sq` must be non-NULL if `labs` are non-NULL!' )
        # better to set up as matrix internally
        # individuals along rows (column matrix for a single level), matching plot_popkin
        if ( !is.matrix( labs ) )
            labs <- cbind( labs )
        # check dimensions now
        if ( nrow( labs ) != n_ind )
            stop( 'Numbers of individuals in `labs` (length of vector or rows of matrix: ', nrow( labs ), ') must match that of genotype data (', n_ind, ')!' )
        n_labs <- ncol( labs )
        if ( length( labs_sigma_sq ) != n_labs )
            stop( 'The length of `labs_sigma_sq` (', length( labs_sigma_sq ), ') must equal the number of columns of `labs` (', n_labs, ')!' )
        # now that all dimensions make sense, make sure variance components are in range and their sum doesn't exceed 1
        if ( any( labs_sigma_sq < 0 ) )
            stop( '`labs_sigma_sq` must be non-negative!')
        if ( any( labs_sigma_sq > 1 ) )
            stop( '`labs_sigma_sq` cannot be greater than 1!')
        # calculate residual variance correctly for these cases
        sigma_sq_residual <- 1 - herit - sum( labs_sigma_sq )
        if ( sigma_sq_residual < 0 )
            stop( 'The sum of `labs_sigma_sq` (', sum( labs_sigma_sq ), ') must not exceed `1 - herit` (', 1 - herit, ')!' )
    }

    # `labs` was edited to always be a matrix, so return it!
    # (`labs_sigma_sq` is unchanged)
    return(
        list(
            labs = labs,
            sigma_sq_residual = sigma_sq_residual
        )
    )
}
