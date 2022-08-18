#' The model covariance matrix of the trait
#'
#' This function returns the expected covariance matrix of a trait vector simulated via `sim_trait`.
#' Below there are `n` individuals.
#'
#' @inheritParams sim_trait
#' @param kinship The `n`-by-`n` kinship matrix of the individuals.
#' These values should be scaled such that an outbred individual has 1/2 self-kinship, the parent-child relationship is 1/4, etc (which is half the values sometimes defined for kinship).
#'
#' @return The `n`-by-`n` trait covariance matrix equal to
#' `sigma_sq * ( herit * 2 * kinship + ( 1 - herit ) * I )`,
#' where `I` is the `n`-by-`n` identity matrix.
#'
#' @examples
#' # create a dummy kinship matrix
#' kinship <- matrix(
#'     data = c(
#'         0.6, 0.1, 0.0,
#'         0.1, 0.6, 0.1,
#'         0.0, 0.1, 0.6
#'     ),
#'     nrow = 3,
#'     byrow = TRUE
#' )
#' # covariance of simulated traits
#' V <- cov_trait(kinship = kinship, herit = 0.8)
#'
#' @seealso
#' [sim_trait()], [sim_trait_mvn()]
#' 
#' @export
cov_trait <- function(
                      kinship,
                      herit,
                      sigma_sq = 1,
                      labs = NULL,
                      labs_sigma_sq = NULL
                      ) {
    # construct V from kinship and herit

    # check for missing parameters
    if (missing(kinship))
        stop('`kinship` matrix is required!')
    if (missing(herit))
        stop('the heritability `herit` is required!')

    # other checks
    if (length(sigma_sq) != 1)
        stop('`sigma_sq` must be a scalar! (input has length ', length(sigma_sq), ')')
    if (length(herit) != 1)
        stop('`herit` must be a scalar! (input has length ', length(herit), ')')
    if (herit < 0)
        stop('`herit` must be non-negative!')
    if (herit > 1)
        stop('`herit` cannot be greater than 1!')
    if (sigma_sq <= 0)
        stop('`sigma_sq` must be positive!')
    if ( !isSymmetric( kinship ) )
        stop( '`kinship` must be a square, symmetric matrix!' )
    
    # if there are group effects, check those
    n_ind <- nrow( kinship )
    labs <- check_labs( labs, labs_sigma_sq, n_ind, herit )
    
    # identity matrix of same dimension as kinship
    I <- diag( nrow( kinship ) )
    # desired covariance matrix (except for overall scale)
    V <- 2 * herit * kinship + (1-herit) * I
    
    # if there are group effects, add them now
    if ( !is.null( labs ) ) {
        n_labs <- ncol( labs )
        for ( i in 1L : n_labs ) {
            # process level i
            labs_i <- labs[ , i ]
            labs_i_sigma_sq <- labs_sigma_sq[ i ]
            # unique groups, implicitly numbered by order of appearance
            groups_i <- unique( labs_i )
            # map individuals to their groups by index
            group_indexes_i <- match( labs_i, groups_i )
            # in this level, we have a simple block matrix reflecting shared group membership
            n_groups_i <- length( groups_i )
            Ci <- matrix( 0, n_ind, n_ind )
            for ( j in 1L : n_groups_i ) {
                Ci <- Ci + tcrossprod( group_indexes_i == j )
            }
            # multiply by variance factor, and add to running covariance sum
            V <- V + Ci * labs_i_sigma_sq
        }
    }
    
    # multiply by overall scale
    V <- V * sigma_sq
    
    # return
    return( V )
}
