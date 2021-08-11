#' Simulate traits from a kinship matrix under the infinitesimal model
#'
#' Simulate matrix of trait replicates given a kinship matrix and model parameters (the desired heritability, total variance scale, and mean).
#' Although these traits have the covariance structure of genetic traits, and have heritabilities that can be estimated, they do not have causal loci (an association test against any locus should fail).
#' Below `n` is the number of individuals.
#'
#' @param rep The number of replicate traits to simulate.
#' Simulating all you need at once is more efficient than simulating each separately (the kinship matrix is eigendecomposed once per run, shared across replicates).
#' @param kinship The `n`-by-`n` kinship matrix of the individuals to simulate from.
#' @param herit The desired heritability (proportion of trait variance due to genetics).
#' @param mu The desired parametric mean value of the trait (default zero).
#' @param sigma_sq The desired parametric variance factor of the trait (default 1).
#' This factor corresponds to the variance of an outbred individual (see [cov_trait()]).
#' @param tol Tolerance factor for an internal test of positive semi-definiteness of the trait covariance matrix.
#' Procedure fails if any eigenvalues are smaller than `-tol` times the absolute value of the largest eigenvalue.
#' Increase this value only if you are getting errors but you're sure your covariance matrix (the output of [cov_trait()]) is positive semi-definite.
#'
#' @return A `rep`-by-`n` matrix containing the simulated traits along the rows, individuals along the columns.
#'
#' @examples
#' # create a dummy kinship matrix
#' # make sure it is positive definite!
#' kinship <- matrix(
#'     data = c(
#'         0.6, 0.1, 0.0,
#'         0.1, 0.5, 0.0,
#'         0.0, 0.0, 0.5
#'     ),
#'     nrow = 3
#' )
#' # draw simulated traits (matrix)
#' traits <- sim_trait_mvn( rep = 10, kinship = kinship, herit = 0.8 )
#' traits
#' 
#' @export
sim_trait_mvn <- function(
                          rep,
                          kinship,
                          herit,
                          mu = 0,
                          sigma_sq = 1,
                          tol = 1e-06
                          ) {
    # check for missing parameters
    if (missing(kinship))
        stop('`kinship` is required!')
    if (missing(herit))
        stop('the heritability `herit` is required!')
    
    # other checks
    if (length(mu) != 1)
        stop('`mu` must be a scalar! (input has length ', length(mu), ')')
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

    # construct the total covariance matrix of the trait
    V <- cov_trait(
        kinship = kinship,
        herit = herit,
        sigma_sq = sigma_sq
    )

    # calculate matrix square root, such that
    # V == tcrossprod( V_sqrt )
    V_sqrt <- sqrt_matrix( V, tol = tol )

    # simulate large matrix with IID standard normal values
    n_ind <- nrow( kinship )
    traits <- matrix(
        stats::rnorm( rep * n_ind ),
        nrow = rep,
        ncol = n_ind
    )

    # apply the affine transformation to get desired mean and covariance structure
    traits <- mu + tcrossprod( traits, V_sqrt )
    
    # return traits matrix
    return( traits )
}
