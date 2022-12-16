#' Simulate traits from a kinship matrix under the infinitesimal model
#'
#' Simulate matrix of trait replicates given a kinship matrix and model parameters (the desired heritability, group effects, total variance scale, and mean).
#' Although these traits have the covariance structure of genetic traits, and have heritabilities that can be estimated, they do not have causal loci (an association test against any locus should fail).
#' Below `n` is the number of individuals.
#'
#' @inheritParams sim_trait
#' @param rep The number of replicate traits to simulate.
#' Simulating all you need at once is more efficient than simulating each separately (the kinship matrix is eigendecomposed once per run, shared across replicates).
#' @param kinship The `n`-by-`n` kinship matrix of the individuals to simulate from.
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
#' @seealso
#' [cov_trait()], [sim_trait()]
#' 
#' @export
sim_trait_mvn <- function(
                          rep,
                          kinship,
                          herit,
                          mu = 0,
                          sigma_sq = 1,
                          labs = NULL,
                          labs_sigma_sq = NULL,
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

    # construct the total covariance matrix of the trait
    # NOTE: `herit`, `sigma_sq`, and `labs` get checked inside this function, won't check beforehand
    V <- cov_trait(
        kinship = kinship,
        herit = herit,
        sigma_sq = sigma_sq,
        labs = labs,
        labs_sigma_sq = labs_sigma_sq
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
