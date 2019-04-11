#' The model covariance matrix of the trait
#'
#' This function returns the expected covariance matrix of a trait vector simulated via \code{sim_trait}.
#' Below let there be \eqn{n} individuals.
#'
#' @param kinship The \eqn{n \times n}{n-by-n} kinship matrix \eqn{\Phi} of the individuals.  This may be the true matrix of the genotype simulation or a good estimate of empirical data obtained via the package \code{popkin}.
#' @param herit The heritability \eqn{h^2} (proportion of trait variance due to genetics).
#' @param sigma_sq Overall variance multiplicative factor \eqn{\sigma^2} (default 1).  This factor corresponds to the variance of an outbred individual.
#'
#' @return The \eqn{n \times n}{n-by-n} trait covariance matrix equal to
#' \deqn{\sigma^2 ( 2 h^2 \Phi + (1-h^2) I ),}
#' where \eqn{I} is an \eqn{n \times n}{n-by-n} identity matrix.
#'
#' @examples
#' # create a dummy kinship matrix
#' kinship <- matrix(
#'              data = c(0.6,0.1,0, 0.1,0.6,0.1, 0,0.1,0.6),
#'              nrow = 3,
#'              byrow = TRUE
#'              )
#' # covariance of simulated traits
#' V <- cov_trait(kinship = kinship, herit = 0.8)
#'
#' @export
cov_trait <- function(kinship, herit, sigma_sq = 1) {
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

    # identity matrix of same dimension as kinship
    I <- diag(rep.int(1, nrow(kinship)))
    # desired covariance matrix (except for overall scale)
    V <- 2 * herit * kinship + (1-herit) * I
    # multiply by scale and return
    V * sigma_sq
}
