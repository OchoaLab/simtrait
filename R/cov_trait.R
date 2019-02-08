#' The model covariance matrix of the trait
#'
#' This function returns the expected covariance matrix of a trait vector simulated via \code{sim_trait}.
#' Below let there be \eqn{n} individuals.
#'
#' @param kinship The \eqn{n \times n}{n-by-n} kinship matrix \eqn{\Phi} of the individuals.  This may be the true matrix of the genotype simulation or a good estimate of empirical data obtained via the package \code{popkin}.
#' @param herit The heritability \eqn{h^2} (proportion of trait variance due to genetics).
#'
#' @return The \eqn{n \times n}{n-by-n} trait covariance matrix equal to
#' \deqn{2 h^2 \Phi + (1-h^2) I,}
#' where \eqn{I} is an \eqn{n \times n}{n-by-n} identity matrix.
#'
#' @examples
#' # create a dummy kinship matrix
#' Phi <- matrix(
#'              data=c(0.6,0.1,0, 0.1,0.6,0.1, 0,0.1,0.6),
#'              nrow=3,
#'              byrow=TRUE
#'              )
#' # covariance of simulated traits
#' V <- cov_trait(kinship=Phi, herit=0.8)
#'
#' @export
cov_trait <- function(kinship, herit) {
    # construct V from kinship and herit
    # resulting overall sigma2=1 (could change later, for now not needed)
    I <- diag(rep.int(1, nrow(kinship))) # identity matrix of same dimension as kinship
    V <- 2 * herit * kinship + (1-herit) * I # return this directly!
}
