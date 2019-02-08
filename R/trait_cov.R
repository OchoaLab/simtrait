# the true covariance matrix of the trait
trait_cov <- function(Phi, herit) {
    # construct V from Phi and herit
    # resulting overall sigma2=1 (could change later, for now not needed)
    I <- diag(rep.int(1, nrow(Phi))) # identity matrix of same dimension as Phi
    V <- 2 * herit * Phi + (1-herit) * I # return this directly!
}
