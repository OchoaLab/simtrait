#' Simulate a complex trait from genotypes
#'
#' Simulate a complex trait \eqn{y} given a SNP genotype matrix and model parameters (the desired heritability and either the true ancestral allele frequencies used to generate the genotypes or the mean kinship).  Users can choose the number of causal loci and minimum marginal allele frequency requirements for the causal loci.  The code selects random loci to be causal, draws random Normal effect sizes for these loci (scaled appropriately) and random independent non-genetic effects.
#' Below let there be \eqn{m} loci and \eqn{n} individuals.
#'
#' In order to scale the locus effect size vector correctly to 
#' 
#' @param X The \eqn{m \times n}{m-by-n} genotype matrix.  This is a numeric matrix consisting of reference allele counts (in \code{c(0,1,2,NA)} for a diploid organism).
#' @param m_causal The number of causal loci desired.
#' @param herit The desired heritability (proportion of trait variance due to genetics).
#' @param p_anc The length-\eqn{m} vector of true ancestral allele frequencies.  Recommended way to scale output parameters to achieve the desired heritability.  Either this or \code{mean_kinship} must be specified.
#' @param mean_kinship The mean kinship coefficient across all individuals in the data.  This offers an alternative way to scale output parameters to achieve the desired heritability, as \code{p_anc} is only known for simulated data.  Either this or \code{p_anc} must be specified.
#' @param maf_cut The optional minimum allele frequency threshold (default 5\%).  This prevents rare alleles from being causal in the simulation.
#'
#' @return A list containing the simulated trait \code{y} (length \eqn{n}), the vector of causal locus indeces \code{i} (length \eqn{m_causal}), and the locus effect size vector \code{beta} (length \eqn{m}, equals zero at every non-causal locus).
#'
#' @examples
#' # construct a dummy genotype matrix
#' X <- matrix(
#'            data=c(0,1,2,1,2,1,0,0,1),
#'            nrow=3,
#'            byrow=TRUE
#'            )
#' # made up mean kinship for example
#' mean_kinship <- 0.1
#' # create simulated trait and associated data
#' obj <- sim_trait(X=X, m_causal=1, herit=0.8, mean_kinship=mean_kinship)
#' # trait vector
#' obj$y
#' # randomly-picked causal locus index
#' obj$i
#' # locus effect size vector
#' obj$beta
#' 
#' @export
# TODO: make it work with BEDMatrix?
sim_trait <- function(X, m_causal, herit, p_anc, mean_kinship, maf_cut=0.05) {
    # check for missing parameters
    if (missing(X)) stop('Fatal: genotype matrix `X` must be specified (no default value)')
    if (missing(m_causal)) stop('Fatal: the number of causal loci `m_causal` must be specified (no default value)')
    if (missing(herit)) stop('Fatal: the heritability `herit` must be specified (no default value)')
    if (missing(p_anc) && missing(mean_kinship)) stop('Fatal: either the true ancestral allele frequency vector `p_anc` or `mean_kinship` must be specified (no default values)')

    # other checks
    if (length(herit) > 1) stop('Fatal: `herit` must be a scalar! (input has length ', length(herit), ')')
    if (herit < 0) stop('Fatal: `herit` must be positive!')
    if (herit > 1) stop('Fatal: `herit` cannot be greater than 1!')

    # data dimensions (will be transposed if X is a BEDMatrix object)
    m <- nrow(X)
    
    # compute marginal allele frequencies
    p_anc_hat <- rowMeans(X, na.rm=TRUE)/2
    
    # select random SNPs! this performs the magic...
    # also runs additional checks
    i <- select_loci(p_anc_hat, m_causal, maf_cut)

    # construct beta coefficients that are mostly zero (sparse)
    beta <- rep.int(0, m)
    # draw random SNP coefficients for selected loci
    beta[i] <- rnorm(m_causal, 0, 1)

    # to scale beta to give correct heritability, we need to estimate the pq = p(1-p) vector
    # calculate pq = p_anc * (1 - p_anc) in one of two ways
    if ( !missing(p_anc) ) { # this takes precedence, it should be most accurate
        # direct calculation
        pq <- p_anc * (1 - p_anc)
    } else if ( !missing(mean_kinship) ) {
        # indirect, infer from genotypes and mean kinship
        # recall E[ p_anc_hat * (1 - p_anc_hat) ] = pq * (1 - mean_kinship), so we solve for pq:
        pq <- p_anc_hat * (1 - p_anc_hat) / (1 - mean_kinship)
    } else {
        # a redundant check (if this were so it should have died earlier)
        stop('Fatal: either `p_anc` or `mean_kinship` must be specified!')
    }
    
    # the genetic variance (in a vector of per-loci terms) is
    varXB <- 2 * pq * beta^2
    # that should equal:
    # sum( varXB ) = herit * sigma^2
    # Let's solve for sigma and divide it out:
    # (NOTE only non-zero terms are expicitly summed over for greater speed and to avoid numerical issues)
    sigma <- sqrt( sum( varXB[i] ) / herit )
    # adjust betas so final sigma=1
    beta[i] <- beta[i] / sigma # divide by standard deviation
    
    # construct genotype signal
    Xi <- X[i,] # the subset of causal data
    if (any(is.na(Xi))) {
        # if any of the causal loci are missing, let's treat them as zeroes
        # this isn't perfect but we must do something to apply this to real data
        Xi[is.na(Xi)] <- 0
    }
    # note: though I could have done (beta %*% X) (no subsetting), this is more efficient and numerically accurate, and allows me to easily replace NAs with zeroes
    G <- drop( beta[i] %*% Xi ) # this is a vector
    # NOTE by construction:
    # Cov(G) = 2 * herit * Phi

    # draw noise
    E <- rnorm(G, 0, 1 - herit) # noise has mean zero but variance (1-herit), length matches that of G (genotypes)
    # NOTE by construction:
    # Cov(E) = (1-herit) * I

    # lastly, here's the trait:
    y <- G + E

    # return all these things
    list(y=y, i=i, beta=beta)
}
