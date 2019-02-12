#' Simulate a complex trait from genotypes
#'
#' Simulate a complex trait \eqn{y} given a SNP genotype matrix and model parameters (the desired heritability and the true ancestral allele frequencies used to generate the genotypes).  Users can choose the number of causal loci and minimum marginal allele frequency requirements for the causal loci.  The code selects random loci to be causal, draws random Normal effect sizes for these loci (scaled appropriately) and random independent non-genetic effects.
#' Below let there be \eqn{m} loci and \eqn{n} individuals.
#'
#' In order to center and scale the trait and locus effect size vector correctly to the desired parameters (mean, variance factor, and heritability), the parametric ancestral allele frequencies (\code{p_anc}) must be known.
#' This is necessary since in the context of Heritability the genotypes are themselves random variables (with means given by \code{p_anc} and a covariance structure given by \code{p_anc} and the kinship matrix), so their parameters must be taken into account as well.
#' If \code{p_anc} are indeed known, which they are when the genotypes were simulated, then the trait will have the specified mean and covariance matrix in agreement with \code{\link{cov_trait}}.
#'
#' If the desire is to simulate a trait using real genotypes, where \code{p_anc} is unknown, an imperfect compromise is possible as long as the kinship matrix (\code{kinship}) is known, and still the desired covariance structure is not exactly achieved.
#' The trait is centered using its arithmetic mean, which overfits, which on average results in the desired mean but reduces the overall variance of the trait.
#' The scaling is also performed using sample \code{p_anc} estimates, which causes a bias parametrized by \code{kinship}, so the bias can be partially adjusted for using \code{kinship} (see the \code{simtrait} package vignette).
#' While the true \code{kinship} is unknown in real data, it can be estimated accurately using the \code{popkin} package!
#' 
#' @param X The \eqn{m \times n}{m-by-n} genotype matrix.  This is a numeric matrix consisting of reference allele counts (in \code{c(0,1,2,NA)} for a diploid organism).
#' @param m_causal The number of causal loci desired.
#' @param herit The desired heritability (proportion of trait variance due to genetics).
#' @param p_anc The length-\eqn{m} vector of true ancestral allele frequencies.  Recommended way to adjust the simulated trait to achieve the desired heritability and covariance structure.  Either this or \code{kinship} must be specified.
#' @param kinship The \eqn{n \times n}{n-by-n} kinship matrix of the individuals in the data.  This offers an alternative way to adjust the simulated parameters parameters to partially achieve the desired covariance structure, since \code{p_anc} is only known for simulated data.  Either this or \code{p_anc} must be specified.
#' @param mu The desired parametric mean value of the trait (default zero).  The sample mean of the trait will not be exactly zero, but instead have an expectation of \code{mu} (with potentially large variance depending on the kinship matrix and the heritability).
#' @param sigmaSq The desired parametric variance factor of the trait (default 1).  This factor corresponds to the variance of an outbred individual (see \code{\link{cov_trait}}).
#' @param maf_cut The optional minimum allele frequency threshold (default 5\%).  This prevents rare alleles from being causal in the simulation.  Note that this threshold is applied to the sample allele frequencies and not their true parametric values (\code{p_anc}), even if these are available.
#'
#' @return A list containing the simulated trait \code{y} (length \eqn{n}), the vector of causal locus indeces \code{i} (length \eqn{m_causal}), and the locus effect size vector \code{beta} (length \eqn{m_causal}) at the causal loci.  However, if \code{herit==0} then \code{i} and \code{beta} will have zero length regardless of \code{m_causal}.
#'
#' @examples
#' # construct a dummy genotype matrix
#' X <- matrix(
#'            data=c(0,1,2,1,2,1,0,0,1),
#'            nrow=3,
#'            byrow=TRUE
#'            )
#' # made up ancestral allele frequency vector for example
#' p_anc <- c(0.5, 0.6, 0.2)
#' # create simulated trait and associated data
#' obj <- sim_trait(X=X, m_causal=2, herit=0.8, p_anc=p_anc)
#' # trait vector
#' obj$y
#' # randomly-picked causal locus indeces
#' obj$i
#' # locus effect size vector
#' obj$beta
#' 
#' @export
# TODO:
# - make it work with BEDMatrix?
sim_trait <- function(X, m_causal, herit, p_anc, kinship, mu=0, sigmaSq=1, maf_cut=0.05) {
    # check for missing parameters
    if (missing(X)) stop('Fatal: genotype matrix `X` must be specified (no default value)')
    if (missing(m_causal)) stop('Fatal: the number of causal loci `m_causal` must be specified (no default value)')
    if (missing(herit)) stop('Fatal: the heritability `herit` must be specified (no default value)')
    if (missing(p_anc) && missing(kinship)) stop('Fatal: either the true ancestral allele frequency vector `p_anc` or `kinship` must be specified (no default values)')
    
    # other checks
    if (length(mu) != 1) stop('Fatal: `mu` must be a scalar! (input has length ', length(mu), ')')
    if (length(sigmaSq) != 1) stop('Fatal: `sigmaSq` must be a scalar! (input has length ', length(sigmaSq), ')')
    if (length(herit) != 1) stop('Fatal: `herit` must be a scalar! (input has length ', length(herit), ')')
    if (herit < 0) stop('Fatal: `herit` must be non-negative!')
    if (herit > 1) stop('Fatal: `herit` cannot be greater than 1!')
    if (sigmaSq <= 0) stop('Fatal: `sigmaSq` must be positive!')

    # data dimensions (will be transposed if X is a BEDMatrix object)
    m <- nrow(X)
    n <- ncol(X)

    if (herit == 0) {
        # lots of work can be avoided in this edge case
        # the index and coefficients vectors are empty
        i <- c()
        beta <- c()
        G <- 0 # construct a trivial genotype effect of zero (becomes vector automatically later)
    } else {
        
        # compute marginal allele frequencies
        p_anc_hat <- rowMeans(X, na.rm=TRUE)/2
        
        # select random SNPs! this performs the magic...
        # also runs additional checks
        i <- select_loci(p_anc_hat, m_causal, maf_cut)
        
        # draw random SNP coefficients for selected loci
        beta <- stats::rnorm(m_causal, 0, 1)
        
        ###############
        ### KINSHIP ###
        ###############
        
        if (!missing(kinship)) {
            # precompute some things when this is present
            mean_kinship <- mean(kinship)
        }
        
        #############
        ### SCALE ###
        #############
        
        # to scale beta to give correct heritability, we need to estimate the pq = p(1-p) vector
        # calculate pq = p_anc * (1 - p_anc) in one of two ways
        # (NOTE only causal loci are needed, subset for speed and memory)
        if ( !missing(p_anc) ) { # this takes precedence, it should be most accurate
            # direct calculation
            pq <- p_anc[i] * (1 - p_anc[i])
        } else if ( !missing(kinship) ) {
            # indirect, infer from genotypes and mean kinship
            # recall E[ p_anc_hat * (1 - p_anc_hat) ] = pq * (1 - mean_kinship), so we solve for pq:
            pq <- p_anc_hat[i] * (1 - p_anc_hat[i]) / (1 - mean_kinship)
        } else {
            # a redundant check (if this were so, it should have died earlier)
            stop('Fatal: either `p_anc` or `kinship` must be specified!')
        }
        
        # the genetic variance (in a vector of per-loci terms) is
        varXB <- 2 * pq * beta^2
        # that should equal:
        # sum( varXB ) = herit * sigma^2
        # Let's solve for sigma and divide it out:
        sigma <- sqrt( sum( varXB ) / herit )
        # adjust betas so final sigma=sqrt(sigmaSq) as desired!
        beta <- beta * sqrt(sigmaSq) / sigma # scale by standard deviations
        
        # construct genotype signal
        Xi <- X[i,] - 1 # the subset of causal data, "centered" so heterozygotes are zero
        if (any(is.na(Xi))) {
            # if any of the causal loci are missing, let's treat them as zeroes
            # this isn't perfect but we must do something to apply this to real data
            Xi[is.na(Xi)] <- 0
        }
        G <- drop( beta %*% Xi ) # this is a vector
        # NOTE by construction:
        # Cov(G) = 2 * herit * Phi
        
        ##############
        ### CENTER ###
        ##############

        # calculate the mean of the genetic effect
        if ( !missing(p_anc) ) {
            # parametric solution
            muXB <- drop( beta %*% (2*p_anc[i]-1) )
        } else {
            muXB <- 0 # an experiment, assume this is nearly zero if p_anc is symmetric around 1/2 (and 
            ## # all other cases, whether kinship is present or otherwise
            ## # remove sample mean
            ## muXB <- mean(G)
            ## # the prev. step dramatically reduces the covariance of the data
            ## # to make things look a tiny bit better, add back some MVN noise!
            ## # NOTE: this particular version adds a single constant to the muXB vector, which corrects for the mean difference between covariance matrices, but does not correct for the "Standard Kinship distortions".
            ## muXB <- muXB + stats::rnorm(1, 0, sqrt(mean_kinship * 2 * sum(beta^2 * pq)) ) # is this right?
        }
        # in all cases:
        # - remove the mean from the genotypes (muXB)
        # - add the desired mean
        G <- G - muXB + mu
    }

    if (herit == 1) {
        E <- 0 # in this edge case there is no "noise", just genotype effects
    } else {
        # draw noise
        E <- stats::rnorm(n, 0, (1 - herit) * sigmaSq ) # noise has mean zero but variance ((1-herit) * sigmaSq)
        # NOTE by construction:
        # Cov(E) = (1-herit) * sigmaSq * I
    }

    # lastly, here's the trait:
    y <- G + E

    # return all these things
    list(y=y, i=i, beta=beta)
}
