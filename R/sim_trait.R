#' Simulate a complex trait from genotypes
#'
#' Simulate a complex trait given a SNP genotype matrix and model parameters (the desired heritability, and either the true ancestral allele frequencies used to generate the genotypes, or the kinship matrix of the individuals).
#' Users can choose the number of causal loci and minimum marginal allele frequency requirements for the causal loci.
#' The code selects random loci to be causal, draws random Normal effect sizes for these loci (scaled appropriately) and random independent non-genetic effects.
#' Suppose there are `m` loci and `n` individuals.
#'
#' In order to center and scale the trait and locus effect size vector correctly to the desired parameters (mean, variance factor, and heritability), the parametric ancestral allele frequencies (`p_anc`) must be known.
#' This is necessary since in the context of Heritability the genotypes are themselves random variables (with means given by `p_anc` and a covariance structure given by `p_anc` and the kinship matrix), so the parameters of the genotypes must be taken into account.
#' If `p_anc` are indeed known (true for simulated genotypes), then the trait will have the specified mean and covariance matrix in agreement with `\link{cov_trait}`.
#'
#' If the desire is to simulate a trait using real genotypes, where `p_anc` is unknown, a compromise that works well in practice is possible if the `kinship` matrix is known (see package vignette).
#' The kinship matrix can be estimated accurately using the `popkin` package!
#' 
#' @param X The `m`-by-`n` genotype matrix (if `loci_on_cols = FALSE`, transposed otherwise), or a BEDMatrix object.
#' This is a numeric matrix consisting of reference allele counts (in `c(0, 1, 2, NA)` for a diploid organism).
#' @param m_causal The number of causal loci desired.
#' @param herit The desired heritability (proportion of trait variance due to genetics).
#' @param p_anc The length-`m` vector of true ancestral allele frequencies.
#' Recommended way to adjust the simulated trait to achieve the desired heritability and covariance structure.
#' Either this or `kinship` must be specified.
#' @param kinship The mean kinship value of the individuals in the data.
#' The `n`-by-`n` kinship matrix of the individuals in the data is also accepted.
#' This offers an alternative way to adjust the simulated parameters parameters to achieve the desired covariance structure for real genotypes, since `p_anc` is only known for simulated data.
#' Either this or `p_anc` must be specified.
#' @param mu The desired parametric mean value of the trait (default zero).
#' The sample mean of the trait will not be exactly zero, but instead have an expectation of `mu` (with potentially large variance depending on the kinship matrix and the heritability).
#' @param sigma_sq The desired parametric variance factor of the trait (default 1).
#' This factor corresponds to the variance of an outbred individual (see `\link{cov_trait}`).
#' @param maf_cut The optional minimum allele frequency threshold (default `NA`, no threshold).
#' This prevents rare alleles from being causal in the simulation.
#' Note that this threshold is applied to the sample allele frequencies and not their true parametric values (`p_anc`), even if these are available.
#' @param loci_on_cols If `TRUE`, `X` has loci on columns and individuals on rows; if `FALSE` (the default), loci are on rows and individuals on columns.
#' If `X` is a BEDMatrix object, loci are taken to be on the columns (regardless of the value of `loci_on_cols`).
#' @param m_chunk_max BEDMatrix-specific, sets the maximum number of loci to process at the time.
#' If memory usage is excessive, set to a lower value than default (expected only for extremely large numbers of individuals).
#' @param const_herit_loci If `TRUE`, causal coefficients are inversely proportional to the square root of `p*(1-p)`, where `p` is the ancestral allele frequency, which ensures equal per-causal-locus contribution to trait variance.
#' Causal coefficient signs (+/-) are drawn randomly with equal probability.
#' If `FALSE` (the default), draws causal coefficients randomly from a standard normal distribution, rescaled to result in the desired heritability.
#'
#' @return A named list containing:
#'
#' - `trait`: length-`n` vector of the simulated trait
#' - `causal_indexes`: length-`m_causal` vector of causal locus indexes
#' - `causal_coeffs`: length-`m_causal` vector of locus effect sizes at the causal loci
#' 
#' However, if `herit = 0` then `causal_indexes` and `causal_coeffs` will have zero length regardless of `m_causal`.
#'
#' @examples
#' # construct a dummy genotype matrix
#' X <- matrix(
#'            data = c(0,1,2,1,2,1,0,0,1),
#'            nrow = 3,
#'            byrow = TRUE
#'            )
#' # made up ancestral allele frequency vector for example
#' p_anc <- c(0.5, 0.6, 0.2)
#' 
#' # create simulated trait and associated data
#' obj <- sim_trait(X = X, m_causal = 2, herit = 0.8, p_anc = p_anc)
#' 
#' # trait vector
#' obj$trait
#' # randomly-picked causal locus indexes
#' obj$causal_indexes
#' # locus effect size vector
#' obj$causal_coeffs
#' 
#' @export
sim_trait <- function(
                      X,
                      m_causal,
                      herit,
                      p_anc,
                      kinship,
                      mu = 0,
                      sigma_sq = 1,
                      maf_cut = NA,
                      loci_on_cols = FALSE,
                      m_chunk_max = 1000,
                      const_herit_loci = FALSE
                      ) {
    # check for missing parameters
    if (missing(X))
        stop('genotype matrix `X` is required!')
    if (missing(m_causal))
        stop('the number of causal loci `m_causal` is required!')
    if (missing(herit))
        stop('the heritability `herit` is required!')
    if (missing(p_anc) && missing(kinship))
        stop('either the true ancestral allele frequency vector `p_anc` or `kinship` are required!')
    
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

    # simplifies subsetting downstream
    if ('BEDMatrix' %in% class(X))
        loci_on_cols <- TRUE

    # get m_loci to check m_causal before more heavy computations
    m_loci <- if (loci_on_cols) ncol(X) else nrow(X)
    if (m_causal > m_loci)
        stop('m_causal (', m_causal, ') exceeds the number of loci (', m_loci, ')!')
    
    if (herit == 0) {
        # lots of work can be avoided in this edge case
        # the index and coefficients vectors are empty
        causal_indexes <- c()
        causal_coeffs <- c()
        G <- 0 # construct a trivial genotype effect of zero (becomes vector automatically later)
    } else {
        
        ###################################
        ### MARGINAL ALLELE FREQUENCIES ###
        ###################################

        # may not need this
        p_anc_hat <- NULL
        # these are needed here to select loci, if there's a MAF threshold
        # (they are also needed in real datasets, but if that's the only need then let's wait until we've subset the genotype matrix)
        if ( !is.na( maf_cut ) ) {
            # compute marginal allele frequencies
            p_anc_hat <- allele_freqs(
                X,
                loci_on_cols = loci_on_cols,
                m_chunk_max = m_chunk_max
            )
        }
        
        ###################
        ### CAUSAL LOCI ###
        ###################
        
        # select random SNPs! this performs the magic...
        # also runs additional checks
        causal_indexes <- select_loci(
            m_causal = m_causal,
            m_loci = m_loci,
            maf = p_anc_hat, # NULL if is.na( maf_cut )
            maf_cut = maf_cut # may be NA
        )

        # subset data to consider causal loci only
        if ( !is.null( p_anc_hat ) ) # if we had this already
            p_anc_hat <- p_anc_hat[ causal_indexes ]
        if ( !missing( p_anc ) )
            p_anc <- p_anc[ causal_indexes ] # subset if available
        # the subset of causal data
        # (drop = FALSE for keeping as a matrix even if m_causal == 1)
        if (loci_on_cols) {
            # also transpose for consistent behavior downstream
            X <- t( X[, causal_indexes, drop = FALSE] )
        } else{
            X <- X[causal_indexes, , drop = FALSE]
        }
        
        ###################################
        ### MARGINAL ALLELE FREQUENCIES ###
        ###################################

        # do here if we don't already have them and need them
        # this will be faster now, if done on the subset of causal genotypes only

        # these are used to select loci, or to simulate from real genotypes
        if ( missing( p_anc ) && is.null( p_anc_hat ) ) {
            # compute marginal allele frequencies
            p_anc_hat <- allele_freqs(
                X,
#                loci_on_cols = loci_on_cols, # now that genotypes are extracted, they are in default orientation (no need for passing this, plus passing it messes things up)
                m_chunk_max = m_chunk_max
            )
        }
        
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
        
        # to scale causal_coeffs to give correct heritability, we need to estimate the pq = p(1-p) vector
        # calculate pq = p_anc * (1 - p_anc) in one of two ways
        if ( !missing(p_anc) ) { # this takes precedence, it should be most accurate
            # direct calculation
            pq <- p_anc * (1 - p_anc)
        } else if ( !missing(kinship) ) {
            # indirect, infer from genotypes and mean kinship
            # recall E[ p_anc_hat * (1 - p_anc_hat) ] = pq * (1 - mean_kinship), so we solve for pq:
            pq <- p_anc_hat * (1 - p_anc_hat) / (1 - mean_kinship)
        } else {
            # a redundant check (if this were so, it should have died earlier)
            stop('either `p_anc` or `kinship` must be specified!')
        }
        
        # construct SNP coefficients for selected loci
        if ( const_herit_loci ) {
            # make them inverse to the pq
            # in this case no scale corrections are needed, direct formula works!
            causal_coeffs <- sqrt( herit * sigma_sq / ( 2 * pq * m_causal ) )
            # best results are obtained when signs are random too
            causal_coeffs <- causal_coeffs * sample( c(1, -1), m_causal, replace = TRUE )
        } else {
            # draw them randomly (standard normal)
            causal_coeffs <- stats::rnorm(m_causal, 0, 1)
            
            # the initial genetic variance is
            sigma_0 <- sqrt( 2 * sum( pq * causal_coeffs^2 ) )
            # adjust causal_coeffs so final variance is sigma_sq*herit as desired!
            causal_coeffs <- causal_coeffs * sqrt( sigma_sq * herit ) / sigma_0 # scale by standard deviations
        }
        
        # construct genotype signal
        if ( anyNA( X ) ) {
            # if any of the causal loci are missing, let's treat them as zeroes
            # this isn't perfect but we must do something to apply this to real data
            X[ is.na( X ) ] <- 0
        }
        G <- drop( causal_coeffs %*% X ) # this is a vector
        # NOTE by construction:
        # Cov(G) = 2 * herit * kinship
        
        ##############
        ### CENTER ###
        ##############

        # calculate the mean of the genetic effect
        if ( !missing( p_anc ) ) {
            # parametric solution
            muXB <- 2 * drop( causal_coeffs %*% p_anc )
        } else {
            # works very well assuming causal_coeffs and p_anc are uncorrelated!
            muXB <- 2 * sum( causal_coeffs ) * mean( p_anc_hat )
        }
        # in all cases:
        # - remove the mean from the genotypes (muXB)
        # - add the desired mean
        G <- G - muXB + mu
    }

    if (herit == 1) {
        E <- 0 # in this edge case there is no "noise", just genotype effects
    } else {
        # length of E
        n_ind <- ncol(X)
        # draw noise
        E <- stats::rnorm(n_ind, 0, (1 - herit) * sigma_sq ) # noise has mean zero but variance ((1-herit) * sigma_sq)
        # NOTE by construction:
        # Cov(E) = (1-herit) * sigma_sq * I
    }

    # lastly, here's the trait:
    trait <- G + E

    # return all these things
    list(
        trait = trait,
        causal_indexes = causal_indexes,
        causal_coeffs = causal_coeffs
    )
}
