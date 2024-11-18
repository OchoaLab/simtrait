#' Simulate a complex trait from genotypes
#'
#' Simulate a complex trait given a SNP genotype matrix and model parameters, which are minimally: the number of causal loci, the heritability, and either the true ancestral allele frequencies used to generate the genotypes or the mean kinship of all individuals.
#' An optional minimum marginal allele frequency for the causal loci can be set.
#' The output traits have by default a zero mean and unit variance (for outbred individuals), but those parameters can be modified.
#' The code selects random loci to be causal, constructs coefficients for these loci (scaled appropriately) and random Normal independent non-genetic effects and random environment group effects if specified.
#' There are two models for constructing causal coefficients: random coefficients (RC; default) and fixed effect sizes (FES; i.e., coefficients roughly inversely proportional to allele frequency; use `fes = TRUE`).
#' Suppose there are `m` loci and `n` individuals.
#'
#' To center and scale the trait and locus coefficients vector correctly to the desired parameters (mean, variance, heritability), the parametric ancestral allele frequencies (`p_anc`) must be known.
#' This is necessary since in the heritability model the genotypes are random variables (with means given by `p_anc` and a covariance structure given by `p_anc` and the kinship matrix), so these genotype distribution parameters are required.
#' If `p_anc` are known (true for simulated genotypes), then the trait will have the specified mean and covariance matrix in agreement with [cov_trait()].
#' To simulate traits using real genotypes, where `p_anc` is unknown, a compromise that works well in practice is possible if the mean `kinship` is known (see package vignette).
#' We recommend estimating the mean kinship using the `popkin` package!
#' 
#' @param X The `m`-by-`n` genotype matrix (if `loci_on_cols = FALSE`, transposed otherwise), or a BEDMatrix object.
#' This is a numeric matrix consisting of reference allele counts (in `c(0, 1, 2, NA)` for a diploid organism).
#' @param m_causal The desired number of causal loci.
#' Ignored if `causal_indexes` is provided.
#' @param herit The desired heritability (proportion of trait variance due to genetics).
#' @param p_anc The length-`m` vector of true ancestral allele frequencies.
#' Optional but recommended for simulations.
#' Either this or `kinship` must be specified.
#' @param kinship The mean kinship value of the individuals in the data.
#' The `n`-by-`n` kinship matrix of the individuals in the data is also accepted.
#' Optional but recommended for real data.
#' Either this or `p_anc` must be specified.
#' @param mu The desired parametric mean value of the trait (scalar, default 0).
#' @param sigma_sq The desired parametric variance factor of the trait (scalar, default 1).
#' Corresponds to the variance of an outbred individual.
#' @param labs Optional labels assigning individuals to groups, to simulate environment group effects.
#' Values can be numeric or strings, simply assigning the same values to individuals in the same group.
#' If vector (single environment), length must be number of individuals.
#' If matrix (multiple environments), individuals must be along rows, and environments along columns.
#' The environments are not required to be nested.
#' If this is non-`NULL`, then `labs_sigma_sq` must also be given!
#' @param labs_sigma_sq Optional vector of environment effect variance proportions, one value for each environment given in `labs` (a scalar if `labs` is a vector, otherwise its length should be the number of columns of `labs`).
#' Ignored unless `labs` is also given.
#' As these are variance proportions, each value must be non-negative and `sum(labs_sigma_sq) + herit <= 1` is required so residual variance is non-negative.
#' @param maf_cut The optional minimum allele frequency threshold (default `NA`, no threshold).
#' This prevents rare alleles from being causal in the simulation.
#' Threshold is applied to the *sample* allele frequencies and not their true parametric values (`p_anc`), even if these are available.
#' Ignored if `causal_indexes` is provided.
#' @param loci_on_cols If `TRUE`, `X` has loci on columns and individuals on rows; if `FALSE` (the default), loci are on rows and individuals on columns.
#' If `X` is a BEDMatrix object, loci are always on the columns (`loci_on_cols` is ignored).
#' @param m_chunk_max BEDMatrix-specific, sets the maximum number of loci to process at the time.
#' If memory usage is excessive, set to a lower value than default (expected only for extremely large numbers of individuals).
#' @param fes If `TRUE`, causal coefficients are inversely proportional to the square root of `p_anc * ( 1 - p_anc )` (estimated when `p_anc` is unavailable), which ensures *fixed effect sizes* (FES) per causal locus.
#' Signs (+/-) are drawn randomly with equal probability.
#' If `FALSE` (the default), *random coefficients* (RC) are drawn from a standard Normal distribution.
#' In both cases coefficients are rescaled to result in the desired heritability.
#' @param fes_kinship_method String specifying the bias correction method for the case when `fes = TRUE` and `kinship` is provided while `p_anc` is absent (option is ignored otherwise).
#' `original` corresponds to a simple formula with substantial bias.
#' `mle` replaces sample allele frequencies with maximum likelihood estimates assuming a Beta distribution with a known variance scale given by `kinship` (see [p_anc_est_beta_mle()]).
#' `bayesian` calculates the expectation over the posterior of a desired inverse variance term, assuming the same Beta distribution as the `mle` option (see [inv_var_est_bayesian()]).
#' @param causal_indexes If provided, will use the loci at these indexes as causal loci.
#' When `NULL` (default), causal indexes are drawn randomly.
#' Thus, parameters `m_loci` and `maf_cut` are ignored and have no effect if this parameter is set.
#'
#' @return A named list containing:
#'
#' - `trait`: length-`n` vector of the simulated trait
#' - `causal_indexes`: length-`m_causal` vector of causal locus indexes
#' - `causal_coeffs`: length-`m_causal` vector of coefficients at the causal loci
#' - `group_effects`: length-`n` vector of simulated environment group effects, or 0 (scalar) if not simulated
#' 
#' However, if `herit = 0` then `causal_indexes` and `causal_coeffs` will have zero length regardless of `m_causal`.
#'
#' @examples
#' # construct a dummy genotype matrix
#' X <- matrix(
#'     data = c(
#'         0, 1, 2,
#'         1, 2, 1,
#'         0, 0, 1
#'     ),
#'     nrow = 3,
#'     byrow = TRUE
#' )
#' # made up ancestral allele frequency vector for example
#' p_anc <- c(0.5, 0.6, 0.2)
#' # made up mean kinship
#' kinship <- 0.2
#' # desired heritability
#' herit <- 0.8
#' 
#' # create simulated trait and associated data
#' # default is *random coefficients* (RC) model
#' obj <- sim_trait(X = X, m_causal = 2, herit = herit, p_anc = p_anc)
#' 
#' # trait vector
#' obj$trait
#' # randomly-picked causal locus indexes
#' obj$causal_indexes
#' # regression coefficients vector
#' obj$causal_coeffs
#' 
#' # *fixed effect sizes* (FES) model
#' obj <- sim_trait(X = X, m_causal = 2, herit = herit, p_anc = p_anc, fes = TRUE)
#'
#' # either model, can apply to real data by replacing `p_anc` with `kinship`
#' obj <- sim_trait(X = X, m_causal = 2, herit = herit, kinship = kinship)
#'
#' @seealso
#' [cov_trait()], [sim_trait_mvn()]
#' 
#' @export
sim_trait <- function(
                      X,
                      m_causal,
                      herit,
                      p_anc = NULL,
                      kinship = NULL,
                      mu = 0,
                      sigma_sq = 1,
                      labs = NULL,
                      labs_sigma_sq = NULL,
                      maf_cut = NA,
                      loci_on_cols = FALSE,
                      m_chunk_max = 1000,
                      fes = FALSE,
                      fes_kinship_method = c('original', 'mle', 'bayesian'),
                      causal_indexes = NULL
                      ) {
    # check for missing parameters
    if (missing(X))
        stop('genotype matrix `X` is required!')
    if (missing(herit))
        stop('the heritability `herit` is required!')
    if (is.null(p_anc) && is.null(kinship))
        stop('either the true ancestral allele frequency vector `p_anc` or `kinship` are required!')
    # if we provided causal indexes...
    if ( !is.null( causal_indexes ) ) {
        # force to ignore this parameter if we don't need it (to avoid calculating p_anc_est)
        maf_cut <- NA
        # infer m_causal from this
        m_causal <- length( causal_indexes )
    } else if ( missing( m_causal ) )
        # require m_causal otherwise
        stop('the number of causal loci `m_causal` is required when `causal_indexes` is not provided!')
    
    # other checks
    if (length(mu) != 1)
        stop('`mu` must be a scalar! (input has length ', length(mu), ')')
    if (length(sigma_sq) != 1)
        stop('`sigma_sq` must be a scalar! (input has length ', length(sigma_sq), ')')
    if (sigma_sq <= 0)
        stop('`sigma_sq` must be positive!')
    fes_kinship_method <- match.arg( fes_kinship_method )
    
    # simplifies subsetting downstream
    if ('BEDMatrix' %in% class(X))
        loci_on_cols <- TRUE

    # get m_loci to check m_causal before more heavy computations
    m_loci <- if (loci_on_cols) ncol(X) else nrow(X)
    n_ind <- if (loci_on_cols) nrow(X) else ncol(X)
    if (m_causal > m_loci)
        stop('m_causal (', m_causal, ') exceeds the number of loci (', m_loci, ')!')
    # compare to p_anc too if it was provided
    if ( !is.null( p_anc ) && length( p_anc ) != m_loci )
        stop( '`p_anc` length (', length( p_anc ) , ') does not equal `m_loci` (', m_loci, ')' )

    # check herit, labs and calculate residual variance with this shared function
    obj <- check_herit_labs( herit, labs, labs_sigma_sq, n_ind )
    labs <- obj$labs
    sigma_sq_residual <- obj$sigma_sq_residual

    if (herit == 0) {
        # lots of work can be avoided in this edge case
        # the index and coefficients vectors are empty
        #causal_indexes <- c() # already null by default
        causal_coeffs <- c()
        G <- 0 # construct a trivial genotype effect of zero (becomes vector automatically later)
    } else {
        
        ###################################
        ### MARGINAL ALLELE FREQUENCIES ###
        ###################################

        # may not need this
        p_anc_est <- NULL
        # these are needed here to select loci, if there's a MAF threshold
        # (they are also needed in real datasets, but if that's the only need then let's wait until we've subset the genotype matrix)
        if ( !is.na( maf_cut ) ) {
            # compute marginal allele frequencies
            p_anc_est <- allele_freqs(
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
        # NOTE: if causal_indexes are provided, they aren't overwritten!
        if ( is.null( causal_indexes ) )
            causal_indexes <- select_loci(
                m_causal = m_causal,
                m_loci = m_loci,
                maf = p_anc_est, # NULL if is.na( maf_cut )
                maf_cut = maf_cut # may be NA
            )
        
        # subset data to consider causal loci only
        if ( !is.null( p_anc_est ) ) # if we had this already
            p_anc_est <- p_anc_est[ causal_indexes ]
        if ( !is.null( p_anc ) )
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
        if ( is.null( p_anc ) && is.null( p_anc_est ) ) {
            # compute marginal allele frequencies
            p_anc_est <- allele_freqs(
                X,
#                loci_on_cols = loci_on_cols, # now that genotypes are extracted, they are in default orientation (no need for passing this, plus passing it messes things up)
                m_chunk_max = m_chunk_max
            )
        }
        
        ###############
        ### KINSHIP ###
        ###############
        
        if (!is.null(kinship)) {
            # precompute some things when this is present
            mean_kinship <- mean(kinship)
        }
        
        #############
        ### SCALE ###
        #############
        
        # to scale causal_coeffs to give correct heritability, we need to estimate the pq = p(1-p) vector
        # calculate pq = p_anc * (1 - p_anc) in various ways
        if ( !is.null( p_anc ) ) { # this takes precedence, it should be most accurate
            # direct calculation
            pq <- p_anc * ( 1 - p_anc )
        } else if ( !is.null( kinship ) ) {
            # for RC this fine, and for FES this otherwise corresponds to `fes_kinship_method == 'original'`
            # indirect, infer from genotypes and mean kinship
            # recall E[ p_anc_est * (1 - p_anc_est) ] = pq * (1 - mean_kinship), so we solve for pq:
            pq <- p_anc_est * ( 1 - p_anc_est ) / ( 1 - mean_kinship )
            if ( fes && fes_kinship_method == 'mle' ) {
                # use this other way to estimate a better value considering the variance of the estimate and assuming Beta
                p_anc_mle <- p_anc_est_beta_mle( p_anc_est, mean_kinship )
                # now use these values directly
                pq <- p_anc_mle * ( 1 - p_anc_mle )
            } else if ( fes && fes_kinship_method == 'bayesian' ) {
                # a Bayesian version for estimating the factor we want (1/sqrt(pq))
                inv_var_est <- inv_var_est_bayesian( p_anc_est, mean_kinship, g = 0.5 )
                # to have the below code work, deconstruct back into the quantity we need
                pq <- 1 / inv_var_est^2
            }
        } else {
            # a redundant check (if this were so, it should have died earlier)
            stop('either `p_anc` or `kinship` must be specified!')
        }
        
        # construct SNP coefficients for selected loci
        if ( fes ) {
            # make them inverse to the pq
            # in this case no scale corrections are needed, direct formula works!
            causal_coeffs <- sqrt( herit * sigma_sq / ( 2 * pq * m_causal ) )
            # best results are obtained when signs are random too
            causal_coeffs <- causal_coeffs * sample( c(1, -1), m_causal, replace = TRUE )
        } else {
            # draw them randomly (standard normal)
            causal_coeffs <- stats::rnorm(m_causal, 0, 1)
            
            # the initial genetic standard deviation is
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
        if ( !is.null( p_anc ) ) {
            # parametric solution
            muXB <- 2 * drop( causal_coeffs %*% p_anc )
        } else {
            # works very well assuming causal_coeffs and p_anc are uncorrelated!
            muXB <- 2 * sum( causal_coeffs ) * mean( p_anc_est )
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
        E <- stats::rnorm(n_ind, 0, sqrt( sigma_sq_residual * sigma_sq ) ) # noise has mean zero but variance (sigma_sq_residual * sigma_sq)
        # NOTE by construction:
        # Cov(E) = sigma_sq_residual * sigma_sq * I
    }

    # environment group effects
    group_effects <- 0
    if ( !is.null( labs ) ) {
        n_labs <- ncol( labs )
        for ( i in 1L : n_labs ) {
            # process environment i
            labs_i <- labs[ , i ]
            labs_i_sigma_sq <- labs_sigma_sq[ i ]
            # unique groups, implicitly numbered by order of appearance
            groups_i <- unique( labs_i )
            # map individuals to their groups by index
            group_indexes_i <- match( labs_i, groups_i )
            # draw their random effects
            group_eff_i <- stats::rnorm( length( groups_i ), 0, sqrt( labs_i_sigma_sq * sigma_sq ) )
            # distribute the effects from groups to individuals
            group_effects_i <- group_eff_i[ group_indexes_i ]
            # add to running sum
            group_effects <- group_effects + group_effects_i
        }
    }
    
    # lastly, here's the trait:
    trait <- G + E + group_effects

    # return all these things
    list(
        trait = trait,
        causal_indexes = causal_indexes,
        causal_coeffs = causal_coeffs,
        group_effects = group_effects
    )
}
