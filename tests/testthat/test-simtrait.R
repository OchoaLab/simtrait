context("test-simtrait")

test_that( "sim_trait_model works", {
    # start with some toy data, separate from `sim_trait`
    # (it will be re-tested later through `sim_trait`)

    # make the data to test with
    m_loci <- 100
    n_ind <- 7
    X <- matrix( rbinom( m_loci * n_ind, 2, 0.5 ), m_loci, n_ind )

    # start with a non-genetic model with zero heritability, which is supposed to work
    model <- list(
        sigma_sq = 1.7,
        herit = 0,
        sigma_sq_residual = 1,
        alpha = -3.7
    )
    expect_silent( obj <- sim_trait_model( model, X ) )
    expect_true( is.numeric( obj$trait ) )
    expect_equal( length( obj$trait ), n_ind )
    expect_true( !anyNA( obj$trait ) )
    expect_equal( obj$group_effects, 0 )
    # these subsets are just copied over
    expect_equal( obj[ names( model ) ], model[ names( model ) ] )
    
    # add genetic params, make it more complicated
    model$herit <- 0.7
    model$sigma_sq_residual <- 1 - model$herit
    # first expect error because we're missing causal info despite non-zero herit
    expect_error( sim_trait_model( model, X ) )
    # now add causal info
    m_causal <- 11
    model$causal_indexes <- sample.int( m_loci, m_causal )
    model$causal_coeffs <- rnorm( m_causal )
    # this should now be successful
    expect_silent( obj <- sim_trait_model( model, X ) )
    expect_true( is.numeric( obj$trait ) )
    expect_equal( length( obj$trait ), n_ind )
    expect_true( !anyNA( obj$trait ) )
    expect_equal( obj$group_effects, 0 )
    expect_equal( obj[ names( model ) ], model[ names( model ) ] )

    # add missingness to genotypes
    p_miss <- 0.01
    X[ sample.int( length( X ), length( X ) * p_miss ) ] <- NA
    # repeat test
    expect_silent( obj <- sim_trait_model( model, X ) )
    expect_true( is.numeric( obj$trait ) )
    expect_equal( length( obj$trait ), n_ind )
    expect_true( !anyNA( obj$trait ) )
    expect_equal( obj$group_effects, 0 )
    expect_equal( obj[ names( model ) ], model[ names( model ) ] )

    # push the last model to have zero residual variance, another edge case
    model2 <- model
    model2$herit <- 1
    model2$sigma_sq_residual <- 0
    # repeat test
    expect_silent( obj <- sim_trait_model( model2, X ) )
    expect_true( is.numeric( obj$trait ) )
    expect_equal( length( obj$trait ), n_ind )
    expect_true( !anyNA( obj$trait ) )
    expect_equal( obj$group_effects, 0 )
    expect_equal( obj[ names( model2 ) ], model2[ names( model2 ) ] )
    
    # last tests are for group effects
    # try single level first
    labs <- sample( c('a', 'b', 'c'), n_ind, replace = TRUE )
    model$labs_sigma_sq <- 0.2
    model$sigma_sq_residual <- 1 - model$herit - sum( model$labs_sigma_sq )
    # repeat test
    expect_silent( obj <- sim_trait_model( model, X, labs ) )
    expect_true( is.numeric( obj$trait ) )
    expect_equal( length( obj$trait ), n_ind )
    expect_true( !anyNA( obj$trait ) )
    expect_true( is.numeric( obj$group_effects ) )
    expect_equal( length( obj$group_effects ), n_ind )
    expect_true( !anyNA( obj$group_effects ) )
    expect_equal( obj[ names( model ) ], model[ names( model ) ] )

    # now two levels
    labs2 <- sample( c('x', 'y', 'z'), n_ind, replace = TRUE )
    labs <- cbind( labs, labs2 )
    model$labs_sigma_sq <- c( 0.1, 0.15 )
    model$sigma_sq_residual <- 1 - model$herit - sum( model$labs_sigma_sq )
    # repeat test
    expect_silent( obj <- sim_trait_model( model, X, labs ) )
    expect_true( is.numeric( obj$trait ) )
    expect_equal( length( obj$trait ), n_ind )
    expect_true( !anyNA( obj$trait ) )
    expect_true( is.numeric( obj$group_effects ) )
    expect_equal( length( obj$group_effects ), n_ind )
    expect_true( !anyNA( obj$group_effects ) )
    expect_equal( obj[ names( model ) ], model[ names( model ) ] )
})

test_that( "fold_allele_freqs works", {
    # dies without arguments
    expect_error( fold_allele_freqs() )

    # create some data for test
    # first create a case where there's no flipping at all because everything is already small enough
    m_loci <- 10
    ps <- runif( m_loci, max = 0.5 )
    expect_equal( ps, fold_allele_freqs( ps ) )
    # conversely, we know what to expect when everything was flipped to be on the other side
    expect_equal( ps, fold_allele_freqs( 1 - ps ) )

    # toy case, fully manually constructed
    expect_equal( c(0.1, 0.5, 0.2), fold_allele_freqs( c(0.1, 0.5, 0.8) ) )
})

test_that("allele_freqs works", {
    # Construct toy data
    X <- matrix(
        c(0, 1, 2,
          1, 0, 1,
          1, NA, 2),
        nrow = 3,
        byrow = TRUE
    )
    # known values
    maf_rows <- c(1/2, 1/3, 3/4)
    maf_cols <- c(1/3, 1/4, 5/6)
    counts_rows <- cbind( c(3, 2, 3), c(3, 4, 1) )
    counts_cols <- cbind( c(2, 1, 5), c(4, 3, 1) )
    obj_rows <- list( p_anc_est = maf_rows, counts = counts_rows )
    obj_cols <- list( p_anc_est = maf_cols, counts = counts_cols )
    
    # row means
    expect_equal(
        allele_freqs( X ),
        maf_rows
    )
    
    # col means
    expect_equal(
        allele_freqs( X, loci_on_cols = TRUE ),
        maf_cols
    )

    # repeat with full count matrix
    expect_equal(
        allele_freqs( X, want_counts = TRUE ),
        obj_rows
    )
    expect_equal(
        allele_freqs( X, want_counts = TRUE, loci_on_cols = TRUE ),
        obj_cols
    )

    # repeat with folded allele frequencies option
    # known values (modified to be *minor* allele frequencies)
    maf_rows <- c(1/2, 1/3, 1/4)
    maf_cols <- c(1/3, 1/4, 1/6)
    
    # row means
    expect_equal(
        allele_freqs(X, fold = TRUE),
        maf_rows
    )
    
    # col means
    expect_equal(
        allele_freqs(X, fold = TRUE, loci_on_cols = TRUE),
        maf_cols
    )

    # test subsetting
    for ( subset_ind in 1L:3L ) {
        # test a somewhat ridiculous case with just one individual, which also tests an edge case
        # row means
        expect_equal(
            allele_freqs( X, subset_ind = subset_ind ),
            X[ , subset_ind ] / 2
        )
        
        # col means
        expect_equal(
            allele_freqs( X, subset_ind = subset_ind, loci_on_cols = TRUE ),
            X[ subset_ind, ] / 2
        )
        
        # now test opposite cases where one is excluded
        # also tests negative index case
        # row means
        expect_equal(
            allele_freqs( X, subset_ind = -subset_ind ),
            rowMeans( X[ , -subset_ind ], na.rm = TRUE ) / 2
        )
        
        # col means
        expect_equal(
            allele_freqs( X, subset_ind = -subset_ind, loci_on_cols = TRUE ),
            colMeans( X[ -subset_ind, ], na.rm = TRUE ) / 2
        )
    }
})

validate_indexes <- function( causal_indexes, m_causal, m_loci, maf_cut = NA, maf = NULL, mac_cut = NA, mac = NULL ) {
    # the length of the index vector equals desired m_causal
    expect_equal( length( causal_indexes ), m_causal )
    # all indexes are equal or smaller than m_loci
    expect_true( all( causal_indexes <= m_loci ) )
    # all indexes are equal or larger than 1
    expect_true( all( causal_indexes >= 1 ) )
    # test that MAF and MAC filters work
    # assumes `maf` and `mac` are already *minor*!
    if ( !is.na( maf_cut ) )
        expect_true( min( maf[ causal_indexes ] ) >= maf_cut )
    if ( !is.na( mac_cut ) )
        expect_true( min( mac[ causal_indexes ] ) >= mac_cut )
}

test_that("select_loci works", {
    # cause errors with missing args
    expect_error( select_loci() )
    
    m_loci <- 1000
    m_causal <- 50
    # test simple version without MAF thresholds
    causal_indexes <- select_loci( m_causal = m_causal, m_loci = m_loci )
    validate_indexes( causal_indexes, m_causal, m_loci )
    
    # construct some simple data for test
    # function needs allele counts, let's make something not too trivial
    n_ind <- 100
    # draw counts of first allele
    counts1 <- rbinom( m_loci, n_ind, 0.3 )
    # usually the second allele is the rest, but this way we allow for missingness
    counts2 <- rbinom( m_loci, n_ind - counts1, 0.9 )
    # input to function
    counts <- cbind( counts1, counts2 )
    # and actual vectors to threshold
    mac <- pmin( counts1, counts2 )
    n_obs <- counts1 + counts2
    maf <- mac / n_obs
    
    # first cause an error on purpose
    # (ask for more causal loci than there are loci)
    expect_error( select_loci( m_causal = 10000, counts = counts ) )
    
    # now proper runs
    maf_cut <- 0.05
    mac_cut <- 5
    # MAF only
    causal_indexes <- select_loci( m_causal = m_causal, counts = counts, maf_cut = maf_cut )
    validate_indexes( causal_indexes, m_causal, m_loci, maf_cut, maf )
    # MAC only
    causal_indexes <- select_loci( m_causal = m_causal, counts = counts, mac_cut = mac_cut )
    validate_indexes( causal_indexes, m_causal, m_loci, mac_cut = mac_cut, mac = mac )
    # both MAF and MAC, which are non-redundant under missingness
    causal_indexes <- select_loci( m_causal = m_causal, counts = counts, maf_cut = maf_cut, mac_cut = mac_cut )
    validate_indexes( causal_indexes, m_causal, m_loci, maf_cut, maf, mac_cut, mac )
})

test_that("check_herit_labs works", {
    # test trivial case without labels, just herit
    herit <- 0.8
    labs <- NULL
    labs_sigma_sq <- NULL
    n_ind <- 1 # not used in this case

    # first errors on purpose, all of which have to do with bad heritabilities
    expect_error( check_herit_labs( -0.1, labs, labs_sigma_sq, n_ind ) )
    expect_error( check_herit_labs( 1.1, labs, labs_sigma_sq, n_ind ) )

    # a successful run
    expect_silent(
        obj <- check_herit_labs( herit, labs, labs_sigma_sq, n_ind )
    )
    expect_true( is.list( obj ) )
    expect_equal( names( obj ), c( 'labs', 'sigma_sq_residual' ) )
    expect_true( is.null( obj$labs ) ) # in this case only
    expect_equal( obj$sigma_sq_residual, 1 - herit ) # in this simple case

    # test non-trivial case where there are labs
    # here one layer only
    labs <- c( rep.int( 'a', 5L ), rep.int( 'b', 5L ) )
    labs_sigma_sq <- 0.3
    n_ind <- length( labs )
    # reduce heritability to accomodate effects
    herit <- 0.3
    # cause errors due to mismatches between variables
    expect_error( check_herit_labs( herit, labs, -0.1, n_ind ) )
    expect_error( check_herit_labs( herit, labs, 0.8, n_ind ) )
    expect_error( check_herit_labs( herit, labs, labs_sigma_sq, n_ind + 1L ) )
    expect_error( check_herit_labs( herit, labs, c( labs_sigma_sq, 0.2 ), n_ind ) )

    # a successful run
    expect_silent(
        obj <- check_herit_labs( herit, labs, labs_sigma_sq, n_ind )
    )
    expect_true( is.list( obj ) )
    expect_equal( names( obj ), c( 'labs', 'sigma_sq_residual' ) )
    expect_equal( obj$labs, cbind( labs ) ) # non-NULL case
    expect_equal( obj$sigma_sq_residual, 1 - herit - sum( labs_sigma_sq ) ) # ditto

    # run with two layers of labs
    # here labs are not hierarchical but should work anyway
    labs <- cbind( labs, c( rep.int( 'x', 3L ), rep.int( 'y', 3L ), rep.int( 'z', 4L ) ) )
    labs_sigma_sq <- c( 0.3, 0.2 )
    expect_silent(
        obj <- check_herit_labs( herit, labs, labs_sigma_sq, n_ind )
    )
    expect_true( is.list( obj ) )
    expect_equal( names( obj ), c( 'labs', 'sigma_sq_residual' ) )
    expect_equal( obj$labs, cbind( labs ) ) # non-NULL case
    expect_equal( obj$sigma_sq_residual, 1 - herit - sum( labs_sigma_sq ) ) # ditto
})

test_that("cov_trait works", {
    # when kinship is trivial no-structure, V is identity (any herit)
    n <- 100 # for experiments
    I <- diag(rep.int(1,n)) # identity matrix
    kinship <- I/2 # trivial kinship
    V <- cov_trait(kinship = kinship, herit = 0.6) # should work for any value of herit
    expect_equal(V, I)

    # zero heritability also leads to identity
    kinship <- matrix(runif(n^2), nrow = n) # complete noise matrix
    kinship <- crossprod(kinship) # this makes this random kinship a proper positive-definite matrix!
    kinship <- kinship / max(diag(kinship)) # hacky way to force kinship to be have a max of one (happens along diagonal)
    expect_true( all(kinship >= 0) ) # sanity check, kinship is all above or equal to zero
    expect_true( all(kinship <= 1) ) # sanity check, kinship is all below or equal to one
    V <- cov_trait(kinship = kinship, herit = 0)
    expect_equal(V, I)

    # repeat with same noise kiship and now random heritability
    herit <- runif(1)
    V <- cov_trait(kinship = kinship, herit = 0)
    # all values should be positive
    expect_true( all(V >= 0) )
    # and bounded above by 2 (due to kinship scaling)
    expect_true( all(V <= 2) )
    # V should be symmetric (since kinship is)
    expect_equal( V, t(V) )

    # repeat with group effects
    # labels for group effects
    labs1 <- sample( c('a', 'b', 'c'), n, replace = TRUE )
    labs2 <- sample( c('a', 'z'), n, replace = TRUE )
    labs <- cbind( labs1, labs2 )
    # reduce heritability to have not tiny group variances
    herit <- 0.3
    labs_sigma_sq1 <- 0.3
    labs_sigma_sq2 <- 0.2
    labs_sigma_sq <- c( labs_sigma_sq1, labs_sigma_sq2 )
    
    # test single level version
    expect_silent(
        V <- cov_trait( kinship = kinship, herit = herit, labs = labs1, labs_sigma_sq = labs_sigma_sq1 )
    )
    # all values should be positive
    expect_true( all(V >= 0) )
    # and bounded above by 2 (due to kinship scaling)
    expect_true( all(V <= 2) )
    # V should be symmetric (since kinship is)
    expect_equal( V, t(V) )
    
    # test multi-level version
    expect_silent(
        V <- cov_trait( kinship = kinship, herit = herit, labs = labs, labs_sigma_sq = labs_sigma_sq )
    )
    # all values should be positive
    expect_true( all(V >= 0) )
    # and bounded above by 2 (due to kinship scaling)
    expect_true( all(V <= 2) )
    # V should be symmetric (since kinship is)
    expect_equal( V, t(V) )
    
})

test_that( "herit_loci works", {
    # create toy random data
    m_loci <- 10
    # ancestral allele frequencies
    p_anc <- runif( m_loci )
    # causal loci
    causal_coeffs <- rnorm( m_loci ) / m_loci
    # default value for now
    sigma_sq <- 1
        
    # cause errors on purpose
    # missing mandatory arguments first
    expect_error( herit_loci( ) )
    expect_error( herit_loci( p_anc = p_anc ) )
    expect_error( herit_loci( causal_coeffs = causal_coeffs ) )

    # now a successful case
    expect_silent(
        herit_vec <- herit_loci( p_anc, causal_coeffs )
    )
    # compare to direct calculation
    expect_equal(
        herit_vec,
        2 * p_anc * ( 1 - p_anc ) * causal_coeffs^2 / sigma_sq
    )

    # repeat with a non-default sigma_sq
    sigma_sq <- 9
    expect_silent(
        herit_vec <- herit_loci( p_anc, causal_coeffs, sigma_sq = sigma_sq )
    )
    # compare to direct calculation
    expect_equal(
        herit_vec,
        2 * p_anc * ( 1 - p_anc ) * causal_coeffs^2 / sigma_sq
    )

    # revert to default variance
    sigma_sq <- 1
    # test trivial causal_indexes that is just everything
    causal_indexes <- 1 : m_loci
    expect_silent(
        herit_vec <- herit_loci( p_anc, causal_coeffs, causal_indexes = causal_indexes )
    )
    # compare to direct calculation
    expect_equal(
        herit_vec,
        2 * p_anc * ( 1 - p_anc ) * causal_coeffs^2 / sigma_sq
    )

    # now actually take a subset
    causal_indexes <- 1 : (m_loci/2)
    expect_silent(
        herit_vec <- herit_loci( p_anc, causal_coeffs, causal_indexes = causal_indexes )
    )
    # compare to direct calculation, more cumbersome in this case but meh
    expect_equal(
        herit_vec,
        2 * p_anc[ causal_indexes ] * ( 1 - p_anc[ causal_indexes ] ) * causal_coeffs[ causal_indexes ]^2 / sigma_sq
    )
    
})

test_that( 'p_anc_est_beta_mle works', {
    kinship <- 0.1
    # try a grid of values, including edge cases (0,1)
    m_loci <- 7
    p_anc_est <- 0 : (m_loci - 1) / (m_loci - 1)
    
    # this works on vector inputs!
    expect_silent(
        p_anc_mle <- p_anc_est_beta_mle( p_anc_est, kinship )
    )
    expect_true( is.numeric( p_anc_mle ) )
    expect_equal( length( p_anc_mle ), m_loci )
    expect_true( !anyNA( p_anc_mle ) )
    expect_true( min( p_anc_mle ) >= 0 )
    expect_true( max( p_anc_mle ) <= 1 )

    # test expected symmetry (works since inputs are regular grid)
    indexes_half <- 1L : ( m_loci / 2 )
    expect_equal(
        p_anc_mle[ indexes_half ],
        1 - p_anc_mle[ m_loci + 1L - indexes_half ]
    )

    # make sure that NAs work without errors
    # however, it won't work if all are NAs, so add a dummy and fixed value, test that too
    # in fact all of these inputs are fixed values, confirm that
    p_anc_est <- c( NA, 0, 0.5, 1 )
    expect_silent(
        p_anc_mle <- p_anc_est_beta_mle( p_anc_est, kinship )
    )
    expect_equal( p_anc_mle, p_anc_est )
})

test_that( 'inv_var_est_bayesian works', {
    kinship <- 0.1
    g <- 0.5
    # try a grid of values, excluding for now edge cases (0,1) that result in infinities
    m_loci <- 7
    p_anc_est <- 1 : m_loci / ( m_loci + 1 )

    # a successful example
    expect_silent(
        inv_var_est <- inv_var_est_bayesian( p_anc_est, kinship, g )
    )
    expect_true( is.numeric( inv_var_est ) )
    expect_equal( length( inv_var_est ), m_loci )
    expect_true( !anyNA( inv_var_est ) )
    expect_true( min( inv_var_est ) > 0 )

    # test expected symmetry (works since inputs are regular grid)
    indexes_half <- 1L : ( m_loci / 2 )
    expect_equal(
        inv_var_est[ indexes_half ],
        inv_var_est[ m_loci + 1L - indexes_half ]
    )

    # test special edge cases with known values
    p_anc_est <- c( NA, 0, 1 )
    expect_silent(
        inv_var_est <- inv_var_est_bayesian( p_anc_est, kinship, g )
    )
    expect_equal( inv_var_est, c(NA, Inf, Inf) )
})

validate_sim_trait <- function( obj, herit, n_ind, m_causal, m_loci, p_anc, sigma_sq = 1, maf_cut = NA, maf = NULL, mac_cut = NA, mac = NULL, fes = FALSE, check_herit = TRUE ) {
    trait <- obj$trait # trait vector
    causal_indexes <- obj$causal_indexes # causal locus indeces
    causal_coeffs <- obj$causal_coeffs # locus effect size vector
    # test trait
    expect_equal( length( trait ), n_ind ) # length as expected
    # test causal locus indeces
    validate_indexes( causal_indexes, m_causal, m_loci, maf_cut, maf, mac_cut, mac )
    # test effect sizes
    expect_equal( length( causal_coeffs ), m_causal ) # length as expected
    # these checks only work if p_anc is known, or under the old unbiasing tricks
    # (want to skip for new MLE and Bayesian strategies)
    if ( check_herit ) {
        # verify heritability, exactly given when p_anc is known
        expect_equal(
            herit,
            sum( herit_loci( p_anc[ causal_indexes ], causal_coeffs, sigma_sq = sigma_sq ) )
        )
        # in this case, the per-locus heritabilities are all equal, test that further
        if ( fes ) 
            expect_equal(
                herit_loci( p_anc[ causal_indexes ], obj$causal_coeffs, sigma_sq = sigma_sq ),
                rep.int( herit / m_causal, m_causal )
            )
    }
    # test separate group effects, if present
    lg <- length( obj$group_effects )
    if ( lg == 1 ) {
        expect_equal( obj$group_effects, 0 )
    } else {
        expect_equal( lg, n_ind )
    }
}

test_that("sim_trait works", {
    # in these tests maf_cut stays with default value
    # create unstructured data for test
    n_ind <- 10
    m_loci <- 100
    mn <- n_ind * m_loci # product recurs
    p_anc <- runif( m_loci ) # random ancestral allele frequencies
    X <- matrix( data = rbinom( mn , 2, p_anc ), nrow = m_loci )
    # sample MAFs
    p_anc_hat <- rowMeans( X ) / 2
    expect_true( all(X %in% c(0,1,2)) ) # sanity check for genotypes
    m_causal <- 5
    herit <- 0.8
    # true kinship for unstructured data
    kinship <- diag( n_ind ) / 2
    
    # first cause an error on purpose
    # (ask for more causal loci than there are loci)
    expect_error( sim_trait(X = X, m_causal = 1000, herit = herit, p_anc = p_anc) )
    # provide p_anc too short compared to X's loci
    expect_error( sim_trait( X = X, m_causal = m_causal, herit = herit, p_anc = p_anc[ 1 : (m_loci/2) ] ) )
    # test scalar p_anc in particular, that was a problem before
    expect_error( sim_trait( X = X, m_causal = m_causal, herit = herit, p_anc = p_anc[1] ) )
    
    # test p_anc version
    obj <- sim_trait(X = X, m_causal = m_causal, herit = herit, p_anc = p_anc)
    validate_sim_trait( obj, herit, n_ind, m_causal, m_loci, p_anc )
    
    # test kinship version
    obj <- sim_trait(X = X, m_causal = m_causal, herit = herit, kinship = kinship)
    validate_sim_trait( obj, herit, n_ind, m_causal, m_loci, p_anc_hat, sigma_sq = 1 - mean(kinship) )

    # throw in random missing values in X, repeat all tests
    missingness <- 0.01 # simulate a reasonably low proportion of missingness
    iM <- sample.int( mn, mn * missingness) # random loci to set to NA
    X[iM] <- NA # introduce random missing values
    # recalculate sample MAFs under missingness
    p_anc_hat <- rowMeans( X, na.rm = TRUE ) / 2
    # also calculate MACs and proper MAFs (*minor*)
    counts1 <- rowSums( X, na.rm = TRUE )
    counts2 <- rowSums( 2 - X, na.rm = TRUE )
    mac <- pmin( counts1, counts2 )
    maf <- mac / ( counts1 + counts2 )
    
    # test p_anc version
    obj <- sim_trait(X = X, m_causal = m_causal, herit = herit, p_anc = p_anc)
    validate_sim_trait( obj, herit, n_ind, m_causal, m_loci, p_anc )
    
    # test kinship version
    obj <- sim_trait(X = X, m_causal = m_causal, herit = herit, kinship = kinship)
    validate_sim_trait( obj, herit, n_ind, m_causal, m_loci, p_anc_hat, sigma_sq = 1 - mean(kinship) )

    # try version that centerscales the old way
    # first version uses true p_anc
    obj <- sim_trait( X = X, m_causal = m_causal, herit = herit, p_anc = p_anc, old_center_scale = TRUE )
    validate_sim_trait( obj, herit, n_ind, m_causal, m_loci, check_herit = FALSE )
    # repeat without true p_anc, ought to estimate them
    obj <- sim_trait( X = X, m_causal = m_causal, herit = herit, old_center_scale = TRUE )
    validate_sim_trait( obj, herit, n_ind, m_causal, m_loci, check_herit = FALSE )
    # and try an alternative way to calculate the variance to standardize
    obj <- sim_trait( X = X, m_causal = m_causal, herit = herit, old_center_scale = TRUE, old_sample_var = TRUE )
    validate_sim_trait( obj, herit, n_ind, m_causal, m_loci, check_herit = FALSE )

    # test now a case with specified causal indexes
    # just use the ones from the previous run
    # must work without m_causal
    # test p_anc version
    expect_silent(
        obj2 <- sim_trait( X = X, causal_indexes = obj$causal_indexes, herit = herit, p_anc = p_anc)
    )
    validate_sim_trait( obj2, herit, n_ind, m_causal, m_loci, p_anc )
    # confirm that they stayed the same
    expect_equal( obj2$causal_indexes, obj$causal_indexes )
    # test kinship version
    expect_silent(
        obj2 <- sim_trait( X = X, causal_indexes = obj$causal_indexes, herit = herit, kinship = kinship )
    )
    validate_sim_trait( obj2, herit, n_ind, m_causal, m_loci, p_anc_hat, sigma_sq = 1 - mean(kinship) )
    expect_equal( obj2$causal_indexes, obj$causal_indexes )
    
    # test with MAF threshold
    maf_cut <- 0.05
    
    # test p_anc version
    obj <- sim_trait(X = X, m_causal = m_causal, herit = herit, p_anc = p_anc, maf_cut = maf_cut)
    validate_sim_trait( obj, herit, n_ind, m_causal, m_loci, p_anc, maf_cut = maf_cut, maf = maf )
    
    # test kinship version
    obj <- sim_trait(X = X, m_causal = m_causal, herit = herit, kinship = kinship, maf_cut = maf_cut)
    validate_sim_trait( obj, herit, n_ind, m_causal, m_loci, p_anc_hat, sigma_sq = 1 - mean(kinship), maf_cut = maf_cut, maf = maf )

    # MAC threshold now
    mac_cut <- 5
    obj <- sim_trait( X = X, m_causal = m_causal, herit = herit, p_anc = p_anc, mac_cut = mac_cut )
    validate_sim_trait( obj, herit, n_ind, m_causal, m_loci, p_anc, mac_cut = mac_cut, mac = mac )
    obj <- sim_trait( X = X, m_causal = m_causal, herit = herit, kinship = kinship, mac_cut = mac_cut )
    validate_sim_trait( obj, herit, n_ind, m_causal, m_loci, p_anc_hat, sigma_sq = 1 - mean(kinship), mac_cut = mac_cut, mac = mac )

    # and both MAF and MAC thresholds
    obj <- sim_trait( X = X, m_causal = m_causal, herit = herit, p_anc = p_anc, maf_cut = maf_cut, mac_cut = mac_cut )
    validate_sim_trait( obj, herit, n_ind, m_causal, m_loci, p_anc, maf_cut = maf_cut, maf = maf, mac_cut = mac_cut, mac = mac )
    obj <- sim_trait( X = X, m_causal = m_causal, herit = herit, kinship = kinship, maf_cut = maf_cut, mac_cut = mac_cut )
    validate_sim_trait( obj, herit, n_ind, m_causal, m_loci, p_anc_hat, sigma_sq = 1 - mean(kinship), maf_cut = maf_cut, maf = maf, mac_cut = mac_cut, mac = mac )

    # test fes version
    # test p_anc version
    obj <- sim_trait( X = X, m_causal = m_causal, herit = herit, p_anc = p_anc, fes = TRUE )
    validate_sim_trait( obj, herit, n_ind, m_causal, m_loci, p_anc, fes = TRUE )
    # consider new tricks for unbiasing this case when there's no p_anc
    # combine all of these with mac_cut because there are fixed SNPs otherwise!  (maf_cut = 0 doesn't remove anything, matching is greater or equal!)
    mac_cut <- 1
    # original!
    obj <- sim_trait( X = X, m_causal = m_causal, herit = herit, kinship = kinship, fes = TRUE, mac_cut = mac_cut )
    validate_sim_trait( obj, herit, n_ind, m_causal, m_loci, p_anc_hat, sigma_sq = 1 - mean(kinship), fes = TRUE, mac_cut = mac_cut, mac = mac )
    # MLE!
    obj <- sim_trait( X = X, m_causal = m_causal, herit = herit, kinship = kinship, fes = TRUE, mac_cut = mac_cut, fes_kinship_method = 'mle' )
    validate_sim_trait( obj, herit, n_ind, m_causal, m_loci, fes = TRUE, mac_cut = mac_cut, mac = mac, check_herit = FALSE )
    # Bayesian!
    obj <- sim_trait( X = X, m_causal = m_causal, herit = herit, kinship = kinship, fes = TRUE, mac_cut = mac_cut, fes_kinship_method = 'bayesian' )
    validate_sim_trait( obj, herit, n_ind, m_causal, m_loci, fes = TRUE, mac_cut = mac_cut, mac = mac, check_herit = FALSE )
    
    # labels for group effects
    labs1 <- sample( c('a', 'b', 'c'), n_ind, replace = TRUE )
    labs2 <- sample( c('a', 'z'), n_ind, replace = TRUE )
    labs <- cbind( labs1, labs2 )
    # reduce heritability to have not tiny group variances
    herit <- 0.3
    labs_sigma_sq1 <- 0.3
    labs_sigma_sq2 <- 0.2
    labs_sigma_sq <- c( labs_sigma_sq1, labs_sigma_sq2 )
    
    # version with group effects
    # for simplicity use p_anc version
    # pass labs without their variances
    expect_error( sim_trait( X = X, m_causal = m_causal, herit = herit, p_anc = p_anc, labs = labs1 ) )
    # variances of wrong length
    expect_error( sim_trait( X = X, m_causal = m_causal, herit = herit, p_anc = p_anc, labs = labs1, labs_sigma_sq = labs_sigma_sq ) )
    # excessive variance value ( herit + labs_sigma_sq > 1 )
    expect_error( sim_trait( X = X, m_causal = m_causal, herit = herit, p_anc = p_anc, labs = labs1, labs_sigma_sq = 0.8 ) )

    # a complete run, one level
    expect_silent(
        obj <- sim_trait( X = X, m_causal = m_causal, herit = herit, p_anc = p_anc, labs = labs1, labs_sigma_sq = labs_sigma_sq1 )
    )
    validate_sim_trait( obj, herit, n_ind, m_causal, m_loci, p_anc )
    
    # a complete run, two levels
    expect_silent(
        obj <- sim_trait( X = X, m_causal = m_causal, herit = herit, p_anc = p_anc, labs = labs, labs_sigma_sq = labs_sigma_sq )
    )
    validate_sim_trait( obj, herit, n_ind, m_causal, m_loci, p_anc )
})

test_that( "sqrt_matrix works", {
    # create a random positive semidefinite matrix
    m <- 20
    n <- 10
    X <- matrix(
        rnorm( m*n ),
        nrow = m,
        ncol = n
    )
    V <- tcrossprod( X )
    expect_equal( ncol( V ), m )
    expect_equal( nrow( V ), m )

    # now use function
    expect_silent( 
        V_sqrt <- sqrt_matrix( V )
    )
    # this is the important property we require of the matrix square root
    expect_equal( V, tcrossprod( V_sqrt ) )
})

test_that("sim_trait_mvn works", {
    # create unstructured data for test
    n <- 10
    herit <- 0.8
    # true kinship for unstructured data
    kinship <- diag(n) / 2
    # trait replicates
    rep <- 10
    
    # make sure function stops if crucial parameters are missing
    expect_error( sim_trait_mvn() )
    expect_error( sim_trait_mvn( rep = rep ) )
    expect_error( sim_trait_mvn( kinship = kinship ) )
    expect_error( sim_trait_mvn( herit = herit ) )
    expect_error( sim_trait_mvn( rep = rep, kinship = kinship ) )
    expect_error( sim_trait_mvn( rep = rep, herit = herit ) )
    expect_error( sim_trait_mvn( kinship = kinship, herit = herit ) )

    # complete invocation
    expect_silent(
        traits <- sim_trait_mvn(
            rep = rep,
            kinship = kinship,
            herit = herit
        )
    )
    expect_true( is.matrix( traits ) )
    expect_true( is.numeric( traits ) )
    expect_equal( nrow( traits ), rep )
    expect_equal( ncol( traits ), n )

    # labels for group effects
    labs1 <- sample( c('a', 'b', 'c'), n, replace = TRUE )
    labs2 <- sample( c('a', 'z'), n, replace = TRUE )
    labs <- cbind( labs1, labs2 )
    # reduce heritability to have not tiny group variances
    herit <- 0.3
    labs_sigma_sq1 <- 0.3
    labs_sigma_sq2 <- 0.2
    labs_sigma_sq <- c( labs_sigma_sq1, labs_sigma_sq2 )

    # test one-level version
    expect_silent(
        traits <- sim_trait_mvn(
            rep = rep,
            kinship = kinship,
            herit = herit,
            labs = labs1,
            labs_sigma_sq = labs_sigma_sq1
        )
    )
    expect_true( is.matrix( traits ) )
    expect_true( is.numeric( traits ) )
    expect_equal( nrow( traits ), rep )
    expect_equal( ncol( traits ), n )

    # test two-level version
    expect_silent(
        traits <- sim_trait_mvn(
            rep = rep,
            kinship = kinship,
            herit = herit,
            labs = labs,
            labs_sigma_sq = labs_sigma_sq
        )
    )
    expect_true( is.matrix( traits ) )
    expect_true( is.numeric( traits ) )
    expect_equal( nrow( traits ), rep )
    expect_equal( ncol( traits ), n )
})

test_that( "rmsd works", {
    # random data to test on
    x <- rnorm(10)
    y <- rnorm(10)
    # introduce a little missingess
    x[ 5 ] <- NA
    y[ 8 ] <- NA
    
    # trigger failures on purpose
    # (missing arguments)
    expect_error( rmsd( ) )
    expect_error( rmsd( x ) )
    expect_error( rmsd( y = y ) )
    # arguments of different lengths
    expect_error( rmsd( x, 1 ) )
    expect_error( rmsd( x, y[ 1:5 ] ) )

    # now the proper run
    expect_silent( d <- rmsd( x, y ) )
    expect_equal( length( d ), 1 )
    expect_true( !is.na( d ) )
    expect_true( d > 0 )

    # other obvious checks
    expect_equal( rmsd( x, x ), 0 )
    expect_equal( rmsd( y, y ), 0 )
    expect_equal( rmsd( x, y ), rmsd( y, x ) ) # transitivity
})

test_that( "pval_srmsd, pval_type_1_err, pval_infl work", {
    # random data to test on
    n <- 10
    n_null <- 8 # must be strictly smaller than n for things to work
    p_miss <- 0.2 # add missingness too
    # all p-values are uniform
    pvals <- runif( n )
    # select a random few to be NA
    pvals[ sample.int( n, n * p_miss ) ] <- NA
    # pick a few random ones to be causal (to be removed inside)
    causal_indexes <- sample.int( n, n - n_null )
    
    ### pval_srmsd
    
    # trigger failures on purpose
    # (missing arguments)
    expect_error( pval_srmsd( ) )
    expect_error( pval_srmsd( pvals ) )
    expect_error( pval_srmsd( causal_indexes = causal_indexes ) )
    # p-values out of range cause errors
    expect_error( pval_srmsd( c(pvals, -1), causal_indexes ) )
    expect_error( pval_srmsd( c(pvals, 10), causal_indexes ) )
    # empty causal_indexes trigger a specific error
    expect_error( pval_srmsd( pvals, numeric(0) ) ) # NOTE: c() is NULL!
    # and all loci being causal (no nulls) also triggers errors
    expect_error( pval_srmsd( pvals, 1:length(pvals) ) )

    # now the successful run, simple version
    expect_silent( srmsd <- pval_srmsd( pvals, causal_indexes ) )
    expect_equal( length( srmsd ), 1 )
    expect_true( !is.na( srmsd ) )
    
    # now the successful run, detailed version
    expect_silent( data <- pval_srmsd( pvals, causal_indexes, detailed = TRUE ) )
    expect_equal( class(data), 'list' )
    expect_equal( length(data), 3 )
    expect_equal( names(data), c('srmsd', 'pvals_null', 'pvals_unif') )
    expect_equal( length( data$srmsd ), 1 )
    expect_true( !is.na( data$srmsd ) )
    # equal lengths
    expect_equal( length( data$pvals_null ), length( data$pvals_unif ) )
    # inequality because NAs get things removed too
    expect_true( length( data$pvals_null ) <= n_null )
    expect_true( length( data$pvals_unif ) <= n_null )
    # all p-values should be between zero and one (NAs were removed by sorting!)
    expect_true( !anyNA( data$pvals_null ) )
    expect_true( !anyNA( data$pvals_unif ) )
    expect_true( all( data$pvals_unif >= 0 ) )
    expect_true( all( data$pvals_null >= 0 ) )
    expect_true( all( data$pvals_unif <= 1 ) )
    expect_true( all( data$pvals_null <= 1 ) )

    # work with NULL version
    expect_silent( srmsd <- pval_srmsd( pvals, causal_indexes = NULL ) )
    expect_equal( length( srmsd ), 1 )
    expect_true( !is.na( srmsd ) )

    ### pval_type_1_err
    
    # trigger failures on purpose
    # (missing arguments)
    expect_error( pval_type_1_err( ) )
    expect_error( pval_type_1_err( pvals ) )
    expect_error( pval_type_1_err( causal_indexes = causal_indexes ) )
    # p-values out of range cause errors
    expect_error( pval_type_1_err( c(pvals, -1), causal_indexes ) )
    expect_error( pval_type_1_err( c(pvals, 10), causal_indexes ) )
    # empty causal_indexes trigger a specific error
    expect_error( pval_type_1_err( pvals, numeric(0) ) ) # NOTE: c() is NULL!
    # and all loci being causal (no nulls) also triggers errors
    expect_error( pval_type_1_err( pvals, 1:length(pvals) ) )

    # now the successful run
    expect_silent( type_1_err <- pval_type_1_err( pvals, causal_indexes ) )
    expect_equal( length( type_1_err ), 1 )
    expect_true( !is.na( type_1_err ) )
    expect_true( type_1_err >= 0 )
    expect_true( type_1_err <= 1 )
    
    # work with NULL version
    expect_silent( type_1_err <- pval_type_1_err( pvals, causal_indexes = NULL ) )
    expect_equal( length( type_1_err ), 1 )
    expect_true( !is.na( type_1_err ) )
    expect_true( type_1_err >= 0 )
    expect_true( type_1_err <= 1 )

    # a vector of alphas!
    alpha <- c(1e-1, 1e-2, 1e-3)
    expect_silent( type_1_err <- pval_type_1_err( pvals, causal_indexes, alpha = alpha ) )
    expect_equal( length( type_1_err ), length( alpha ) )
    expect_true( !anyNA( type_1_err ) )
    expect_true( all( type_1_err >= 0 ) )
    expect_true( all( type_1_err <= 1 ) )
    # gives same answer
    pvals_null_cum <- ecdf( pvals[ -causal_indexes ] )
    type_1_err_ecdf <- pvals_null_cum( alpha )
    expect_equal( type_1_err, type_1_err_ecdf )
    # if p-values are NULL, return should be a vector of NAs
    expect_silent( type_1_err <- pval_type_1_err( NULL, causal_indexes, alpha = alpha ) )
    expect_equal( type_1_err, rep.int( NA, length( alpha ) ) )
    
    ### pval_infl
    
    # trigger failures on purpose
    # (missing arguments)
    expect_error( pval_infl( ) )
    # p-values out of range cause errors
    expect_error( pval_infl( c(pvals, -1) ) )
    expect_error( pval_infl( c(pvals, 10) ) )

    # now the successful run, simple version
    expect_silent( infl <- pval_infl( pvals ) )
    expect_equal( length( infl ), 1 )
    expect_true( !is.na( infl ) )

    # other sanity checks
    # perfect non-inflation
    expect_equal( pval_infl( 0.5 ), 1 )
    # inflation
    expect_true( pval_infl( 0.4 ) > 1 )
    # deflation
    expect_true( pval_infl( 0.6 ) < 1 )
    # NULL input
    expect_true( is.na( pval_infl( NULL ) ) )
})

test_that( "pval_gc works", {
    # simulate some non-uniform p-values, so we know that inflation factor is not 1
    m <- 97
    pvals <- rbeta( m, 1, 10 )
    # sprinkle missingness in
    p_miss <- 0.2 # add missingness too
    pvals[ sample.int( m, m * p_miss ) ] <- NA

    # error only if input is missing
    expect_error( pval_gc() )
    
    # confirm that inflation factor is not 1 (it should be much higher) the way data was simulated
    expect_silent( 
        lambda <- pval_infl( pvals )
    )
    expect_true( lambda > 1 )

    # now apply GC correction
    expect_silent( 
        obj <- pval_gc( pvals )
    )
    # confirm that the inflation factor matches what we calculated before
    expect_equal( obj$lambda, lambda )
    # test p-values
    expect_true( is.numeric( obj$pvals ) )
    expect_equal( length( obj$pvals ), m )
    expect_true( min( obj$pvals, na.rm = TRUE ) >= 0 )
    expect_true( max( obj$pvals, na.rm = TRUE ) <= 1 )
    # confirm that GC's inflation is actually exactly 1 now
    expect_silent(
        lambda_gc <- pval_infl( obj$pvals )
    )
    # precision here is lower than it is elsewhere
    expect_equal( lambda_gc, 1, tolerance = 1e-2 )
    # a more direct test, median should be 0.5
    expect_equal( median( obj$pvals, na.rm = TRUE ), 0.5, tolerance = 1e-3 )
})

test_that( "pval_aucpr, pval_power_calib work", {
    # random data to test on
    n <- 10
    n_null <- 5 # must be strictly smaller than n for things to work
    # all p-values are uniform
    pvals <- runif( n )
    # pick a few random ones to be causal (to be removed inside)
    causal_indexes <- sample.int( n, n - n_null )
    # pick one of each class to be NA
    pvals[ causal_indexes ][1] <- NA
    pvals[ -causal_indexes ][1] <- NA

    ### pval_aucpr
    
    # trigger failures on purpose
    # (missing arguments)
    expect_error( pval_aucpr( ) )
    expect_error( pval_aucpr( pvals ) )
    expect_error( pval_aucpr( causal_indexes = causal_indexes ) )
    # p-values out of range cause errors
    expect_error( pval_aucpr( c(pvals, -1), causal_indexes ) )
    expect_error( pval_aucpr( c(pvals, 10), causal_indexes ) )
    # empty causal_indexes trigger a specific error
    expect_error( pval_aucpr( pvals, c() ) )
    # and all loci being causal (no nulls) also triggers errors
    expect_error( pval_aucpr( pvals, 1:length(pvals) ) )

    # now a successful case, simple version
    expect_silent( auc <- pval_aucpr( pvals, causal_indexes ) )
    expect_equal( length(auc), 1 )
    expect_true( !is.na(auc) )
    expect_true( auc >= 0 )
    expect_true( auc <= 1 )

    # now a successful case, curve version
    expect_silent( obj <- pval_aucpr( pvals, causal_indexes, curve = TRUE ) )
    expect_equal( class( obj ), 'PRROC' )
    auc <- obj$auc.integral
    expect_equal( length(auc), 1 )
    expect_true( !is.na(auc) )
    expect_true( auc >= 0 )
    expect_true( auc <= 1 )

    ### pval_power_calib

    # trigger failures on purpose
    # (missing arguments)
    expect_error( pval_power_calib( ) )
    expect_error( pval_power_calib( pvals ) )
    expect_error( pval_power_calib( causal_indexes = causal_indexes ) )
    # p-values out of range cause errors
    expect_error( pval_power_calib( c(pvals, -1), causal_indexes ) )
    expect_error( pval_power_calib( c(pvals, 10), causal_indexes ) )
    # empty causal_indexes trigger a specific error
    expect_error( pval_power_calib( pvals, c() ) )
    # and all loci being causal (no nulls) also triggers errors
    expect_error( pval_power_calib( pvals, 1:length(pvals) ) )

    # now a successful case
    expect_silent( power <- pval_power_calib( pvals, causal_indexes ) )
    expect_equal( length(power), 1 )
    expect_true( !is.na(power) )
    expect_true( power >= 0 )
    expect_true( power <= 1 )

    # a vector of alphas!
    alpha <- c(1e-1, 1e-2, 1e-3)
    expect_silent( power <- pval_power_calib( pvals, causal_indexes, alpha ) )
    expect_equal( length(power), length( alpha ) )
    expect_true( !anyNA(power) )
    expect_true( all( power >= 0 ) )
    expect_true( all( power <= 1 ) )
    # if p-values are NULL, return should be a vector of NAs
    expect_silent( power <- pval_power_calib( NULL, causal_indexes, alpha = alpha ) )
    expect_equal( power, rep.int( NA, length( alpha ) ) )
})
