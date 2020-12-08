context("test-simtrait")

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
    
    # row means
    expect_equal(
        allele_freqs(X),
        maf_rows
    )
    
    # col means
    expect_equal(
        allele_freqs(X, loci_on_cols = TRUE),
        maf_cols
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

})

test_that("select_loci works", {
    # cause errors with missing args
    expect_error( select_loci() )
    
    m_loci <- 1000
    m_causal <- 50
    # test simple version without MAF thresholds
    indexes <- select_loci( m_causal = m_causal, m_loci = m_loci )
    # the length of the index vector equals desired m_causal
    expect_equal( length( indexes ), m_causal )
    # all indexes are equal or smaller than m_loci
    expect_true( all( indexes <= m_loci ) )
    # all indexes are equal or larger than 1
    expect_true( all( indexes >= 1 ) )
    
    # construct some simple data for test
    maf <- ( 1 : m_loci ) / m_loci
    maf_cut <- 0.05
    # first cause an error on purpose
    # (ask for more causal loci than there are loci)
    expect_error( select_loci( m_causal = 10000, maf = maf ) )
    # now a proper run
    indexes <- select_loci( m_causal = m_causal, maf = maf, maf_cut = maf_cut )
    # the length of the index vector equals desired m_causal
    expect_equal( length( indexes ), m_causal )
    # all indexes are equal or smaller than m=length(maf)
    expect_true( all( indexes <= length( maf ) ) )
    # all indexes are equal or larger than 1
    expect_true( all( indexes >= 1 ) )
    # test that MAF filter works
    expect_true( all( maf[ indexes ] >= maf_cut ) ) # test bottom of range
    expect_true( all( maf[ indexes ] <= 1 - maf_cut ) ) # test top of range
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

test_that("sim_trait works", {
    # in these tests maf_cut stays with default value
    # create unstructured data for test
    n <- 10
    m <- 100
    mn <- n*m # product recurs
    p_anc <- runif(m) # random ancestral allele frequencies
    X <- matrix(data = rbinom(mn, 2, p_anc), nrow = m)
    expect_true( all(X %in% c(0,1,2)) ) # sanity check for genotypes
    m_causal <- 5
    herit <- 0.8
    # true kinship for unstructured data
    kinship <- diag(n) / 2
    
    # first cause an error on purpose
    # (ask for more causal loci than there are loci)
    expect_error( sim_trait(X = X, m_causal = 1000, herit = herit, p_anc = p_anc) )
    
    # test p_anc version
    obj <- sim_trait(X = X, m_causal = m_causal, herit = herit, p_anc = p_anc)
    trait <- obj$trait # trait vector
    causal_indexes <- obj$causal_indexes # causal locus indeces
    causal_coeffs <- obj$causal_coeffs # locus effect size vector
    # test trait
    expect_equal( length(trait), n) # length as expected
    # test causal locus indeces
    expect_equal( length(causal_indexes), m_causal ) # length as expected
    expect_true( all(causal_indexes <= m) ) # range as expected
    expect_true( all(causal_indexes >= 1) ) # range as expected
    # test effect sizes
    expect_equal( length(causal_coeffs), m_causal) # length as expected
    # verify heritability, exactly given when p_anc is known
    expect_equal(
        herit,
        sum( herit_loci( p_anc[ causal_indexes ], causal_coeffs ) )
    )
    
    # test kinship version
    obj <- sim_trait(X = X, m_causal = m_causal, herit = herit, kinship = kinship)
    trait <- obj$trait # trait vector
    causal_indexes <- obj$causal_indexes # causal locus indeces
    causal_coeffs <- obj$causal_coeffs # locus effect size vector
    # test trait
    expect_equal( length(trait), n) # length as expected
    # test causal locus indeces
    expect_equal( length(causal_indexes), m_causal ) # length as expected
    expect_true( all(causal_indexes <= m) ) # range as expected
    expect_true( all(causal_indexes >= 1) ) # range as expected
    # test effect sizes
    expect_equal( length(causal_coeffs), m_causal) # length as expected
    # verify heritability, version for unknown p_anc
    p_anc_hat <- rowMeans(X)/2
    expect_equal(
        herit,
        sum( herit_loci( p_anc_hat[ causal_indexes ], causal_coeffs, sigma_sq = 1 - mean(kinship) ) )
    )

    # throw in random missing values in X, repeat all tests
    missingness <- 0.01 # simulate a reasonably low proportion of missingness
    iM <- sample(1:mn, mn * missingness) # random loci to set to NA
    X[iM] <- NA # introduce random missing values
    
    # test p_anc version
    obj <- sim_trait(X = X, m_causal = m_causal, herit = herit, p_anc = p_anc)
    trait <- obj$trait # trait vector
    causal_indexes <- obj$causal_indexes # causal locus indeces
    causal_coeffs <- obj$causal_coeffs # locus effect size vector
    # test trait
    expect_equal( length(trait), n) # length as expected
    # test causal locus indeces
    expect_equal( length(causal_indexes), m_causal ) # length as expected
    expect_true( all(causal_indexes <= m) ) # range as expected
    expect_true( all(causal_indexes >= 1) ) # range as expected
    # test effect sizes
    expect_equal( length(causal_coeffs), m_causal) # length as expected
    # verify heritability, exactly given when p_anc is known
    expect_equal(
        herit,
        sum( herit_loci( p_anc[ causal_indexes ], causal_coeffs ) )
    )
    
    # test kinship version
    obj <- sim_trait(X = X, m_causal = m_causal, herit = herit, kinship = kinship)
    trait <- obj$trait # trait vector
    causal_indexes <- obj$causal_indexes # causal locus indeces
    causal_coeffs <- obj$causal_coeffs # locus effect size vector
    # test trait
    expect_equal( length(trait), n) # length as expected
    # test causal locus indeces
    expect_equal( length(causal_indexes), m_causal ) # length as expected
    expect_true( all(causal_indexes <= m) ) # range as expected
    expect_true( all(causal_indexes >= 1) ) # range as expected
    # test effect sizes
    expect_equal( length(causal_coeffs), m_causal) # length as expected
    # verify heritability, version for unknown p_anc
    p_anc_hat <- rowMeans( X, na.rm = TRUE ) / 2 # repeat for new data with missingness
    expect_equal(
        herit,
        sum( herit_loci( p_anc_hat[ causal_indexes ], causal_coeffs, sigma_sq = 1 - mean(kinship) ) )
    )

    # test with MAF threshold
    maf_cut <- 0.05
    
    # test p_anc version
    obj <- sim_trait(X = X, m_causal = m_causal, herit = herit, p_anc = p_anc, maf_cut = maf_cut)
    trait <- obj$trait # trait vector
    causal_indexes <- obj$causal_indexes # causal locus indeces
    causal_coeffs <- obj$causal_coeffs # locus effect size vector
    # test trait
    expect_equal( length(trait), n) # length as expected
    # test causal locus indeces
    expect_equal( length(causal_indexes), m_causal ) # length as expected
    expect_true( all(causal_indexes <= m) ) # range as expected
    expect_true( all(causal_indexes >= 1) ) # range as expected
    # test effect sizes
    expect_equal( length(causal_coeffs), m_causal) # length as expected
    # verify heritability, exactly given when p_anc is known
    expect_equal(
        herit,
        sum( herit_loci( p_anc[ causal_indexes ], causal_coeffs ) )
    )
    
    # test kinship version
    obj <- sim_trait(X = X, m_causal = m_causal, herit = herit, kinship = kinship, maf_cut = maf_cut)
    trait <- obj$trait # trait vector
    causal_indexes <- obj$causal_indexes # causal locus indeces
    causal_coeffs <- obj$causal_coeffs # locus effect size vector
    # test trait
    expect_equal( length(trait), n) # length as expected
    # test causal locus indeces
    expect_equal( length(causal_indexes), m_causal ) # length as expected
    expect_true( all(causal_indexes <= m) ) # range as expected
    expect_true( all(causal_indexes >= 1) ) # range as expected
    # test effect sizes
    expect_equal( length(causal_coeffs), m_causal) # length as expected
    # verify heritability, version for unknown p_anc
    p_anc_hat <- rowMeans( X, na.rm = TRUE ) / 2 # repeat for new data with missingness
    expect_equal(
        herit,
        sum( herit_loci( p_anc_hat[ causal_indexes ], causal_coeffs, sigma_sq = 1 - mean(kinship) ) )
    )

    # test const_herit_loci version
    # suffices to test p_anc version
    obj <- sim_trait(X = X, m_causal = m_causal, herit = herit, p_anc = p_anc, const_herit_loci = TRUE)
    trait <- obj$trait # trait vector
    causal_indexes <- obj$causal_indexes # causal locus indeces
    causal_coeffs <- obj$causal_coeffs # locus effect size vector
    # test trait
    expect_equal( length(trait), n) # length as expected
    # test causal locus indeces
    expect_equal( length(causal_indexes), m_causal ) # length as expected
    expect_true( all(causal_indexes <= m) ) # range as expected
    expect_true( all(causal_indexes >= 1) ) # range as expected
    # test effect sizes
    expect_equal( length(causal_coeffs), m_causal) # length as expected
    # verify heritability, exactly given when p_anc is known
    expect_equal(
        herit,
        sum( herit_loci( p_anc[ causal_indexes ], causal_coeffs ) )
    )
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

test_that( "pval_srmsd works", {
    # random data to test on
    n <- 10
    n_null <- 8 # must be strictly smaller than n for things to work
    p_miss <- 0.2 # add missingness too
    # all p-values are uniform
    pvals <- runif( n )
    # select a random few to be NA
    pvals[ sample.int( n, n * p_miss ) ] <- NA
    # pick a few random ones to be causal (to be removed inside)
    causal_loci <- sample.int( n, n - n_null )
    
    # trigger failures on purpose
    # (missing arguments)
    expect_error( pval_srmsd( ) )
    expect_error( pval_srmsd( pvals ) )
    expect_error( pval_srmsd( causal_loci = causal_loci ) )
    # p-values out of range cause errors
    expect_error( pval_srmsd( c(pvals, -1), causal_loci ) )
    expect_error( pval_srmsd( c(pvals, 10), causal_loci ) )
    # empty causal_loci trigger a specific error
    expect_error( pval_srmsd( pvals, c() ) )
    # and all loci being causal (no nulls) also triggers errors
    expect_error( pval_srmsd( pvals, 1:length(pvals) ) )

    # now the successful run, simple version
    expect_silent( srmsd <- pval_srmsd( pvals, causal_loci ) )
    expect_equal( length( srmsd ), 1 )
    expect_true( !is.na( srmsd ) )
    
    # now the successful run, detailed version
    expect_silent( data <- pval_srmsd( pvals, causal_loci, detailed = TRUE ) )
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
})

test_that( "pval_infl works", {
    # random data to test on
    n <- 10
    p_miss <- 0.2 # add missingness too
    # all p-values are uniform
    pvals <- runif( n )
    # select a random few to be NA
    pvals[ sample.int( n, n * p_miss ) ] <- NA
    
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
})

test_that( "pval_aucpr works", {
    # random data to test on
    n <- 10
    n_null <- 5 # must be strictly smaller than n for things to work
    # all p-values are uniform
    pvals <- runif( n )
    # pick a few random ones to be causal (to be removed inside)
    causal_loci <- sample.int( n, n - n_null )
    # pick one of each class to be NA
    pvals[ causal_loci ][1] <- NA
    pvals[ -causal_loci ][1] <- NA
    
    # trigger failures on purpose
    # (missing arguments)
    expect_error( pval_aucpr( ) )
    expect_error( pval_aucpr( pvals ) )
    expect_error( pval_aucpr( causal_loci = causal_loci ) )
    # p-values out of range cause errors
    expect_error( pval_aucpr( c(pvals, -1), causal_loci ) )
    expect_error( pval_aucpr( c(pvals, 10), causal_loci ) )
    # empty causal_loci trigger a specific error
    expect_error( pval_aucpr( pvals, c() ) )
    # and all loci being causal (no nulls) also triggers errors
    expect_error( pval_aucpr( pvals, 1:length(pvals) ) )

    # now a successful case, simple version
    expect_silent( auc <- pval_aucpr( pvals, causal_loci ) )
    expect_equal( length(auc), 1 )
    expect_true( !is.na(auc) )
    expect_true( auc >= 0 )
    expect_true( auc <= 1 )

    # now a successful case, curve version
    expect_silent( obj <- pval_aucpr( pvals, causal_loci, curve = TRUE ) )
    expect_equal( class( obj ), 'PRROC' )
    auc <- obj$auc.integral
    expect_equal( length(auc), 1 )
    expect_true( !is.na(auc) )
    expect_true( auc >= 0 )
    expect_true( auc <= 1 )
})
