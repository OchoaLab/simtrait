context("test-simtrait")

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
})

test_that("select_loci works", {
    
    # construct some simple data for test
    maf <- (1:999)/1000
    maf_cut <- 0.05
    m_causal <- 50
    # first cause an error on purpose
    # (ask for more causal loci than there are loci)
    expect_error( select_loci(maf = maf, m_causal = 10000) )
    # now a proper runx
    i <- select_loci(maf = maf, m_causal = m_causal, maf_cut = maf_cut)
    # the length of the index vector equals desired m_causal
    expect_equal( length(i), m_causal )
    # all indexes are equal or smaller than m=length(maf)
    expect_true( all(i <= length(maf)) )
    # all indexes are equal or larger than 1
    expect_true( all(i >= 1) )
    # test that MAF filter works
    expect_true( all(maf[i] >= maf_cut) ) # test bottom of range
    expect_true( all(maf[i] <= 1-maf_cut) ) # test top of range
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
})
