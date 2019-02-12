context("test-simtrait")

test_that("select_loci works", {
    # construct some simple data for test
    maf <- (1:999)/1000
    maf_cut <- 0.05
    m_causal <- 50
    i <- select_loci(maf=maf, m_causal=m_causal, maf_cut=maf_cut)
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
    Phi <- I/2 # trivial kinship
    V <- cov_trait(kinship=Phi, herit=0.6) # should work for any value of herit
    expect_equal(V, I)

    # zero heritability also leads to identity
    Phi <- matrix(runif(n^2), nrow=n) # complete noise matrix
    Phi <- crossprod(Phi) # this makes this random Phi a proper positive-definite matrix!
    Phi <- Phi / max(diag(Phi)) # hacky way to force kinship to be have a max of one (happens along diagonal)
    expect_true( all(Phi >= 0) ) # sanity check, kinship is all above or equal to zero
    expect_true( all(Phi <= 1) ) # sanity check, kinship is all below or equal to one
    V <- cov_trait(kinship=Phi, herit=0)
    expect_equal(V, I)

    # repeat with same noise kiship and now random heritability
    herit <- runif(1)
    V <- cov_trait(kinship=Phi, herit=0)
    # all values should be positive
    expect_true( all(V >= 0) )
    # and bounded above by 2 (due to kinship scaling)
    expect_true( all(V <= 2) )
    # V should be symmetric (since Phi is)
    expect_equal( V, t(V) )
})

test_that("sim_trait works", {
    # in these tests maf_cut stays with default value
    # create unstructured data for test
    n <- 10
    m <- 100
    mn <- n*m # product recurs
    p_anc <- runif(m) # random ancestral allele frequencies
    X <- matrix(data=rbinom(mn, 2, p_anc), nrow=m)
    expect_true( all(X %in% c(0,1,2)) ) # sanity check for genotypes
    m_causal <- 5
    herit <- 0.8
    Phi <- diag(1/2, n) # true kinship for unstructured data
    # test p_anc version
    obj <- sim_trait(X=X, m_causal=m_causal, herit=herit, p_anc=p_anc)
    y <- obj$y # trait vector
    i <- obj$i # causal locus indeces
    beta <- obj$beta # locus effect size vector
    # test trait
    expect_equal( length(y), n) # length as expected
    # test causal locus indeces
    expect_equal( length(i), m_causal ) # length as expected
    expect_true( all(i <= m) ) # range as expected
    expect_true( all(i >= 1) ) # range as expected
    # test effect sizes
    expect_equal( length(beta), m_causal) # length as expected
    
    # test kinship version
    obj <- sim_trait(X=X, m_causal=m_causal, herit=herit, kinship=Phi)
    y <- obj$y # trait vector
    i <- obj$i # causal locus indeces
    beta <- obj$beta # locus effect size vector
    # test trait
    expect_equal( length(y), n) # length as expected
    # test causal locus indeces
    expect_equal( length(i), m_causal ) # length as expected
    expect_true( all(i <= m) ) # range as expected
    expect_true( all(i >= 1) ) # range as expected
    # test effect sizes
    expect_equal( length(beta), m_causal) # length as expected

    # throw in random missing values in X, repeat all tests
    missingness <- 0.01 # simulate a reasonably low proportion of missingness
    iM <- sample(1:mn, mn*missingness) # random loci to set to NA
    X[iM] <- NA # introduce random missing values
    
    # test p_anc version
    obj <- sim_trait(X=X, m_causal=m_causal, herit=herit, p_anc=p_anc)
    y <- obj$y # trait vector
    i <- obj$i # causal locus indeces
    beta <- obj$beta # locus effect size vector
    # test trait
    expect_equal( length(y), n) # length as expected
    # test causal locus indeces
    expect_equal( length(i), m_causal ) # length as expected
    expect_true( all(i <= m) ) # range as expected
    expect_true( all(i >= 1) ) # range as expected
    # test effect sizes
    expect_equal( length(beta), m_causal) # length as expected
    
    # test kinship version
    obj <- sim_trait(X=X, m_causal=m_causal, herit=herit, kinship=Phi)
    y <- obj$y # trait vector
    i <- obj$i # causal locus indeces
    beta <- obj$beta # locus effect size vector
    # test trait
    expect_equal( length(y), n) # length as expected
    # test causal locus indeces
    expect_equal( length(i), m_causal ) # length as expected
    expect_true( all(i <= m) ) # range as expected
    expect_true( all(i >= 1) ) # range as expected
    # test effect sizes
    expect_equal( length(beta), m_causal) # length as expected
})
