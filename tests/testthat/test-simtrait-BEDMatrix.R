context("test-simtrait-BEDMatrix")

if (suppressMessages(suppressWarnings(require(BEDMatrix)))) {

    # load BED genotypes using BEDMatrix
    X <- suppressMessages(suppressWarnings(BEDMatrix('dummy-33-101-0.1')))
    # convert to my usual format for comparison
    # (regular R matrix, transposed)
    X_R <- t( X[] )
    # dimensions for validations below
    # (transposed from usual, for BEDMatrix)
    n <- nrow(X)
    m <- ncol(X)
    
    test_that("allele_freqs works with BEDMatrix", {
        expect_equal(
            allele_freqs(X_R),
            allele_freqs(X)
        )
    })
    
    test_that("sim_trait works with BEDMatrix", {
        # in these tests maf_cut stays with default value
        m_causal <- 5
        herit <- 0.8
        # trivial kinship, ok for these tests
        kinship <- diag(n) / 2
        # treat sample allele frequencies as truth, just to have something in these tests
        p_anc <- allele_freqs(X)
        
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
}
