context("test-simtrait-BEDMatrix")

if (suppressMessages(suppressWarnings(require(BEDMatrix)))) {

    # load BED genotypes using BEDMatrix
    X_BM <- suppressMessages(suppressWarnings(BEDMatrix('dummy-33-101-0.1')))
    # convert to my usual format for comparison
    X_R <- t( X_BM[] )
    

    test_that("allele_freqs works with BEDMatrix", {
        expect_equal(
            allele_freqs(X_R),
            allele_freqs(X_BM)
        )
    })
}
