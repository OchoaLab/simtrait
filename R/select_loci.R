# internal function
select_loci <- function(
                        m_causal,
                        m_loci = NA,
                        maf = NULL,
                        maf_cut = NA
                        ) {
    
    # check for missing parameters
    if ( missing( m_causal ) )
        stop('the number of causal loci `m_causal` is required!')

    # get number of loci, one way or another
    if ( is.na( m_loci ) ) {
        # get it from maf, which must be given then
        if ( is.null( maf ) )
            stop('either `m_loci` or the marginal allele frequency vector `maf` are required!')
        
        # data dimensions
        m_loci <- length( maf )
    }

    # other checks
    if (m_causal > m_loci)
        stop('the number of causal loci cannot be larger than the total number of loci (', m_causal, ' > ', m_loci, ')')

    # select random loci!
    if ( !is.na( maf_cut ) && !is.null( maf ) ) {
        # we might not want to pick extremely rare alleles, so set MAF thresholds
        indexes <- which( maf_cut <= maf & maf <= 1 - maf_cut ) # candidate locus indexes
        sample( indexes, m_causal ) # these are the chosen locus indeces!
    } else {
        # in extremely large settings, let's just pick completely at random
        sample.int( m_loci, m_causal )
    }
}
