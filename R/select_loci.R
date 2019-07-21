# internal function
select_loci <- function(maf, m_causal, maf_cut = 0.05) {
    # check for missing parameters
    if (missing(maf))
        stop('marginal allele frequency vector `maf` is required!')
    if (missing(m_causal))
        stop('the number of causal loci `m_causal` is required!')
    
    # data dimensions
    m_loci <- length(maf)
    # other checks
    if (m_causal > m_loci)
        stop('the number of causal loci cannot be larger than the total number of loci (', m_causal, ' > ', m_loci, ')')
    
    # select random loci!
    # we might not want to pick extremely rare alleles, so set MAF thresholds
    i <- which(maf_cut <= maf & maf <= 1 - maf_cut) # candidate locus indexes
    sample(i, m_causal) # these are the chosen locus indeces!
}
