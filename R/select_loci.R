# internal function
select_loci <- function(maf, m_causal, maf_cut=0.05) {
    # check for missing parameters
    if (missing(maf)) stop('Fatal: marginal allele frequency vector `maf` must be specified (no default value)')
    if (missing(m_causal)) stop('Fatal: the number of causal loci `m_causal` must be specified (no default value)')
    
    # data dimensions
    m <- length(maf)
    # other checks
    if (m_causal > m) stop('Fatal: the number of causal loci selected cannot be larger than the total number of loci (', m_causal, ' > ', m, ')')
    
    # select random loci!
    # we might not want to pick extremely rare alleles, so set MAF thresholds
    is <- which(maf_cut <= maf & maf <= 1 - maf_cut) # candidate locus indexes
    sample(is, m_causal) # these are the chosen locus indeces!
}
