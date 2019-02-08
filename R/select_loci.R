select_loci <- function(X, mafCut=0.05, nSnps=1) {
    # select random SNPs!
    # one caveat is we don't want to pick extremely rare alleles for this round, so let's pick something that's at least >5% MAF...
    pHat <- rowMeans(X)/2 # not quite "minor" AF, but it's easy to go from here...
    is <- which(mafCut < pHat & pHat < 1 - mafCut) # candidate SNP indexes
    sample(is, nSnps) # these are the chosen SNP indeces!
}
