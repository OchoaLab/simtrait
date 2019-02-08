# a package for simulating traits

# NOTE: variance model requires true pAnc to work, technically
# can we somehow make it work with just observed genotypes?

# this is a complex trait simulation
# in order for heritability interpretation to be correct, we need the pAnc model parameters to set the scale of the estimates.
sim_trait <- function(X, pAnc, herit, mafCut=0.05, nSnps=1) {
    # parameters:
    # mafCut: 5% is the threshold we want to consider for minor allele frequencies
    # nSnps: number of causative SNPs, to have polygenic traits!

    # data dimensions (will be transposed if X is a BEDMatrix object)
    m <- nrow(X)
    
    # select random SNPs! this performs the magic...
    i <- select_loci(X, mafCut, nSnps)

    # construct beta coefficients that are mostly zero (sparse)
    beta <- rep.int(0, m)
    # draw random SNP coefficients for selected loci
    beta[i] <- rnorm(nSnps, 0, 1)
    # the genetic variance (in a vector of per-loci terms) is
    varXB <- 4 * pAnc * (1 - pAnc) * beta^2
    # that should equal the heritability is
    # sum( varXB ) = 2 * herit * s
    # so the scale factor we want is
    # (NOTE only non-zero terms are expicitly summed over for greater speed and to avoid numerical issues)
    s <- sum( varXB[i] ) / (2 * herit)
    # adjust betas so total variance agrees with heritability
    beta[i] <- beta[i] / sqrt(s) # divide by standard deviation
    
    # construct genotype signal
    G <- drop( beta %*% X[i,] ) # this is a vector
    # NOTE by construction
    # Cov(G) = 2 * herit * Phi

    # old empirical normalization, probably failed...
    ## varG <- var(G) # sample variance of genotype contributions
    ## s <- sqrt( herit / varG ) # rescale to have desired signal to noise
    ## G <- G * s # rescale to have desired signal to noise
    ## beta <- beta * s # we care about this!

    # draw noise
    E <- rnorm(G, 0, 1 - herit) # noise has mean zero but variance (1-herit), length matches that of G (genotypes)
    # By construction
    # Cov(E) = (1-herit) * I

    # lastly, here's the trait:
    y <- G + E

    # return all these things
    list(y=y, i=i, beta=beta) #, G=G, E=E)
}



