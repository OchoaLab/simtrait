# internal function
select_loci <- function(
                        m_causal,
                        m_loci = NA,
                        counts = NULL,
                        maf_cut = NA,
                        mac_cut = NA
                        ) {
    
    # check for missing parameters
    if ( missing( m_causal ) )
        stop('the number of causal loci `m_causal` is required!')
    # validate counts if provided
    if ( !is.null( counts ) ) {
        if ( !is.matrix( counts ) )
            stop( '`counts` must be a matrix!' )
        if ( ncol( counts ) != 2 )
            stop( '`counts` must have two columns!' )
        if ( !is.numeric( counts ) )
            stop( '`counts` must be numeric!' )
        if ( anyNA( counts ) )
            stop( '`counts` cannot have missing values!' )
    }
    
    # get number of loci, one way or another, and perform loads of checks
    if ( is.na( m_loci ) ) {
        # get it from counts matrix
        if ( !is.null( counts ) ) {
            m_loci <- nrow( counts )
        } else 
            stop('`m_loci` or `counts` are required!')
    } else {
        # if `m_loci` was given, check the matrix too if given
        if ( !is.null( counts ) && m_loci != nrow( counts ) )
            stop( 'The number of rows of `counts` must equal `m_loci`!' )
    }
    
    # other checks
    if (m_causal > m_loci)
        stop('the number of causal loci cannot be larger than the total number of loci (', m_causal, ' > ', m_loci, ')')

    # make sure counts are available for these cases, then calculate the items I need
    if ( !is.na( maf_cut ) || !is.na( mac_cut ) ) {
        if ( is.null( counts ) )
            stop( '`counts` is required if either `maf_cut` or `mac_cut` are provided!' )
        # getting minor counts is always a necessity
        mac <- pmin( counts[,1], counts[,2] )
        # calculate MAF (*minor* AF) now if we want an MAF threshold
        if ( !is.na( maf_cut ) )
            maf <- mac / rowSums( counts )
    }

    # we might not want to pick extremely rare alleles, so set MAF thresholds
    # define narrower candidate locus indexes
    indexes <- NULL
    if ( !is.na( maf_cut ) && !is.na( mac_cut ) ) {
        # if both MAF and MAC filters were defined, use them both simultaneously
        # this is the smartest way to do it
        indexes <- which( maf_cut <= maf & mac_cut <= mac )
    } else if ( !is.na( maf_cut ) ) {
        indexes <- which( maf_cut <= maf )
    } else if ( !is.na( mac_cut ) ) {
        indexes <- which( mac_cut <= mac )
    }
    
    # select random loci!
    if ( !is.null( indexes ) ) {
        causal_indexes <- sample( indexes, m_causal ) # these are the chosen locus indeces!
    } else {
        # in extremely large settings, let's just pick completely at random
        causal_indexes <- sample.int( m_loci, m_causal )
    }
    return( causal_indexes )
}
