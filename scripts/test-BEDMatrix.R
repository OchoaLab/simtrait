library(optparse)  # for terminal options
library(BEDMatrix) # to load real genotypes with low memory usage
library(simtrait)  # to simulate complex trait

# define options
option_list = list(
    make_option(c("-i", "--input"), type = 'character', default = NA, 
                help = "Path for input genotypes (BED/BIM/FAM files, excluding extension)"),
    make_option("--herit", type = "double", default = 0.8, 
                help = "heritability", metavar = "double"),
    make_option("--m_causal", type = "integer", default = 100, 
                help = "num causal loci", metavar = "int"),
    ## make_option(c("-s", "--sort_causal"), action = "store_true", default = FALSE, 
    ##             help = "num causal loci")
)

opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

# get values
name_in <- opt$input
if (is.na(name_in))
    stop('`-i/--input` is mandatory!')
m_causal <- opt$m_causal
herit <- opt$herit
#sort_causal <- opt$sort_causal

# load genotypes
message('BEDMatrix')
X <- BEDMatrix(name_in)

# HACK
kinship_estimate <- 0

# simulate trait!
# can be done inside an inner loop to simulate thousands of traits per genotype dataset
message('simtrait')
# NOTE: using kinship version, only choice for real data
obj_trait <- sim_trait(
    X = X,
    m_causal = m_causal,
    herit = herit,
    kinship = kinship_estimate
#    sort_causal = sort_causal
)
#trait <- obj_trait$trait
#causal_indexes <- obj_trait$causal_indexes
