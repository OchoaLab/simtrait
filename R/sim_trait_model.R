#' Simulate a complex trait from a fixed trait model
#'
#' This function requires a trait `model`, which contains values for all of its variance components and intercept coeficient, as well as pre-determined causal loci indexes along with their coefficients, and simply simulates the traits with genetic effects and random non-genetic effects, for the individuals provided.
#' This is useful to simulate traits for new individuals from a pre-existing model that was constructed using some previous data, so that the trait is fundamentally the same in all individuals (rather than a new random trait as `sim_trait` does for each call).
#'
#' @inheritParams sim_trait
#' @param model A list returned by a previous `sim_trait` call.
#' This list must have scalar values with names `sigma_sq`, `herit`, `sigma_sq_residual`, `alpha`;
#' if `herit > 0`, then `causal_indexes` and `causal_coeffs` vectors are also required in this list.
#' Lastly, either `labs` and `model$labs_sigma_sq` must be both missing, or both must be present.
#' @param X_subset If `TRUE`, `X` is treated as if it were already subset to the desired causal loci.  Otherwise (default) `X` will be subset internally using `model$causal_indexes`.
#' @param skip_checks If `TRUE`, skips validations of input model and other parameters.
#' Not recommended for regular users; this is provided for [sim_trait()] internal use only, which calls [sim_trait_model()], since these checks are redundant with [sim_trait()]'s own validations.
#'
#' @return The same type of named list returned by [sim_trait()], please see that for more details.  This includes the desired `trait` values, as well as copies of the input model parameters.
#'
#' @examples
#' # make two genotype matrices, one to construct the initial trait,
#' # the second one is a new dataset from which we want to simulate the same trait.
#' # They should have the same loci (number of rows), but individuals should be different!
#' X <- matrix(
#'     data = c(
#'         0, 1, 2,
#'         1, 2, 1,
#'         0, 0, 1
#'     ),
#'     nrow = 3,
#'     byrow = TRUE
#' )
#' X2 <- matrix(
#'     data = c(
#'         2, 1,
#'         2, 0,
#'         1, 1
#'     ),
#'     nrow = 3,
#'     byrow = TRUE
#' )
#' 
#' # here's a minimal `sim_trait` call that simulates a trait of interest from the first dataset
#' model <- sim_trait( X = X, m_causal = 2, herit = 0.7, kinship = 0.2 )
#'
#' # use the same model to simulate traits for the new individuals!
#' obj <- sim_trait_model( model, X2 )
#' # trait vector for new individuals:
#' obj$trait
#'
#' @seealso
#' [sim_trait()] for simulating traits from scratch, picking random causal variants and their coefficients from various models.
#'
#' @export
sim_trait_model <- function(
                            model,
                            X,
                            labs = NULL,
                            loci_on_cols = FALSE,
                            X_subset = FALSE,
                            skip_checks = FALSE
                            ) {
    if ( !skip_checks ) {
        # check for missing parameters
        if ( missing( model ) )
            stop( '`model` is required!' )
        if ( missing( X ) )
            stop( 'genotype matrix `X` is required!' )
        # simplifies subsetting downstream
        if ( 'BEDMatrix' %in% class( X ) )
            loci_on_cols <- TRUE
        
        # check elements of `model` further
        if ( !is.list( model ) )
            stop( '`model` must be a list!labs_sigma_sq' )
        keys <- c('sigma_sq', 'herit', 'sigma_sq_residual', 'alpha')
        keys <- keys[ !( keys %in% names( model ) ) ]
        if ( length( keys ) > 0 )
            stop( '`model` is missing these keys: ', toString( keys ) )
        # NOTE: `causal_coeffs` can be missing when heritability is zero, though that'd be weird.  Otherwise both must be present.
        if ( model$herit > 0 ) {
            if ( is.null( model$causal_indexes ) )
                stop( '`model` must contain `causal_indexes` if `model$herit > 0`!' )
            if ( is.null( model$causal_coeffs ) )
                stop( '`model` must contain `causal_coeffs` if `model$herit > 0`!' )
            if ( length( model$causal_indexes ) != length( model$causal_coeffs ) )
                stop( '`model$causal_indexes` and `model$causal_coeffs` must have the same length!' )
        }
        # `labs_sigma_sq` must be present if we passed `labs`, and vice versa.  One direction is tested with `check_herit_labs`, here we test both directions anyway because that's easier
        if ( is.null( labs ) != is.null( model$labs_sigma_sq ) )
            stop( '`labs` and `model$labs_sigma_sq` must both or neither be present!' )
    }
    # get dimensions, never skip since this is always needed
    n_ind <- if (loci_on_cols) nrow(X) else ncol(X)
    if ( !skip_checks )
        # check `labs` against `n_ind` and `labs_sigma_sq`, normalize into matrix form if needed
        # (also performs checks on herit and labs_sigma_sq that are not necessary for a model that already exists and was already checked, but meh; labs is for new individuals so that must be checked anew)
        labs <- check_herit_labs( model$herit, labs, model$labs_sigma_sq, n_ind )$labs
    
    # genetic component
    G <- model$alpha # construct a trivial genotype effect that is desired shift only
    if ( model$herit > 0 ) {
        # subset to causal data, if it's not already subset
        if ( !X_subset ) {
            if (loci_on_cols) {
                # also transpose for consistent behavior downstream
                X <- t( X[ , model$causal_indexes, drop = FALSE ] )
            } else{
                X <- X[ model$causal_indexes, , drop = FALSE ]
            }
        } else if (loci_on_cols)
            # if matrix is already subset, we still have to pull out (for BEDMatrix) and transpose whole matrix in this case, or the code below won't work
            X <- t( X[ , ] )
        
        # "impute" missing values if needed, replace by mean value per locus
        if ( anyNA( X ) )
            X <- t( apply( X, 1, function( xi ) {
                # skip entirely if there are no missing values in this row
                if ( anyNA( X ) )
                    xi[ is.na( xi ) ] <- mean( xi, na.rm = TRUE )
                return( xi )
            } ) )
        
        # construct genotype signal
        G <- G + drop( model$causal_coeffs %*% X ) # this is a vector
        # NOTE by construction:
        # Cov(G) = 2 * herit * kinship
    }
    
    # residual
    if ( model$sigma_sq_residual == 0 ) {
        E <- 0 # in this edge case there is no "noise", just genotype effects
    } else {
        # draw noise
        E <- stats::rnorm(n_ind, 0, sqrt( model$sigma_sq_residual * model$sigma_sq ) ) # noise has mean zero but variance (sigma_sq_residual * sigma_sq)
        # NOTE by construction:
        # Cov(E) = sigma_sq_residual * sigma_sq * I
    }
    
    # environment group effects
    group_effects <- 0
    if ( !is.null( labs ) ) {
        labs_sigma_sq <- model$labs_sigma_sq
        n_labs <- ncol( labs )
        for ( i in 1L : n_labs ) {
            # process environment i
            labs_i <- labs[ , i ]
            labs_i_sigma_sq <- labs_sigma_sq[ i ]
            # unique groups, implicitly numbered by order of appearance
            groups_i <- unique( labs_i )
            # map individuals to their groups by index
            group_indexes_i <- match( labs_i, groups_i )
            # draw their random effects
            group_eff_i <- stats::rnorm( length( groups_i ), 0, sqrt( labs_i_sigma_sq * model$sigma_sq ) )
            # distribute the effects from groups to individuals
            group_effects_i <- group_eff_i[ group_indexes_i ]
            # add to running sum
            group_effects <- group_effects + group_effects_i
        }
    }
    
    # lastly, here's the trait:
    trait <- G + E + group_effects

    # a bit of cleanup, so return object doesn't return old trait and group effect values
    if ( !is.null( model$trait ) )
        model$trait <- NULL
    if ( !is.null( model$group_effects ) )
        model$group_effects <- NULL
    
    # return all these things
    c( 
        list(
            trait = trait,
            group_effects = group_effects
        ),
        model
    )
}
