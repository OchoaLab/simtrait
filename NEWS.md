# simtrait 1.0.0.9000 (2019-02-13)

* Public GitHub release!

# simtrait 1.0.1.9000 (2019-03-26)

* Fixed `m_causal=1` bug (used to die, now works correctly).

# simtrait 1.0.2.9000 (2019-04-10)

* Renamed argument `sigmaSq` to `sigma_sq` (parameter of `sim_trait` and `cov_trait`).
* `sim_trait` return list has more descriptive names:
  * `trait` (old `y`)
  * `causal_indexes` (old `i`)
  * `causal_coeffs` (old `beta`)
* Overall coding style changes.
* Vignette changes in parallel with `popkin` 1.2.0.9000 changes.

# simtrait 1.0.3.9000 (2019-04-16)

* Vignette and `README` changes in parallel with `bnpsd` 1.1.0.9000 changes.
* 2019-05-13: added ORCID to author info

# simtrait 1.0.4.9000 (2019-07-19)

* Added `allele_freqs`, which handles `BEDMatrix` objects correctly.
  Rest of `sim_trait` doesn't handle `BEDMatrix` objects yet.

# simtrait 1.0.5.9000 (2019-07-20)

* Now `sim_trait`, and consequently the whole package, supports `BEDMatrix` objects.
  This is ideal for simulating traits from real and potentially very large genotype data.

# simtrait 1.0.6.9000 (2019-08-02)

* Improved `allele_freqs` to load genotype chunks from a `BEDMatrix` object more efficiently, using as much available memory as possible.
  Requires an updated `popkin` package (>= 1.2.6.9000).

# simtrait 1.0.7.9000 (2019-12-17)

* Preemptively updated `class` usage now that matrices return a two-element array in R-devel
* Moved logo to `man/figures/`
* Minor Roxygen-related updates

# simtrait 1.0.8.9000 (2020-07-07)

* Function `sim_trait` changed default `maf_cut` from 0.05 to `NA`, and in that case code no longer computes marginal allele frequencies unless they are needed because `p_anc` is missing (and `kinship` is provided).
  This is much faster and memory efficient in extremely large simulations, and just makes more sense to me generally (the MAF threshold was a quirky feature, still available but non-default now).

# simtrait 1.0.9.9000 (2020-07-20)

* Function `sim_trait` now delays calculating allele frequencies if `maf_cut = NA` and `p_anc` is missing (and `kinship` is provided), so that these frequencies are only calculated on the small subset of loci selected to be causal (rather than the whole genome, which was the original behavior).
  This is expected to speed-up trait simulations from real genotypes.
* Vignette had minor corrections and edits (focused on model and algorithm description).

# simtrait 1.0.10.9000 (2020-08-05)

* Added function `sim_trait_mvn`, which draws traits from the multivariate normal distribution with covariance structure matching that of genetic traits (also called "infinitesimal" model).
  These traits are useful for testing heritability estimation.
  A visual validation of these simulated traits is available in the vignette.

# simtrait 1.0.11.9000 (2020-10-14)

* Function `alelle_freq` was modified to reduce memory usage for the BEDMatrix case.
  - For both this function and `sim_trait`, the BEDMatrix-specific options `mem_factor` and `mem_lim` were replaced by the option `m_chunk_max`.
  - As a consequence of the changes, `popkin` is no longer a dependency of this `simtrait` package.
    However, `popkin` is still recommended for estimating the kinship matrices required by some of the functions of this `simtrait` package.

# simtrait 1.0.12.9000 (2020-11-05)

* Added functions used for benchmarking genetic association methods based on traits simulated with `simtrait`:
  - `rmsd`: General root-mean-square deviation (RMSD) function
  - `pval_srmsd`: Signed RMSD between observed and expected (uniform) p-values
  - `pval_aucpr`: Area under the curve (AUC) of precision-recall (PR) curve
    - Added `PRROC` package dependency

# simtrait 1.0.13.9000 (2020-11-06)

* Added function `pval_infl` (classic inflation factor, but from p-values) to complement yesterday's new functions

# simtrait 1.0.14.9000 (2020-12-02)

* Function `allele_freqs` added option `fold`, which if `TRUE` returns *folded* i.e. *minor* allele frequencies.
  Default is `FALSE`, to return allele frequencies for the alternative allele (the allele counted as it given in the genotype matrix, whether it is the minor or major allele).

# simtrait 1.0.15.9000 (2020-12-07)

* Function `sim_trait` added option `const_herit_loci`, which when `TRUE` constructs causal coefficients as inversely proportional to the square root of `p*(1-p)`, where `p` is the ancestral allele frequency.
  This ensures equal per-causal-locus contribution to trait variance.
  Default draws causal coefficients randomly from a standard normal distribution, rescaled to result in the desired heritability, for unequal per-locus contribution to trait variance.
* Added function `herit_loci`, which calculates per-locus heritabilities based on variance formula (in terms of allele frequencies, coefficients, and overall trait variance factor).
  * Used to validate function `sim_trait` in unit tests.

# simtrait 1.0.16.9000 (2020-12-16)

* Function `sim_trait` with option `const_herit_loci = TRUE` now adds random signs (+/-) to the causal coefficients.
  * Added usage example for this option to vignette.

# simtrait 1.0.17.9000 (2021-01-21)

* Function `sim_trait` now requires that `p_anc` have the same length as the number of loci in `X` (stops with an error otherwise).  Previously this was not checked and could return traits that were `NA` for all individuals without clear indications that anything was wrong.

# simtrait 1.0.18.9000 (2021-02-16)

* Documentation updates:
  - Fixed links to functions, in many cases these were broken because of incompatible mixed Rd and markdown syntax (now markdown is used more fully).

# simtrait 1.0.19.9000 (2021-05-27)

* Function `pval_srmsd` now accepts `causal_indexes = NULL` to handle cases where all p-values are null (before a `NULL` input would cause an error).

# simtrait 1.0.20 (2021-08-10)

- First CRAN submission!
- Major changes
  - Function `sim_trait` renamed option `const_herit_loci` to `fes` (fixed effect sizes), the language used in the paper.
  - Updated documentation generally to agree with paper and to be more concise.
    - "betas" are exclusively referred to as regression coefficients (used to be referred to incorrectly as effect sizes)
	- Package-wide doc now includes all key functions in examples.
- Minor changes
  - Reformatted this `NEWS.md` slightly to improve its automatic parsing.
  - DESCRIPTION
    - Removed `LazyData: true` (to avoid a new "NOTE" on CRAN).
    - Added GitHub links (`URL` and `BugReports`)
  - Fixed spelling in documentation.

# simtrait 1.0.21 (2021-08-12)

- CRAN resubmission
- Added bioRxiv paper reference to description.
- Reset `par()` in vignette examples.

# simtrait 1.0.21.9000 (2022-08-15)

- Function `sim_trait` fixed an important bug resulting in misspecified heritability!
  - The previous buggy version:
    - Non-genetic variance was misspecified, accidentally passing to `rnorm` the desired variance `(1 - herit) * sigma_sq` where the standard deviation (its square root) was required!
    - The resulting effective heritability was given by the requested `herit` value by `herit / ( herit + (1-herit)^2 * sigma_sq )`.  If `sigma_sq = 1` (default), the effective heritability was always larger than desired!
- Vignette updates
  - Added a more detailed comparison between theoretical and empirical covariance matrices from simulations, to help assess the effect of the previous `sim_trait` bug and its fix.
  - Removed `inbr_diag` in `plot_popkin` calls, which in this case made diagonal values larger (as they were larger than 1), among other minor adjustments.
- Updated reference DOI to newest preprint
- README now includes CRAN installation instructions alongside GitHub version.

# simtrait 1.0.22.9000 (2022-08-17)

- Functions `sim_trait`, `sim_trait_mvn`, and `cov_trait`: added parameters `labs` and `labs_sigma_sq` to simulate/model non-genetic group effects.
  - Updated vignette to illustrate this new feature.
  - Function `sim_trait` also now returns `group_effects` as one of the named elements of the return list.
- Function `sim_trait` parameters `p_anc` and `kinship` now have default `NULL` values (used to not have default values), to facilitate scripting in cases that can be either simulated or real genotypes.

# simtrait 1.0.23.9000 (2022-08-18)

- Added functions `pval_type_1_err` and `pval_power_calib` for calculating type I error rates and calibrated power, respectively.

# simtrait 1.0.24.9000 (2022-08-19)

- Functions `pval_type_1_err` and `pval_power_calib`:
  - Now option `alpha` can be a vector, and return value is a vector of estimates, one for each `alpha`.
  - Switched documentation wording from "calculate" to "estimate", and acknowledges that estimates are accurate only when the number of p-values is much larger than `1/alpha`.

# simtrait 1.1.0 (2022-08-22)

- Corrected spelling
- Updated CRAN comments
