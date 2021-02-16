# 2019-02-13 - simtrait 1.0.0.9000

* Public GitHub release!

# 2019-03-26 - simtrait 1.0.1.9000

* Fixed `m_causal=1` bug (used to die, now works correctly).

# 2019-04-10 - simtrait 1.0.2.9000

* Renamed argument `sigmaSq` to `sigma_sq` (parameter of `sim_trait` and `cov_trait`).
* `sim_trait` return list has more descriptive names:
  * `trait` (old `y`)
  * `causal_indexes` (old `i`)
  * `causal_coeffs` (old `beta`)
* Overall coding style changes.
* Vignette changes in parallel with popkin 1.2.0.9000 changes.

# 2019-04-16 - simtrait 1.0.3.9000

* Vignette and README changes in parallel with bnpsd 1.1.0.9000 changes.
* 2019-05-13: added ORCID to author info

# 2019-07-19 - simtrait 1.0.4.9000

* Added `allele_freqs`, which handles `BEDMatrix` objects correctly.
  Rest of `sim_trait` doesn't handle `BEDMatrix` objects yet.

# 2019-07-20 - simtrait 1.0.5.9000

* Now `sim_trait`, and consequently the whole package, supports `BEDMatrix` objects.
  This is ideal for simulating traits from real and potentially very large genotype data.

# 2019-08-02 - simtrait 1.0.6.9000

* Improved `allele_freqs` to load genotype chunks from a `BEDMatrix` object more efficiently, using as much available memory as possible.
  Requires an updated `popkin` package (>= 1.2.6.9000).

# 2019-12-17 - simtrait 1.0.7.9000

* Preemptively updated `class` usage now that matrices return a two-element array in R-devel
* Moved logo to `man/figures/`
* Minor Roxygen-related updates

# 2020-07-07 - simtrait 1.0.8.9000

* Function `sim_trait` changed default `maf_cut` from 0.05 to `NA`, and in that case code no longer computes marginal allele frequencies unless they are needed because `p_anc` is missing (and `kinship` is provided).
  This is much faster and memory efficient in extremely large simulations, and just makes more sense to me generally (the MAF threshold was a quirky feature, still available but non-default now).

# 2020-07-20 - simtrait 1.0.9.9000

* Function `sim_trait` now delays calculating allele frequencies if `maf_cut = NA` and `p_anc` is missing (and `kinship` is provided), so that these frequencies are only calculated on the small subset of loci selected to be causal (rather than the whole genome, which was the original behavior).
  This is expected to speed-up trait simulations from real genotypes.
* Vignette had minor corrections and edits (focused on model and algorithm description).

# 2020-08-05 - simtrait 1.0.10.9000

* Added function `sim_trait_mvn`, which draws traits from the multivariate normal distribution with covariance structure matching that of genetic traits (also called "infinitesimal" model).
  These traits are useful for testing heritability estimation.
  A visual validation of these simulated traits is available in the vignette.

# 2020-10-14 - simtrait 1.0.11.9000

* Function `alelle_freq` was modified to reduce memory usage for the BEDMatrix case.
  - For both this function and `sim_trait`, the BEDMatrix-specific options `mem_factor` and `mem_lim` were replaced by the option `m_chunk_max`.
  - As a consequence of the changes, `popkin` is no longer a dependency of this `simtrait` package.
    However, `popkin` is still recommended for estimating the kinship matrices required by some of the functions of this `simtrait` package.

# 2020-11-05 - simtrait 1.0.12.9000

* Added functions used for benchmarking genetic association methods based on traits simulated with `simtrait`:
  - `rmsd`: General root-mean-square deviation (RMSD) function
  - `pval_srmsd`: Signed RMSD between observed and expected (uniform) p-values
  - `pval_aucpr`: Area under the curve (AUC) of precision-recall (PR) curve
    - Added `PRROC` package dependency

# 2020-11-06 - simtrait 1.0.13.9000

* Added function `pval_infl` (classic inflation factor, but from p-values) to complement yesterday's new functions

# 2020-12-02 - simtrait 1.0.14.9000

* Function `allele_freqs` added option `fold`, which if `TRUE` returns *folded* i.e. *minor* allele frequencies.
  Default is `FALSE`, to return allele frequencies for the alternative allele (the allele counted as it given in the genotype matrix, whether it is the minor or major allele).

# 2020-12-07 - simtrait 1.0.15.9000

* Function `sim_trait` added option `const_herit_loci`, which when `TRUE` constructs causal coefficients as inversely proportional to the square root of `p*(1-p)`, where `p` is the ancestral allele frequency.
  This ensures equal per-causal-locus contribution to trait variance.
  Default draws causal coefficients randomly from a standard normal distribution, rescaled to result in the desired heritability, for unequal per-locus contribution to trait variance.
* Added function `herit_loci`, which calculates per-locus heritabilities based on variance formula (in terms of allele frequencies, coefficients, and overall trait variance factor).
  * Used to validate function `sim_trait` in unit tests.

# 2020-12-16 - simtrait 1.0.16.9000

* Function `sim_trait` with option `const_herit_loci = TRUE` now adds random signs (+/-) to the causal coefficients.
  * Added usage example for this option to vignette.

# 2021-01-21 - simtrait 1.0.17.9000

* Function `sim_trait` now requires that `p_anc` have the same length as the number of loci in `X` (stops with an error otherwise).  Previously this was not checked and could return traits that were `NA` for all individuals without clear indications that anything was wrong.

# 2021-02-16 - simtrait 1.0.18.9000

* Documentation updates:
  - Fixed links to functions, in many cases these were broken because of incompatible mixed Rd and markdown syntax (now markdown is used more fully).
