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
