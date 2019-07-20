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

