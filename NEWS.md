# TiDEomics 0.99.5

SIGNIFICANT USER-VISIBLE CHANGES

* `enrich_p_threshold` renamed to `enrich_threshold` in `plot_modules_h()`.
* `assay` parameter in `plot_modules_h()` and `plot_modules_v()` now
  must be specified explicitly (no default).
* `prepare_WGCNA()`: `assay` parameter must be specified explicitly
  (no default).
* `plot_trend()`: `ylab` parameter renamed to `ylabel`, default changed from
  `"Abundance"` to `"Log2 abundance"`.
* `plot_modules_v()`: default `ylabel` changed from
  `"Log2 abundance normalised to Time 0"` to `"Log2 abundance"`.

# TiDEomics 0.99.4

* Accepted by Bioconductor.
