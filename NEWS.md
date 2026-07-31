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
* `plot_cor_matrix()`: default `title` changed from
  `"Correlation between samples"` to `NULL`; the title is now auto-generated
  as `"Sample correlation (<Method>)"` with the method capitalised (e.g.
  "Sample correlation (Spearman)").
* `plot_pca()` and `plot_umap()`: when `circle = TRUE`, auto-computed axis
  limits now use 0.35x padding around the data range (via `coord_cartesian()`)
  instead of 1.5x the min/max (via `xlim()`/`ylim()`).

# TiDEomics 0.99.4

* Accepted by Bioconductor.
