# TiDEomics 0.99.7

SIGNIFICANT USER-VISIBLE CHANGES

* `plot_scatter()` now shows the dashed `y = x` line by default only
  when comparing same-named properties from two data frames. This can
  be overridden by `show_diagonal`.
* `DE_timecourse()` uses clearer unified column names in result tables, the
  changes are summarised as below.
  `fits` and `design` retain the original model coefficient names.

  | Test | 0.99.6 result-column name | 0.99.7 result-column name |
  | --- | --- | --- |
  | Factor time effect | `TimeFactor<time>`, or `logFC` when two time points | `logFC_TimeFactor<time>` |
  | Numeric time effect | `TimeNumeric` in `global_time_effect`, or `logFC` in `time_effect` and `reference_time_effect` | `logFC_TimeSlope` |
  | Spline time effect | `splines::ns(...)<basis>` terms converted by `make.names()` | `SplineBasis_1`, `SplineBasis_2`, etc. |
  | Factor group-time interaction | `GroupX.TimeFactor<time>`, or `logFC` when one interaction | `delta_logFC_GroupX_TimeFactor<time>` |
  | Numeric group-time interaction | `GroupX.TimeNumeric`, or `logFC` when two groups | `delta_logFC_GroupX_TimeSlope` |
  | Spline group-time interaction | `GroupX:splines::ns(...)<basis>` terms converted by `make.names()` | `delta_SplineBasis_GroupX_1`, etc. |
  | Group contrast at reference time | `GroupX_at_reference_time`, or `logFC` when two groups | `logFC_GroupX_at_reference_time` |

  `X` represents a non-reference group; `<time>` represents a non-reference
  time level.

# TiDEomics 0.99.6

NEW FEATURES

* TiDEomics now supports multi-dataset / multi-omics analysis:
  * `combine_tide()` combines multiple `prepare_tide()` outputs, identifies
    features shared by the datasets, and returns all-feature filtered data, shared-feature subsets, and filtering summaries.
  * `run_consensus_WGCNA()` performs consensus WGCNA across datasets
    using the shared features.
  * `extract_consensus_hubs()` returns consensus and dataset-specific
    features with highest kME to their modules.
  * Examples and documentations now include a second simulated dataset to
    demonstrate the multi-dataset / multi-omics workflows.

* Some existing functions now accept multi-dataset / multi-omics input:
  * `plot_trend()` supports multiple layouts for multi-dataset input,
    with shared or independent y-axis scales.
  * `plot_pca()` and `plot_umap()` calculate and plot per dataset
    using all available complete features or shared features, and
    combine the plots with shared legends.
  * `plot_cor_matrix()` plots correlation heatmaps per dataset side by side
    with shared legends, and reports pairwise correlations of cophenetic
    distances when samples are shared across datasets.
  * `plot_modules_v()` vertically stacks module heatmaps and mean
    profile views for multi-dataset input.
  * `plot_modules_h()` horizontally aligns module heatmaps and mean
    profile views for multi-dataset input, with the optional enrichment
    annotation on the right and shared legends.
  * `summarise_module_metrics()` accepts both individual and consensus
    WGCNA results, reporting metrics including kME distribution summaries.

* New analysis and visualisation functions:
  * `DE_timecourse()` performs limma-based temporal differential-expression
    analysis within selected groups or jointly across groups.
  * `plot_scatter()` compares per-feature properties between datasets or
    within one table, with quadrant labelling and correlation summaries.
  * `plot_upset_h()` provides Upset plot of feature set intersections with
    optional horizontally aligned per-intersection enrichment textboxes.
  * `plot_upset_v()` provides Upset plot of feature set intersections with
    optional vertically aligned enrichment dot matrices with configurable dot
    size and colour, term ordering, and custom term colours.

SIGNIFICANT USER-VISIBLE CHANGES

* Enrichment workflows now provide more consistent result structures.
  `enrichR_list()` reports feature counts. `enrichGO_rank()` follows
  `enrichGO_list`, returns full results under `all` and optionally adds
  simplified terms under `simplified`.

* WGCNA workflows are improved. In addition to the multi-dataset / multi-omics
  compatibility, `prepare_WGCNA()` now reports a custom power estimate, and
  `plot_modules_h()` has improved display controls, including feature marking
  with custom colours, enrichment term wrapping, and selected module profiles.

* `prepare_tide()` has filtering mode `below_quantile` renamed to `quantile`,
  and reports residual cutoffs.

* `decomp_variance()` now reports numbers of evaluated features, including
  successes, errors, warnings, and skippings.

* `calc_feature_property()` adds `Net_FC`, reporting signed
  abundance difference between last and first time points.

* `theme_custom()` supports font family control.

BUG FIXES

* `plot_WGCNA()` now plots dendrogram for every feature block.
* `prepare_tide()` now accepts `feature_threshold = NULL`.
* Enrichment legends in `plot_modules_h()` and `plot_upset_h()` now
  omit unused colours.
* `enrichGO_list()` now omits universe when no background is supplied,
  rather than passing an explicit NULL to `clusterProfiler::enrichGO()`.

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
