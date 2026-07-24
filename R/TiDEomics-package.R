#' TiDEomics: Time-course Differential Expression analysis of omics data
#'
#' @description
#' TiDEomics provides a workflow for **multi-group
#' time-course** omics data analysis, analysing
#' **time-dominant**, **group-dominant**, and **group-specific
#' temporal** effects through pairwise differential expression,
#' variance decomposition, and co-expression module analysis (WGCNA).
#' A **residual variance** filtering strategy prioritises features
#' with structured differential expression.
#' It supports datasets with missing values (e.g. mass
#' spectrometry-based proteomics) and operates on
#' `SummarizedExperiment` objects for
#' compatibility with the Bioconductor ecosystem.
#'
#' @details
#' Key functionality:
#' - **Data preparation & normalisation:** [create_input()],
#'   [prepare_tide()], [split_groups()], [merge_replicates()],
#'   [merge_groups()], [normalise_to_start()], [impute_groups()]
#' - **Quality control & exploration:** [plot_missing()],
#'   [plot_distribution()], [plot_ID()], [plot_cv()],
#'   [plot_cor_matrix()], [plot_pca()], [plot_pca_3D()],
#'   [plot_pca_arrows()], [plot_pca_by_group()],
#'   [plot_umap()], [plot_umap_by_group()]
#' - **Feature properties & variance decomposition:** [calc_feature_property()],
#'   [summarise_feature_property()], [group_specific_features()],
#'   [decomp_variance()], [plot_variance()], [plot_trend()]
#' - **Pairwise differential expression:** [DE_between_group()],
#'   [DE_between_time()], [plot_DE_between_group()],
#'   [plot_DE_between_time()], [plot_volcano()]
#' - **Segmentation regression with Trendy:** [run_Trendy()],
#'   [summarise_Trendy()], [plot_segments()], [plot_breakpoints()]
#' - **Co-expression module identification with WGCNA:** [prepare_WGCNA()],
#'   [run_WGCNA()], [WGCNA_module()], [extract_hubs()], [plot_WGCNA()],
#'   [plot_modules_h()], [plot_modules_v()],
#'   [summarise_module_pattern()], [summarise_module_metrics()]
#' - **Functional enrichment:** [enrichGO_list()], [enrichGO_rank()],
#'   [enrichR_list()], [enrich_msigdb()], [plot_GO()]
#' - **Interoperability & settings:** [flatten_DE()],
#'   [flatten_enrich()], [set_custom_palette()],
#'   [get_custom_palette()], [theme_custom()]
#'
#' Two experimental designs are auto-detected from the sample
#' annotation: independent samples and repeated measures.
#'
#' See `vignette("TiDEomics")` for a step-by-step tutorial.
#'
#' @seealso
#' Useful links:
#' - \url{https://hte123.github.io/TiDEomics}
#' - \url{https://github.com/hte123/TiDEomics}
#' - Report bugs at \url{https://github.com/hte123/TiDEomics/issues}
#'
#' @import SummarizedExperiment
#' @import ggplot2
#' @import patchwork
#' @keywords internal
"_PACKAGE"
