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
#' @importFrom circlize colorRamp2
#' @importFrom clusterProfiler bitr cnetplot dotplot emapplot enrichGO gseGO
#' @importFrom clusterProfiler merge_result simplify
#' @importFrom ComplexHeatmap Heatmap HeatmapAnnotation Legend anno_block
#' @importFrom ComplexHeatmap anno_link anno_mark anno_zoom draw packLegend
#' @importFrom ComplexHeatmap rowAnnotation
#' @importFrom dplyr across all_of anti_join arrange bind_rows count desc
#' @importFrom dplyr distinct everything filter full_join group_by mutate n
#' @importFrom dplyr n_distinct pull rename row_number select slice slice_head
#' @importFrom dplyr slice_min summarise ungroup
#' @importFrom enrichplot pairwise_termsim
#' @importFrom ggforce geom_mark_ellipse
#' @importFrom ggh4x elem_list_rect facet_grid2 strip_themed
#' @importFrom ggplotify as.ggplot
#' @importFrom ggpubr annotate_figure ggarrange stat_compare_means text_grob
#' @importFrom ggridges geom_density_ridges
#' @importFrom ggrepel geom_text_repel
#' @importFrom ggsci pal_iterm pal_jco pal_simpsons
#' @importFrom limma contrasts.fit duplicateCorrelation eBayes lmFit
#' @importFrom limma makeContrasts topTable
#' @importFrom lme4 VarCorr lmer lmerControl
#' @importFrom pbapply pblapply pboptions
#' @importFrom PCAtools getComponents pairsplot pca plotloadings screeplot
#' @importFrom randtests bartels.rank.test
#' @importFrom scales pal_hue pal_viridis
#' @importFrom tibble column_to_rownames rownames_to_column
#' @importFrom tidyr complete pivot_longer pivot_wider unite
#' @importFrom Trendy breakpointDist formatResults plotFeature results
#' @importFrom Trendy topTrendy trendy
#' @importFrom umap umap
#' @importFrom WGCNA allowWGCNAThreads binarizeCategoricalColumns
#' @importFrom WGCNA blockwiseModules cor corPvalueStudent goodGenes
#' @importFrom WGCNA goodSamples labeledHeatmap labels2colors moduleEigengenes
#' @importFrom WGCNA orderMEs pickSoftThreshold plotDendroAndColors plotMEpairs
#' @importFrom WGCNA signedKME
#' @keywords internal
"_PACKAGE"
