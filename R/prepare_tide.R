#' Prepare TiDEomics input
#'
#' One-step data preparation wrapper for TiDEomics: creates the
#'   SummarizedExperiment object, normalises to starting time point,
#'   filters out group-specific features, runs variance
#'   decomposition, and filters out noisy features with high residual
#'   variance. Returns both unmerged and merged
#'   (mean of replicates) outputs, with all-feature and filtered
#'   versions of each, ready for downstream analysis.
#'
#' The pipeline:
#'   1. `create_input()`: validate and build SE
#'   2. `normalise_to_start()`: add time-0-normalised assay
#'   3. `split_groups()` -> `merge_replicates()` -> `calc_feature_property()`
#'      -> `summarise_feature_property()`: per-feature expression statistics
#'   4. `group_specific_features()`: filter out group-specific features
#'   5. `decomp_variance()` : LMM variance decomposition on both `"orig"`
#'      and `"norm"` assays
#'   6. Residual filter: filter out features with high residual
#'      variance, calculated from `decomp_variance()` on each assay
#'   7. Produce filtered-unmerged and filtered-merged SEs for each assay
#'
#' @param data A data frame with rows as features (e.g., genes, proteins)
#' and columns as samples. The first column should contain feature identifiers
#' (e.g., gene symbols) and be named 'Feature'. The data should have been
#' normalised and log transformed. Missing values are allowed.
#' @param sample_ann A data frame containing sample annotations with required
#' columns: 'Sample', 'Group', and 'Time'. 'Replicate' and 'Batch' are optional
#' columns and will be auto-generated if not provided. 'Time',
#' 'Replicate' and 'Batch' should be numeric.
#' @param subject_col Optional: name of a column in `sample_ann` identifying
#' biological subjects measured repeatedly across time points (e.g.,
#' `"PatientID"`). If provided, the column is renamed to 'Subject' and used
#' for repeated-measures analyses downstream. If NULL (default), all samples
#' are treated as independent, e.g. for cell culture experiments.
#' @param replicate_col Optional: name of a column in `sample_ann` identifying
#' replicate IDs. If provided, the column is renamed to 'Replicate'. If NULL
#' (default), the function looks for a column named 'Replicate'; if absent,
#' replicate IDs are auto-generated within each Group and Time (and Subject,
#' if provided).
#' @param batch_col Optional: name of a column in `sample_ann` identifying
#' batch information. If provided, the column is renamed to 'Batch'. If NULL
#' (default), the function looks for a column named 'Batch'; if absent, all
#' samples are assigned to batch 1.
#' @param feature_threshold Threshold for `calc_feature_property()`. Feature
#'   values > feature_threshold are considered expressed (default: 0).
#' @param filter_ratio Minimum ratio of time points a feature must be expressed
#'   in, to be considered expressed in a group (default: 0.5, meaning that a
#'   feature must be > feature_threshold in >=50% time points in a group to be
#'   considered expressed). Passed to `group_specific_features()`.
#' @param min_groups Minimum number of groups a feature must be expressed in,
#'   with the above filter_ratio, to pass the cross group filter
#'   (default: 2). Passed to `group_specific_features()`.
#' @param keep How to select features after variance decomposition:
#'   `"below_quantile"` (default): keep features with residual variance below
#'   `residual_threshold` quantile (0-1 scale).
#'   `"top_n"`: keep the top `n_keep` lowest-residual features.
#'   `"threshold"`: keep features with residual percentage below
#'   `residual_threshold` (0-100 scale; use 100 to skip filtering).
#' @param residual_threshold Numeric: when `keep = "below_quantile"`,
#'   the quantile of residual variance to use as cutoff (0-1 scale).
#'   When `keep = "threshold"`, the percentage cutoff
#'   (0-100 scale, 100 = keep all features).
#'   (default: NULL, auto-assigned as 0.75 for "below_quantile",
#'   100 for "threshold")
#' @param n_keep Number of features to keep when `keep = "top_n"`
#'   (default: 5000).
#' @param interaction Passed to `decomp_variance()`: if TRUE, adds
#'   `(1|Group:Time)` (or
#'   `(1|Subject:Time)` when Subject present) interaction term (default: FALSE).
#' @param n_cores Number of cores for variance decomposition (default: 1).
#'
#' @returns A named list with elements. The elements suffixed with
#'   "_filtered" have passed both the cross-group and residual filters,
#'   and are split into $orig and $norm sub-lists,
#'   with variance decomposition and residual filtering done per assay.
#'   1. `se` (SE with replicates, all features, both assays);
#'   2. `se_filtered` (a list of two SEs with replicates,
#'     filtered features by `orig` / `norm`, use for WGCNA);
#'   3. `merged_list` (per-group merged SEs, all features, use for Trendy);
#'   4. `merged_list_filtered` (a list of two per-group merged SEs,
#'     filtered features by `orig` / `norm`);
#'   5. `merged_se` (single merged SE, all groups combined, all features,
#'     use for plot_modules_v/h);
#'   6. `merged_se_filtered` (a list of two single merged SEs, filtered features
#'     by `orig` / `norm`, use for plot_modules_v/h);
#'   7. `variance` (a list of two variance decomposition data.frames,
#'     computed from `orig` / `norm` assay);
#'   8. `feature_property` (feature property summary, all features);
#'   9. `filter_summary` (per-stage filtering statistics);
#'   10. `group_specific_filter` (a vector of group-specific features,
#'     removed by the cross-group filter);
#'   11-13. `DE`, `enrichment`, `WGCNA` (empty on initialization).
#'
#' @export
#' @examples
#' data(tutorial_data)
#' data(tutorial_sample_info)
#' tide <- prepare_tide(tutorial_data, tutorial_sample_info,
#'     keep = "threshold", residual_threshold = 100)
#' tide$filter_summary
prepare_tide <- function(
    data,
    sample_ann,
    subject_col = NULL,
    replicate_col = NULL,
    batch_col = NULL,
    feature_threshold = 0,
    filter_ratio = 0.5,
    min_groups = 2,
    keep = c("below_quantile", "top_n", "threshold"),
    residual_threshold = NULL,
    n_keep = 5000,
    interaction = FALSE,
    n_cores = 1
) {
    keep <- match.arg(keep)
    if (is.null(residual_threshold)) {
        residual_threshold <- if (keep == "below_quantile") 0.75 else 100
    }

    .check_numeric(feature_threshold, "feature_threshold")
    .check_df(data, "data")
    .check_df(sample_ann, "sample_ann")
    .check_pval(filter_ratio, "filter_ratio")
    .check_positive_int(min_groups, "min_groups")
    if (keep == "below_quantile") {
        .check_pval(residual_threshold, "residual_threshold")
    } else if (keep == "threshold") {
        if (residual_threshold < 0 || residual_threshold > 100) {
            stop("'residual_threshold' must be 0-100 when keep = 'threshold'.")
        }
    }
    .check_positive_int(n_keep, "n_keep")
    .check_logical(interaction, "interaction")
    .check_positive_int(n_cores, "n_cores")

    # Build SE and normalise
    se <- create_input(data, sample_ann,
        subject_col = subject_col,
        replicate_col = replicate_col,
        batch_col = batch_col
    )

    se <- normalise_to_start(se)

    n_total <- nrow(se)
    n_groups <- length(unique(se$Group))
    if (min_groups > n_groups) {
        stop("min_groups (", min_groups, ") exceeds the number of groups (",
            n_groups, ").")
    }
    message("Preparing TiDEomics input: ", n_total, " features, ",
            ncol(se), " samples, ", n_groups, " groups")

    # Feature detection statistics
    se_list <- split_groups(se)
    merged_list_all <- merge_replicates(se_list)
    merged_list_all <- calc_feature_property(merged_list_all,
        threshold = feature_threshold)
    prop <- summarise_feature_property(merged_list_all)

    # Cross-group filter via group_specific_features
    cross_pass <- group_specific_features(prop,
        groups = NULL,
        filter_ratio = filter_ratio,
        group_pct = min_groups / n_groups,
        genename = FALSE,
        GO = FALSE
    )$features

    n_cross_pass <- length(cross_pass)
    message("Cross-group filter (Exp_ratio >= ", filter_ratio,
            " in >= ", min_groups, " groups): kept ",
            n_cross_pass, " of ", n_total, " features (",
            round(100 * n_cross_pass / n_total, 1), "%)")

    if (is.null(cross_pass) || n_cross_pass == 0) {
        stop("No features pass the cross-group filter. ",
            "Relax filter_ratio or min_groups.")
    }

    se_cross <- se[cross_pass, ]

    # Variance decomposition and residual filter on both assays
    .run_one_assay <- function(assay_name, se_cross) {
        message("--- assay: ", assay_name, " ---")
        var_decomp <- decomp_variance(se_cross,
            assay = assay_name,
            interaction = interaction,
            core = n_cores
        )

        res_var <- var_decomp$Residual
        names(res_var) <- var_decomp$Feature
        res_var <- res_var[!is.na(res_var)]

        if (keep == "below_quantile") {
            cutoff <- stats::quantile(res_var, probs = residual_threshold,
                na.rm = TRUE)
            keep_features <- names(res_var)[res_var <= cutoff]
            residual_cutoff <- as.numeric(cutoff)
        } else if (keep == "threshold") {
            keep_features <- names(res_var)[res_var < residual_threshold]
            residual_cutoff <- residual_threshold
        } else {
            keep_features <- names(sort(res_var))[seq_len(
                min(n_keep, length(res_var)))]
            residual_cutoff <- max(res_var[keep_features], na.rm = TRUE)
        }

        n_resid_pass <- length(keep_features)
        message("Residual filter (", keep, "): kept ", n_resid_pass,
                " of ", nrow(se_cross), " features (",
                round(100 * n_resid_pass / nrow(se_cross), 1), "%)")

        se_f <- se[keep_features, ]
        merged_f <- merge_replicates(split_groups(se_f))
        list(
            variance = var_decomp,
            se_filtered = se_f,
            merged_list_filtered = merged_f,
            merged_se_filtered = merge_groups(merged_f),
            n_resid_pass = n_resid_pass,
            residual_cutoff = residual_cutoff
        )
    }

    res_orig <- .run_one_assay("orig", se_cross)
    res_norm <- .run_one_assay("norm", se_cross)

    # Shared outputs
    merged_se_all <- merge_groups(merged_list_all)

    .format_cutoff <- function(cutoff) {
        if (keep == "threshold" && residual_threshold >= 100) {
            "skipped (threshold = 100)"
        } else {
            paste0(keep, " = ", round(cutoff, 3))
        }
    }

    filter_summary <- data.frame(
        stage = c("input", "cross_group_filter",
            "residual_filter", "residual_filter"),
        assay = c(NA_character_, NA_character_, "orig", "norm"),
        n_features = c(n_total, n_cross_pass,
            res_orig$n_resid_pass, res_norm$n_resid_pass),
        pct_kept = c(100,
            round(100 * n_cross_pass / n_total, 1),
            round(100 * res_orig$n_resid_pass / n_total, 1),
            round(100 * res_norm$n_resid_pass / n_total, 1)),
        threshold = c(NA_character_,
            paste0("Exp_ratio >= ", filter_ratio,
                " in >= ", min_groups, " groups"),
            .format_cutoff(res_orig$residual_cutoff),
            .format_cutoff(res_norm$residual_cutoff)),
        stringsAsFactors = FALSE
    )

    list(
        se                   = se,
        se_filtered          = list(orig = res_orig$se_filtered,
                                    norm = res_norm$se_filtered),
        merged_list          = merged_list_all,
        merged_list_filtered = list(orig = res_orig$merged_list_filtered,
                                    norm = res_norm$merged_list_filtered),
        merged_se            = merged_se_all,
        merged_se_filtered   = list(orig = res_orig$merged_se_filtered,
                                    norm = res_norm$merged_se_filtered),
        variance             = list(orig = res_orig$variance,
                                    norm = res_norm$variance),
        feature_property     = prop,
        filter_summary       = filter_summary,
        group_specific_filter = setdiff(rownames(se), cross_pass),
        DE                   = list(),
        enrichment           = list(),
        WGCNA                = list()
    )
}
