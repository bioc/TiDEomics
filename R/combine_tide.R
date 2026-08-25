#' Combine multiple `prepare_tide()` outputs
#'
#' Combine a list of `prepare_tide()` output objects for multi-dataset or
#'   multi-omics analysis.
#'
#' @param tide_list A non-empty named list of results from `prepare_tide()`.
#' @param assay The assay to use: a numeric index or character name, e.g. 1 or
#'   "orig" for original data, 2 or "norm" for time 0 normalised data.
#'
#' @returns A list containing the filtered and merged SE objects (all features
#'   and shared features), variance decomposition, shared feature names,
#'   the resolved assay name, and a per-dataset filtering summary.
#' @export
#' @examples
#' data(tutorial_data)
#' data(tutorial_sample_info)
#' tide <- prepare_tide(tutorial_data, tutorial_sample_info,
#'     keep = "threshold", residual_threshold = 100)
#' combined <- combine_tide(list(RNA = tide, Protein = tide), assay = "norm")
#' combined$filter_summary
combine_tide <- function(tide_list, assay = 1) {
    if (!is.list(tide_list) || length(tide_list) == 0L) {
        stop("'tide_list' must be a non-empty list of prepare_tide() results.")
    }
    if (is.null(names(tide_list)) || anyNA(names(tide_list)) ||
            any(names(tide_list) == "") || anyDuplicated(names(tide_list))) {
        stop("'tide_list' must have unique, non-empty names.")
    }

    required <- c("se", "se_filtered", "merged_se_filtered",
        "variance", "filter_summary")
    for (nm in names(tide_list)) {
        tide <- tide_list[[nm]]
        if (!is.list(tide) || !all(required %in% names(tide))) {
            stop("'tide_list[[", nm,"]] must be a result from prepare_tide().")
        }
        .check_se(tide$se, paste0("tide_list[[", nm,"]]$se"))
        .check_se_list(tide$se_filtered,
            paste0("tide_list[[", nm,"]]$se_filtered"))
        .check_se_list(tide$merged_se_filtered,
            paste0("tide_list[[", nm,"]]$merged_se_filtered"))
        if (!is.list(tide$variance) || length(tide$variance) == 0L) {
            stop("'tide_list[[", nm,"]]$variance' must be a non-empty list.")
        }
        if (!is.data.frame(tide$filter_summary)) {
            stop("'tide_list[[", nm,"]]$filter_summary' must be a data.frame.")
        }
    }

    assay_names <- names(tide_list[[1L]]$se_filtered)
    if (is.null(assay_names) || anyNA(assay_names) ||
            any(assay_names == "") || anyDuplicated(assay_names)) {
        stop("Selected assay lists must have unique, non-empty names.")
    }
    for (nm in names(tide_list)) {
        if (!identical(names(tide_list[[nm]]$se_filtered), assay_names) ||
            !identical(names(tide_list[[nm]]$merged_se_filtered),
                assay_names) ||
            !identical(names(tide_list[[nm]]$variance), assay_names)) {
            stop("All prepare_tide results must have the same assay names.")
        }
    }

    if (is.character(assay)) {
        if (length(assay) != 1L || is.na(assay) || !nzchar(assay) ||
                !assay %in% assay_names) {
            stop("'assay' must be one of: ",
                paste(assay_names, collapse = ", "), ".")
        }
        assay_name <- assay
    } else if (is.numeric(assay)) {
        if (length(assay) != 1L || is.na(assay) || !is.finite(assay) ||
                assay != as.integer(assay) || assay < 1L ||
                assay > length(assay_names)) {
            stop("Numeric 'assay' must be an integer from 1 to ",
                length(assay_names), ".")
        }
        assay_name <- assay_names[as.integer(assay)]
    } else {
        stop("'assay' must be a single assay name or numeric position.")
    }

    se_list_filtered <- lapply(tide_list,
        function(tide) tide$se_filtered[[assay_name]])
    se_list_merged <- lapply(tide_list,
        function(tide) tide$merged_se_filtered[[assay_name]])
    variance_decomp <- lapply(tide_list,
        function(tide) tide$variance[[assay_name]])

    feature_sets <- lapply(se_list_filtered, rownames)
    if (any(vapply(feature_sets, is.null, logical(1)))) {
        stop("Filtered SummarizedExperiment objects must have row names.")
    }
    shared_features <- Reduce(intersect, feature_sets)
    if (length(shared_features) == 0L) {
        stop("No features are shared across the selected filtered datasets.")
    }

    se_list_filtered_shared <- lapply(se_list_filtered,
        function(se) se[shared_features, , drop = FALSE])
    se_list_merged_shared <- lapply(se_list_merged,
        function(se) se[shared_features, , drop = FALSE])

    filter_summary <- do.call(rbind, lapply(names(tide_list), function(nm) {
        tide <- tide_list[[nm]]
        selected <- tide$filter_summary[
            tide$filter_summary$stage == "residual_filter" &
                tide$filter_summary$assay == assay_name, , drop = FALSE]
        n_kept <- nrow(se_list_filtered[[nm]])
        data.frame(
            Dataset = nm,
            n_total = nrow(tide$se),
            n_kept = n_kept,
            pct_kept = if (nrow(selected) == 1L) selected$pct_kept[1] else
                round(100 * n_kept / nrow(tide$se), 1),
            residual_cutoff = if (nrow(selected) == 1L)
                selected$residual_cutoff[1] else NA_real_,
            stringsAsFactors = FALSE)
    }))
    row.names(filter_summary) <- NULL

    message("Combining ", length(tide_list),
        " prepare_tide() results using assay '", assay_name, "'")
    message("Shared filtered features: ", length(shared_features),
        " (", paste(names(feature_sets), lengths(feature_sets),
            sep = "=", collapse = ", "), ")")

    list(
        se_list_filtered = se_list_filtered,
        se_list_merged = se_list_merged,
        se_list_filtered_shared = se_list_filtered_shared,
        se_list_merged_shared = se_list_merged_shared,
        variance_decomp = variance_decomp,
        shared_features = shared_features,
        assay = assay_name,
        filter_summary = filter_summary
    )
}
