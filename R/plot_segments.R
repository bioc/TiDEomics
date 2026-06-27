#' Plot segmented regression
#'
#' @description Plot segmented regression results for specified features
#' in each group, using the `Trendy::plotFeature()` function. The input is
#' a list of SummarizedExperiment objects with no missing values and a list
#' of Trendy analysis results for each group.
#'
#' @param se_obj_imp_list A list of SummarizedExperiment objects with no
#' missing values
#' @param res_list A list of the Trendy analysis results, output of
#' `run_Trendy()`
#' @param feature A character vector of feature names to be plotted
#' @param group A character vector of group names to be analysed. If NULL,
#' all groups in the input list will be used (default is NULL)
#' @param nrow Number of rows in the plot layout. If NULL and multiple features
#' are provided, it will be set to half the number of features (rounded up)
#' (default is NULL)
#' @param ... Additional arguments to be passed to the `Trendy::plotFeature()`
#'
#' @import SummarizedExperiment
#'
#' @returns Plots of the segmented regression fits for the specified features
#' in each group.
#' @export
#' @examples
#' data(example_obj)
#' example_obj <- normalise_to_start(example_obj)
#' example_obj_list <- split_groups(example_obj)
#' example_obj_merged_list <- merge_replicates(example_obj_list)
#' example_obj_merged_list <-
#'     calc_feature_property(example_obj_merged_list, threshold = 0)
#' # no missing value in the example dataset, so imputation is not necessary
#' example_obj_merged_imp_list <- impute_groups(example_obj_merged_list)
#'
#' data("example_res_list")
#' plot_segments(example_obj_merged_imp_list, example_res_list,
#'     feature = c("Mctp1"))
plot_segments <- function(se_obj_imp_list, res_list, feature,
    group = NULL, nrow = NULL, ...) {
    .check_character(feature, "feature")
    .check_se_list(se_obj_imp_list, "se_obj_imp_list")
    .check_se_list_merged(se_obj_imp_list, "se_obj_imp_list")
    .check_se_list_no_na(se_obj_imp_list, "se_obj_imp_list")
    .check_list(res_list, "res_list", "run_Trendy")
    if (is.null(group)) {
        group <- names(se_obj_imp_list) |>
            intersect(names(res_list))
    } else if (!all(group %in% names(se_obj_imp_list))) {
        stop("At least one of the specified groups is not found in the ",
        "input object list.")
    }

    for (i in group) {
        message("Plotting segmented regression for group: ", i)
        .plot_segments_one_group(se_obj_imp_list[[i]], res_list[[i]],
            feature, group = i, nrow = nrow, ...)
    }
}


#' Plot segmented regression (one group)
#'
#' @description Plot segmented regression results for specified features
#'
#' @param se_obj_imp A SummarizedExperiment object with no missing values
#' @param res Result of the Trendy analysis, output of `run_Trendy()`
#' @param feature A character vector of feature names to be plotted
#' @param group A character string specifying the group name (for title
#' purposes)
#' @param nrow Number of rows in the plot layout. If NULL and multiple features
#' are provided, it will be set to half the number of features (rounded up)
#' (default is NULL)
#' @param ... Additional arguments to be passed to the `Trendy::plotFeature()`
#'
#' @import SummarizedExperiment
#'
#' @returns A plot of the segmented regression fits for the specified features.
#' @keywords internal
.plot_segments_one_group <- function(se_obj_imp, res, feature,
    group, nrow = NULL, ...) {
    if (length(feature) > 1) {
        if (is.null(nrow)) {
            nrow <- ceiling(length(feature) / 2)
            ncol <- 2
        } else {
            ncol <- ceiling(length(feature) / nrow)
        }
    } else {
        nrow <- 1
        ncol <- 1
    }

    M_imp <- assays(se_obj_imp)[[1]]
    time.vector <- colData(se_obj_imp)$Time

    graphics::par(mfrow = c(nrow, ncol))
    Trendy::plotFeature(
        Data = M_imp[feature, ],
        tVectIn = time.vector,
        featureNames = feature,
        trendyOutData = res,
        ...
    )
}
