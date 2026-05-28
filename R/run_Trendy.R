#' Segmented regression analysis
#'
#' @description Run segmented regression analysis with Trendy on imputed
#' data for multiple groups of samples in a list of SummarizedExperiment
#' objects.
#' The results will be returned in a list format. Each element in the list
#' corresponds to one group of samples and contains the Trendy analysis results
#' for that group.
#'
#' If the `feature` parameter is not specified, all features expressed in at
#' least `minExp` time points in each group will be included in the analysis
#' for the
#' group. `calc_feature_property()` should be run before imputation to
#' calculate the expression ratio for each feature in each group.
#' Missing values can be
#' imputed after `calc_feature_property()` using the `impute_groups()` function.
#'
#' If a group has less than 4 time points are available, Trendy analysis cannot
#' be performed and NULL will be returned for the group.
#'
#' @param se_obj_imp_list A list of SummarizedExperiment objects with no
#' missing values
#' @param group A character vector of group names to be analysed. If NULL, all
#' groups in the input list will be used (default is NULL)
#' @param feature A character vector of feature names to be included in the
#' analysis. If NULL, all features with expression ratio (Exp_ratio) >= minExp
#' will be
#' used, based on results of `calc_feature_property()` before imputation.
#' (default is NULL)
#' @param minExp Minimum expression ratio for a feature to be included in
#' the analysis (default is 0.5)
#' @param maxK Parameter of `Trendy::trendy()`, maximum number of breakpoints
#' allowed in the segmented regression model (default is 1)
#' @param meanCut Parameter of `Trendy::trendy()`, minimum mean expression
#' required for a feature to be included in the analysis (default is 0)
#' @param minNumInSeg Parameter of `Trendy::trendy()`, minimum number of
#' samples required in each segment (default is 3)
#' @param NCores Number of cores to use for parallel processing (default is 1)
#' @param ... Additional arguments to be passed to the `Trendy::trendy()`
#'
#' @import SummarizedExperiment
#'
#' @returns A list of Trendy analysis results, including the fitted model
#' parameters and statistics for each feature in each group.
#' @export
#' @examples
#' data("example")
#' example_obj <- normalise_to_start(example_obj)
#' example_obj_list <- split_groups(example_obj)
#' example_obj_merged_list <- merge_replicates(example_obj_list)
#' example_obj_merged_list <- calc_feature_property(example_obj_merged_list,
#'     threshold = 0)
#' # no missing value in the example dataset, so imputation is not necessary
#' example_obj_merged_imp_list <- impute_groups(example_obj_merged_list)
#'
#' # "untreated" group has only 3 time points, so Trendy analysis will not be
#' # performed for this group
#' example_res_list <- run_Trendy(example_obj_merged_imp_list, maxK = 1,
#'     minNumInSeg = 2, meanCut = 0)
#' # usethis::use_data(example_res_list)
#'
#' # plot_segments(example_obj_merged_imp_list, example_res_list,
#' #     feature = c("Mctp1"))
#' # plot_breakpoints(example_res_list)
#' # trendy_summary <- summarise_Trendy(example_res_list)
#' # trendy_list <- extract_segment_trends(trendy_summary)
run_Trendy <- function(se_obj_imp_list, group = NULL,
    feature = NULL, minExp = 0.5, maxK = 1,
    meanCut = 0, minNumInSeg = 3, NCores = 1, ...) {
    if (is.null(group)) {
        group <- names(se_obj_imp_list)
    } else if (!all(group %in% names(se_obj_imp_list))) {
        stop("At least one of the specified groups is not found in the ",
        "input object list.")
    }

    if (is.null(feature)) {
        for (i in group) {
            if (!"Exp_ratio" %in% colnames(rowData(se_obj_imp_list[[i]]))) {
                stop("Feature not specified and Exp_ratio haven't been ",
                "calculated. To run the function on all features filtered by ",
                "Exp_ratio, please run `calc_feature_property()` after ",
                "`merge_replicates` and before `impute_groups` on the input ",
                "object list to calculate Exp_ratio for each feature in each ",
                "group.")
            }
        }
    }

    message("Max number of breakpoints: ", maxK)
    message("Min mean expression: ", meanCut)
    message("Min number of samples in each segment: ", minNumInSeg)

    res_list <- list()
    for (i in group) {
        message("Running Trendy for group: ", i)
        se_obj_imp <- se_obj_imp_list[[i]]
        res <- .run_Trendy_one_group(se_obj_imp,
            minExp = minExp, feature = feature,
            maxK = maxK, meanCut = meanCut,
            minNumInSeg = minNumInSeg, NCores = NCores, ...
        )
        res_list[[i]] <- res
    }

    return(res_list)
}


#' Segmented regression analysis (for one group)
#'
#' @description Run segmented regression analysis with Trendy on imputed data
#' for one group of samples in a SummarizedExperiment object. For use in the
#' `run_Trendy()` function.
#'
#' If less than 2 * minNumInSeg time points are available, Trendy analysis will
#' not be performed and NULL will be returned.
#'
#' @param se_obj_imp A SummarizedExperiment object with no missing values
#' @param minExp Minimum expression ratio for a feature to be included in
#' the analysis (default is 0.5)
#' @param feature A character vector of feature names to be included in the
#' analysis. If NULL, all features with expression ratio (Exp_ratio) >= minExp
#' will be
#' used, based on results of `calc_feature_property()` before imputation.
#' (default is NULL)
#' @param maxK Parameter of `Trendy::trendy()`, maximum number of breakpoints
#' allowed in the segmented regression model (default is 1)
#' @param meanCut Parameter of `Trendy::trendy()`, minimum mean expression
#' required for a feature to be included in the analysis (default is 0)
#' @param minNumInSeg Parameter of `Trendy::trendy()`, minimum number of
#' samples required in each segment (default is 3)
#' @param NCores Number of cores to use for parallel processing (default is 1)
#' @param ... Additional arguments to be passed to the `Trendy::trendy()`
#'
#' @import SummarizedExperiment
#' @import magrittr
#'
#' @returns Trendy analysis result, including the fitted model parameters
#' and statistics for each feature.
#' @keywords internal
.run_Trendy_one_group <- function(se_obj_imp,
    minExp = 0.5, feature = NULL, maxK = 1,
    meanCut = 0, minNumInSeg = 3, NCores = 1, ...) {

    n_time_points <- length(unique(colData(se_obj_imp)$Time))
    if (n_time_points < ((maxK + 1) * minNumInSeg)) {
        group <- se_obj_imp$Group %>% unique()
        message("Trendy analysis is not performed for group ", group,
            ": number of time points (", n_time_points,
            ") less than required ((maxK + 1) * minNumInSeg = ",
            ((maxK + 1) * minNumInSeg))
        return(NULL)
    }

    if (is.null(feature)) {
        feature <- rowData(se_obj_imp) %>%
            as.data.frame() %>%
            dplyr::filter(Exp_ratio >= minExp) %>%
            dplyr::pull(Feature)
        message("Feature not specified. Using ", length(feature),
            " features expressed in >=", minExp * 100, "% time points.")
    } else {
        feature <- intersect(feature, row.names(se_obj_imp))
        if (length(feature) == 0) {
            stop("None of the specified features are found in the input ",
            "object.")
        }
        message("Using ", length(feature),
            " specified features present in the data.")
    }

    input <- se_obj_imp[feature, ]
    M_imp <- assays(input)[[1]]
    time.vector <- colData(input)$Time
    trendy_out <- Trendy::trendy(
        Data = M_imp, tVectIn = time.vector, maxK = maxK,
        meanCut = meanCut, minNumInSeg = minNumInSeg, NCores = NCores, ...
    )
    res <- Trendy::results(trendy_out)

    return(res)
}
