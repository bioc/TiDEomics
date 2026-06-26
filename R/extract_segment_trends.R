#' Extract feature trends
#' @description Extract features with different segment trends from trendy
#' segmented regression results summary
#'
#' @param trendy.summary A data frame, output of `summarise_Trendy()`,
#' containing the segmented regression results for all features in all groups.
#'
#'
#' @returns A nested list of features grouped by their segment trend patterns
#' within each group.
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
#' data("example_res_list")
#' trendy_summary <- summarise_Trendy(example_res_list)
#' trendy_list <- extract_segment_trends(trendy_summary)
extract_segment_trends <- function(trendy.summary) {
    .check_df_trendy_summary(trendy.summary, "trendy.summary")
    trendy_list <- list()

    for (i in unique(trendy.summary$Group)) {
        trendy_list[[i]] <- list()

        trend_comb <- unique(trendy.summary |> dplyr::filter(Group == i) |>
            dplyr::pull(Pattern))

        for (j in trend_comb) {
            trendy_list[[i]][[j]] <- trendy.summary |>
                dplyr::filter(Group == i) |>
                dplyr::filter(Pattern == j) |>
                dplyr::pull(Feature)
        }
    }

    return(trendy_list)
}
