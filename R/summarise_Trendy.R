#' Summarise Trendy results
#' @description Summarise Trendy results into data frame
#'
#' @param res_list A list of Trendy analysis results, output of `run_Trendy()`
#' @param ... Additional arguments to be passed to the `Trendy::topTrendy()`

#' @import magrittr
#' @importFrom dplyr select mutate
#'
#' @returns A data frame of summary results of Trendy, including breakpoints,
#' segment slopes and p-values for each fitted feature in each group. A
#' "Pattern" column is added to show the combination of segment trends (up,
#' down, stable) for each feature.
#' @export
#' @examples
#' data("example")
#' example_obj <- normalise_to_start(example_obj)
#' example_obj_list <- split_groups(example_obj)
#' example_obj_merged_list <- merge_replicates(example_obj_list)
#' example_obj_merged_list <-
#'     calc_feature_property(example_obj_merged_list, threshold = 0)
#' # no missing value in the example dataset, so imputation is not necessary
#' example_obj_merged_imp_list <- impute_groups(example_obj_merged_list)
#'
#' # "untreated" group has only 3 time points, so Trendy analysis will not be
#' # performed for this group
#' # example_res_list <- run_Trendy(example_obj_merged_imp_list, maxK = 1,
#' #     minNumInSeg = 2, meanCut = 0)
#' data("example_res_list")
#'
#' plot_segments(example_obj_merged_imp_list, example_res_list,
#'     feature = c("Mctp1"))
#' plot_breakpoints(example_res_list)
#' trendy_summary <- summarise_Trendy(example_res_list)
#' trendy_list <- extract_segment_trends(trendy_summary)
summarise_Trendy <- function(res_list, ...) {
    res_summary <- list()
    for (i in names(res_list)) {
        trendy.summary <- .summarise_Trendy_one_group(res_list[[i]])
        trendy.summary$Group <- i
        trendy.summary <- trendy.summary %>%
            dplyr::select(.data$Group, dplyr::everything())
        res_summary[[i]] <- trendy.summary
    }
    res_summary_df <- do.call(rbind, res_summary) %>% as.data.frame()

    # convert segment trends from 1, -1, 0 to up, down, stable
    df <- res_summary_df

    seg_cols <- grep("^Segment[0-9]+\\.Trend$", colnames(df), value = TRUE)
    df[seg_cols] <- df[seg_cols] %>% lapply(function(col) {
        plyr::mapvalues(col, c(1, 0, -1), c("up", "stable", "down"))
    })

    # all combinations of segment trends
    patterns <- apply(df[seg_cols], 1, function(row) {
        prs <- which(!is.na(row))
        paste(row[prs], collapse = "_")
    })
    df_patterns <- df %>%
        mutate(Pattern = patterns)

    return(df_patterns)
}

#' Summarise Trendy results (one group)
#' @description Summarise Trendy results into data frame (one group)
#'
#' @param res Result of the Trendy analysis, output of `run_Trendy()`
#' @param ... Additional arguments to be passed to the `Trendy::topTrendy`
#'
#' @importFrom dplyr rename
#' @import magrittr
#'
#' @returns A data frame of summary results of Trendy, including breakpoints,
#' segment slopes and p-values for each fitted feature.
#' @keywords internal
.summarise_Trendy_one_group <- function(res, ...) {
    res.top <- Trendy::topTrendy(res, ...)
    trendy.summary <- Trendy::formatResults(res.top)

    # bug in colname when maxK = 1
    if ("topTrendyData.Breakpoints.featureNames..." %in%
        colnames(trendy.summary)) {
        trendy.summary <- trendy.summary %>%
            dplyr::rename(Breakpoint =
                .data$`topTrendyData.Breakpoints.featureNames...`)
    }

    return(trendy.summary)
}
