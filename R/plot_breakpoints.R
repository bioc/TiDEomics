#' Plot breakpoint distribution
#' @description Plot breakpoint distribution among time points for each group
#'
#' @param res_list A list of the Trendy analysis results, output of 
#' `run_Trendy()`
#' @param group (Optional) A character vector of group names to be plotted. 
#' If NULL, all groups in the input list will be used (default is NULL)
#' @param fontsize (Optional) Font size for the plot (default is 8)
#' @param ... Additional arguments to be passed to `Trendy::topTrendy()`
#'
#' @import SummarizedExperiment
#' @import magrittr
#' @import ggplot2
#' @importFrom dplyr mutate rename
#'
#' @returns A plot showing the distribution of breakpoints over time for 
#' each specified group
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
plot_breakpoints <- function(res_list, group = NULL, fontsize = 8, ...) {
    if (is.null(group)) {
        group <- names(res_list)
    } else if (!all(group %in% names(res_list))) {
        stop("At least one of the specified groups is not found in the ", 
        "input object list.")
    }

    p_list <- list()
    for (i in group) {
        res <- res_list[[i]]
        res.top <- Trendy::topTrendy(res, ...)

        # Breakpoint distribution over the time course
        res.bp <- Trendy::breakpointDist(res.top)

        res.bp.df <- res.bp %>%
            as.data.frame() %>%
            tibble::rownames_to_column("Day") %>%
            dplyr::rename("Count" = ".") %>%
            mutate(
                Day = as.numeric(Day) %>% factor(),
                Group = i
            )

        if (i == group[1]) {
            res.bp.df.all <- res.bp.df
        } else {
            res.bp.df.all <- rbind(res.bp.df.all, res.bp.df)
        }
    }

    res.bp.df.all <- res.bp.df.all %>%
        as.data.frame() %>%
        mutate(Group = factor(Group, levels = group))

    p <- ggplot(res.bp.df.all, aes(x = Day, y = Count, color = Group)) +
        geom_point(position = position_dodge(width = 0.2)) +
        geom_line(aes(group = Group), position = position_dodge(width = 0.2)) +
        labs(
            title = "Breakpoint distribution",
            x = "Day",
            y = "Count of breakpoints"
        ) +
        theme_custom(base_size = fontsize) +
        scale_color_manual(values = get_custom_palette(group))

    print(p)
}
