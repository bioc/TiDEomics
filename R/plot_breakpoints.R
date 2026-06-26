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
#' @import ggplot2
#'
#' @returns A plot showing the distribution of breakpoints over time for 
#' each specified group
#' @export
#' @examples
#' data("example_res_list")
#' plot_breakpoints(example_res_list)
plot_breakpoints <- function(res_list, group = NULL, fontsize = 8, ...) {
    .check_positive(fontsize, "fontsize")
    .check_list(res_list, "res_list", "run_Trendy")
    if (!is.null(group)) .check_character(group, "group")
    if (is.null(group)) {
        group <- names(res_list)
    } else if (!all(group %in% names(res_list))) {
        stop("At least one of the specified groups is not found in the ", 
        "input object list.")
    }

    bp_list <- list()
    for (i in group) {
        res <- res_list[[i]]
        if (is.null(res)) {
            message("Skipping group '", i,
                "': Trendy analysis was not performed.")
            next
        }
        res.top <- Trendy::topTrendy(res, ...)

        # Breakpoint distribution over the time course
        res.bp <- Trendy::breakpointDist(res.top)

        res.bp.df <- as.data.frame(res.bp)
        res.bp.df <- tibble::rownames_to_column(res.bp.df, "Day")
        colnames(res.bp.df)[ncol(res.bp.df)] <- "Count"
        res.bp.df <- res.bp.df |>
            dplyr::mutate(
                Day = as.numeric(Day) |> factor(),
                Group = i
            )

        bp_list[[i]] <- res.bp.df
    }

    res.bp.df.all <- do.call(rbind, bp_list) |>
        as.data.frame() |>
        dplyr::mutate(Group = factor(Group, levels = group))

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
    return(invisible(p))
}
