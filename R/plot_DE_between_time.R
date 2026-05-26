#' DE number between time points
#' @description Plot the number of differentially expressed features between 
#' time points within each group with heatmaps
#'
#' @param se_obj A SummarizedExperiment object containing the data and 
#' metadata, with a "Time" column in the colData indicating the time points of 
#' the samples
#' @param de_list Output of `DE_between_time()`, a nested list of DE results 
#' for each group and time point comparison
#' @param value Whether to display the actual number of DE features in each 
#' heatmap cell (default is TRUE)
#' @param fontsize Font size for the heatmap displaying the number of DE 
#' features between time points (default is 8)
#' @param nrow Number of rows for arranging the heatmaps (default is 1)
#' @param heatmap_width Width of each heatmap (default is 4)
#' @param heatmap_unit Unit for the heatmap width (default is "cm")
#'
#' @import SummarizedExperiment
#' @import magrittr
#'
#' @returns One heatmap per group showing the number of DE features between 
#' time points for each group, with the same scale across groups for easy 
#' comparison. The heatmap cells are annotated with the actual number of DE 
#' features.
#' @export
#' @examples
#' data("example")
#' example_obj <- normalise_to_start(example_obj)
#'
#' DE_between_time_out <- DE_between_time(example_obj, assay = 1)
#' plot_DE_between_time(example_obj,
#'     de_list = DE_between_time_out$de_list,
#'     fontsize = 8, value = TRUE, nrow = 1, heatmap_width = 3)
plot_DE_between_time <- function(
    se_obj, de_list,
    value = TRUE,
    fontsize = 8,
    nrow = 1,
    heatmap_width = 4,
    heatmap_unit = "cm"
) {
    # heatmap of DE numbers
    de_num_list <- list()
    de_num_max <- 0
    time_series <- sort(unique(colData(se_obj)$Time)) %>% as.character() 
    # groups may have different time points

    for (i in names(de_list)) {
        de_num <- matrix(data = NA, nrow = length(time_series),
            ncol = length(time_series) - 1)

        colnames(de_num) <- time_series[-1]
        rownames(de_num) <- time_series

        for (d1 in time_series) {
            for (d2 in time_series[-1]) {
                if (as.numeric(d1) < as.numeric(d2)) {
                    label <- paste0("t", d2, "-t", d1)
                    if (!(label %in% names(de_list[[i]]))) {
                        next
                    }
                    de_num[d1, d2] <- de_list[[i]][[label]] %>%
                        dim() %>%
                        magrittr::extract2(1)
                }
            }
        }

        de_num_list[[i]] <- de_num
        de_num_max <- max(de_num_max, max(de_num, na.rm = TRUE))
    }

    col_fun <- circlize::colorRamp2(
        c(0, de_num_max / 2, de_num_max),
        c("#4575B4", "#FFFFBF", "#D73027")
    )

    # plot heatmap with same color scale
    p_list <- list()
    for (i in names(de_list)) {
        de_num <- de_num_list[[i]]

        legend <- FALSE

        if (value) {
            p <- ComplexHeatmap::Heatmap(de_num,
                name = "DE number",
                cluster_rows = FALSE, cluster_columns = FALSE,
                show_row_names = TRUE, show_column_names = TRUE,
                show_heatmap_legend = legend,
                column_names_rot = 0,
                column_names_centered = TRUE,
                cell_fun = function(j, i, x, y, width, height, fill) {
                    if (!is.na(de_num[i, j])) {
                        grid::grid.text(sprintf("%.0f", de_num[i, j]), x, y,
                            gp = grid::gpar(fontsize = fontsize, 
                                color = "black")
                        )
                    }
                },
                heatmap_width = unit(heatmap_width, heatmap_unit),
                col = col_fun,
                na_col = "white"
            )
        } else {
            p <- ComplexHeatmap::Heatmap(de_num,
                name = "DE number",
                cluster_rows = FALSE, cluster_columns = FALSE,
                show_row_names = TRUE, show_column_names = TRUE,
                show_heatmap_legend = legend,
                column_names_rot = 0,
                column_names_centered = TRUE,
                heatmap_width = unit(heatmap_width, heatmap_unit),
                col = col_fun,
                na_col = "white"
            )
        }

        p_list[[i]] <- grid::grid.grabExpr(ComplexHeatmap::draw(p)) %>% 
            ggplotify::as.ggplot()
    }

    # heatmap legend as a separate ggplot object
    lgd <- ComplexHeatmap::Legend(
        col_fun = col_fun,
        title = "DE number"
    )

    legend_grob <- grid::grid.grabExpr(
        grid::grid.draw(ComplexHeatmap::packLegend(lgd))
    ) %>%
        ggplotify::as.ggplot()

    # arrange heatmaps
    ncol <- ceiling(length(p_list) / nrow)

    heatmap_grid <- ggpubr::ggarrange(
        plotlist = p_list,
        labels = names(p_list),
        nrow = nrow,
        ncol = ncol,
        hjust = 0, vjust = 0.5,
        font.label = list(size = fontsize + 2)
    )

    # combine heatmaps and legend, to ensure legend is on the side
    final_plot <- ggpubr::ggarrange(
        heatmap_grid,
        legend_grob,
        ncol = 2,
        widths = c(length(p_list), 0.8)
    )

    print(final_plot %>%
        ggpubr::annotate_figure(top = ggpubr::text_grob(paste0(
            "Number of DE features between time points\n" 
        ), face = "bold", size = fontsize + 4))
    )
}
