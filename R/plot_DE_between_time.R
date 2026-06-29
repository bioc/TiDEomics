#' DE number between time points
#' @description Plot the number of differentially expressed features between
#' time points within each group with heatmaps
#'
#' DE features were pre-filtered with `DE_between_time()`, but by specifying
#' `adjP_thres` and `logFC_thres`, users can re-filter the DE features for
#' plotting.
#'
#' @param DE_between_time_out Output of `DE_between_time()`, a list with
#' `all_list`, `de_list`, `time_series` elements.
#' @param value Whether to display the actual number of DE features in each
#' heatmap cell (default is TRUE)
#' @param fontsize Font size for the heatmap displaying the number of DE
#' features between time points (default is 8)
#' @param nrow Number of rows for arranging the heatmaps (default is 1)
#' @param heatmap_width Width of each heatmap (default is 4)
#' @param heatmap_unit Unit for the heatmap width (default is "cm")
#' @param adjP_thres (Optional) Threshold for adjusted p-value to consider a
#' feature as differentially expressed, for re-filtering
#' the DE features. (default is NULL, no re-filtering)
#' @param logFC_thres (Optional) Threshold for log2 fold change to consider a
#' feature as differentially expressed, for re-filtering
#' the DE features. (default is NULL, no re-filtering)
#'
#' @import SummarizedExperiment
#'
#' @returns One heatmap per group showing the number of DE features between
#' time points for each group, with the same scale across groups for easy
#' comparison. The heatmap cells are optionally annotated with the number of DE
#' features.
#' @export
#' @examples
#' data(example_obj)
#' example_obj <- normalise_to_start(example_obj)
#'
#' DE_between_time_out <- DE_between_time(example_obj, assay = 1)
#' plot_DE_between_time(DE_between_time_out,
#'     fontsize = 8, value = TRUE, nrow = 1, heatmap_width = 3)
plot_DE_between_time <- function(
    DE_between_time_out,
    value = TRUE,
    fontsize = 8,
    nrow = 1,
    heatmap_width = 4,
    heatmap_unit = "cm",
    adjP_thres = NULL, logFC_thres = NULL
) {
    .check_character(heatmap_unit, "heatmap_unit")
    .check_positive(fontsize, "fontsize")
    .check_positive_int(nrow, "nrow")
    .check_positive(heatmap_width, "heatmap_width")
    .check_logical(value, "value")

    if (!is.list(DE_between_time_out) ||
        !all(c("all_list", "de_list") %in% names(DE_between_time_out))) {
        stop("'DE_between_time_out' must be the output of DE_between_time().")
    }

    all_list <- DE_between_time_out$all_list

    if (!is.null(adjP_thres) && !is.null(logFC_thres)) {
        .check_pval(adjP_thres, "adjP_thres")
        .check_nonneg(logFC_thres, "logFC_thres")
        message("Re-filtering DE features with adjP_thres = ", adjP_thres,
            " and logFC_thres = ", logFC_thres, ".")

        # Extract DE features with new thresholds
        de_list <- list()

        for (i in names(all_list)) {
            de_list[[i]] <- list()

            for (tb_name in names(all_list[[i]])) {
                tb <- all_list[[i]][[tb_name]]
                tb_de <- tb |>
                    dplyr::select(Feature, Comparison, Group, Cond1, Cond2, 
                        logFC, adj.P.Val) |>
                    dplyr::filter(adj.P.Val < adjP_thres) |>
                    dplyr::filter(logFC > logFC_thres | logFC < -logFC_thres)
                de_list[[i]][[tb_name]] <- tb_de
            }
        }
    } else {
        de_list <- DE_between_time_out$de_list
    }

    # heatmap of DE numbers
    de_num_list <- list()
    de_num_max <- 0
    time_series <- DE_between_time_out$time_series
    if (length(time_series) < 2) {
        stop("Need at least 2 time points for DE-between-time plot, ",
            "found ", length(time_series), ".")
    }
    # groups may have different time points

    # Last time point has no later points to compare to - exclude it
    time_rows <- time_series[-length(time_series)]

    for (i in names(de_list)) {
        de_num <- matrix(data = NA, nrow = length(time_rows),
            ncol = length(time_series) - 1)

        colnames(de_num) <- time_series[-1]
        rownames(de_num) <- time_rows

        for (d1 in time_rows) {
            for (d2 in time_series[-1]) {
                if (as.numeric(d1) < as.numeric(d2)) {
                    label <- paste0("t", d2, "-t", d1)
                    if (!(label %in% names(de_list[[i]]))) {
                        next
                    }
                    de_num[d1, d2] <- nrow(de_list[[i]][[label]])
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

        p_list[[i]] <- grid::grid.grabExpr(ComplexHeatmap::draw(p)) |>
            ggplotify::as.ggplot()
    }

    # heatmap legend as a separate ggplot object
    lgd <- ComplexHeatmap::Legend(
        col_fun = col_fun,
        title = "DE number"
    )

    legend_grob <- grid::grid.grabExpr(
        grid::grid.draw(ComplexHeatmap::packLegend(lgd))
    ) |>
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

    final_plot <- final_plot |>
        ggpubr::annotate_figure(top = ggpubr::text_grob(paste0(
            "Number of DE features between time points\n"
        ), face = "bold", size = fontsize + 4))
    return(final_plot)
}
