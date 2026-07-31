#' Plot correlation matrix
#'
#' @description Plot correlation matrix between samples as a heatmap
#'
#' @param se_obj A SummarizedExperiment object created by `create_input()`
#' @param assay The assay to use in the SummarizedExperiment object: a numeric
#'   index or character name, e.g. 1 or "orig" for original data, 2 or "norm"
#'   for time 0 normalised data. (Default: 1.)
#' @param use Parameter of `stats::cor()` (default is "pairwise.complete.obs")
#' @param method Parameter of `stats::cor()` (default is "spearman")
#' @param label_group Whether to label Group (default is TRUE)
#' @param label_time Whether to label Time (default is TRUE)
#' @param label_rep Whether to label Replicate (default is TRUE)
#' @param label_batch Whether to label Batch (default is TRUE)
#' @param show_rownames Whether to show row names in the heatmap (default is
#' FALSE)
#' @param show_colnames Whether to show column names in the heatmap (default is
#' FALSE)
#' @param fontsize Font size for the plot (default is 8)
#' @param cellwidth Cell width for the heatmap (default is 1)
#' @param cellheight Cell height for the heatmap (default is 1)
#' @param title Title of the heatmap. If NULL, auto-generated as
#'   "Sample correlation (<Method>)" with the method capitalised
#'   (e.g. "Sample correlation (Spearman)") (default: NULL)
#' @param ... Additional arguments to be passed to `ComplexHeatmap::Heatmap()`
#'
#' @import SummarizedExperiment
#'
#' @returns A heatmap showing the correlation between samples.
#' @export
#' @examples
#' data(example_obj)
#' plot_cor_matrix(example_obj)
plot_cor_matrix <- function(
    se_obj,
    assay = 1,
    use = "pairwise.complete.obs",
    method = c("spearman", "pearson", "kendall"),
    label_group = TRUE, label_time = TRUE,
    label_rep = TRUE, label_batch = TRUE,
    show_rownames = FALSE, show_colnames = FALSE,
    fontsize = 8, cellwidth = 1, cellheight = 1,
    title = NULL,
    ...
) {
    .check_se(se_obj)
    assay <- .match_assay(assay, se_obj)
    .check_logical(label_group, "label_group")
    .check_logical(label_time, "label_time")
    .check_logical(label_rep, "label_rep")
    .check_logical(label_batch, "label_batch")
    .check_logical(show_rownames, "show_rownames")
    .check_logical(show_colnames, "show_colnames")
    .check_positive(fontsize, "fontsize")
    .check_positive(cellwidth, "cellwidth")
    .check_positive(cellheight, "cellheight")
    .check_character(use, "use")
    if (is.null(method) || length(method) != 1) {
        method <- "spearman"
    }
    method <- match.arg(method, c("spearman", "pearson", "kendall"))
    if (is.null(title)) {
        method_display <- paste0(
            toupper(substr(method, 1, 1)),
            substr(method, 2, nchar(method))
        )
        title <- paste0("Sample correlation (", method_display, ")")
    }
    .check_character(title, "title")
    cor_table <- stats::cor(assay(se_obj, assay),
        use = use,
        method = method
    )

    ann_colors <- list(
        Group = get_custom_palette(unique(colData(se_obj)$Group)),
        Time = scales::pal_viridis()(length(unique(colData(se_obj)$Time))) |>
            stats::setNames(nm = colData(se_obj)$Time |> unique() |> sort()),
        Replicate =
            ggsci::pal_iterm()(length(unique(colData(se_obj)$Replicate))) |>
            stats::setNames(nm =
                colData(se_obj)$Replicate |> unique() |> sort()),
        Batch = ggsci::pal_simpsons()(length(unique(colData(se_obj)$Batch))) |>
            stats::setNames(nm = colData(se_obj)$Batch |> unique() |> sort())
    )

    if (label_group) {
        ann <- c("Group")
    } else {
        ann <- c()
    }
    if (label_time) {
        ann <- c(ann, "Time")
    }
    if (label_rep) {
        ann <- c(ann, "Replicate")
    }
    if (label_batch) {
        ann <- c(ann, "Batch")
    }

    ht <- ComplexHeatmap::Heatmap(cor_table,
        name = "Correlation",
        clustering_distance_rows = stats::as.dist(1 - cor_table),
        clustering_distance_columns = stats::as.dist(1 - cor_table),
        col = grDevices::colorRampPalette(c("#3C5488FF", "white",
            "#E64B35FF"))(100),
        show_row_names = show_rownames,
        show_column_names = show_colnames,
        column_title = title,
        row_title = NULL,
        column_title_gp = grid::gpar(fontsize = fontsize + 2,
            fontface = "bold"),
        top_annotation = if (length(ann) > 0) {
            ComplexHeatmap::HeatmapAnnotation(
                df = colData(se_obj) |> as.data.frame() |>
                    dplyr::mutate(Time = factor(Time,
                        levels = as.character(colData(se_obj)$Time |>
                            unique() |> sort())
                    )) |>
                    dplyr::select(dplyr::all_of(ann)),
                col = ann_colors[ann],
                annotation_name_side = "left",
                annotation_legend_param = list(
                    title_gp = grid::gpar(fontsize = fontsize),
                    labels_gp = grid::gpar(fontsize = fontsize)
                ),
                annotation_name_gp = grid::gpar(fontsize = fontsize),
                border = TRUE
            )
        } else {
            NULL
        },
        width = unit(ncol(cor_table) * cellwidth / 1.5, "mm"),
        height = unit(nrow(cor_table) * cellheight / 1.5, "mm"),
        border = TRUE,
        heatmap_legend_param = list(
            title = "r", title_gp = grid::gpar(fontsize = fontsize),
            labels_gp = grid::gpar(fontsize = fontsize)
        ),
        ...
    )
    return(ht)
}
