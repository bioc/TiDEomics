#' Plot UMAP by group
#'
#' @description Plot UMAP for each group separately. Accepts either a
#' SummarizedExperiment object or a named list of them (from
#' `split_groups()`).
#'
#' @param se_obj A SummarizedExperiment object, or a named list of them
#'   (e.g. from `split_groups()`).
#' @param seed Random seed for UMAP (default: 1234).
#' @param nrow Number of rows for arranging the plots (default: 1).
#' @param umap_neighbors UMAP n_neighbors parameter (default: auto-selected
#'   based on sample count).
#' @param fontsize Base font size (default: 8).
#' @param assay The assay to use in the SummarizedExperiment object: a numeric
#'   index or character name, e.g. 1 or "orig" for original data, 2 or "norm"
#'   for time 0 normalised data. (Default: 1.)
#' @param legend_pos Legend position (default: "right").
#'
#' @import ggplot2
#' @import SummarizedExperiment
#'
#' @returns A series of UMAP plots showing the distribution of samples in
#' each group, coloured by Time.
#' @export
#' @examples
#' data(example_obj)
#' plot_umap_by_group(example_obj)
#'
#' # Also accepts a list from split_groups()
#' example_obj_list <- split_groups(example_obj)
#' plot_umap_by_group(example_obj_list)
plot_umap_by_group <- function(se_obj, seed = 1234, nrow = 1,
    umap_neighbors = NULL,
    fontsize = 8, assay = 1, legend_pos = "right") {
    if (is.list(se_obj) && !methods::is(se_obj, "SummarizedExperiment")) {
        .check_se_list(se_obj)
        se_list <- se_obj
    } else {
        .check_se(se_obj)
        se_list <- split_groups(se_obj)
    }
    .check_positive_int(seed, "seed")
    .check_positive_int(nrow, "nrow")
    .check_character(legend_pos, "legend_pos")
    .check_positive(fontsize, "fontsize")
    if (!is.null(umap_neighbors))
        .check_positive_int(umap_neighbors, "umap_neighbors")
    assay <- .match_assay(assay, se_list[[1]])

    umap_list <- list()
    for (group in names(se_list)) {
        umap_layout <- plot_umap(se_list[[group]],
            seed = seed,
            assay = assay,
            umap_neighbors = umap_neighbors
        )
        umap_list[[group]] <- ggplot(umap_layout$umap_layout,
            aes(x = V1, y = V2, color = Time)) +
            geom_point(size = 2) +
            scale_color_viridis_c() +
            theme_custom(base_size = fontsize)
    }

    ncol <- ceiling(length(umap_list) / nrow)

    p <- ggpubr::ggarrange(
        plotlist = umap_list, labels = names(umap_list),
        font.label = list(size = fontsize + 2),
        hjust = 0, vjust = 0.5,
        nrow = nrow, ncol = ncol,
        common.legend = TRUE,
        legend = legend_pos
    ) |>
        ggpubr::annotate_figure(top =
        ggpubr::text_grob("UMAP - by group (features without missing values)\n",
            face = "bold", size = fontsize + 4))
    return(p)
}
