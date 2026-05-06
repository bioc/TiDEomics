#' Plot UMAP by group (one object)
#'
#' @description Plot UMAP by group, input is a SummarizedExperiment object
#' with a "Group" column in the colData.
#'
#' @param se_obj A SummarizedExperiment object created by `create_input()`
#' @param seed Random seed for UMAP (default is 1234)
#' @param nrow Number of rows for arranging the UMAP plots (default is 1)
#' @param umap_neighbors UMAP n_neighbors parameter (default is selected by
#' `umap_n_neighbors()` function based on the number of samples)
#' @param fontsize Font size for the plot (default is 8)
#' @param assay Assay index to use, where 1 is the original data and 2 is
#' normalised to time 0 (if available) (default is 1)
#' @param legend_pos Legend position for the PCA plots (default is "right")
#'
#' @import ggplot2
#' @import SummarizedExperiment
#' @import magrittr
#'
#' @returns A series of UMAP plots showing the distribution of samples in
#' each group, coloured by Time.
#' @export
#' @examples
#' data("example")
#' plot_umap_by_group(example_obj)
plot_umap_by_group <- function(se_obj, seed = 1234, nrow = 1,
    umap_neighbors = NULL,
    fontsize = 8, assay = 1, legend_pos = "right") {
    umap_list <- list()
    for (group in unique(se_obj$Group)) {
        umap_layout <- plot_umap(se_obj[, se_obj$Group == group],
            plot = FALSE,
            seed = seed,
            assay = assay,
            umap_neighbors = umap_neighbors
        )
        umap_list[[group]] <- ggplot(umap_layout,
            aes(x = .data$V1, y = .data$V2, color = .data$Time)) +
            geom_point(size = 2) +
            scale_color_viridis_c() +
            theme_custom(base_size = fontsize)
    }

    ncol <- ceiling(length(umap_list) / nrow)

    print(ggpubr::ggarrange(
        plotlist = umap_list, labels = names(umap_list),
        font.label = list(size = fontsize + 2),
        hjust = 0, vjust = 0.5,
        nrow = nrow, ncol = ncol,
        common.legend = TRUE,
        legend = legend_pos
    ) %>%
        ggpubr::annotate_figure(top =
        ggpubr::text_grob("UMAP - by group (features without missing values)\n",
            face = "bold", size = fontsize + 4)))
}


#' Plot UMAP by group (list of objects)
#'
#' @description Plot UMAP by group, input is a list of SummarizedExperiment
#' objects, with each object corresponding to a group
#'
#' @param se_obj_list A list of SummarizedExperiment objects, such as output
#' of `split_groups()`
#' @param seed Random seed for UMAP (default is 1234)
#' @param nrow Number of rows for arranging the UMAP plots (default is 1)
#' @param umap_neighbors UMAP n_neighbors parameter (default is selected by
#' `umap_n_neighbors()` function based on the number of samples)
#' @param fontsize Font size for the plot (default is 8)
#' @param assay Assay index to use, where 1 is the original data and 2 is
#' normalised to time 0 (if available) (default is 1)
#' @param legend_pos Legend position for the PCA plots (default is "right")
#'
#' @import ggplot2
#' @import SummarizedExperiment
#' @import magrittr
#'
#' @returns A series of UMAP plots showing the distribution of samples in
#' each group, coloured by Time.
#' @export
#' @examples
#' data("example")
#' example_obj_list <- split_groups(example_obj)
#' plot_umap_by_group_list(example_obj_list)
plot_umap_by_group_list <- function(se_obj_list, seed = 1234, nrow = 1,
    umap_neighbors = NULL,
    fontsize = 8, assay = 1, legend_pos = "right") {
    umap_list <- list()
    for (group in names(se_obj_list)) {
        umap_layout <- plot_umap(se_obj_list[[group]],
            plot = FALSE,
            seed = seed,
            assay = assay,
            umap_neighbors = umap_neighbors
        )
        umap_list[[group]] <- ggplot(umap_layout,
            aes(x = .data$V1, y = .data$V2, color = .data$Time)) +
            geom_point(size = 2) +
            scale_color_viridis_c() +
            theme_custom(base_size = fontsize)
    }

    ncol <- ceiling(length(umap_list) / nrow)

    print(ggpubr::ggarrange(
        plotlist = umap_list, labels = names(umap_list),
        font.label = list(size = fontsize + 2),
        hjust = 0, vjust = 0.5,
        nrow = nrow, ncol = ncol,
        common.legend = TRUE,
        legend = legend_pos
    ) %>%
        ggpubr::annotate_figure(top =
        ggpubr::text_grob("UMAP - by group (features without missing values)\n",
            face = "bold", size = fontsize + 4)))
}
