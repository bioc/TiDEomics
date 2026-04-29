#' UMAP
#'
#' @description Plot UMAP of samples, using features without missing values
#'
#' The UMAP plots are labelled by Group, Time, Replicate (if more than 1),
#' Batch (if more than 1), and number of identified features (non-missing
#' values).
#'
#' @param se_obj A SummarizedExperiment object created by `create_input()`
#' @param seed Random seed for UMAP (default is 1234)
#' @param plot Whether to plot the figures (default is TRUE)
#' @param plot_ID Whether to include a UMAP plot coloured by number of
#' identified features (default is FALSE)
#' @param circle Logical, whether to draw circles (ellipses) around samples of
#' each group (default is FALSE)
#' @param xlim_min Minimum x-axis limit when drawing ellipses (default 1.5*min
#' UMAP x)
#' @param xlim_max Maximum x-axis limit when drawing ellipses (default 1.5*max
#' UMAP x)
#' @param ylim_min Minimum y-axis limit when drawing ellipses (default 1.5*min
#' UMAP y)
#' @param ylim_max Maximum y-axis limit when drawing ellipses (default 1.5*max
#' UMAP y)
#' @param umap_neighbors UMAP n_neighbors parameter (default is selected by
#' `umap_n_neighbors()` function based on the number of samples)
#' @param fontsize Font size for the plot (default is 8)
#' @param assay Assay index to use, where 1 is the original data and 2 is
#' normalised to time 0 (if available) (default is 1)
#' @param ... Additional arguments passed to `ggforce::geom_mark_ellipse()`
#' for customizing the ellipses
#'
#' @import SummarizedExperiment
#' @import magrittr
#' @import ggplot2
#' @importFrom dplyr mutate
#' @importFrom stats complete.cases
#'
#' @returns A UMAP plot showing the distribution of samples. And a data
#' frame containing UMAP coordinates and sample annotations for custom plotting.
#' @export
#' @examples
#' data("example")
#' umap_layout <- plot_umap(example_obj)
plot_umap <- function(
    se_obj, seed = 1234, plot = TRUE, plot_ID = FALSE,
    circle = FALSE,
    xlim_min = NULL, xlim_max = NULL,
    ylim_min = NULL, ylim_max = NULL,
    umap_neighbors = NULL,
    fontsize = 8,
    assay = 1,
    ...
) {
    M_0 <- se_obj@assays@data[[assay]]
    if (sum(is.na(M_0)) > 0) {
        message("Input data contains missing values. Only complete rows ",
        "will be used for UMAP.")
    }
    M <- M_0[complete.cases(M_0), ]
    sp_info <- as.data.frame(colData(se_obj))

    n_neighbors <- ifelse(is.null(umap_neighbors),
        umap_n_neighbors(ncol(M)),
        umap_neighbors)

    message("Using n_neighbors = ", n_neighbors)

    umap <- umap::umap(M %>% t(), random_state = seed,
        n_neighbors = n_neighbors)
    umap_layout <- umap$layout %>% as.data.frame()
    umap_layout$Sample <- rownames(umap_layout)
    umap_layout <- merge(umap_layout, sp_info, by = "Sample")
    umap_title <- paste0("UMAP (", nrow(M), " features without missing values)")

    # number of identified proteins
    id_count <- apply(M_0, 2, function(x) sum(!is.na(x))) %>% as.data.frame()
    colnames(id_count) <- "ID"
    id_count$Sample <- rownames(id_count)
    umap_layout <- merge(umap_layout, id_count, by = "Sample")

    if (!plot) {
        return(umap_layout)
    }

    umap_list <- list()
    umap_list[["Group"]] <- ggplot(umap_layout,
        aes(x = V1, y = V2, color = Group)) +
        geom_point(size = 2) +
        scale_color_manual(values =
            get_custom_palette(unique(colData(se_obj)$Group))) +
        theme_custom(base_size = fontsize)
    umap_list[["Time"]] <- ggplot(umap_layout,
        aes(x = V1, y = V2, color = Time)) +
        geom_point(size = 2) +
        scale_color_viridis_c() +
        theme_custom(base_size = fontsize)
    if (plot_ID) {
        umap_list[["ID"]] <- ggplot(umap_layout,
            aes(x = V1, y = V2, color = ID)) +
            geom_point(size = 2) +
            scale_color_viridis_c() +
            theme_custom(base_size = fontsize)
    }

    if (length(unique(umap_layout$Replicate)) < 2) {
        umap_list[["Replicate"]] <- NULL
    } else {
        umap_list[["Replicate"]] <- ggplot(umap_layout,
            aes(x = V1, y = V2, color = Replicate)) +
            geom_point(size = 2) +
            scale_color_hue() +
            theme_custom(base_size = fontsize)
    }
    if (length(unique(umap_layout$Batch)) < 2) {
        umap_list[["Batch"]] <- NULL
    } else {
        umap_list[["Batch"]] <- ggplot(umap_layout,
            aes(x = V1, y = V2, color = Batch)) +
            geom_point(size = 2) +
            scale_color_hue() +
            theme_custom(base_size = fontsize)
    }

    print(ggpubr::ggarrange(
        plotlist = umap_list,
        nrow = 2, ncol = ceiling(length(umap_list) / 2),
        common.legend = FALSE
    ) %>%
        ggpubr::annotate_figure(top = ggpubr::text_grob(
            paste0(umap_title, "\n"),
            face = "bold", size = fontsize + 2)))

    p1 <- umap_layout %>%
        mutate(Time = Time %>% factor()) %>%
        ggplot(aes(x = V1, y = V2, fill = Group)) +
        geom_point(aes(size = Time), alpha = 0.8, shape = 21) +
        scale_fill_manual(values =
            get_custom_palette(unique(colData(se_obj)$Group))) +
        guides(fill = guide_legend(override.aes = list(size = 4))) +
        theme_custom(base_size = fontsize) +
        ggtitle(umap_title)
    print(p1)

    if (is.null(xlim_min)) {
        xlim_min <- 1.5 * min(umap_layout$V1)
    }
    if (is.null(xlim_max)) {
        xlim_max <- 1.5 * max(umap_layout$V1)
    }
    if (is.null(ylim_min)) {
        ylim_min <- 1.5 * min(umap_layout$V2)
    }
    if (is.null(ylim_max)) {
        ylim_max <- 1.5 * max(umap_layout$V2)
    }

    if (circle) {
        p2 <- p1 + ggforce::geom_mark_ellipse(aes(fill = Group, label = Group),
            con.cap = 0, alpha = 0.1, ...
        ) +
            xlim(xlim_min, xlim_max) + ylim(ylim_min, ylim_max)
        print(p2)
    }

    return(umap_layout)
}

#' UMAP neighbors number
#'
#' @description Select UMAP n_neighbors parameter based on the number of
#' samples
#' For sample number > 5, use sample number / 5, with a min of 5 and max of 15.
#' For sample number <= 5, use sample number - 1.
#'
#' @param sample_n Number of samples
#'
#' @returns UMAP n_neighbors parameter
#' @keywords internal
umap_n_neighbors <- function (sample_n) {
    if (sample_n <= 5) {
        return(sample_n - 1)
    } else {
        return(floor(min(15, max(sample_n / 5, 5))))
    }
}
