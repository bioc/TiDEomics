#' Plot PCA by group
#'
#' @description Plot PCA for each group separately, with optional circles
#' around time points and arrows indicating trajectory over time. Accepts
#' either a SummarizedExperiment object or a list of them (from
#' `split_groups()`).
#'
#' @param se_obj A SummarizedExperiment object, or a named list of them
#'   (e.g. from `split_groups()`).
#' @param circle Logical, whether to draw ellipses around samples of each
#'   time point (default: TRUE).
#' @param arrow Logical, whether to draw arrows indicating the trajectory
#'   over time (default: TRUE).
#' @param pc1 Principal component for the x-axis (default: 1).
#' @param pc2 Principal component for the y-axis (default: 2).
#' @param nrow Number of rows for arranging the plots (default: 1).
#' @param fontsize Base font size (default: 8).
#' @param assay Assay index to use, where 1 is the original data and 2 is
#'   normalised to time 0 (if available) (default: 1).
#' @param legend_pos Legend position (default: "right").
#'
#' @import ggplot2
#' @import SummarizedExperiment
#'
#' @returns A series of PCA plots showing the distribution of samples in
#' each group, coloured by Time.
#' @export
#' @examples
#' data(example_obj)
#' plot_pca_by_group(example_obj, circle = TRUE, arrow = TRUE)
#'
#' # Also accepts a list from split_groups()
#' example_obj_list <- split_groups(example_obj)
#' plot_pca_by_group(example_obj_list)
plot_pca_by_group <- function(se_obj, circle = TRUE, arrow = TRUE,
    pc1 = 1, pc2 = 2, nrow = 1, fontsize = 8, assay = 1,
    legend_pos = "right") {
    if (is.list(se_obj) && !methods::is(se_obj, "SummarizedExperiment")) {
        .check_se_list(se_obj)
        se_list <- se_obj
    } else {
        .check_se(se_obj)
        se_list <- split_groups(se_obj)
    }
    .check_logical(circle, "circle")
    .check_logical(arrow, "arrow")
    .check_positive_int(pc1, "pc1")
    .check_positive_int(pc2, "pc2")
    .check_positive_int(nrow, "nrow")
    .check_character(legend_pos, "legend_pos")
    .check_positive(fontsize, "fontsize")
    assay <- .match_assay(assay, se_list[[1]])

    pca_list <- list()
    for (group in names(se_list)) {
        p1 <- plot_pca_arrows(se_list[[group]], circle = circle,
            arrow = arrow, pc1 = pc1, pc2 = pc2,
            fontsize = fontsize, assay = assay)
        pca_list[[group]] <- p1
    }

    ncol <- ceiling(length(pca_list) / nrow)

    p <- ggpubr::ggarrange(
        plotlist = pca_list, labels = names(pca_list),
        font.label = list(size = fontsize + 2),
        hjust = 0, vjust = 0.5,
        nrow = nrow, ncol = ncol,
        common.legend = TRUE,
        legend = legend_pos
    ) |>
        ggpubr::annotate_figure(top =
        ggpubr::text_grob(
            "PCA by group (features without missing values)\n",
            face = "bold", size = fontsize + 4))

    return(p)
}
#' Plot PCA with arrows
#'
#' @description Plot PCA with arrows indicating trajectory over time, for
#' a single group.
#'
#' @param se_obj A SummarizedExperiment object (single group).
#' @param circle Logical, whether to draw ellipses around samples of each
#'   time point (default: TRUE).
#' @param arrow Logical, whether to draw arrows indicating the trajectory
#'   over time (default: TRUE).
#' @param pc1 Principal component for the x-axis (default: 1).
#' @param pc2 Principal component for the y-axis (default: 2).
#' @param fontsize Base font size (default: 8).
#' @param assay Assay index to use, where 1 is the original data and 2 is
#'   normalised to time 0 (if available) (default: 1).
#'
#' @import ggplot2
#' @import SummarizedExperiment
#'
#' @returns A PCA plot coloured by Time, with optional ellipses and
#' trajectory arrows.
#' @export
#' @examples
#' data(example_obj)
#' plot_pca_arrows(example_obj[, example_obj$Group == "IFNbeta"])
plot_pca_arrows <- function(se_obj, circle = TRUE, arrow = TRUE,
    pc1 = 1, pc2 = 2, fontsize = 8, assay = 1) {
    .check_se(se_obj)
    .check_logical(circle, "circle")
    .check_logical(arrow, "arrow")
    .check_positive_int(pc1, "pc1")
    .check_positive_int(pc2, "pc2")
    .check_positive(fontsize, "fontsize")
    assay <- .match_assay(assay, se_obj)

    if (length(unique(se_obj$Group)) > 1) {
        warning("Input object contains multiple groups, circles / arrows ",
        "will be calculated across all samples. Please use ",
        "`plot_pca_by_group()` or subset the object to one group ",
        "before plotting.")
    }

    pca2 <- plot_pca(se_obj, assay = assay,
        plot_screeplot = FALSE, plot_loadings = FALSE,
        plot_morepc = FALSE)$pca
    pc <- as.data.frame(pca2$rotated)
    pc$Sample <- rownames(pc)
    pc <- merge(pc, as.data.frame(colData(se_obj)), by = "Sample")

    xlab <- paste0("PC", pc1, ": ",
        round(pca2$variance[[paste0("PC", pc1)]], 2), "%")
    ylab <- paste0("PC", pc2, ": ",
        round(pca2$variance[[paste0("PC", pc2)]], 2), "%")

    pc1_col <- paste0("PC", pc1); pc2_col <- paste0("PC", pc2)
    colnames(pc)[colnames(pc) == pc1_col] <- "PC1"
    colnames(pc)[colnames(pc) == pc2_col] <- "PC2"
    arrows_in <- pc[, c("PC1", "PC2", "Sample", "Time")] |>
        dplyr::arrange(Time)
    x_arrows <- stats::aggregate(PC1 ~ Time, data = arrows_in, FUN = mean)
    y_arrows <- stats::aggregate(PC2 ~ Time, data = arrows_in, FUN = mean)
    arrows <- merge(x_arrows, y_arrows, by = "Time") |>
        dplyr::rename(x_start = PC1, y_start = PC2)

    time_series <- sort(unique(x_arrows$Time))

    for (i in seq(1, length(time_series))) {
        if (i == length(time_series)) next
        arrows[which(arrows$Time == time_series[i]), "x_end"] <-
            arrows[which(arrows$Time == time_series[i + 1]), "x_start"]
        arrows[which(arrows$Time == time_series[i]), "y_end"] <-
            arrows[which(arrows$Time == time_series[i + 1]), "y_start"]
    }

    arrows <- arrows[!is.na(arrows$x_end), ]
    if (circle && arrow) {
        p1 <- pc |>
            dplyr::mutate(Time = factor(Time,
                levels = unique(pc$Time) |> sort())) |>
            ggplot(aes(x = PC1, y = PC2, color = Time)) +
            geom_segment(
                data = arrows, aes(
                    x = x_start, y = y_start,
                    xend = x_end, yend = y_end
                ),
                arrow = arrow(length = unit(0.25, "cm")),
                alpha = 0.75, color = "black"
            ) +
            geom_point(size = 3) +
            ggforce::geom_mark_ellipse(aes(fill = Time),
                con.cap = 0, expand = unit(2, "mm")
            ) +
            scale_color_viridis_d(option = "D") +
            scale_fill_viridis_d(option = "D") +
            theme_custom(base_size = fontsize) +
            xlim(range(pc$PC1) * 1.2) +
            ylim(range(pc$PC2) * 1.2) +
            xlab(xlab) +
            ylab(ylab)
    } else if (arrow) {
        p1 <- pc |>
            dplyr::mutate(Time = factor(Time,
                levels = unique(pc$Time) |> sort())) |>
            ggplot(aes(x = PC1, y = PC2, color = Time)) +
            geom_segment(
                data = arrows, aes(
                    x = x_start, y = y_start,
                    xend = x_end, yend = y_end
                ),
                arrow = arrow(length = unit(0.25, "cm")), alpha = 0.75,
                color = "black"
            ) +
            geom_point(size = 3) +
            scale_color_viridis_d(option = "D") +
            scale_fill_viridis_d(option = "D") +
            theme_custom(base_size = fontsize) +
            xlab(xlab) +
            ylab(ylab)
    } else if (circle) {
        p1 <- pc |>
            dplyr::mutate(Time = factor(Time,
                levels = unique(pc$Time) |> sort())) |>
            ggplot(aes(x = PC1, y = PC2, color = Time)) +
            geom_point(size = 3) +
            ggforce::geom_mark_ellipse(aes(fill = Time),
                con.cap = 0, expand = unit(2, "mm")
            ) +
            scale_color_viridis_d(option = "D") +
            scale_fill_viridis_d(option = "D") +
            theme_custom(base_size = fontsize) +
            xlim(range(pc$PC1) * 1.2) +
            ylim(range(pc$PC2) * 1.2) +
            xlab(xlab) +
            ylab(ylab)
    } else {
        p1 <- pc |>
            dplyr::mutate(Time = factor(Time,
                levels = unique(pc$Time) |> sort())) |>
            ggplot(aes(x = PC1, y = PC2, color = Time)) +
            geom_point(size = 3) +
            scale_color_viridis_d(option = "D") +
            scale_fill_viridis_d(option = "D") +
            theme_custom(base_size = fontsize) +
            xlab(xlab) +
            ylab(ylab)
    }

    return(p1)
}
