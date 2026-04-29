#' Plot PCA by group (one object)
#'
#' @description Plot PCA for each group separately, with optional circles 
#' around time points and arrows indicating trajectory over time.
#'
#' @param se_obj A SummarizedExperiment object created by `create_input()`
#' @param circle Logical, whether to draw circles (ellipses) around samples 
#' of each time point (default is TRUE)
#' @param arrow Logical, whether to draw arrows indicating the trajectory 
#' over time (default is TRUE)
#' @param nrow Number of rows for arranging the PCA plots (default is 1)
#' @param fontsize Font size for the PCA plots (default is 8)
#' @param assay Assay index to use, where 1 is the original data and 2 is 
#' normalised to time 0 (if available) (default is 1)
#'
#' @import ggplot2
#' @import SummarizedExperiment
#' @import magrittr
#' @importFrom dplyr arrange
#'
#' @returns A series of PCA plots showing the distribution of samples in 
#' each group, coloured by Time.
#' @export
#' @examples
#' data("example")
#' plot_pca_by_group(example_obj, circle = TRUE, arrow = TRUE)
plot_pca_by_group <- function(se_obj, circle = TRUE, arrow = TRUE, 
    nrow = 1, fontsize = 8, assay = 1) {
    pca_list <- list()
    for (group in unique(se_obj$Group)) {
        p1 <- plot_pca_arrows(se_obj[, se_obj$Group == group],
            circle = circle,
            arrow = arrow,
            fontsize = fontsize,
            assay = assay
        )
        pca_list[[group]] <- p1
    }

    ncol <- ceiling(length(pca_list) / nrow)

    print(ggpubr::ggarrange(
        plotlist = pca_list, labels = names(pca_list),
        font.label = list(size = fontsize + 2),
        hjust = 0, vjust = 0.5,
        nrow = nrow, ncol = ncol,
        common.legend = TRUE
    ) %>%
        ggpubr::annotate_figure(top = 
            ggpubr::text_grob("PCA by group (features without missing values)", 
            face = "bold", size = fontsize + 2)))
}


#' Plot PCA by group (list of objects)
#'
#' @description Plot PCA by group, input is a list of SummarizedExperiment 
#' objects, with each object corresponding to a group
#'
#' @param se_obj_list A list of SummarizedExperiment objects, such as output 
#' of `split_groups()`
#' @param circle Logical, whether to draw circles (ellipses) around samples 
#' of each time point (default is TRUE)
#' @param arrow Logical, whether to draw arrows indicating the trajectory 
#' over time (default is TRUE)
#' @param nrow Number of rows for arranging the PCA plots (default is 1)
#' @param fontsize Font size for the PCA plots (default is 8)
#' @param assay Assay index to use, where 1 is the original data and 2 is 
#' normalised to time 0 (if available) (default is 1)
#'
#' @import ggplot2
#' @import SummarizedExperiment
#' @import magrittr
#' @importFrom dplyr arrange
#' @returns A series of PCA plots showing the distribution of samples in 
#' each group, coloured by Time.
#' @export
#' @examples
#' data("example")
#' example_obj_list <- split_groups(example_obj)
#' plot_pca_by_group_list(example_obj_list)
plot_pca_by_group_list <- function(se_obj_list, circle = TRUE, arrow = TRUE, 
    nrow = 1, fontsize = 8, assay = 1) {
    pca_list <- list()
    for (group in names(se_obj_list)) {
        p1 <- plot_pca_arrows(se_obj_list[[group]], circle = circle,
            arrow = arrow, fontsize = fontsize, assay = assay)
        pca_list[[group]] <- p1
    }

    ncol <- ceiling(length(pca_list) / nrow)

    print(ggpubr::ggarrange(
        plotlist = pca_list, labels = names(pca_list),
        font.label = list(size = fontsize + 2),
        hjust = 0, vjust = 0.5,
        nrow = nrow, ncol = ncol,
        common.legend = TRUE
    ) %>%
        ggpubr::annotate_figure(top =
            ggpubr::text_grob("PCA by group (features without missing values)", 
            face = "bold", size = fontsize + 2)))
}


#' Plot PCA with arrows
#'
#' @description Plot PCA with arrows indicating trajectory over time
#'
#' @param se_obj A SummarizedExperiment object
#' @param circle Logical, whether to draw circles (ellipses) around samples 
#' of each time point (default is TRUE)
#' @param arrow Logical, whether to draw arrows indicating the trajectory 
#' over time (default is TRUE)
#' @param fontsize Font size for the PCA plot (default is 8)
#' @param assay Assay index to use, where 1 is the original data and 2 is 
#' normalised to time 0 (if available) (default is 1)
#'
#' @import ggplot2
#' @import SummarizedExperiment
#' @import magrittr
#' @importFrom dplyr arrange
#' @importFrom stats aggregate
#'
#' @returns A PCA plot showing the distribution of samples, coloured by 
#' Time, with arrows indicating the trajectory over time.
#' @export
#' @examples
#' data("example")
#' plot_pca_arrows(example_obj[, example_obj$Group == "IFNbeta"])
plot_pca_arrows <- function(se_obj, circle = TRUE, arrow = TRUE, 
    fontsize = 8, assay = 1) {

    if (length(unique(se_obj$Group)) > 1) {
        warning("Input object contains multiple groups, circles / arrows ", 
        "will be calculated across all samples. Please use ", 
        "`plot_pca_by_group()` or subset the object to one group ", 
        "before plotting.")
    }

    pca2 <- plot_pca(se_obj, plot = FALSE, assay = assay)
    pc <- as.data.frame(pca2$rotated)
    pc$Sample <- rownames(pc)
    pc <- merge(pc, as.data.frame(colData(se_obj)), by = "Sample")

    xlab <- paste0("PC1: ", round(pca2$variance[["PC1"]], 2), "%")
    ylab <- paste0("PC2: ", round(pca2$variance[["PC2"]], 2), "%")

    arrows_in <- pc[, c("PC1", "PC2", "Sample", "Time")] %>% arrange(Time)
    x_arrows <- aggregate(PC1 ~ Time, data = arrows_in, FUN = mean)
    y_arrows <- aggregate(PC2 ~ Time, data = arrows_in, FUN = mean)
    arrows <- merge(x_arrows, y_arrows, by = "Time") %>%
        dplyr::rename(x_start = PC1, y_start = PC2)

    time_series <- sort(unique(x_arrows$Time)) # not necessarily 1,2,3...

    for (i in seq(1, length(time_series))) {
        if (i == length(time_series)) next
        arrows[which(arrows$Time == time_series[i]), "x_end"] <- 
            arrows[which(arrows$Time == time_series[i + 1]), "x_start"]
        arrows[which(arrows$Time == time_series[i]), "y_end"] <- 
            arrows[which(arrows$Time == time_series[i + 1]), "y_start"]
    }

    if (circle & arrow) {
        p1 <- pc %>%
            mutate(Time = factor(Time, levels = unique(pc$Time) %>% sort())) %>%
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
            ggforce::geom_mark_ellipse(aes(fill = Time),
                con.cap = 0, expand = unit(2, "mm")
            ) +
            scale_color_viridis_d(option = "D") +
            scale_fill_viridis_d(option = "D") +
            theme_custom(base_size = fontsize) +
            xlim(range(pc$PC1) * 1.2) + # extra space for circles
            ylim(range(pc$PC2) * 1.2) +
            xlab(xlab) +
            ylab(ylab)
    } else if (arrow) {
        p1 <- pc %>%
            mutate(Time = factor(Time, levels = unique(pc$Time) %>% sort())) %>%
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
        p1 <- pc %>%
            mutate(Time = factor(Time, levels = unique(pc$Time) %>% sort())) %>%
            ggplot(aes(x = PC1, y = PC2, color = Time)) +
            geom_point(size = 3) +
            ggforce::geom_mark_ellipse(aes(fill = Time),
                con.cap = 0, expand = unit(2, "mm")
            ) +
            scale_color_viridis_d(option = "D") +
            scale_fill_viridis_d(option = "D") +
            theme_custom(base_size = fontsize) +
            xlim(range(pc$PC1) * 1.2) + # extra space for circles
            ylim(range(pc$PC2) * 1.2) +
            xlab(xlab) +
            ylab(ylab)
    } else {
        p1 <- pc %>%
            mutate(Time = factor(Time, levels = unique(pc$Time) %>% sort())) %>%
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
