#' Plot PCA
#'
#' @description Plot principal component analysis (PCA), using features
#' without missing values
#'
#' The samples are coloured by Group and sized by Time. Ellipses are drawn
#' for each Group.
#'
#' @param se_obj A SummarizedExperiment object created by `create_input()`
#' @param plot Logical, whether to plot PCA and other plots (default is TRUE)
#' @param plot_screeplot Logical, whether to plot screeplot with
#' `PCAtools::screeplot()` (default is TRUE)
#' @param plot_loadings Logical, whether to plot loadings with
#' `PCAtools::plotloadings()` (default is TRUE)
#' @param plot_morepc Logical, whether to plot pairs of more PCs with
#' `PCAtools::pairsplot()` (default is TRUE)
#' @param morepc A numeric vector of PCs to plot in the pairs plot
#' (default is 1:5)
#' @param pc1 Numeric, which PC to use for the x-axis (default is 1)
#' @param pc2 Numeric, which PC to use for the y-axis (default is 2)
#' @param circle Logical, whether to draw circles (ellipses) around samples
#' of each group (default is FALSE)
#' @param xlim_min Minimum x-axis limit when drawing ellipses
#' (default 1.5*min PC1)
#' @param xlim_max Maximum x-axis limit when drawing ellipses
#' (default 1.5*max PC1)
#' @param ylim_min Minimum y-axis limit when drawing ellipses
#' (default 1.5*min PC2)
#' @param ylim_max Maximum y-axis limit when drawing ellipses
#' (default 1.5*max PC2)
#' @param fontsize Font size for the PCA plot (default is 8)
#' @param assay Assay index to use, where 1 is the original data and 2 is
#' normalised to time 0 (if available) (default is 1)
#'
#' @import ggplot2
#' @import SummarizedExperiment
#'
#' @returns A PCA plot showing the distribution of samples, other plots
#' provided by PCAtools package, and PCAtools output object for custom plotting.
#' @export
#' @examples
#' data("example")
#' PC = plot_pca(example_obj, morepc = seq(1, 3))
plot_pca <- function(
    se_obj,
    plot = TRUE,
    plot_screeplot = TRUE,
    plot_loadings = TRUE,
    plot_morepc = TRUE,
    morepc = seq(1, 5),
    pc1 = 1, pc2 = 2,
    circle = FALSE,
    xlim_min = NULL, xlim_max = NULL,
    ylim_min = NULL, ylim_max = NULL,
    fontsize = 8,
    assay = 1
) {
    .check_se(se_obj)
    .check_logical(plot, "plot")
    .check_logical(circle, "circle")
    .check_logical(plot_screeplot, "plot_screeplot")
    .check_logical(plot_loadings, "plot_loadings")
    .check_logical(plot_morepc, "plot_morepc")
    if (!missing(pc1)) .check_positive_int(pc1, "pc1")
    if (!missing(pc2)) .check_positive_int(pc2, "pc2")
    .check_positive(fontsize, "fontsize")
    assay <- .match_assay(assay, se_obj)
    M_0 <- assay(se_obj, assay)
    if (sum(is.na(M_0)) > 0) {
        message("Input data contains missing values. Only complete rows ",
        "will be used for PCA.")
    }
    M <- M_0[stats::complete.cases(M_0), ]
    sp_info <- as.data.frame(colData(se_obj))

    pca2 <- PCAtools::pca(M, metadata = sp_info, center = TRUE)
    if (!plot) {
        return(pca2)
    }

    if (plot_screeplot) {
        print(PCAtools::screeplot(pca2,
            components = PCAtools::getComponents(pca2)[seq(1, 10)],
            titleLabSize = fontsize,
            axisLabSize = fontsize
        ))
    }
    if (plot_loadings) {
        print(PCAtools::plotloadings(pca2,
            titleLabSize = fontsize,
            legendLabSize = fontsize,
            # labSize = fontsize - 4,
            axisLabSize = fontsize
        ))
    }

    pc <- as.data.frame(pca2$rotated)
    pc$Sample <- rownames(pc)
    pc <- merge(pc, sp_info, by = "Sample")

    pc1_name <- paste0("PC", pc1)
    pc2_name <- paste0("PC", pc2)

    xlab <- paste0(pc1_name, ": ", round(pca2$variance[[pc1_name]], 2), "%")
    ylab <- paste0(pc2_name, ": ", round(pca2$variance[[pc2_name]], 2), "%")

    title <- paste0("PCA (", nrow(M), " features without missing values)")
    p1 <- pc |>
        dplyr::mutate(Time = factor(Time)) |>
        ggplot(aes(x = !!sym(pc1_name), y = !!sym(pc2_name),
            fill = Group, label = Sample)) +
        geom_point(aes(size = Time), alpha = 0.8, shape = 21) +
        scale_fill_manual(values =
            get_custom_palette(unique(colData(se_obj)$Group))) +
        scale_size_discrete() +
        theme_custom(base_size = fontsize) +
        guides(fill = guide_legend(override.aes = list(size = 4))) +
        xlab(xlab) +
        ylab(ylab) +
        ggtitle(title)
    print(p1)

    if (is.null(xlim_min)) {
        xlim_min <- 1.5 * min(pc[[pc1_name]])
    }
    if (is.null(xlim_max)) {
        xlim_max <- 1.5 * max(pc[[pc1_name]])
    }
    if (is.null(ylim_min)) {
        ylim_min <- 1.5 * min(pc[[pc2_name]])
    }
    if (is.null(ylim_max)) {
        ylim_max <- 1.5 * max(pc[[pc2_name]])
    }

    if (circle) {
        p1$layers <- c(p1$layers, list(ggforce::geom_mark_ellipse(
            aes(fill = Group, label = Group),
            con.cap = 0, alpha = 0.1, label.fontsize = fontsize,
            label.buffer = unit(0, "mm")
        )))
        p2 <- p1 +
            xlim(xlim_min, xlim_max) + ylim(ylim_min, ylim_max)
        print(p2)
    }

    # plot more PCs
    if (plot_morepc) {
        if (!all(morepc %in% seq(1, ncol(pca2$rotated)))) {
            warning("Some specified PCs for the pairs plot are out of range. ",
            "Please check the number of PCs available.")
            morepc <- morepc[morepc %in% seq(1, ncol(pca2$rotated))]
        }
        if (length(morepc) < 2) {
            warning("At least two valid PCs are required for the pairs ",
            "plot. Please specify more PCs.")
        } else {
            print(PCAtools::pairsplot(pca2,
                colby = "Group", components = paste0("PC", morepc),
                hline = 0, vline = 0, gridlines.major = FALSE,
                gridlines.minor = FALSE,
                plotaxes = FALSE,
                trianglelabSize = fontsize,
                colkey = get_custom_palette(colData(se_obj)$Group |> unique())
            ))
        }
    }
    return(pca2)
}
