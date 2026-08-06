#' Plot PCA
#'
#' @description Plot principal component analysis (PCA), using features
#' without missing values
#'
#' The samples are coloured by Group and sized by Time. Ellipses are drawn
#' for each Group.
#'
#' @param se_obj A SummarizedExperiment object created by `create_input()`,
#'   or a named list of SummarizedExperiment objects (one per dataset)
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
#' (default: NULL, auto-computed with 0.35x padding around the PC1 range)
#' @param xlim_max Maximum x-axis limit when drawing ellipses
#' (default: NULL, auto-computed with 0.35x padding around the PC1 range)
#' @param ylim_min Minimum y-axis limit when drawing ellipses
#' (default: NULL, auto-computed with 0.35x padding around the PC2 range)
#' @param ylim_max Maximum y-axis limit when drawing ellipses
#' (default: NULL, auto-computed with 0.35x padding around the PC2 range)
#' @param fontsize Font size for the PCA plot (default is 8)
#' @param assay The assay to use in the SummarizedExperiment object: a numeric
#'   index or character name, e.g. 1 or "orig" for original data, 2 or "norm"
#'   for time 0 normalised data. (Default: 1.)
#' @param shared_only Logical, whether to subset each dataset to the features
#'   shared by all datasets before PCA (default is FALSE; only applies
#'   when `se_obj` is a list of multiple datasets)
#' @param names Optional character vector of dataset names for the plot
#'   titles. Default: NULL, i.e. the names of `se_obj` (when a list).
#'
#' @import ggplot2
#' @import SummarizedExperiment
#'
#' @returns A list with `p_list` (a named list of ggplot objects) and `pca`
#' (a PCAtools output object for custom plotting). When `se_obj` is a named
#' list of datasets, `p_list` and `pca` are named lists (one per dataset) and
#' the combined PCA scatter plots are returned as `plot`.
#' @export
#' @examples
#' data(example_obj)
#' data(example_obj_2)
#' PC = plot_pca(example_obj, morepc = seq(1, 3))
#' PC$p_list
#' PC_multi = plot_pca(list(Exp1 = example_obj, Exp2 = example_obj_2))
#' PC_multi$plot
plot_pca <- function(
    se_obj,
    plot_screeplot = TRUE,
    plot_loadings = TRUE,
    plot_morepc = TRUE,
    morepc = seq(1, 5),
    pc1 = 1, pc2 = 2,
    circle = FALSE,
    xlim_min = NULL, xlim_max = NULL,
    ylim_min = NULL, ylim_max = NULL,
    fontsize = 8,
    assay = 1,
    shared_only = FALSE,
    names = NULL
) {
    .check_logical(shared_only, "shared_only")
    .check_logical(circle, "circle")
    .check_logical(plot_screeplot, "plot_screeplot")
    .check_logical(plot_loadings, "plot_loadings")
    .check_logical(plot_morepc, "plot_morepc")
    if (!missing(pc1)) .check_positive_int(pc1, "pc1")
    if (!missing(pc2)) .check_positive_int(pc2, "pc2")
    .check_positive(fontsize, "fontsize")
    .check_numeric(morepc, "morepc")
    if (!is.null(xlim_min)) .check_numeric(xlim_min, "xlim_min")
    if (!is.null(xlim_max)) .check_numeric(xlim_max, "xlim_max")
    if (!is.null(ylim_min)) .check_numeric(ylim_min, "ylim_min")
    if (!is.null(ylim_max)) .check_numeric(ylim_max, "ylim_max")
    se_obj_list <- .as_se_list(se_obj)
    nms <- base::names(se_obj_list)
    multi <- length(se_obj_list) > 1L
    multi_style <- multi || !is.null(names)
    if (!is.null(names)) {
        .check_character(names, "names")
        if (length(names) != length(se_obj_list)) {
            stop("'names' must have one entry per dataset.")
        }
    } else if (multi_style) {
        names <- nms
    }
    resolved <- vapply(nms, function(nm) {
        .match_assay(assay, se_obj_list[[nm]])
    }, character(1))
    if (!all(resolved == resolved[1])) {
        stop("The 'assay' must be the same across all datasets.")
    }
    assay <- unname(resolved[1])
    if (multi && shared_only) {
        common <- Reduce(intersect, lapply(se_obj_list, rownames))
        if (length(common) == 0L) {
            stop("No features are shared by all datasets.")
        }
        se_obj_list <- lapply(se_obj_list, function(se) {
            se[intersect(rownames(se), common), , drop = FALSE]
        })
        message("Subsetting features: ",
            paste(nms, vapply(se_obj_list, nrow, integer(1)),
            sep = "=", collapse = ", "))
    }
    results <- lapply(seq_along(se_obj_list), function(i) {
        se <- se_obj_list[[i]]
        M_0 <- assay(se, assay)
        if (sum(is.na(M_0)) > 0) {
            if (multi) {
                message("Input data contains missing values in '", nms[[i]],
                    "'. Only complete rows will be used for PCA.")
            } else {
                message("Input data contains missing values. Only complete",
                    " rows will be used for PCA.")
            }
        }
        M <- M_0[stats::complete.cases(M_0), ]
        sp_info <- as.data.frame(colData(se))

        pca2 <- PCAtools::pca(M, metadata = sp_info, center = TRUE)

        p_list <- list()

        if (plot_screeplot) {
            p3 <- PCAtools::screeplot(pca2,
                components = PCAtools::getComponents(pca2)[seq(1, 10)],
                titleLabSize = fontsize,
                axisLabSize = fontsize
            )
            p_list$p3 <- p3
        }
        if (plot_loadings) {
            p4 <- PCAtools::plotloadings(pca2,
                titleLabSize = fontsize,
                legendLabSize = fontsize,
                # labSize = fontsize - 4,
                axisLabSize = fontsize
            )
            p_list$p4 <- p4
        }

        pc <- as.data.frame(pca2$rotated)
        pc$Sample <- rownames(pc)
        pc <- merge(pc, sp_info, by = "Sample")

        pc1_name <- paste0("PC", pc1)
        pc2_name <- paste0("PC", pc2)

        n_pcs <- ncol(pca2$rotated)
        if (pc1 > n_pcs || pc2 > n_pcs) {
            stop("'pc1' and 'pc2' must be within the available PCs (",
                n_pcs, ").")
        }

        xlab <- paste0(pc1_name, ": ", round(pca2$variance[[pc1_name]], 2), "%")
        ylab <- paste0(pc2_name, ": ", round(pca2$variance[[pc2_name]], 2), "%")

        if (multi_style) {
            title <- paste0(names[[i]], " (", nrow(M), " features)")
        } else {
            title <- paste0("PCA (", nrow(M),
                " features without missing values)")
        }
        p1 <- pc |>
            dplyr::mutate(Time = factor(Time)) |>
            ggplot(aes(x = !!sym(pc1_name), y = !!sym(pc2_name),
                fill = Group, label = Sample)) +
            geom_point(aes(size = Time), alpha = 0.8, shape = 21) +
            scale_fill_manual(values =
                get_custom_palette(unique(colData(se)$Group))) +
            scale_size_discrete() +
            theme_custom(base_size = fontsize) +
            guides(fill = guide_legend(override.aes = list(size = 4))) +
            xlab(xlab) +
            ylab(ylab) +
            ggtitle(title)
        p_list$p1 <- p1

        if (is.null(xlim_min) || is.null(xlim_max)) {
            x_rng <- range(pc[[pc1_name]], na.rm = TRUE)
            x_pad <- 0.35 * diff(x_rng)
            if (is.null(xlim_min)) {
                xlim_min <- x_rng[1] - x_pad
            }
            if (is.null(xlim_max)) {
                xlim_max <- x_rng[2] + x_pad
            }
        }
        if (is.null(ylim_min) || is.null(ylim_max)) {
            y_rng <- range(pc[[pc2_name]], na.rm = TRUE)
            y_pad <- 0.35 * diff(y_rng)
            if (is.null(ylim_min)) {
                ylim_min <- y_rng[1] - y_pad
            }
            if (is.null(ylim_max)) {
                ylim_max <- y_rng[2] + y_pad
            }
        }

        if (circle) {
            p1$layers <- c(p1$layers, list(ggforce::geom_mark_ellipse(
                aes(fill = Group, label = Group),
                con.cap = 0, alpha = 0.1, label.fontsize = fontsize,
                label.buffer = unit(0, "mm")
            )))
            p2 <- p1 +
                coord_cartesian(
                    xlim = c(xlim_min, xlim_max),
                    ylim = c(ylim_min, ylim_max)
                )
            p_list$p2 <- p2
        }

        # plot more PCs
        if (plot_morepc) {
            morepc_i <- morepc[morepc %in% seq(1, ncol(pca2$rotated))]
            if (length(morepc_i) != length(morepc)) {
                warning("Some specified PCs for the pairs plot are out of ",
                    "range. Please check the number of PCs available.")
            }
            if (length(morepc_i) < 2) {
                warning("At least two valid PCs are required for the pairs ",
                    "plot. Please specify more PCs.")
            } else {
                p5 <- PCAtools::pairsplot(pca2,
                    colby = "Group", components = paste0("PC", morepc_i),
                    hline = 0, vline = 0, gridlines.major = FALSE,
                    gridlines.minor = FALSE,
                    plotaxes = FALSE,
                    trianglelabSize = fontsize,
                    colkey = get_custom_palette(colData(se)$Group |> unique())
                )
                p_list$p5 <- p5
            }
        }
        list(p_list = p_list, pca = pca2)
    })
    if (!multi) {
        return(results[[1]])
    }
    scatter <- patchwork::wrap_plots(
        lapply(results, function(r) {
            if (circle) r$p_list[["p2"]] else r$p_list[["p1"]]
        }),
        guides = "collect"
    ) + patchwork::plot_annotation(
        title = paste0("PCA across datasets",
            if (shared_only) " (shared features without missing values)"
            else " (all features without missing values)"),
        theme = ggplot2::theme(legend.position = "right")
    )
    return(list(
        p_list = stats::setNames(lapply(results, function(r) r$p_list), nms),
        pca = stats::setNames(lapply(results, function(r) r$pca), nms),
        plot = scatter
    ))
}
