#' Plot UMAP
#'
#' @description Plot UMAP of samples, using features without missing values
#'
#' The UMAP plots are labelled by Group, Time, Replicate (if more than 1),
#' Batch (if more than 1), and number of identified features (non-missing
#' values).
#'
#' @param se_obj A SummarizedExperiment object created by `create_input()`,
#'   or a named list of SummarizedExperiment objects (one per dataset)
#' @param seed Random seed for UMAP (default is 1234)
#' @param plot_ID Whether to include a UMAP plot coloured by number of
#' identified features (default is FALSE)
#' @param circle Logical, whether to draw circles (ellipses) around samples of
#' each group (default is FALSE)
#' @param xlim_min Minimum x-axis limit when drawing ellipses
#' (default: NULL, auto-computed with 0.35x padding around the UMAP x range)
#' @param xlim_max Maximum x-axis limit when drawing ellipses
#' (default: NULL, auto-computed with 0.35x padding around the UMAP x range)
#' @param ylim_min Minimum y-axis limit when drawing ellipses
#' (default: NULL, auto-computed with 0.35x padding around the UMAP y range)
#' @param ylim_max Maximum y-axis limit when drawing ellipses
#' (default: NULL, auto-computed with 0.35x padding around the UMAP y range)
#' @param umap_neighbors UMAP n_neighbors parameter (default is selected by
#' `.umap_n_neighbors()` function based on the number of samples)
#' @param fontsize Font size for the plot (default is 8)
#' @param assay The assay to use in the SummarizedExperiment object: a numeric
#'   index or character name, e.g. 1 or "orig" for original data, 2 or "norm"
#'   for time 0 normalised data. (Default: 1.)
#' @param shared_only Logical, whether to subset each dataset to the features
#'   shared by all datasets before UMAP (default is FALSE; only applies
#'   when `se_obj` is a list of multiple datasets)
#' @param names Optional character vector of dataset names for the plot
#'   titles. Default: NULL, i.e. the names of `se_obj` (when a list).
#'
#' @import SummarizedExperiment
#' @import ggplot2
#'
#' @returns A list with `p_list` (a named list of ggplot objects) and
#' `umap_layout` (a data frame containing UMAP coordinates and sample
#' annotations for custom plotting). When `se_obj` is a named list of datasets,
#' both are named lists (one per dataset) and the combined UMAP plots are
#' returned as `plot`.
#' @export
#' @examples
#' data(example_obj)
#' data(example_obj_2)
#' umap_out <- plot_umap(example_obj)
#' umap_out$p_list
#' umap_multi <- plot_umap(list(Exp1 = example_obj, Exp2 = example_obj_2))
#' umap_multi$plot
plot_umap <- function(
    se_obj, seed = 1234, plot_ID = FALSE,
    circle = FALSE,
    xlim_min = NULL, xlim_max = NULL,
    ylim_min = NULL, ylim_max = NULL,
    umap_neighbors = NULL,
    fontsize = 8,
    assay = 1,
    shared_only = FALSE,
    names = NULL
) {
    .check_logical(shared_only, "shared_only")
    .check_logical(plot_ID, "plot_ID")
    .check_logical(circle, "circle")
    .check_positive_int(seed, "seed")
    .check_positive(fontsize, "fontsize")
    if (!is.null(xlim_min)) .check_numeric(xlim_min, "xlim_min")
    if (!is.null(xlim_max)) .check_numeric(xlim_max, "xlim_max")
    if (!is.null(ylim_min)) .check_numeric(ylim_min, "ylim_min")
    if (!is.null(ylim_max)) .check_numeric(ylim_max, "ylim_max")
    if (!is.null(umap_neighbors))
        .check_positive_int(umap_neighbors, "umap_neighbors")
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
                    "'. Only complete rows will be used for UMAP.")
            } else {
                message("Input data contains missing values. Only complete ",
                    "rows will be used for UMAP.")
            }
        }
        M <- M_0[stats::complete.cases(M_0), ]
        sp_info <- as.data.frame(colData(se))

        n_neighbors <- ifelse(is.null(umap_neighbors),
            .umap_n_neighbors(ncol(M)),
            umap_neighbors)

        if (multi) {
            message(nms[[i]], ": using n_neighbors = ", n_neighbors)
        } else {
            message("Using n_neighbors = ", n_neighbors)
        }

        umap <- umap::umap(M |> t(), random_state = seed,
            n_neighbors = n_neighbors)
        umap_layout <- umap$layout |> as.data.frame()
        umap_layout$Sample <- rownames(umap_layout)
        umap_layout <- merge(umap_layout, sp_info, by = "Sample")
        if (multi_style) {
            umap_title <- paste0(names[[i]], " (", nrow(M), " features)")
        } else {
            umap_title <- paste0("UMAP (", nrow(M),
                " features without missing values)")
        }

        # number of identified proteins
        id_count <- apply(M_0, 2, function(x) sum(!is.na(x))) |> as.data.frame()
        colnames(id_count) <- "ID"
        id_count$Sample <- rownames(id_count)
        umap_layout <- merge(umap_layout, id_count, by = "Sample")

        umap_list <- list()
        umap_list[["Group"]] <- ggplot(umap_layout,
            aes(x = V1, y = V2, color = Group)) +
            geom_point(size = 2) +
            scale_color_manual(values =
                get_custom_palette(unique(colData(se)$Group))) +
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

        p_list <- list()

        p_list$p3 <- ggpubr::ggarrange(
            plotlist = umap_list,
            nrow = 2, ncol = ceiling(length(umap_list) / 2),
            common.legend = FALSE
        ) |>
            ggpubr::annotate_figure(top = ggpubr::text_grob(
                paste0(umap_title, "\n"),
                face = "bold", size = fontsize + 2))

        p1 <- umap_layout |>
            dplyr::mutate(Time = Time |> factor()) |>
            ggplot(aes(x = V1, y = V2, fill = Group)) +
            geom_point(aes(size = Time), alpha = 0.8, shape = 21) +
            scale_fill_manual(values =
                get_custom_palette(unique(colData(se)$Group))) +
            scale_size_discrete() +
            guides(fill = guide_legend(override.aes = list(size = 4))) +
            theme_custom(base_size = fontsize) +
            ggtitle(umap_title)
        p_list$p1 <- p1

        if (is.null(xlim_min) || is.null(xlim_max)) {
            x_rng <- range(umap_layout$V1, na.rm = TRUE)
            x_pad <- 0.35 * diff(x_rng)
            if (is.null(xlim_min)) {
                xlim_min <- x_rng[1] - x_pad
            }
            if (is.null(xlim_max)) {
                xlim_max <- x_rng[2] + x_pad
            }
        }
        if (is.null(ylim_min) || is.null(ylim_max)) {
            y_rng <- range(umap_layout$V2, na.rm = TRUE)
            y_pad <- 0.35 * diff(y_rng)
            if (is.null(ylim_min)) {
                ylim_min <- y_rng[1] - y_pad
            }
            if (is.null(ylim_max)) {
                ylim_max <- y_rng[2] + y_pad
            }
        }

        if (circle) {
            p1$layers <- c(p1$layers, list(
                ggforce::geom_mark_ellipse(
                    aes(fill = Group, label = Group),
                    con.cap = 0, alpha = 0.1, label.fontsize = fontsize,
                    label.buffer = unit(0, "mm")
                )))
            p_list$p2 <- p1 +
                coord_cartesian(
                    xlim = c(xlim_min, xlim_max),
                    ylim = c(ylim_min, ylim_max)
                )
        }

        list(p_list = p_list, umap_layout = umap_layout)
    })
    if (!multi) {
        return(results[[1]])
    }
    umap_all <- patchwork::wrap_plots(
        lapply(results, function(r) {
            if (circle) r$p_list[["p2"]] else r$p_list[["p1"]]
        }),
        guides = "collect"
    ) + patchwork::plot_annotation(
        title = paste0("UMAP across datasets",
            if (shared_only) " (shared features without missing values)"
            else " (all features without missing values)"),
        theme = ggplot2::theme(legend.position = "right")
    )
    return(list(
        p_list = stats::setNames(lapply(results, function(r) r$p_list), nms),
        umap_layout = stats::setNames(
            lapply(results, function(r) r$umap_layout), nms),
        plot = umap_all
    ))
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
.umap_n_neighbors <- function (sample_n) {
    if (sample_n <= 5) {
        return(sample_n - 1)
    } else {
        return(floor(min(15, max(sample_n / 5, 5))))
    }
}
