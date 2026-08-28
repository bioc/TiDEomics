#' Plot correlation matrix
#'
#' @description Plot correlation matrix between samples as a heatmap
#'
#' @details With a named list of datasets that share the same samples,
#'   the per-dataset dendrograms are compared pairwise with the
#'   cophenetic correlation.
#'
#' @param se_obj A SummarizedExperiment object created by `create_input()`,
#'   or a named list of SummarizedExperiment objects (one per dataset).
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
#'   FALSE)
#' @param show_colnames Whether to show column names in the heatmap (default is
#'   FALSE)
#' @param fontsize Font size for the plot (default is 8)
#' @param cellwidth Cell width for the heatmap (default is 1)
#' @param cellheight Cell height for the heatmap (default is 1)
#' @param title Title of the heatmap. If NULL, auto-generated as
#'   "Sample correlation (<Method>)" with the method capitalised
#'   (e.g. "Sample correlation (Spearman)") (default: NULL)
#' @param shared_only Logical, whether to subset each dataset to the features
#'   shared by all datasets before computing correlations (default is FALSE;
#'   only applies when `se_obj` is a list of multiple datasets)
#' @param names Optional character vector of dataset names for the per-heatmap
#'   titles. Default: NULL, i.e. the names of `se_obj` (when a list).
#'   Only applies when `se_obj` is a list of multiple datasets.
#' @param plot_cophenetic Logical, whether to include the cophenetic
#'   correlation in the combined plot title (default is FALSE) in addition
#'   to reporting by text.
#' @param ... Additional arguments to be passed to `ComplexHeatmap::Heatmap()`
#'
#' @import SummarizedExperiment
#'
#' @returns A heatmap showing the correlation between samples. With a named
#'   list of datasets, all heatmaps share the same colour scale,
#'   include a cophenetic-correlation subtitle when `plot_cophenetic` is TRUE,
#'   and are returned invisibly as a named list of Heatmap objects.
#' @export
#' @examples
#' data(example_obj)
#' data(example_obj_2)
#' plot_cor_matrix(example_obj)
#' plot_cor_matrix(list(Exp1 = example_obj, Exp2 = example_obj_2))
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
    shared_only = FALSE,
    names = NULL,
    plot_cophenetic = FALSE,
    ...
) {
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
    .check_logical(shared_only, "shared_only")
    .check_logical(plot_cophenetic, "plot_cophenetic")
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

    se_obj_list <- .as_se_list(se_obj)
    nms <- base::names(se_obj_list)
    multi <- length(se_obj_list) > 1L
    if (multi && !is.null(names)) {
        .check_character(names, "names")
        if (length(names) != length(se_obj_list)) {
            stop("'names' must have one entry per dataset.")
        }
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
    disp_names <- if (is.null(names)) nms else names

    cor_mats <- lapply(se_obj_list, function(se) {
        stats::cor(assay(se, assay), use = use, method = method)
    })

    all_cors <- unlist(lapply(cor_mats, function(cm) {
        as.vector(cm[upper.tri(cm)])
    }))
    min_val <- min(all_cors, na.rm = TRUE)
    max_val <- max(all_cors, na.rm = TRUE)
    col_fun <- circlize::colorRamp2(
        c(min_val, (min_val + max_val) / 2, max_val),
        c("#3C5488FF", "white", "#E64B35FF")
    )

    sim_text <- NULL
    if (multi) {
        colnames_list <- lapply(se_obj_list, function(se) {
            colnames(assay(se, assay))
        })
        no_dup <- all(vapply(colnames_list, function(cn) {
            !anyDuplicated(cn)
        }, logical(1)))
        matched <- no_dup && all(vapply(colnames_list, function(cn) {
            setequal(cn, colnames_list[[1]])
        }, logical(1)))
        if (matched) {
            first_cn <- colnames_list[[1]]
            cophs <- lapply(seq_along(cor_mats), function(i) {
                cm <- cor_mats[[i]][first_cn, first_cn]
                as.vector(stats::cophenetic(
                    stats::hclust(stats::as.dist(1 - cm))))
            })
            pair_names <- utils::combn(disp_names, 2, paste, collapse = "-")
            r_vals <- utils::combn(seq_along(nms), 2, function(idx) {
                stats::cor(cophs[[idx[1]]], cophs[[idx[2]]])
            })
            r_str <- sprintf("%.2f", r_vals)
            if (length(r_vals) == 1) {
                sim_text <- paste0("Cophenetic r = ", r_str)
            } else {
                pair_text <- paste0(pair_names, " = ", r_str)
                chunks <- split(pair_text,
                    ceiling(seq_along(pair_text) / 3))
                lines <- vapply(chunks, paste, character(1),
                    collapse = ", ")
                sim_text <- paste0("Cophenetic r: ",
                    paste(lines, collapse = "\n"))
            }
            message(sim_text)
        } else {
            message("Skipping similarity statistics: datasets have ",
                "different samples.")
        }
    }

    legend_i <- length(se_obj_list)
    heatmaps <- lapply(seq_along(se_obj_list), function(i) {
        se <- se_obj_list[[i]]
        nm <- nms[i]
        cor_table <- cor_mats[[i]]

        cd <- as.data.frame(colData(se))

        ann_df <- cd[, character(0), drop = FALSE]
        ann_colors <- list()
        if (label_group && "Group" %in% colnames(cd)) {
            ann_df$Group <- cd$Group
            ann_colors$Group <- get_custom_palette(unique(cd$Group))
        }
        if (label_time && "Time" %in% colnames(cd)) {
            time_levels <- unique(cd$Time) |> sort()
            ann_df$Time <- factor(cd$Time, levels = time_levels)
            ann_colors$Time <- scales::pal_viridis()(
                length(time_levels)) |>
                stats::setNames(nm = time_levels)
        }
        if (label_rep && "Replicate" %in% colnames(cd)) {
            ann_df$Replicate <- cd$Replicate
            rep_levels <- unique(cd$Replicate) |> sort()
            ann_colors$Replicate <- ggsci::pal_iterm()(
                length(rep_levels)) |>
                stats::setNames(nm = rep_levels)
        }
        if (label_batch && "Batch" %in% colnames(cd)) {
            ann_df$Batch <- cd$Batch
            batch_levels <- unique(cd$Batch) |> sort()
            ann_colors$Batch <- ggsci::pal_simpsons()(
                length(batch_levels)) |>
                stats::setNames(nm = batch_levels)
        }

        ha <- if (ncol(ann_df) > 0) {
            ComplexHeatmap::HeatmapAnnotation(
                df = ann_df,
                col = ann_colors,
                annotation_name_side = "left",
                show_annotation_name = TRUE,
                show_legend = i == legend_i,
                annotation_legend_param = list(
                    title_gp = grid::gpar(fontsize = fontsize),
                    labels_gp = grid::gpar(fontsize = fontsize)
                ),
                annotation_name_gp = grid::gpar(fontsize = fontsize),
                border = TRUE
            )
        } else {
            NULL
        }

        cor_cluster <- stats::hclust(stats::as.dist(1 - cor_table))

        ComplexHeatmap::Heatmap(cor_table,
            name = nm,
            col = col_fun,
            cluster_rows = cor_cluster,
            cluster_columns = cor_cluster,
            show_row_names = show_rownames,
            show_column_names = show_colnames,
            show_row_dend = TRUE,
            show_column_dend = TRUE,
            show_heatmap_legend = i == legend_i,
            column_title = if (multi) {
                paste0(disp_names[i], " (", nrow(se), ")")
            } else {
                title
            },
            column_title_gp = grid::gpar(fontsize = fontsize + 2,
                fontface = "bold"),
            top_annotation = ha,
            width = unit(ncol(cor_table) * cellwidth / 1.5, "mm"),
            height = unit(nrow(cor_table) * cellheight / 1.5, "mm"),
            heatmap_legend_param = list(
                title = "r",
                title_gp = grid::gpar(fontsize = fontsize),
                labels_gp = grid::gpar(fontsize = fontsize)
            ),
            ...
        )
    })

    if (!multi) {
        return(heatmaps[[1]])
    }
    heatmaps <- stats::setNames(heatmaps, nms)

    title_final <- if (!plot_cophenetic || is.null(sim_text)) {
        title
    } else {
        paste0(title, "\n", sim_text)
    }
    title_lines <- length(strsplit(title_final, "\n",
        fixed = TRUE)[[1]]) + 1L

    # Capture each heatmap as a grob and
    # arrange side-by-side sized to exact content width
    # (eliminates white gaps from equal-width columns).
    grobs <- lapply(heatmaps, function(ht) {
        grid::grid.grabExpr(ComplexHeatmap::draw(ht,
            padding = unit(c(0, 0, 0, 0), "mm")
        ))
    })

    # Compute exact content dimensions
    content_dims <- lapply(heatmaps, function(ht) {
        show_al <- !is.null(ht@top_annotation)
        ht_list <- methods::new("HeatmapList")
        ht_list <- ComplexHeatmap::add_heatmap(ht_list, ht)
        ht_list <- ComplexHeatmap::make_layout(ht_list,
            padding = unit(c(0, 0, 0, 0), "mm"),
            show_heatmap_legend = ht@heatmap_param$show_heatmap_legend,
            show_annotation_legend = show_al
        )
        cw <- ComplexHeatmap::component_width(ht_list)
        ch <- ComplexHeatmap::component_height(ht_list)
        c(
            w = grid::convertWidth(sum(cw), "mm", valueOnly = TRUE),
            h = grid::convertHeight(sum(ch), "mm", valueOnly = TRUE)
        )
    })
    content_w <- vapply(content_dims, `[[`, numeric(1), "w")
    content_h <- max(vapply(content_dims, `[[`, numeric(1), "h"))

    grid::grid.newpage()
    grid::pushViewport(grid::viewport(
        layout = grid::grid.layout(2, 1,
            heights = unit(c(title_lines, content_h),
                c("lines", "mm"))
        )
    ))
    grid::pushViewport(grid::viewport(layout.pos.row = 1))
    grid::grid.text(title_final,
        gp = grid::gpar(fontsize = fontsize + 4, fontface = "bold")
    )
    grid::popViewport()
    grid::pushViewport(grid::viewport(layout.pos.row = 2))
    grid::pushViewport(grid::viewport(
        layout = grid::grid.layout(1, length(heatmaps),
            widths = unit(content_w, "mm")
        )
    ))
    for (i in seq_along(grobs)) {
        grid::pushViewport(grid::viewport(layout.pos.col = i))
        grid::grid.draw(grobs[[i]])
        grid::popViewport()
    }
    grid::popViewport()
    grid::popViewport()
    grid::popViewport()
    invisible(heatmaps)
}
