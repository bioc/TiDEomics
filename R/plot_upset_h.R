#' UpSet plot of feature sets with horizontally aligned enrichment terms
#'
#' Combines an UpSet plot with row annotations for intersection
#' size and optional textboxes of top enriched functional terms.
#'
#' The binary combination pattern and displayed intersection label are
#' printed when the function is called. These mappings can be used as keys
#' for `intersection_labels` and `intersection_colours`.
#'
#' @param set_list A named list of character vectors of feature
#'   sets (e.g., `list(Exp1_Up = ..., Exp1_Down = ..., Exp2_Up = ...,
#'   Exp2_Down = ...)`). Empty input sets and zero-size intersections are
#'   excluded from the plot.
#' @param enrich_list A named enrichment list for intersections, such as the
#'   output of `enrichGO_list()$all` or `enrichR_list()`. Each named element
#'   should be a data.frame with `Cluster`, `Description`, and the column
#'   named by `enrich_rank_by`. If NULL
#'   (default), no enrichment terms are displayed.
#'   The `Cluster` column may contain binary combination names or displayed
#'   UpSet labels. Displayed labels are auto-generated from `names(set_list)`
#'   by joining participating set names with `_`, and may be overridden with
#'   `intersection_labels`.
#'   See `plot_modules_h()` for the same per-category enrichment parameter
#'   convention.
#'   When `enrich_category` is a vector, each category
#'   is shown as a separate column in the annotation textbox.
#' @param enrich_category Category name(s) of enrichment to display: a single
#'   string (e.g. `"BP"`), or a vector for multiple categories (e.g.
#'   `c("BP", "CC")`). Categories must be present in `enrich_list`.
#' @param enrich_top_n Top N enrichment terms per intersection. A single
#'   integer (default: 3) is recycled for all categories in
#'   `enrich_category`. Supply an integer vector of matching length for
#'   per-category control.
#' @param enrich_rank_by Numeric column used to rank enrichment terms within
#'   each intersection. A single string (default: `"p.adjust"`) is recycled
#'   for all categories in `enrich_category`. Supply a character vector of
#'   matching length for per-category control. Also used for colour-coding:
#'   terms with
#'   values below `enrich_threshold` are shown in black, terms above in grey.
#' @param enrich_threshold Threshold for colour-coding enrichment terms
#'   (default: 0.05). Terms where `enrich_rank_by` column value falls
#'   below this threshold are drawn in black, terms above in grey.
#'   Supply a vector matching `enrich_category` for per-category control,
#'   or set to `NULL` to skip colour-coding (all terms black).
#' @param enrich_wrap_width Optional maximum enrichment term width in cm.
#'   The default `NULL` disables wrapping.
#' @param set_order Optional order of the sets. When
#'   `NULL` (default), sets are ordered by decreasing size. Pass a character
#'   vector of set names to control display order.
#' @param intersection_order Optional order of intersection rows.
#'   When `NULL` (default), intersections are ordered by decreasing
#'   combination size. Pass a character vector of combination names
#'   (as returned by `ComplexHeatmap::comb_name()`), or a numeric
#'   permutation of their indices, to specify a custom order.
#' @param intersection_labels Optional named character vector mapping
#'   `comb_name` binary patterns to custom intersection names,
#'   overriding the auto-generated defaults. E.g.,
#'   `c("1010" = "Concordant_Up", "0101" = "Concordant_Down")`.
#'   Entries for patterns not present in the combination matrix cause an
#'   error.
#' @param intersection_colours Optional named colour vector for intersections.
#'   Names must match the binary combination patterns used by
#'   `intersection_labels` (e.g. `"1010"`). Colours are used for the
#'   intersection bars and combination marks.
#' @param fontsize Base font size (default: 8).
#' @param ... Additional arguments passed to `ComplexHeatmap::UpSet()`.
#'
#' @returns A ComplexHeatmap UpSet plot, with optional enrichment
#'   term annotations.
#' @export
#' @examples
#' # Mock feature sets with reproducible intersections
#' genes <- paste0("G", 1:200)
#' set_list <- list(Exp1_Up = genes[1:60], Exp1_Down = genes[61:100],
#'     Exp2_Up = genes[c(1:40, 101:115)], Exp2_Down = genes[c(61:80, 116:130)])
#' enrich_list <- list(BP = data.frame(
#'     Cluster = c("Exp1_Up_Exp2_Up", "Exp1_Down_Exp2_Down"),
#'     Description = c("term A", "term B"), p.adjust = c(0.001, 0.04)))
#' intersection_colours <- c("1010" = "#0072B5", "0101" = "#D55E00")
#' plot_upset_h(set_list, enrich_list = enrich_list,
#'     enrich_category = "BP", intersection_colours = intersection_colours)
#' # Custom intersection order uses the binary combination patterns.
#' plot_upset_h(set_list,
#'     intersection_order = c("1010", "1000", "0101", "0100", "0010", "0001"))
plot_upset_h <- function(
    set_list,
    enrich_list = NULL,
    enrich_category = NULL,
    enrich_top_n = 3,
    enrich_rank_by = "p.adjust",
    enrich_threshold = 0.05,
    enrich_wrap_width = NULL,
    set_order = NULL,
    intersection_order = NULL,
    intersection_labels = NULL,
    intersection_colours = NULL,
    fontsize = 8,
    ...
) {
    .check_list(set_list, "set_list")
    if (is.null(names(set_list)) || anyNA(names(set_list)) ||
        any(names(set_list) == "")) {
        stop("'set_list' must be a named list.")
    }
    if (anyDuplicated(names(set_list))) {
        stop("'set_list' names must be unique.")
    }
    if (!all(vapply(set_list, is.character, logical(1)))) {
        stop("'set_list' must contain character vectors.")
    }
    if (any(vapply(set_list, function(x) {
        anyNA(x) || any(x == "") || anyDuplicated(x) > 0
    }, logical(1)))) {
        stop("Each set must contain non-missing, non-empty, unique feature ",
            "identifiers.")
    }
    .check_positive(fontsize, "fontsize")
    fontfamily <- "sans"
    if (!is.null(set_order)) .check_character(set_order, "set_order")
    if (!is.null(intersection_labels)) {
        if (!is.character(intersection_labels) ||
            is.null(names(intersection_labels)) ||
            anyNA(names(intersection_labels)) ||
            any(names(intersection_labels) == "") ||
            anyNA(intersection_labels) || any(intersection_labels == "")) {
            stop("'intersection_labels' must be a named character vector.")
        }
        if (anyDuplicated(names(intersection_labels))) {
            stop("'intersection_labels' names must be unique.")
        }
    }
    if (!is.null(intersection_colours)) {
        if (!is.character(intersection_colours) ||
            is.null(names(intersection_colours)) ||
            anyNA(names(intersection_colours)) ||
            any(names(intersection_colours) == "")) {
            stop("'intersection_colours' must be a named colour vector.")
        }
        if (anyDuplicated(names(intersection_colours))) {
            stop("'intersection_colours' names must be unique.")
        }
        invisible(grDevices::col2rgb(unname(intersection_colours)))
    }

    # Remove empty sets (ComplexHeatmap errors on zero-size sets)
    set_sizes <- lengths(set_list)
    empty_sets <- names(set_list)[set_sizes == 0]
    if (length(empty_sets) > 0) {
        warning("Removing empty set(s): ", paste(empty_sets, collapse = ", "))
        set_list <- set_list[set_sizes > 0]
    }
    if (length(set_list) < 1) {
        stop("Need at least 1 non-empty set for an UpSet plot.")
    }

    # ---- Build combination matrix ----
    # (intersections on rows for right_annotation)
    m <- ComplexHeatmap::make_comb_mat(set_list, set_on_rows = FALSE)
    comb_names <- ComplexHeatmap::comb_name(m)

    # Auto-generate human-readable intersection names from set names
    # by parsing comb_name binary patterns (e.g. "1010" with set names
    # c("Exp1_Up","Exp1_Down","Exp2_Up","Exp2_Down") -> "Exp1_Up_Exp2_Up")
    set_nms <- names(set_list)
    comb_labels <- vapply(comb_names, function(cn) {
        bits <- strsplit(cn, "")[[1]]
        paste(set_nms[bits == "1"], collapse = "_")
    }, character(1))
    # intersection_labels overrides auto-generated names when provided
    if (!is.null(intersection_labels)) {
        unmatched <- setdiff(names(intersection_labels), comb_names)
        if (length(unmatched) > 0L) {
            stop("'intersection_labels' contains unknown combination ",
                "pattern(s): ",
                paste(unmatched, collapse = ", "))
        }
        matched <- intersect(names(intersection_labels), comb_names)
        comb_labels[match(matched, comb_names)] <- intersection_labels[matched]
    }
    if (anyDuplicated(comb_labels)) {
        stop("Displayed intersection labels must be unique.")
    }
    message("Intersection pattern to label mapping:\n",
        paste0("  ", comb_names, ": ", comb_labels, collapse = "\n"))

    if (is.null(set_order)) {
        set_order <- names(set_list)[order(
            lengths(set_list), decreasing = TRUE
        )]
    }
    if (!setequal(set_order, names(set_list)) ||
        length(set_order) != length(set_list)) {
        stop("'set_order' must contain each set name exactly once.")
    }

    if (is.null(intersection_order)) {
        intersection_order <- order(
            ComplexHeatmap::comb_size(m), decreasing = TRUE
        )
    } else if (is.character(intersection_order)) {
        if (!setequal(intersection_order, comb_names) ||
            length(intersection_order) != length(comb_names)) {
            stop("Character 'intersection_order' must contain each ",
                "combination ",
                "name exactly once.")
        }
        intersection_order <- match(intersection_order, comb_names)
    } else if (is.numeric(intersection_order)) {
        if (length(intersection_order) != length(comb_names) ||
            !setequal(intersection_order, seq_along(comb_names))) {
            stop("Numeric 'intersection_order' must be a permutation of ",
                "combination indices.")
        }
    } else {
        stop("'intersection_order' must be NULL, numeric, or character.")
    }

    upset_bar_col <- "#8C8C8C"
    comb_scaffold_col <- "#666666"
    comb_col <- rep(comb_scaffold_col, length(comb_names))
    names(comb_col) <- comb_names
    if (!is.null(intersection_colours)) {
        unknown <- setdiff(names(intersection_colours), comb_names)
        if (length(unknown) > 0L) {
            stop("'intersection_colours' contains unknown combination ",
                "pattern(s): ",
                paste(unknown, collapse = ", "))
        }
        comb_col[names(intersection_colours)] <- unname(intersection_colours)
    }

    extra_upset_args <- list(...)
    if (length(extra_upset_args) > 0L) {
        reserved_args <- c(
            "comb_col", "set_order", "comb_order", "row_names_side",
            "show_column_names", "row_names_gp",
            "top_annotation", "right_annotation"
        )
        conflicts <- intersect(names(extra_upset_args), reserved_args)
        if (length(conflicts) > 0L) {
            stop("Reserved UpSet argument(s) cannot be overridden: ",
                paste(conflicts, collapse = ", "))
        }
    }
    top_ann <- ComplexHeatmap::upset_top_annotation(
        m,
        gp = grid::gpar(fill = upset_bar_col),
        annotation_name_gp = grid::gpar(
            fontsize = fontsize, fontface = "bold",
            fontfamily = fontfamily
        ),
        axis_param = list(gp = grid::gpar(
            fontsize = fontsize, fontfamily = fontfamily
        ))
    )
    right_ann <- ComplexHeatmap::upset_right_annotation(
        m,
        gp = grid::gpar(fill = upset_bar_col),
        annotation_name_side = "bottom",
        annotation_name_gp = grid::gpar(
            fontsize = fontsize, fontface = "bold",
            fontfamily = fontfamily
        ),
        axis_param = list(gp = grid::gpar(
            fontsize = fontsize, fontfamily = fontfamily
        ))
    )
    enrichment_requested <- !is.null(enrich_list) &&
        !is.null(enrich_category) && length(enrich_category) > 0L
    has_enrichment_terms <- FALSE
    enrich_legend_colours <- stats::setNames(
        vector("list", length(enrich_category)), enrich_category
    )
    if (enrichment_requested) {
        if (!is.null(enrich_wrap_width) &&
            (!is.numeric(enrich_wrap_width) ||
                length(enrich_wrap_width) != 1L ||
                !is.finite(enrich_wrap_width) || enrich_wrap_width <= 0)) {
            stop("'enrich_wrap_width' must be one positive number or NULL.")
        }
        .check_character(enrich_category, "enrich_category")
        if (anyNA(enrich_category) || any(enrich_category == "")) {
            stop("'enrich_category' must not contain missing or empty strings.")
        }
        if (anyDuplicated(enrich_category)) {
            stop("'enrich_category' values must be unique.")
        }
        .check_list(enrich_list, "enrich_list")
        if (is.null(names(enrich_list)) || anyNA(names(enrich_list)) ||
            any(names(enrich_list) == "")) {
            stop("'enrich_list' must be a named list.")
        }
        if (anyDuplicated(names(enrich_list))) {
            stop("'enrich_list' names must be unique.")
        }
        missing_categories <- setdiff(enrich_category, names(enrich_list))
        if (length(missing_categories) > 0L) {
            stop("'enrich_category' contains categories not present in ",
                "'enrich_list': ", paste(missing_categories, collapse = ", "))
        }
        .check_character(enrich_rank_by, "enrich_rank_by")

        # Recycle per-category parameters
        n_categories <- length(enrich_category)
        enrich_rank_by <- .recycle_enrich_param(
            enrich_rank_by, "enrich_rank_by", n_categories)
        if (anyNA(enrich_rank_by) || any(enrich_rank_by == "")) {
            stop("'enrich_rank_by' must not contain missing or empty strings.")
        }
        enrich_top_n <- .recycle_enrich_param(
            enrich_top_n, "enrich_top_n", n_categories)
        if (!is.numeric(enrich_top_n) || any(!is.finite(enrich_top_n))) {
            stop("'enrich_top_n' must contain positive integers.")
        }
        for (v in enrich_top_n) .check_positive_int(v, "enrich_top_n")
        if (!is.null(enrich_threshold)) {
            enrich_threshold <- .recycle_enrich_param(
                enrich_threshold, "enrich_threshold", n_categories)
            for (v in enrich_threshold[!is.na(enrich_threshold)]) {
                .check_pval(v, "enrich_threshold")
            }
        } else {
            enrich_threshold <- rep(NA_real_, n_categories)
        }

        # Build cat_terms: per-category -> per-intersection term data.frames
        cat_terms <- list()
        for (ci in seq_along(enrich_category)) {
            cate <- enrich_category[ci]
            if (is.null(enrich_list[[cate]])) {
                message("enrich_list[['", cate, "']] is NULL. Skipping.")
                next
            }
            enr_df <- as.data.frame(enrich_list[[cate]])
            rank_by <- enrich_rank_by[ci]
            pt <- enrich_threshold[ci]
            top_n <- enrich_top_n[ci]
            required <- c("Cluster", "Description", rank_by)
            missing_cols <- setdiff(required, colnames(enr_df))
            if (length(missing_cols) > 0L) {
                stop("enrich_list[['", cate,
                    "']] is missing required column(s): ",
                    paste(missing_cols, collapse = ", "))
            }
            if (!is.numeric(enr_df[[rank_by]])) {
                stop("enrich_list[['", cate, "']][['", rank_by,
                    "']] must be numeric.")
            }

            # Filter to intersections present in this UpSet
            cluster_values <- as.character(enr_df$Cluster)
            original_idx <- match(cluster_values, comb_names)
            display_idx <- match(cluster_values, comb_labels)
            matched_idx <- original_idx
            use_display <- is.na(matched_idx)
            matched_idx[use_display] <- display_idx[use_display]
            matched <- !is.na(matched_idx)
            unmatched <- unique(cluster_values[!matched])
            if (length(unmatched) > 0L) {
                message("In enrich_list[['", cate, "']], ",
                    length(unmatched), " Cluster value(s) were not matched: ",
                    paste(unmatched, collapse = ", "))
            }
            if (!any(matched)) {
                message("No Cluster values in enrich_list[['", cate,
                    "']] matched an UpSet intersection.")
                next
            }
            enr_df$Cluster <- cluster_values
            enr_df$Cluster[matched] <- comb_labels[matched_idx[matched]]
            enr_df <- enr_df[matched, , drop = FALSE]
            enr_df <- enr_df[is.finite(enr_df[[rank_by]]), , drop = FALSE]
            if (nrow(enr_df) == 0) next

            # Per-intersection top N terms
            term_list <- lapply(comb_labels, function(lbl) {
                tmp <- enr_df[enr_df$Cluster == lbl, , drop = FALSE]
                if (nrow(tmp) == 0) {
                    return(data.frame(
                        text = character(0), col = character(0),
                        fontsize = numeric(0),
                        fontfamily = character(0),
                        stringsAsFactors = FALSE
                    ))
                }
                tmp <- tmp[order(tmp[[rank_by]]), , drop = FALSE]
                tmp <- utils::head(tmp, top_n)
                cols <- if (is.na(pt)) {
                    rep("black", nrow(tmp))
                } else {
                    ifelse(
                        !is.na(tmp[[rank_by]]) &
                        tmp[[rank_by]] < pt,
                        "black", "grey60"
                    )
                }
                data.frame(
                    text = as.character(tmp$Description),
                    col = cols,
                    fontsize = fontsize,
                    fontfamily = fontfamily,
                    stringsAsFactors = FALSE
                )
            })
            names(term_list) <- comb_labels
            cat_terms[[cate]] <- term_list
            enrich_legend_colours[[cate]] <- intersect(
                c("black", "grey60"), unique(unlist(
                    lapply(term_list, function(x) x$col),
                    use.names = FALSE
                ))
            )
        }

        if (length(cat_terms) > 0) {
            # Reformat for multi-column textbox
            text <- .reformat_cat_terms(cat_terms)

            # text is keyed by comb_labels (e.g. "Concordant_Up");
            # remap to comb_names for alignment with combination matrix rows
            text_for_anno <- text[comb_labels]
            names(text_for_anno) <- comb_names

            # Include only intersections with >= 1 enrichment term.
            has_terms <- vapply(text_for_anno, function(x) {
                any(vapply(x, function(df) {
                    is.data.frame(df) && nrow(df) > 0
                }, logical(1)))
            }, logical(1))
            if (any(has_terms)) {
                has_enrichment_terms <- TRUE
                text_for_anno <- text_for_anno[has_terms]

                align_to <- factor(comb_names, levels = comb_names)
                textbox_args <- list(
                    word_wrap = !is.null(enrich_wrap_width),
                    max_width = if (is.null(enrich_wrap_width)) {
                        NULL
                    } else {
                        grid::unit(enrich_wrap_width, "cm")
                    },
                    line_space = grid::unit(0.15 * fontsize, "pt"),
                    term_space = grid::unit(0.75 * fontsize, "pt"),
                    text_padding = grid::unit(0.65 * fontsize, "pt")
                )
                enrich_list_anno <- do.call(
                    .anno_multicol_textbox,
                    c(
                        list(
                            align_to = align_to,
                            text = text_for_anno,
                            background_gp = grid::gpar(
                                fill = "grey98", col = "grey50"
                            ),
                            header_gp = grid::gpar(
                                fontsize = fontsize, fontface = "bold",
                                fontfamily = fontfamily
                            ),
                            by = "anno_link", side = "right"
                        ),
                        if (is.null(enrich_wrap_width)) {
                            list()
                        } else {
                            textbox_args
                        }
                    )
                )
                enrich_ra <- ComplexHeatmap::rowAnnotation(
                    Enrichment = enrich_list_anno)
                right_ann <- c(right_ann, enrich_ra)
            }
        }
    }
    if (enrichment_requested && !has_enrichment_terms) {
        message("No matching enrichment terms, ",
            "returning the UpSet plot only.")
    }

    # ---- Build UpSet ----
    p <- ComplexHeatmap::UpSet(m,
        comb_col = unname(comb_col),
        top_annotation = top_ann,
        right_annotation = right_ann,
        set_order = set_order,
        comb_order = intersection_order,
        row_names_side = "left",
        show_column_names = TRUE,
        column_names_gp = grid::gpar(
            fontsize = fontsize, fontfamily = fontfamily
        ),
        row_names_gp = grid::gpar(
            fontsize = fontsize, fontfamily = fontfamily
        ),
        ...
    )

    # ---- Enrichment colour legend ----
    lgd_list <- list()
    if (enrichment_requested && has_enrichment_terms) {
        for (ci in seq_along(enrich_category)) {
            cate <- enrich_category[ci]
            pt <- enrich_threshold[ci]
            if (is.null(enrich_list[[cate]])) next
            if (is.na(pt)) next
            legend_colours <- enrich_legend_colours[[cate]]
            if (length(legend_colours) == 0L) next
            go_lgd <- ComplexHeatmap::Legend(
                title = cate,
                labels = c(
                    black = sprintf("%s < %.2f", enrich_rank_by[ci], pt),
                    grey60 = sprintf("%s >= %.2f", enrich_rank_by[ci], pt)
                )[legend_colours],
                type = "lines",
                legend_gp = grid::gpar(col = legend_colours, lwd = 2),
                labels_gp = grid::gpar(
                    fontsize = fontsize, fontfamily = fontfamily
                ),
                title_gp = grid::gpar(
                    fontsize = fontsize, fontfamily = fontfamily
                )
            )
            lgd_list <- c(lgd_list, list(go_lgd))
        }
    }
    if (length(lgd_list) > 0) {
        ComplexHeatmap::draw(p,
            heatmap_legend_list = lgd_list,
            heatmap_legend_side = "right",
            merge_legends = TRUE)
    } else {
        ComplexHeatmap::draw(p)
    }
    invisible(p)
}
