#' UpSet plot of feature sets across datasets with functional enrichment
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
#' @param enrich_list A named enrichment list for modules,
#'   as produced by `enrichGO_list()$all`. Each named element
#'   should be a data.frame with columns `Cluster`, `Description`,
#'   and the rank column. Accepts any enrichment result with compatible format.
#'   The `Cluster` column must contain intersection names matching the
#'   displayed UpSet labels, which are auto-generated from `names(set_list)`
#'   by joining participating set names with `_`.
#'   See `plot_modules_h()` for the same per-category enrichment parameter
#'   convention.
#'   When `enrich_category` is a vector, each category
#'   is shown as a separate column in the annotation textbox.
#'   If NULL (default), no enrichment terms will be displayed.
#' @param enrich_category Category name(s) of enrichment to display:
#'   a single string (e.g. `"BP"`), or a vector for multiple
#'   categories (e.g. `c("BP", "CC")`).
#' @param enrich_rank_by Numeric column used to rank enrichment terms within
#'   each module. A single string (default: `"p.adjust"`) is recycled for all
#'   categories in `enrich_category`. Supply a character vector of matching
#'   length for per-category control. Also used for colour-coding: terms with
#'   values below `enrich_threshold` are shown in black, terms above in grey.
#' @param enrich_top_n Top N enrichment terms per module. A single integer
#'   (default: 3) is recycled for all categories in `enrich_category`.
#'   Supply an integer vector of matching length for per-category control.
#' @param enrich_threshold Threshold for colour-coding enrichment terms
#'   (default: 0.05). Terms where `enrich_rank_by` column value falls
#'   below this threshold are drawn in black, terms above in grey.
#'   Supply a vector matching `enrich_category` for per-category control,
#'   or set to `NULL` to skip colour-coding (all terms black).
#' @param set_order Optional ordering of the set rows (bottom axis). When
#'   `NULL` (default), sets are ordered by decreasing size. Pass a character
#'   vector of set names to control display order.
#' @param intersection_order Optional ordering of the intersection rows
#'   (top axis).
#'   When `NULL` (default), intersections are ordered by decreasing
#'   combination size. Pass a character vector of combination names
#'   (as returned by `ComplexHeatmap::comb_name()`) for custom order.
#' @param intersection_labels Optional named character vector mapping
#'   `comb_name` binary patterns to custom intersection names,
#'   overriding the auto-generated defaults. E.g.,
#'   `c("1010" = "Concordant_Up", "0101" = "Concordant_Down")`.
#'   Entries for patterns not present in the combination matrix cause an
#'   error.
#' @param intersection_colours Optional named colour vector for intersections.
#'   Names must match the binary combination patterns used by
#'   `intersection_labels` (e.g. `"1010"`). Colours are used for the
#'   intersection bars and combination marks. Unknown patterns cause an error.
#' @param fontsize Base font size (default: 8).
#' @param ... Additional arguments passed to `ComplexHeatmap::UpSet()`.
#'
#' @returns A ComplexHeatmap Upset plot with optional enrichment annotation.
#' @export
#'
#' @examples
#' # Mock feature sets with reproducible intersections
#' genes <- paste0("G", 1:200)
#' set_list <- list(Exp1_Up = genes[1:60], Exp1_Down = genes[61:100],
#'     Exp2_Up = genes[c(1:40, 101:115)], Exp2_Down = genes[c(61:80, 116:130)])
#' enrich_list <- list(BP = data.frame(
#'     Cluster = c("Exp1_Up_Exp2_Up", "Exp1_Down_Exp2_Down"),
#'     Description = c("term A", "term B"), p.adjust = c(0.001, 0.04)))
#' intersection_colours <- c("1010" = "#0072B5", "0101" = "#D55E00")
#' plot_upset(set_list, enrich_list = enrich_list,
#'     enrich_category = "BP", intersection_colours = intersection_colours)
#' # Custom intersection order uses the binary combination patterns.
#' plot_upset(set_list,
#'     intersection_order = c("1010", "1000", "0101", "0100", "0010", "0001"))
plot_upset <- function(
    set_list,
    enrich_list = NULL,
    enrich_category = NULL,
    enrich_top_n = 3,
    enrich_rank_by = "p.adjust",
    enrich_threshold = 0.05,
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
    if (!is.null(enrich_category)) {
        .check_character(enrich_category, "enrich_category")
        if (anyNA(enrich_category) || any(enrich_category == "")) {
            stop("'enrich_category' must not contain missing or empty strings.")
        }
    }
    if (!is.null(set_order)) .check_character(set_order, "set_order")
    if (!is.null(intersection_order) &&
        !is.numeric(intersection_order) &&
        !is.character(intersection_order)) {
        stop("'intersection_order' must be numeric or character.")
    }
    if (!is.null(intersection_labels)) {
        if (!is.character(intersection_labels) ||
            is.null(names(intersection_labels)) ||
            anyNA(names(intersection_labels)) ||
            any(names(intersection_labels) == "") ||
            anyNA(intersection_labels) ||
            any(intersection_labels == "")) {
            stop("'intersection_labels' must be a named character vector",
                " with non-empty names and values.")
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
        # Convert colours to RGB to validate names, hex codes, and R colours.
        invisible(grDevices::col2rgb(unname(intersection_colours)))
    }

    # ---- Remove empty sets (ComplexHeatmap errors on zero-size sets) ----
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
    row_labels <- vapply(comb_names, function(cn) {
        bits <- strsplit(cn, "")[[1]]
        paste(set_nms[bits == "1"], collapse = "_")
    }, character(1))
    # intersection_labels overrides auto-generated names when provided
    if (!is.null(intersection_labels)) {
        overrides <- intersect(names(intersection_labels), comb_names)
        unmatched <- setdiff(names(intersection_labels), comb_names)
        if (length(unmatched) > 0L) {
            stop(
                "'intersection_labels' contains unknown combination ",
                "pattern(s): ",
                paste(unmatched, collapse = ", "))
        }
        if (length(overrides) > 0) {
            row_labels[match(overrides, comb_names)] <-
                intersection_labels[overrides]
        }
    }
    if (anyDuplicated(row_labels)) {
        stop("Intersection labels must be unique.")
    }
    message("Intersection pattern to label mapping:\n",
        paste0("  ", comb_names, ": ", row_labels, collapse = "\n"))

    # Enrichment textbox (if provided)
    if (!is.null(enrich_list) && !is.null(enrich_category)) {
        .check_list(enrich_list, "enrich_list")
        .check_character(enrich_rank_by, "enrich_rank_by")

        # Recycle per-category parameters
        n_categories <- length(enrich_category)
        enrich_rank_by <- .recycle_enrich_param(
            enrich_rank_by, "enrich_rank_by", n_categories)
        if (!is.numeric(enrich_top_n) ||
            any(!is.finite(enrich_top_n))) {
            stop("'enrich_top_n' must contain positive integers.")
        }
        enrich_top_n <- .recycle_enrich_param(
            enrich_top_n, "enrich_top_n", n_categories)
        for (v in enrich_top_n) {
            .check_positive_int(v, "enrich_top_n")
        }
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
            if (!cate %in% names(enrich_list)) {
                warning("enrich_category '", cate,
                    "' not found in enrich_list. Skipping.")
                next
            }
            if (is.null(enrich_list[[cate]])) next
            enr_df <- as.data.frame(enrich_list[[cate]])
            if (!"Cluster" %in% colnames(enr_df)) next
            rank_by <- enrich_rank_by[ci]
            pt <- enrich_threshold[ci]
            top_n <- enrich_top_n[ci]
            if (!rank_by %in% colnames(enr_df)) {
                warning("enrich_rank_by column '", rank_by,
                    "' not found in enrich_list '",
                    cate, "'. Skipping.")
                next
            }

            # Filter to intersections present in this UpSet
            enr_clusters <- unique(enr_df$Cluster)
            matched <- enr_clusters %in% row_labels
            if (any(!matched)) {
                warning("In enrich_list[['", cate, "']], ",
                    sum(!matched), " Cluster value(s) not found in ",
                    "UpSet intersections: ",
                    paste(enr_clusters[!matched], collapse = ", "),
                    ". Available intersections: ",
                    paste(row_labels, collapse = ", "))
            }
            enr_df <- enr_df[enr_df$Cluster %in% row_labels, , drop = FALSE]
            if (nrow(enr_df) == 0) next

            # Per-intersection top N terms
            term_list <- lapply(row_labels, function(lbl) {
                tmp <- enr_df[enr_df$Cluster == lbl, , drop = FALSE]
                if (nrow(tmp) == 0) {
                    return(data.frame(
                        text = character(0), col = character(0),
                        fontsize = numeric(0), stringsAsFactors = FALSE
                    ))
                }
                tmp <- tmp[order(tmp[[rank_by]]), , drop = FALSE]
                tmp <- utils::head(tmp, top_n)
                desc_col <- if ("Description" %in% colnames(tmp)) {
                    "Description"
                } else {
                    colnames(tmp)[1]
                }
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
                    text = as.character(tmp[[desc_col]]),
                    col = cols,
                    fontsize = fontsize,
                    stringsAsFactors = FALSE
                )
            })
            names(term_list) <- row_labels
            cat_terms[[cate]] <- term_list
        }

        if (length(cat_terms) > 0) {
            # Reformat for multi-column textbox
            text <- .reformat_cat_terms(cat_terms)

            # text is keyed by row_labels (e.g. "Concordant_Up");
            # remap to comb_names for alignment with combination matrix rows
            text_for_anno <- text[row_labels]
            names(text_for_anno) <- comb_names

            # Include only intersections with >= 1 enrichment term.
            has_terms <- vapply(text_for_anno, function(x) {
                any(vapply(x, function(df) {
                    is.data.frame(df) && nrow(df) > 0
                }, logical(1)))
            }, logical(1))
            text_for_anno <- text_for_anno[has_terms]

            align_to <- factor(comb_names, levels = comb_names)
            enrich_list_anno <- .anno_multicol_textbox(
                align_to = align_to,
                text = text_for_anno,
                background_gp = grid::gpar(fill = "grey95", col = "grey50"),
                by = "anno_link", side = "right"
            )
            enrich_ra <- ComplexHeatmap::rowAnnotation(
                Enrichment = enrich_list_anno)
            right_ann <- c(
                ComplexHeatmap::upset_right_annotation(m,
                    gp = grid::gpar(fill = "grey30", col = NA),
                    annotation_name_side = "bottom",
                    axis_param = list(gp = grid::gpar(fontsize = fontsize))),
                enrich_ra
            )
        } else {
            right_ann <- ComplexHeatmap::upset_right_annotation(m,
                gp = grid::gpar(fill = "grey30", col = NA),
                annotation_name_side = "bottom",
                axis_param = list(gp = grid::gpar(fontsize = fontsize)))
        }
    } else {
        right_ann <- ComplexHeatmap::upset_right_annotation(m,
            gp = grid::gpar(fill = "grey30", col = NA),
            annotation_name_side = "bottom",
            axis_param = list(gp = grid::gpar(fontsize = fontsize)))
    }

    # ---- Mode colouring ----
    comb_col <- "black"
    if (!is.null(intersection_colours)) {
        unknown_colours <- setdiff(names(intersection_colours), comb_names)
        if (length(unknown_colours) > 0L) {
            stop(
                "'intersection_colours' contains unknown combination ",
                "pattern(s): ",
                paste(unknown_colours, collapse = ", "))
        }
        # Map colours to the final intersection labels.
        named_colours <- setNames(
            intersection_colours[comb_names], comb_names)
        # Override default black for intersections with colors
        if (any(!is.na(named_colours))) {
            comb_col <- named_colours
            comb_col[is.na(comb_col)] <- "black"
        }
    }

    # ---- Build UpSet ----
    if (is.null(set_order)) {
        set_order <- names(set_list)[order(
            lengths(set_list), decreasing = TRUE)]
    }
    if (!is.character(set_order) ||
        !setequal(set_order, names(set_list)) ||
        length(set_order) != length(set_list)) {
        stop("'set_order' must contain each set name exactly once.")
    }
    if (is.null(intersection_order)) {
        intersection_order <- order(
            ComplexHeatmap::comb_size(m), decreasing = TRUE)
    }
    if (is.character(intersection_order)) {
        if (!setequal(intersection_order, comb_names) ||
            length(intersection_order) != length(comb_names)) {
            stop("Character 'intersection_order' must contain each ",
                "combination name exactly once.")
        }
        intersection_order <- match(intersection_order, comb_names)
    } else if (is.numeric(intersection_order)) {
        if (length(intersection_order) != length(comb_names) ||
            !setequal(intersection_order, seq_along(comb_names))) {
            stop("Numeric 'intersection_order' must be a permutation of ",
                "the combination indices.")
        }
    } else {
        stop("'intersection_order' must be NULL, numeric, or character.")
    }
    p <- ComplexHeatmap::UpSet(m,
        comb_col = comb_col,
        right_annotation = right_ann,
        set_order = set_order,
        comb_order = intersection_order,
        row_names_side = "left",
        ...
    )

    # ---- Enrichment colour legend ----
    lgd_list <- list()
    if (!is.null(enrich_list) && !is.null(enrich_category)) {
        for (ci in seq_along(enrich_category)) {
            cate <- enrich_category[ci]
            pt <- enrich_threshold[ci]
            if (is.null(enrich_list[[cate]])) next
            if (is.na(pt)) next
            go_lgd <- ComplexHeatmap::Legend(
                title = cate,
                labels = c(
                    sprintf("%s < %.2f", enrich_rank_by[ci], pt),
                    sprintf("%s >= %.2f", enrich_rank_by[ci], pt)
                ),
                type = "lines",
                legend_gp = grid::gpar(
                    col = c("black", "grey60"), lwd = 2),
                labels_gp = grid::gpar(fontsize = fontsize),
                title_gp = grid::gpar(fontsize = fontsize)
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
