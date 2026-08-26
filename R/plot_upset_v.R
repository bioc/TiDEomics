#' UpSet plot with vertically aligned enrichment dot matrix
#'
#' Draws a ComplexHeatmap UpSet plot with intersections as columns and a
#' clusterProfiler-style enrichment dot matrix underneath. The two panels use
#' the same intersection order, so each enrichment column is vertically
#' aligned with its UpSet intersection.
#'
#' @param set_list A named list of character vectors of feature sets.
#' @param enrich_list A named enrichment list for intersections, such as the
#'   output of `enrichGO_list()$all` or `enrichR_list()`. Each named element
#'   should be a data.frame with `Cluster`, `Description`, and the columns
#'   named by `enrich_rank_by`, `dot_size_col`, and `dot_colour_col`. If NULL
#'   (default), no enrichment terms are displayed.
#' @param enrich_category Category name(s) of enrichment to display: a single
#'   string (e.g. `"BP"`), or a vector for multiple categories (e.g.
#'   `c("BP", "CC")`). Categories must be present in `enrich_list`.
#' @param enrich_top_n Top N enrichment terms per intersection. A single
#'   integer (default: 3) is recycled for all categories in
#'   `enrich_category`. Supply an integer vector of matching length for
#'   per-category control.
#' @param include_all Whether to include every occurrence of a term selected
#'   in any intersection, matching clusterProfiler's `includeAll`.
#'   (default: `TRUE`).
#' @param enrich_rank_by Numeric column used to rank enrichment terms within
#'   each intersection. A single string (default: `"p.adjust"`) is recycled
#'   for all categories in `enrich_category`. Supply a character vector of
#'   matching length for per-category control.
#' @param set_order Optional order of set rows.
#' @param intersection_order Optional order of intersections, supplied as
#'   combination names or a permutation of combination indices.
#' @param intersection_labels Optional named mapping from combination patterns
#'   to displayed intersection labels. Enrichment `Cluster` values may use
#'   either the original combination patterns or the displayed labels.
#' @param intersection_colours Optional named colours keyed by combination
#'   patterns. Unspecified intersections remain neutral.
#' @param dot_size_col Numeric column used for dot size in the enrichment
#'   annotation. Default is `"Count"`.
#' @param dot_colour_col Numeric column used for dot colour in the enrichment
#'   annotation. Default is `"p.adjust"`.
#' @param fontsize Base font size. Default is 8.
#' @param dot_min_size Minimum dot radius in mm. Default is 1.
#' @param dot_max_size Maximum dot radius in mm. Default is 3.
#' @param dot_outline Outline colour for enrichment dots. Use `NA` to remove
#'   outlines. Default is `"black"`.
#' @param dot_outline_width Outline width in points. Default is 0.5.
#' @param term_wrap_width Maximum characters per displayed enrichment term
#'   before wrapping to a new line. Default is 30.
#' @param row_gap_mm Gap between enrichment categories in mm. Default is 2.
#' @param show_dot_column_names Whether to show intersection labels on the
#'   enrichment dot matrix. Default: `FALSE`.
#' @param ... Additional arguments passed to `ComplexHeatmap::UpSet()`.
#'
#' @returns A ComplexHeatmap UpSet plot, with optional enrichment
#'   dot annotations.
#' @export
#' @examples
#' genes <- paste0("G", 1:200)
#' set_list <- list(
#'     Exp1_Up = genes[1:60],
#'     Exp1_Down = genes[61:100],
#'     Exp2_Up = genes[c(1:40, 101:115)],
#'     Exp2_Down = genes[c(61:80, 116:130)]
#' )
#' enrich_list <- list(BP = data.frame(
#'     Cluster = c("Exp1_Up_Exp2_Up", "Exp1_Down_Exp2_Down"),
#'     Description = c("term A", "term B"), Count = c(20, 10),
#'     p.adjust = c(0.001, 0.02), stringsAsFactors = FALSE))
#' plot_upset_v(set_list, enrich_list = enrich_list,
#'     enrich_category = "BP")
#' @references Dot matrix mimicking `clusterProfiler::dotplot()`
plot_upset_v <- function(
    set_list,
    enrich_list = NULL,
    enrich_category = NULL,
    enrich_top_n = 3,
    include_all = TRUE,
    enrich_rank_by = "p.adjust",
    set_order = NULL,
    intersection_order = NULL,
    intersection_labels = NULL,
    intersection_colours = NULL,
    dot_size_col = "Count",
    dot_colour_col = "p.adjust",
    fontsize = 8,
    dot_min_size = 1,
    dot_max_size = 3,
    dot_outline = "black",
    dot_outline_width = 0.5,
    term_wrap_width = 30,
    row_gap_mm = 2,
    show_dot_column_names = FALSE,
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
            stop("'intersection_labels' must be a named character ",
                "vector.")
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
            stop("'intersection_colours' must be a named colour ",
                "vector.")
        }
        if (anyDuplicated(names(intersection_colours))) {
            stop("'intersection_colours' names must be unique.")
        }
        invisible(grDevices::col2rgb(unname(intersection_colours)))
    }

    # Remove empty sets before constructing the combination matrix.
    set_sizes <- lengths(set_list)
    empty_sets <- names(set_list)[set_sizes == 0]
    if (length(empty_sets) > 0L) {
        warning("Removing empty set(s): ", paste(empty_sets, collapse = ", "))
        set_list <- set_list[set_sizes > 0]
    }
    if (length(set_list) < 1L) {
        stop("Need at least 1 non-empty set for an UpSet plot.")
    }

    # Intersections are columns in the vertical layout.
    m <- ComplexHeatmap::make_comb_mat(set_list, set_on_rows = TRUE)
    comb_names <- ComplexHeatmap::comb_name(m)
    set_nms <- names(set_list)
    comb_labels <- vapply(comb_names, function(cn) {
        bits <- strsplit(cn, "")[[1]]
        paste(set_nms[bits == "1"], collapse = "_")
    }, character(1))

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

    upset_args <- list(
        m,
        comb_col = unname(comb_col),
        pt_size = grid::unit(2.5, "mm"),
        bg_pt_col = "#D9D9D9",
        set_order = set_order,
        comb_order = intersection_order,
        row_names_side = "left",
        show_column_names = TRUE,
        row_names_gp = grid::gpar(
            fontsize = fontsize, fontfamily = fontfamily
        )
    )
    extra_upset_args <- list(...)
    if (length(extra_upset_args) > 0L) {
        reserved_args <- c(
            "comb_col", "set_order", "comb_order", "row_names_side",
            "show_column_names", "row_names_gp", "top_annotation",
            "right_annotation"
        )
        conflicts <- intersect(names(extra_upset_args), reserved_args)
        if (length(conflicts) > 0L) {
            stop("Reserved UpSet argument(s) cannot be overridden: ",
                paste(conflicts, collapse = ", "))
        }
        upset_args[names(extra_upset_args)] <- extra_upset_args
    }
    build_upset <- function(colours) {
        current_args <- upset_args
        current_args$comb_col <- unname(colours)
        ht <- do.call(ComplexHeatmap::UpSet, current_args)
        ht@top_annotation <- ComplexHeatmap::upset_top_annotation(
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
        ht@right_annotation <- ComplexHeatmap::upset_right_annotation(
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
        ht
    }

    if (is.null(enrich_list) || is.null(enrich_category) ||
        length(enrich_category) == 0L) {
        upset_ht <- build_upset(comb_col)
        return(invisible(upset_ht))
    }

    # validate enrichment parameters
    .check_logical(include_all, "include_all")
    if (length(dot_min_size) != 1L || length(dot_max_size) != 1L) {
        stop("'dot_min_size' and 'dot_max_size' must each contain one value.")
    }
    .check_positive(dot_min_size, "dot_min_size")
    .check_positive(dot_max_size, "dot_max_size")
    if (length(dot_outline) != 1L ||
        (!is.character(dot_outline) &&
            !(is.logical(dot_outline) && is.na(dot_outline)))) {
        stop("'dot_outline' must contain one colour or NA.")
    }
    if (is.logical(dot_outline)) dot_outline <- NA_character_
    if (length(dot_outline_width) != 1L || length(row_gap_mm) != 1L) {
        stop("'dot_outline_width' and 'row_gap_mm' must each contain one ",
            "value.")
    }
    .check_positive(dot_outline_width, "dot_outline_width")
    .check_positive_int(term_wrap_width, "term_wrap_width")
    .check_logical(show_dot_column_names, "show_dot_column_names")
    if (!is.na(dot_outline)) {
        invisible(grDevices::col2rgb(dot_outline))
    }
    if (dot_max_size < dot_min_size) {
        stop("'dot_max_size' must be at least 'dot_min_size'.")
    }
    .check_positive(row_gap_mm, "row_gap_mm")
    .check_character(dot_size_col, "dot_size_col")
    .check_character(dot_colour_col, "dot_colour_col")
    if (length(dot_size_col) != 1L || is.na(dot_size_col) ||
        dot_size_col == "" || length(dot_colour_col) != 1L ||
        is.na(dot_colour_col) || dot_colour_col == "") {
        stop("'dot_size_col' and 'dot_colour_col' must each contain ",
            "one value.")
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
    enrich_rank_by <- .recycle_enrich_param(
        enrich_rank_by, "enrich_rank_by", length(enrich_category)
    )
    if (anyNA(enrich_rank_by) || any(enrich_rank_by == "")) {
        stop("'enrich_rank_by' must not contain missing or empty strings.")
    }
    enrich_top_n <- .recycle_enrich_param(
        enrich_top_n, "enrich_top_n", length(enrich_category)
    )
    if (!is.numeric(enrich_top_n) || any(!is.finite(enrich_top_n))) {
        stop("'enrich_top_n' must contain positive integers.")
    }
    for (v in enrich_top_n) .check_positive_int(v, "enrich_top_n")
    dot_records <- list()
    record_i <- 0L
    for (ci in seq_along(enrich_category)) {
        cate <- enrich_category[ci]
        if (is.null(enrich_list[[cate]])) {
            message("enrich_list[['", cate, "']] is NULL. Skipping.")
            next
        }
        enr_df <- as.data.frame(enrich_list[[cate]])
        required <- c("Cluster", "Description", enrich_rank_by[ci],
            dot_size_col, dot_colour_col)
        missing_cols <- setdiff(required, colnames(enr_df))
        if (length(missing_cols) > 0L) {
            stop("enrich_list[['", cate,
                "']] is missing required column(s): ",
                paste(missing_cols, collapse = ", "))
        }
        if (!is.numeric(enr_df[[enrich_rank_by[ci]]])) {
            stop("enrich_list[['", cate, "']][['",
                enrich_rank_by[ci], "']] must be numeric.")
        }
        if (!is.numeric(enr_df[[dot_size_col]])) {
            stop("enrich_list[['", cate, "']][['", dot_size_col,
                "']] must be numeric.")
        }
        if (!is.numeric(enr_df[[dot_colour_col]])) {
            stop("enrich_list[['", cate, "']][['", dot_colour_col,
                "']] must be numeric.")
        }
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
        if (nrow(enr_df) == 0L) next
        rank_values <- enr_df[[enrich_rank_by[ci]]]
        size_values <- enr_df[[dot_size_col]]
        colour_values <- enr_df[[dot_colour_col]]
        keep <- is.finite(rank_values) & is.finite(size_values) &
            size_values > 0 & is.finite(colour_values) &
            colour_values >= 0
        enr_df <- enr_df[keep, , drop = FALSE]
        if (nrow(enr_df) == 0L) next
        rank_values <- rank_values[keep]
        size_values <- size_values[keep]
        colour_values <- colour_values[keep]
        enr_df$.rank <- rank_values
        enr_df$.size <- size_values
        enr_df$.colour <- colour_values
        term_ids <- as.character(enr_df$Description)
        if ("ID" %in% colnames(enr_df)) {
            id_values <- as.character(enr_df$ID)
            use_description <- is.na(id_values) | id_values == ""
            term_ids[!use_description] <- id_values[!use_description]
        }
        enr_df$.term_id <- term_ids
        enr_df <- enr_df[order(enr_df$Cluster, enr_df$.rank), , drop = FALSE]
        enr_df <- enr_df[!duplicated(
            paste(enr_df$Cluster, enr_df$.term_id, sep = "\r")
        ), , drop = FALSE]
        selected_top <- do.call(rbind, lapply(comb_labels, function(cl) {
            utils::head(enr_df[enr_df$Cluster == cl, , drop = FALSE],
                enrich_top_n[ci])
        }))
        if (is.null(selected_top) || nrow(selected_top) == 0L) next
        if (include_all) {
            selected_ids <- unique(selected_top$.term_id)
            selected <- enr_df[enr_df$.term_id %in% selected_ids,
                , drop = FALSE]
        } else {
            selected <- selected_top
        }
        for (ri in seq_len(nrow(selected))) {
            record_i <- record_i + 1L
            dot_records[[record_i]] <- data.frame(
                category = cate,
                term = as.character(selected$Description[ri]),
                cluster = as.character(selected$Cluster[ri]),
                rank = selected$.rank[ri],
                size = selected$.size[ri],
                colour = selected$.colour[ri],
                stringsAsFactors = FALSE
            )
        }
    }

    if (length(dot_records) == 0L) {
        upset_ht <- build_upset(comb_col)
        message("No matching enrichment terms, ",
            "returning the UpSet plot only.")
        return(invisible(upset_ht))
    }
    dots <- do.call(rbind, dot_records)
    upset_ht <- build_upset(comb_col)
    dots$term_key <- paste(dots$category, dots$term, sep = "\r")
    term_info <- unique(dots[, c("term_key", "category", "term")])
    term_info$category <- factor(term_info$category, levels = enrich_category)
    term_order <- unlist(lapply(enrich_category, function(cate) {
        keys <- term_info$term_key[term_info$category == cate]
        keys[order(vapply(keys, function(key) {
            min(dots$rank[dots$term_key == key], na.rm = TRUE)
        }, numeric(1)))]
    }), use.names = FALSE)
    term_info <- term_info[
        match(term_order, term_info$term_key), , drop = FALSE
    ]

    count_mat <- matrix(
        NA_real_, nrow = length(term_order), ncol = length(comb_names),
        dimnames = list(term_order, comb_names)
    )
    colour_mat <- count_mat
    for (ri in seq_len(nrow(dots))) {
        row_i <- match(dots$term_key[ri], term_order)
        col_i <- match(dots$cluster[ri], comb_labels)
        count_mat[row_i, col_i] <- dots$size[ri]
        colour_mat[row_i, col_i] <- dots$colour[ri]
    }

    max_count <- max(count_mat, na.rm = TRUE)
    observed_sizes <- count_mat[
        is.finite(count_mat) & count_mat > 0
    ]
    size_min <- min(observed_sizes)
    size_to_radius <- function(x) {
        if (max_count == size_min) {
            return(rep(mean(c(dot_min_size, dot_max_size)), length(x)))
        }
        # Match clusterProfiler::dotplot()'s area-based size mapping.
        scaled <- (sqrt(x) - sqrt(size_min)) /
            (sqrt(max_count) - sqrt(size_min))
        scaled <- pmax(0, pmin(1, scaled))
        dot_min_size + scaled * (dot_max_size - dot_min_size)
    }
    q_values <- colour_mat[is.finite(colour_mat)]
    positive_q <- q_values[q_values > 0]
    q_floor <- if (length(positive_q) > 0L) {
        min(positive_q)
    } else {
        .Machine$double.xmin
    }
    q_values <- pmax(q_values, q_floor)
    q_min <- min(q_values)
    q_max <- max(q_values)
    if (q_min == q_max) {
        q_min <- max(.Machine$double.xmin, q_min / 10)
        q_max <- q_max * 10
    }
    log_q_min <- log10(q_min)
    log_q_max <- log10(q_max)
    # Match clusterProfiler::dotplot()'s log10 colour transformation.
    q_col_fun <- circlize::colorRamp2(
        c(log_q_min, log_q_max),
        c("#E06663", "#327EBA")
    )

    term_labels <- term_info$term
    term_labels <- vapply(term_labels, function(x) {
        paste(strwrap(x, width = term_wrap_width), collapse = "\n")
    }, character(1))
    row_split <- droplevels(term_info$category)
    dot_cell <- function(j, i, x, y, w, h, fill) {
        count_value <- count_mat[i, j]
        colour_value <- colour_mat[i, j]
        if (!is.finite(count_value) || count_value <= 0 ||
            !is.finite(colour_value)) {
            return()
        }
        radius <- size_to_radius(count_value)
        grid::grid.circle(
            x = x, y = y, r = grid::unit(radius, "mm"),
            gp = grid::gpar(
                fill = q_col_fun(log10(max(colour_value, q_floor))),
                col = dot_outline,
                lwd = dot_outline_width
            )
        )
    }

    dot_ht <- ComplexHeatmap::Heatmap(
        count_mat,
        name = "Enrichment",
        col = circlize::colorRamp2(c(0, max_count), c("white", "white")),
        rect_gp = grid::gpar(type = "none"),
        cell_fun = dot_cell,
        cluster_rows = FALSE,
        cluster_columns = FALSE,
        row_order = seq_len(nrow(count_mat)),
        column_order = intersection_order,
        row_split = row_split,
        row_gap = grid::unit(row_gap_mm, "mm"),
        row_title_rot = 0,
        row_title_gp = grid::gpar(
            fontsize = fontsize, fontface = "bold",
            fontfamily = fontfamily
        ),
        show_row_names = TRUE,
        row_labels = term_labels,
        row_names_side = "left",
        row_names_gp = grid::gpar(
            fontsize = fontsize, fontfamily = fontfamily
        ),
        row_names_max_width = grid::unit(8, "cm"),
        show_column_names = show_dot_column_names,
        column_labels = comb_labels,
        column_names_rot = 30,
        column_names_gp = grid::gpar(
            fontsize = fontsize, fontfamily = fontfamily
        ),
        show_heatmap_legend = FALSE,
        column_gap = grid::unit(1, "mm")
    )

    size_at <- pretty(c(size_min, max_count), n = 4)
    size_at <- size_at[size_at >= size_min & size_at <= max_count]
    size_at <- sort(unique(c(size_min, size_at, max_count)))
    if (length(size_at) > 4L) {
        internal_at <- size_at[seq.int(2L, length(size_at) - 1L)]
        keep_internal <- unique(round(seq(
            1, length(internal_at), length.out = 2
        )))
        size_at <- sort(unique(c(
            size_min, internal_at[keep_internal], max_count
        )))
    }
    size_legend_graphics <- lapply(size_at, function(count_value) {
        radius <- size_to_radius(count_value)
        function(x, y, w, h) {
            grid::grid.circle(
                x = x,
                y = y,
                r = grid::unit(radius, "mm"),
                gp = grid::gpar(
                    fill = "#4D4D4D",
                    col = dot_outline,
                    lwd = dot_outline_width
                )
            )
        }
    })
    size_lgd <- ComplexHeatmap::Legend(
        title = dot_size_col,
        at = size_at,
        labels = format(size_at, trim = TRUE),
        graphics = size_legend_graphics,
        grid_width = grid::unit(2 * dot_max_size, "mm"),
        grid_height = grid::unit(2 * dot_max_size, "mm"),
        title_gp = grid::gpar(
            fontsize = fontsize, fontfamily = fontfamily
        ),
        labels_gp = grid::gpar(
            fontsize = fontsize, fontfamily = fontfamily
        )
    )
    log_colour_at <- pretty(c(log_q_min, log_q_max), n = 3)
    log_colour_at <- log_colour_at[
        log_colour_at >= log_q_min & log_colour_at <= log_q_max
    ]
    colour_lgd <- ComplexHeatmap::Legend(
        title = dot_colour_col,
        at = log_colour_at,
        labels = format(10^log_colour_at, trim = TRUE),
        col_fun = q_col_fun,
        title_gp = grid::gpar(
            fontsize = fontsize, fontfamily = fontfamily
        ),
        labels_gp = grid::gpar(
            fontsize = fontsize, fontfamily = fontfamily
        )
    )

    legend_width_mm <- max(
        grid::convertWidth(
            grid::grobWidth(size_lgd@grob), "mm", valueOnly = TRUE
        ),
        grid::convertWidth(
            grid::grobWidth(colour_lgd@grob), "mm", valueOnly = TRUE
        )
    ) + 2
    dot_ht@right_annotation <- ComplexHeatmap::rowAnnotation(
        enrichment_legends = ComplexHeatmap::anno_empty(
            which = "row",
            width = grid::unit(legend_width_mm, "mm"),
            border = FALSE,
            show_name = FALSE
        )
    )
    legend_top_gap_mm <- 15
    draw_dot_legends <- function(ht) {
        ComplexHeatmap::decorate_annotation(
            "enrichment_legends",
            {
                ComplexHeatmap::draw(
                    size_lgd,
                    x = grid::unit(0, "mm"),
                    y = grid::unit(1, "npc") -
                        grid::unit(legend_top_gap_mm, "mm"),
                    just = c("left", "top")
                )
                size_height_mm <- grid::convertHeight(
                    grid::grobHeight(size_lgd@grob),
                    "mm",
                    valueOnly = TRUE
                )
                ComplexHeatmap::draw(
                    colour_lgd,
                    x = grid::unit(0, "mm"),
                    y = grid::unit(
                        1, "npc"
                    ) - grid::unit(
                        legend_top_gap_mm + size_height_mm + 4,
                        "mm"
                    ),
                    just = c("left", "top")
                )
            },
            slice = 1
        )
    }

    vertical_ht <- ComplexHeatmap::`%v%`(upset_ht, dot_ht)
    ComplexHeatmap::draw(
        vertical_ht,
        show_heatmap_legend = FALSE,
        show_annotation_legend = FALSE,
        post_fun = draw_dot_legends
    )
    invisible(vertical_ht)
}
