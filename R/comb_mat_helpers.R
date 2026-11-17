#' Modified from `ComplexHeatmap::make_comb_mat()` (ComplexHeatmap 2.29.0).
#'
#' @param set_list Named list of non-empty character vectors.
#' @param set_on_rows Whether set membership is stored in rows.
#' @return A `comb_mat` matrix with set and combination sizes.
#' @keywords internal
.make_comb_mat <- function(set_list, set_on_rows = TRUE) {
    universe <- unique(unlist(set_list, use.names = FALSE))
    membership <- matrix(
        unlist(lapply(set_list, function(x) universe %in% x),
            use.names = FALSE),
        nrow = length(universe),
        ncol = length(set_list),
        dimnames = list(NULL, names(set_list))
    )
    patterns <- apply(membership, 1L, function(x) {
        paste0(as.integer(x), collapse = "")
    })
    comb_names <- unique(patterns)
    comb_mat <- matrix(
        as.integer(unlist(strsplit(comb_names, ""),
            use.names = FALSE)),
        nrow = length(comb_names),
        byrow = TRUE,
        dimnames = list(NULL, names(set_list))
    )
    comb_size <- tabulate(
        match(patterns, comb_names),
        nbins = length(comb_names)
    )
    comb_order <- order(rowSums(comb_mat), comb_names, decreasing = TRUE)
    comb_mat <- comb_mat[comb_order, , drop = FALSE]
    comb_size <- as.integer(comb_size[comb_order])

    if (set_on_rows) {
        m <- t(comb_mat)
        rownames(m) <- names(set_list)
    } else {
        m <- comb_mat
        colnames(m) <- names(set_list)
    }
    attr(m, "set_size") <- as.integer(colSums(membership))
    attr(m, "comb_size") <- comb_size
    attr(m, "data") <- set_list
    attr(m, "param") <- list(
        mode = "distinct",
        universal_set = NULL,
        set_on_rows = set_on_rows
    )
    class(m) <- c("comb_mat", "matrix")
    m
}

#' Modified from `ComplexHeatmap::comb_name()` (ComplexHeatmap 2.29.0).
#'
#' @param m A combination matrix from `.make_comb_mat()`.
#' @return Character vector of binary combination patterns.
#' @keywords internal
.comb_name <- function(m) {
    set_on_rows <- attr(m, "param")$set_on_rows
    comb_mat <- m
    class(comb_mat) <- "matrix"
    if (set_on_rows) {
        apply(comb_mat, 2, paste, collapse = "")
    } else {
        apply(comb_mat, 1, paste, collapse = "")
    }
}

#' Modified from `ComplexHeatmap::comb_size()` (ComplexHeatmap 2.29.0).
#'
#' @param m A combination matrix from `.make_comb_mat()`.
#' @return Named integer vector of combination sizes.
#' @keywords internal
.comb_size <- function(m) {
    size <- attr(m, "comb_size")
    names(size) <- .comb_name(m)
    size
}

#' Modified from ComplexHeatmap 2.29.0's `upset_top_annotation()` and
#' `upset_right_annotation()` to avoid the failing combination accessors.
#'
#' @param m A combination matrix from `.make_comb_mat()`.
#' @param position Annotation position, either `"top"` or `"right"`.
#' @param gp Graphics parameters for the bars.
#' @param annotation_name_gp Graphics parameters for the annotation name.
#' @param axis_param Axis parameters for the bar plot.
#' @param annotation_name_side Side on which to display the annotation name.
#' @return A ComplexHeatmap annotation object.
#' @keywords internal
.comb_annotation <- function(
    m,
    position,
    gp,
    annotation_name_gp,
    axis_param,
    annotation_name_side = NULL
) {
    set_on_rows <- attr(m, "param")$set_on_rows
    if (position == "top") {
        if (set_on_rows) {
            values <- .comb_size(m)
            name <- "intersection_size"
            label <- "Intersection\nsize"
            height <- 3
        } else {
            values <- attr(m, "set_size")
            name <- "set_size"
            label <- "Set\nsize"
            height <- 2
        }
        annotation <- list(ComplexHeatmap::anno_barplot(
            values,
            border = FALSE,
            gp = gp,
            height = grid::unit(height, "cm"),
            axis_param = axis_param
        ))
        names(annotation) <- name
        annotation <- c(annotation, list(
            show_annotation_name = TRUE,
            annotation_name_gp = annotation_name_gp,
            annotation_name_side = if (is.null(annotation_name_side)) {
                "left"
            } else {
                annotation_name_side
            },
            annotation_name_rot = 0,
            annotation_label = label
        ))
        return(do.call(ComplexHeatmap::HeatmapAnnotation, annotation))
    }

    if (position == "right") {
        if (set_on_rows) {
            values <- attr(m, "set_size")
            name <- "set_size"
            label <- "Set size"
            width <- 2
        } else {
            values <- .comb_size(m)
            name <- "intersection_size"
            label <- "Intersection\nsize"
            width <- 3
        }
        annotation <- list(ComplexHeatmap::anno_barplot(
            values,
            border = FALSE,
            gp = gp,
            width = grid::unit(width, "cm"),
            which = "row",
            axis_param = axis_param
        ))
        names(annotation) <- name
        annotation <- c(annotation, list(
            show_annotation_name = TRUE,
            annotation_name_gp = annotation_name_gp,
            annotation_name_side = if (is.null(annotation_name_side)) {
                "bottom"
            } else {
                annotation_name_side
            },
            annotation_label = label
        ))
        return(do.call(ComplexHeatmap::rowAnnotation, annotation))
    }

    stop("'position' must be either 'top' or 'right'.")
}
