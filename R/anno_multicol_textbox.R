# ---- Local textbox grob (avoids ComplexHeatmap:::textbox_grob) ----
#
# Adapted from ComplexHeatmap:::textbox_grob, with internal helper functions
# .recycle_gp and .subset_gp
# (ComplexHeatmap:::recycle_gp, ComplexHeatmap:::subset_gp)
# included here to avoid dependency on ComplexHeatmap's internal functions.
#
# Simplified text layout: stacks text strings vertically, each with its
# own gp (color, fontsize, fontfamily, fontface). No word-wrapping or
# random colours - our caller always provides explicit colours and
# controls text length via top_n.

#' Local textbox grob
#'
#' Stacks text strings vertically with per-string graphical parameters.
#' Adapted from `ComplexHeatmap:::textbox_grob` to avoid dependency on
#' unexported ComplexHeatmap internals.
#'
#' @param text Character vector of text strings to display.
#' @param x X position as a `grid::unit` (default: 0.5 npc).
#' @param y Y position as a `grid::unit` (default: 0.5 npc).
#' @param just Justification of the viewport (default: `"centre"`).
#' @param gp A `grid::gpar` object with graphical parameters for text.
#'   Elements `col`, `fontsize`, `fontfamily`, and `fontface` are
#'   recycled to match `length(text)`.
#' @param background_gp A `grid::gpar` object for the background rectangle
#'   (default: transparent fill and border).
#' @param round_corners Logical; whether to draw rounded corners on the
#'   background rectangle (default: `FALSE`).
#' @param r Corner radius as a `grid::unit` when `round_corners = TRUE`
#'   (default: 0.1 snpc).
#' @param line_space Spacing between lines as a `grid::unit`
#'   (default: 4 pt).
#' @param text_space Spacing between inline text segments as a
#'   `grid::unit` (default: 4 pt).
#' @param max_width Maximum width as a `grid::unit`
#'   (default: 100 mm).
#' @param padding Padding around content as a `grid::unit` vector of
#'   length 1, 2, or 4 (default: 4 pt all sides).
#' @param first_text_from `"top"` (default) or `"bottom"`.
#' @param add_new_line Logical; force each text string to a new line
#'   (default: `FALSE`).
#' @param word_wrap Logical; enable word wrapping (default: `FALSE`).
#'
#' @return A `grid::gTree` of class `"textbox"`.
#' @keywords internal
.textbox_grob <- function(
    text, x = unit(0.5, "npc"),
    y = unit(0.5, "npc"), just = "centre",
    gp = grid::gpar(),
    background_gp = grid::gpar(col = "transparent", fill = "transparent"),
    round_corners = FALSE, r = unit(0.1, "snpc"),
    line_space = unit(4, "pt"),
    text_space = unit(4, "pt"),
    max_width = unit(100, "mm"),
    padding = unit(4, "pt"),
    first_text_from = "top",
    add_new_line = FALSE, word_wrap = FALSE
) {
    n_text <- length(text)
    for (nm in c("fontsize", "fontfamily", "fontface")) {
        if (is.null(gp[[nm]])) {
            gp[[nm]] <- rep(grid::get.gpar(nm)[[1]], n_text)
        }
    }
    # From ComplexHeatmap:::recycle_gp
    .recycle_gp <- function(gp, n = 1) {
        for (i in seq_along(gp)) {
            x <- gp[[i]]
            if (n > 0) {
                gp[[i]] <- c(
                    rep(x, floor(n / length(x))),
                    x[seq_len(n %% length(x))]
                )
            } else {
                gp[[i]] <- x[1]
            }
        }
        return(gp)
    }

    # From ComplexHeatmap:::subset_gp
    .subset_gp <- function(gp, i) {
        gp <- lapply(gp, function(x) {
            if (length(x) == 1) {
                x
            } else {
                x[i]
            }
        })
        class(gp) <- "gpar"
        return(gp)
    }

    gp <- .recycle_gp(gp, n_text)
    vp_x <- x
    vp_y <- y
    vp_just <- just
    text_lt <- list(
        text = text, space_after = rep(TRUE, n_text),
        new_line_after = rep(FALSE, n_text)
    )
    if (add_new_line && !word_wrap) {
        text_lt$new_line_after <- rep(TRUE, n_text)
    } else if (word_wrap) {
        if (first_text_from == "bottom") {
            od <- rev(seq_len(n_text))
            text <- text[od]
            gp <- .subset_gp(gp, od)
            first_text_from <- "top"
        }
        words <- strsplit(text, "\\s+")
        words <- lapply(words, function(x) {
            n <- length(x)
            x2 <- rep(" ", 2 * n - 1)
            x2[2 * seq_len(n) - 1] <- x
            x2
        })
        nw <- vapply(words, length, integer(1))
        gp2 <- grid::gpar()
        gp2$col <- rep(gp$col, times = nw)
        gp2$fontsize <- rep(gp$fontsize, times = nw)
        gp2$fontfamily <- rep(gp$fontfamily, times = nw)
        gp2$fontface <- rep(gp$fontface, times = nw)
        text_lt <- list(
            text = unlist(words), space_after = rep(FALSE, sum(nw)),
            new_line_after = rep(FALSE, sum(nw))
        )
        text_lt$space_after[cumsum(nw)] <- TRUE
        if (add_new_line) {
            text_lt$new_line_after[cumsum(nw)] <- TRUE
        }
        gp <- gp2
    }
    n <- length(text_lt$text)
    text_gb_lt <- lapply(seq_len(n), function(i) {
        grid::textGrob(text_lt$text[i],
            gp = .subset_gp(gp, i)
        )
    })
    text_width <- vapply(text_gb_lt, function(gb) {
        grid::convertWidth(grid::grobWidth(gb),
            "mm",
            valueOnly = TRUE
        )
    }, 0)
    text_height <- vapply(text_gb_lt, function(gb) {
        grid::convertHeight(grid::grobHeight(gb),
            "mm",
            valueOnly = TRUE
        )
    }, 0)
    if (grid::is.unit(line_space)) {
        line_space <- grid::convertHeight(line_space, "mm", valueOnly = TRUE)
    }
    if (grid::is.unit(text_space)) {
        text_space <- grid::convertWidth(text_space, "mm", valueOnly = TRUE)
    }
    x <- numeric(n)
    y <- numeric(n)
    if (grid::is.unit(max_width)) {
        max_width <- grid::convertWidth(max_width, "mm", valueOnly = TRUE)
    }
    w <- max(text_width)
    max_width <- max(max_width, max(text_width))
    if (first_text_from == "bottom") {
        current_line_width <- text_width[1]
        x[1] <- 0
        y[1] <- 0
        h <- text_height[1]
        for (i in seq_len(n)[-1]) {
            if (current_line_width + text_width[i] + text_space >
                max_width || text_lt$new_line_after[i - 1]) {
                x[i] <- 0
                y[i] <- h + line_space
                current_line_width <- text_width[i]
                w <- max(w, current_line_width)
                h <- y[i] + text_height[i]
            } else {
                x[i] <- current_line_width +
                    text_space * text_lt$space_after[i - 1]
                y[i] <- y[i - 1]
                current_line_width <- x[i] + text_width[i]
                w <- max(w, current_line_width)
                h <- max(h, y[i] + text_height[i])
            }
        }
        just <- c(0, 0)
    } else if (first_text_from == "top") {
        current_line_width <- text_width[1]
        x[1] <- 0
        y[1] <- 0
        h <- -text_height[1]
        prev_line_ind <- 1
        for (i in seq_len(n)[-1]) {
            if (current_line_width + text_width[i] + text_space >
                max_width || text_lt$new_line_after[i - 1]) {
                y[prev_line_ind] <- h
                prev_line_ind <- i
                x[i] <- 0
                y[i] <- h - line_space
                current_line_width <- text_width[i]
                w <- max(w, current_line_width)
                h <- y[i] - text_height[i]
            } else {
                x[i] <- current_line_width +
                    text_space * text_lt$space_after[i - 1]
                y[i] <- y[i - 1]
                current_line_width <- x[i] + text_width[i]
                w <- max(w, current_line_width)
                h <- min(h, y[i] - text_height[i])
                prev_line_ind <- c(prev_line_ind, i)
            }
        }
        y[prev_line_ind] <- h
        y <- y - h
        h <- -h
        just <- c(0, 0)
    } else {
        stop("`first_text_from` can be 'top' or 'bottom'")
    }
    if (length(padding) == 1) {
        padding <- rep(padding, 4)
    } else if (length(padding) == 2) {
        padding <- grid::unit.c(
            padding[1], padding[2], padding[1],
            padding[2]
        )
    }
    padding <- grid::convertWidth(padding, "mm", valueOnly = TRUE)
    w <- w + padding[2] + padding[4]
    h <- h + padding[1] + padding[3]
    x <- x + padding[2]
    y <- y + padding[1]
    gl <- grid::gList(if (round_corners) {
        grid::roundrectGrob(gp = background_gp, r = r)
    } else {
        grid::rectGrob(gp = background_gp)
    }, grid::textGrob(text_lt$text,
        x = x, y = y, gp = gp, default.units = "mm",
        just = just
    ))
    gb <- grid::gTree(children = gl, cl = "textbox", vp = grid::viewport(
        x = vp_x,
        y = vp_y, just = vp_just, width = unit(w, "mm"), height = unit(
            h,
            "mm"
        )
    ))
    return(gb)
}

# ---- Internal: reformat cat_terms to text ----

#' Reformat enrichment category terms into per-module text lists
#'
#' Transposes `cat_terms` from category-first (`cat_terms[[cate]][[mod]]`)
#' to module-first (`text[[mod]][[cate]]`) layout expected by
#' `.multicol_textbox_grob`. Missing categories get empty data.frames.
#'
#' @param cat_terms A named list of categories, each containing a named
#'   list of modules with term annotation data.frames.
#' @return A named list of modules, each containing a named list of
#'   categories with term data.frames.
#' @keywords internal
.reformat_cat_terms <- function(cat_terms) {
    mods <- unique(unlist(lapply(cat_terms, names)))
    out <- vector("list", length(mods))
    names(out) <- mods
    for (mod in mods) {
        out[[mod]] <- list()
        for (cate in names(cat_terms)) {
            if (mod %in% names(cat_terms[[cate]])) {
                out[[mod]][[cate]] <- cat_terms[[cate]][[mod]]
            } else {
                out[[mod]][[cate]] <- data.frame(
                    text = character(0), col = character(0),
                    fontsize = numeric(0), stringsAsFactors = FALSE
                )
            }
        }
    }
    out
}

# ---- Internal: per-column textbox helpers ----

#' Build a single-column textbox grob
#'
#' Stacks term text strings vertically using the local `.textbox_grob`.
#' Returns an invisible placeholder grob if `df` has no rows.
#'
#' @param df A data.frame with columns `text`, `col`, `fontsize`, and
#'   optionally `fontfamily` and `fontface`.
#' @param ... Passed to `.textbox_grob`.
#' @return A textbox grob.
#' @keywords internal
.multicol_textbox_col_grob <- function(df, ...) {
    if (!is.data.frame(df) || nrow(df) == 0) {
        return(grid::nullGrob())
    }
    gp <- grid::gpar(
        col = as.character(df$col),
        fontsize = df$fontsize,
        fontfamily = if ("fontfamily" %in% names(df)) df$fontfamily,
        fontface = if ("fontface" %in% names(df)) df$fontface
    )
    .textbox_grob(
        text = df$text, gp = gp,
        add_new_line = TRUE, ...
    )
}

#' Measure column widths for multi-column textbox layout
#'
#' Computes the maximum grob width per category in millimetres, used to
#' align columns across modules.
#'
#' @param text A per-module list of category term data.frames (output from
#'   `.reformat_cat_terms`).
#' @param textbox_args Additional arguments for
#'   `.multicol_textbox_col_grob`.
#' @return A named numeric vector of column widths in mm.
#' @keywords internal
.measure_multicol_widths <- function(text, textbox_args = list()) {
    cats <- unique(unlist(lapply(text, names), use.names = FALSE))
    if (length(cats) == 0) {
        return(stats::setNames(numeric(0), character(0)))
    }

    textbox_args <- textbox_args[
        setdiff(
            names(textbox_args),
            c("col_gap", "background_gp", "padding", "column_widths_mm")
        )
    ]
    widths_mm <- stats::setNames(numeric(length(cats)), cats)
    for (cate in cats) {
        widths_mm[cate] <- max(vapply(text, function(x) {
            df <- if (cate %in% names(x)) {
                x[[cate]]
            } else {
                data.frame(
                    text = character(0), col = character(0),
                    fontsize = numeric(0), stringsAsFactors = FALSE
                )
            }
            gb <- do.call(
                .multicol_textbox_col_grob,
                c(list(df = df), textbox_args)
            )
            grid::convertWidth(grid::grobWidth(gb), "mm", valueOnly = TRUE)
        }, numeric(1)))
    }
    widths_mm
}

# ---- Internal: multi-column textbox grob ----

#' Multi-column textbox grob
#'
#' Arranges per-category textbox grobs side-by-side in absolute mm units.
#' Returns a `gTree` of class `"multicol_textbox"` with attributes
#' `column_labels`, `column_widths_mm`, `column_gap_mm`, and `padding_mm`
#' for downstream header decoration.
#'
#' @param x A named list of category term data.frames (one per column).
#' @param col_gap Gap between columns as a `grid::unit` (default: 2 mm).
#' @param background_gp Background graphical parameters.
#' @param padding Padding around content as a `grid::unit` vector
#'   (top, left, bottom, right). Default: 0 mm all sides.
#' @param column_widths_mm Optional numeric vector of pre-computed column
#'   widths in mm. If NULL, auto-computed from grob widths.
#' @param ... Passed to `.multicol_textbox_col_grob`.
#' @return A `gTree` of class `"multicol_textbox"`.
#' @keywords internal
.multicol_textbox_grob <- function(
    x, col_gap = unit(2, "mm"),
    background_gp = grid::gpar(fill = "#F7F7F7", col = "#CCCCCC"),
    padding = unit(c(0, 0, 0, 0), "mm"),
    column_widths_mm = NULL, ...
) {
    cats <- names(x)
    if (length(cats) == 0) {
        return(grid::nullGrob())
    }

    col_grobs <- lapply(cats, function(cate) {
        .multicol_textbox_col_grob(x[[cate]], ...)
    })
    names(col_grobs) <- cats

    widths_mm <- vapply(col_grobs, function(g) {
        grid::convertWidth(grid::grobWidth(g), "mm",
            valueOnly = TRUE
        )
    }, numeric(1))
    slot_widths_mm <- widths_mm
    if (!is.null(column_widths_mm)) {
        slot_widths_mm <- column_widths_mm[cats]
        slot_widths_mm[is.na(slot_widths_mm)] <-
            widths_mm[is.na(slot_widths_mm)]
    }
    heights_mm <- vapply(col_grobs, function(g) {
        grid::convertHeight(grid::grobHeight(g), "mm",
            valueOnly = TRUE
        )
    }, numeric(1))
    gap_mm <- grid::convertWidth(col_gap, "mm", valueOnly = TRUE)
    pad_l <- grid::convertWidth(padding[[2]], "mm", valueOnly = TRUE)
    pad_r <- grid::convertWidth(padding[[4]], "mm", valueOnly = TRUE)
    pad_t <- grid::convertHeight(padding[[1]], "mm", valueOnly = TRUE)
    pad_b <- grid::convertHeight(padding[[3]], "mm", valueOnly = TRUE)

    n <- length(cats)
    total_w_mm <- sum(slot_widths_mm) + gap_mm * max(n - 1L, 0L) 
        + pad_l + pad_r
    total_h_mm <- max(heights_mm) + pad_t + pad_b

    vp <- grid::viewport(
        width  = unit(total_w_mm, "mm"),
        height = unit(total_h_mm, "mm")
    )

    children <- list()
    x_mm <- pad_l
    for (i in seq_len(n)) {
        children <- c(children, list(grid::editGrob(col_grobs[[i]],
            vp = grid::viewport(
                x = unit(x_mm, "mm"), y = unit(total_h_mm / 2, "mm"),
                width = unit(slot_widths_mm[i], "mm"),
                height = unit(heights_mm[i], "mm"),
                just = c(0, 0.5)
            )
        )))

        x_mm <- x_mm + slot_widths_mm[i] + gap_mm
    }
    gb <- grid::gTree(
        children = do.call(grid::gList, children), vp = vp,
        cl = "multicol_textbox"
    )
    attr(gb, "column_labels") <- cats
    attr(gb, "column_widths_mm") <- unname(slot_widths_mm)
    attr(gb, "column_gap_mm") <- gap_mm
    attr(gb, "padding_mm") <- c(
        top = pad_t, left = pad_l,
        bottom = pad_b, right = pad_r
    )
    gb
}

# ---- .anno_multicol_textbox (based on ComplexHeatmap:::anno_textbox) ----

#' Multi-column textbox row annotation for ComplexHeatmap
#'
#' Row annotation that displays multiple enrichment categories side by side
#' in a single annotation block. Follows the same interface convention as
#' `ComplexHeatmap:::anno_textbox`. Column headers are drawn above the
#' first module via `decorate_annotation`-style header placement.
#'
#' @param align_to A vector of group labels, a list of indices, or a
#'   named list mapping groups to row indices (same as
#'   `ComplexHeatmap::anno_textbox`).
#' @param text A list of per-module enrichment terms, where each element
#'   is a named list of category term data.frames (columns `text`, `col`,
#'   `fontsize`). Output from `.reformat_cat_terms`.
#' @param background_gp Background graphical parameters for the annotation
#'   block (fill, col, lty, lwd).
#' @param which `"row"` (only row annotations are supported).
#' @param by Link type: `"anno_link"` (default, aligns to heatmap rows)
#'   or `"anno_block"`.
#' @param side `"right"` (default) or `"left"`.
#' @param ... Additional arguments passed to `.multicol_textbox_grob`.
#'   Notable options: `show_headers` (logical), `header_gp` (gpar for
#'   header text), `header_offset` (unit offset above annotation).
#' @return A ComplexHeatmap annotation object (from `anno_link` or
#'   `anno_block`).
#' @keywords internal
.anno_multicol_textbox <- function(
    align_to, text,
    background_gp = grid::gpar(
        fill = "#DDDDDD",
        col = "#AAAAAA"
    ),
    which = c("row", "column"), by = "anno_link", side = c("right", "left"),
    ...
) {
    if (is.null(background_gp$fill)) background_gp$fill <- "#DDDDDD"
    if (is.null(background_gp$col)) background_gp$col <- "#AAAAAA"
    if (is.null(background_gp$lty)) background_gp$lty <- 1
    if (is.null(background_gp$lwd)) background_gp$lwd <- 1
    which <- match.arg(which)
    if (which == "column") {
        stop(".anno_multicol_textbox() can only be used as row annotation.")
    }

    if (is.numeric(align_to) && (is.character(text) || is.data.frame(text))) {
        align_to <- list(v = align_to)
        text <- list(v = text)
    } else if (is.atomic(align_to) && is.list(text)) {
        align_to <- split(seq_along(align_to), align_to)
        cn <- intersect(names(align_to), names(text))
        if (length(cn) == 0) {
            stop("names of `text` should have overlap to levels in ",
                "`align_to`.")
        }
        align_to <- align_to[cn]
        text <- text[cn]
    } else if (is.list(align_to)) {
        if (!is.list(text)) {
            stop("Since `align_to` is a list, `text` should be a list.")
        }
        if (length(align_to) != length(text)) {
            stop("`align_to` and `text` should have the same length.")
        }
        if (!is.null(names(align_to)) && !is.null(names(text))) {
            if (length(setdiff(names(align_to), names(text))) == 0) {
                text <- text[names(align_to)]
            } else {
                stop("`align_to` and `text` should have the same names.")
            }
        }
    } else {
        stop("Format of `align_to` or `text` is wrong.")
    }

    textbox_args <- list(...)
    show_headers <- isTRUE(textbox_args$show_headers)
    header_gp <- textbox_args$header_gp
    header_offset <- textbox_args$header_offset
    textbox_args$show_headers <- NULL
    textbox_args$header_gp <- NULL
    textbox_args$header_offset <- NULL

    if (is.null(header_gp)) {
        header_gp <- grid::gpar(fontface = "bold")
    }
    if (is.null(header_gp$fontsize)) {
        header_sizes <- unlist(lapply(text, function(cat_list) {
            vapply(cat_list, function(df) {
                if (!is.data.frame(df) || !nrow(df) ||
                    !"fontsize" %in% names(df)) {
                    return(8)
                }
                max(df$fontsize, na.rm = TRUE)
            }, numeric(1))
        }), use.names = FALSE)
        header_gp$fontsize <- if (length(header_sizes)) 
            max(header_sizes) else 8
    }
    if (is.null(header_offset)) {
        header_offset <- unit(1.5, "mm")
    }
    dev_was_open <- length(grDevices::dev.list()) > 0
    if (!dev_was_open) {
        tmp <- tempfile(fileext = ".pdf")
        grDevices::pdf(tmp)
    }

    column_widths_mm <- .measure_multicol_widths(text, textbox_args)
    gbl <- lapply(text, function(x) {
        do.call(.multicol_textbox_grob, c(
            list(x = x, column_widths_mm = column_widths_mm),
            textbox_args
        ))
    })

    if (!dev_was_open) {
        grDevices::dev.off()
        unlink(tmp)
    }
    margin <- unit(1, "mm")
    gbl_h <- lapply(gbl, function(x) {
        grid::convertHeight(x$vp$height, "cm") + margin
    })
    gbl_h <- do.call(grid::unit.c, gbl_h)
    gbl_w <- lapply(gbl, function(x) grid::convertWidth(x$vp$width, "cm"))
    gbl_w <- do.call(grid::unit.c, gbl_w)
    gbl_w <- max(gbl_w) + margin

    side <- match.arg(side)
    first_nm <- names(gbl)[1]

    # Detect empty panels (all categories have no terms) for
    # transparent background
    empty_panels <- vapply(text, function(txt) {
        all(vapply(txt, function(df) {
            !is.data.frame(df) || nrow(df) == 0
        }, logical(1)))
    }, logical(1))

    if (by %in% c("anno_link", "anno_zoom")) {
        panel_fun <- function(index, nm) {
            grid::pushViewport(grid::viewport(clip = "off"))
            is_empty <- isTRUE(empty_panels[[nm]])
            bg_fill <- if (is_empty) "transparent"
                else background_gp$fill
            bg_col  <- if (is_empty) "transparent"
                else background_gp$col

            grid::grid.rect(gp = grid::gpar(
                fill = bg_fill, col = bg_col,
                lty = background_gp$lty, lwd = background_gp$lwd
            ))
            if (!is_empty) {
                if (side == "right") {
                    grid::grid.lines(c(0, 1, 1, 0), c(0, 0, 1, 1),
                        gp = grid::gpar(
                            col = background_gp$col,
                            lty = background_gp$lty,
                            lwd = background_gp$lwd
                        ),
                        default.units = "npc"
                    )
                } else {
                    grid::grid.lines(c(1, 0, 0, 1), c(0, 0, 1, 1),
                        gp = grid::gpar(
                            col = background_gp$col,
                            lty = background_gp$lty,
                            lwd = background_gp$lwd
                        ),
                        default.units = "npc"
                    )
                }
            }
            gb <- gbl[[nm]]
            if (show_headers && identical(nm, first_nm)) {
                col_labels <- attr(gb, "column_labels")
                col_widths <- attr(gb, "column_widths_mm")
                col_gap_mm <- attr(gb, "column_gap_mm")
                padding_mm <- attr(gb, "padding_mm")
                x_left_mm <- padding_mm["left"]
                x_pos_mm <- numeric(length(col_widths))

                cur_x <- x_left_mm

                for (i in seq_along(col_widths)) {
                    x_pos_mm[i] <- cur_x + col_widths[i] / 2
                    cur_x <- cur_x + col_widths[i] + col_gap_mm
                }
                grid::grid.text(
                    label = col_labels,
                    x = unit(x_pos_mm, "mm"),
                    y = unit(1, "npc") + header_offset,
                    just = c("center", "bottom"),
                    gp = header_gp
                )
            }
            grid::pushViewport(grid::viewport(
                width = unit(1, "npc") - margin,
                height = unit(1, "npc") - margin
            ))
            gb$vp$x <- gb$vp$width * 0.5
            gb$vp$y <- gb$vp$height * 0.5
            grid::grid.draw(gb)
            grid::popViewport()
            grid::popViewport()
        }
        ComplexHeatmap::anno_link(
            align_to = align_to, panel_fun = panel_fun,
            which = "row",
            size = gbl_h, gap = unit(2, "mm"),
            width = gbl_w + unit(5, "mm"),
            link_gp = background_gp, internal_line = FALSE,
            side = side
        )
    } else {
        panel_fun <- function(index, nm) {
            grid::pushViewport(grid::viewport())
            is_empty <- isTRUE(empty_panels[[nm]])
            bg_fill <- if (is_empty) "transparent"
                else background_gp$fill
            bg_col  <- if (is_empty) "transparent"
                else background_gp$col
            grid::grid.rect(gp = grid::gpar(
                fill = bg_fill,
                col = bg_col, lty = background_gp$lty,
                lwd = background_gp$lwd
            ))
            grid::pushViewport(grid::viewport(
                width  = unit(1, "npc") - margin,
                height = unit(1, "npc") - margin
            ))
            gb <- gbl[[nm]]
            gb$vp$x <- unit(0, "npc")
            gb$vp$justification[1] <- 0
            gb$vp$valid.just[1] <- 0
            grid::grid.draw(gb)
            grid::popViewport()
            grid::popViewport()
        }
        ComplexHeatmap::anno_block(
            align_to = align_to, which = "row",
            panel_fun = panel_fun, width = gbl_w + unit(5, "mm")
        )
    }
}
