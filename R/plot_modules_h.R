#' Normalise and order module labels
#'
#' Converts numeric or digit-string module labels to a factor with levels
#' in natural numeric order. When `prefix = TRUE` (default), labels are
#' formatted as `"M0", "M1", "M2", ...`. When `prefix = FALSE`, bare
#' digit strings are returned (e.g. `"0", "1", "2"`). WGCNA colour names
#' (e.g. `"turquoise"`, `"blue"`) are unaffected by `prefix` and returned
#' as a factor preserving their original order of first appearance.
#'
#' Callers that need a different ordering (e.g. by module size) should
#' re-factor the result after calling this function, as `plot_modules_h()`
#' does.
#'
#' @param x A vector of module labels (numeric, character digits, or
#'   WGCNA colour names).
#' @param prefix Logical; if `TRUE` (default), numeric labels are prefixed
#'   with `"M"`. Ignored for WGCNA colour names.
#' @return A factor of module labels. Numeric labels are ordered by their
#'   numeric value; non-numeric labels preserve their original order of
#'   first appearance.
#' @keywords internal
.module_labels <- function(x, prefix = TRUE) {
    if (is.numeric(x)) {
        nums <- as.integer(x)
        labels <- if (prefix) paste0("M", nums) else as.character(nums)
        return(factor(labels, levels = unique(labels[order(nums)])))
    }
    labels <- as.character(x)
    # Strip optional "M" prefix and check if numeric
    stripped <- gsub("^M", "", labels)
    if (all(grepl("^[0-9]+$", stripped))) {
        nums <- as.integer(stripped)
        labels <- if (prefix) paste0("M", nums) else as.character(nums)
        return(factor(labels, levels = unique(labels[order(nums)])))
    }
    # WGCNA colour names: keep as-is, preserve input order
    factor(labels, levels = unique(labels))
}

#' Build marked-feature term_list for multi-column textbox
#'
#' Build a per-module term_list from `mark_features` for use as a column
#' in `.anno_multicol_textbox`. Accepts a named list or a character vector.
#' Colours are taken from `mark_state$feat_col`, including any custom
#' `feature_colours` supplied to `plot_modules_h()`.
#'
#' @param mark_features A named list of character vectors, mapping
#'   category names to feature IDs, or a plain character vector of
#'   feature IDs (shown in black by default).
#' @param module_df A data frame with columns `Feature` and `Module`.
#' @param mark_state A list with element `feat_col`, a named character
#'   vector mapping feature IDs to colours.
#' @param fontsize Numeric font size for term text.
#' @param mod_levels Character vector of module levels to include.
#' @return A named list of per-module data.frames (columns `text`,
#'   `col`, `fontsize`), or `NULL` if no marked features are found in
#'   the module data.
#' @keywords internal
.build_marked_features_anno <- function(mark_features, module_df, mark_state,
                                        fontsize, mod_levels) {
    if (!is.list(mark_features)) {
        # Character vector: all black, like anno_mark
        feats <- mark_features
        marked_in_mod <- intersect(feats, module_df$Feature)
    } else {
        # Named list: map each feature to its category
        feats <- stats::setNames(
            rep(names(mark_features), lengths(mark_features)),
            unlist(mark_features, use.names = FALSE)
        )
        marked_in_mod <- intersect(names(feats), module_df$Feature)
    }
    if (length(marked_in_mod) == 0) return(NULL)

    mod_map <- stats::setNames(as.character(module_df$Module),
        module_df$Feature)
    cols <- mark_state$feat_col[marked_in_mod]

    termanno <- data.frame(
        id = mod_map[marked_in_mod],
        term = marked_in_mod,
        col = cols,
        fontsize = fontsize,
        stringsAsFactors = FALSE
    ) |>
        dplyr::filter(id %in% mod_levels) |>
        dplyr::group_by(id) |>
        dplyr::slice_head(n = 3) |>
        dplyr::ungroup() |>
        dplyr::select(id, term, col, fontsize)

    if (nrow(termanno) == 0) return(NULL)

    term_list <- lapply(mod_levels, function(mod) {
        tmp <- termanno[termanno$id == mod, ]
        if (nrow(tmp) == 0) {
            data.frame(text = character(0), col = character(0),
                fontsize = numeric(0), stringsAsFactors = FALSE)
        } else {
            data.frame(text = tmp$term, col = tmp$col,
                fontsize = tmp$fontsize, stringsAsFactors = FALSE)
        }
    })
    names(term_list) <- mod_levels
    term_list
}

#' Build the selected-module profile and textbox annotation
#'
#' Combines aligned profile and enrichment panels inside a row `anno_zoom`.
#'
#' @return A ComplexHeatmap row annotation.
#' @keywords internal
.make_selected_anno_zoom <- function(
    profile_row_groups,
    profile_width_mm,
    link_width_mm,
    heatmap_height_mm = NULL,
    profile_height_mm = NULL,
    profile_panel_fun,
    text_panel_fun = NULL,
    text_row_groups = NULL,
    text_heights_mm = numeric(0),
    textbox_width_mm = 0,
    textbox_gap_mm = 2,
    link_gp = grid::gpar(
        fill = "grey98", col = "grey50", lwd = 1
    )
) {
    has_text <- !is.null(text_panel_fun) &&
        length(text_row_groups) > 0L
    valid_text_heights <- text_heights_mm[
        is.finite(text_heights_mm)
    ]
    max_text_height_mm <- if (
        has_text && length(valid_text_heights) > 0L
    ) {
        max(valid_text_heights)
    } else {
        0
    }
    n_panels <- length(profile_row_groups)
    if (is.null(heatmap_height_mm)) {
        panel_size <- rep(1, n_panels)
        min_panel_gap_mm <- 0
        target_panel_height_mm <- max(
            profile_width_mm, max_text_height_mm
        )
        target_gap <- (
            grid::unit(1, "npc") -
                grid::unit(n_panels * target_panel_height_mm, "mm")
        ) / (n_panels + 1)
        panel_gap <- grid::unit.pmax(
            grid::unit(min_panel_gap_mm, "mm"), target_gap
        )
        profile_panel_height <- grid::unit.pmin(
            grid::unit(1, "npc"),
            grid::unit(profile_width_mm, "mm")
        )
    } else {
        panel_height_mm <- max(
            profile_height_mm, max_text_height_mm
        )
        panel_gap_mm <- (
            heatmap_height_mm - n_panels * panel_height_mm
        ) / (n_panels + 1)
        if (panel_gap_mm < 0) {
            stop(
                "The heatmap is too short for selected anno_zoom panels."
            )
        }
        panel_size <- grid::unit(
            rep(panel_height_mm, n_panels), "mm"
        )
        panel_gap <- grid::unit(panel_gap_mm, "mm")
        profile_panel_height <- grid::unit(
            profile_height_mm, "mm"
        )
    }
    panel_extend <- if (n_panels == 1L) {
        grid::unit(c(0, 0), "mm")
    } else {
        -grid::unit.c(panel_gap, panel_gap)
    }
    anno_zoom_fun <- ComplexHeatmap::anno_zoom
    if (n_panels == 1L) {
        replace_mean_pos <- function(x) {
            if (!is.call(x)) {
                return(x)
            }
            if (length(x) > 2L &&
                (identical(x[[1]], as.name("=")) ||
                    identical(x[[1]], as.name("<-"))) &&
                identical(x[[2]], as.name("mean_pos"))) {
                original_mean_pos <- x[[3]]
                x[[3]] <- bquote(
                    if (which == "row" && n_boxes == 1) {
                        0.5
                    } else {
                        .(original_mean_pos)
                    }
                )
                return(x)
            }
            as.call(lapply(x, replace_mean_pos))
        }
        body(anno_zoom_fun) <- replace_mean_pos(
            body(anno_zoom_fun)
        )
    }
    content_width_mm <- profile_width_mm +
        if (has_text) textbox_gap_mm + textbox_width_mm else 0

    panel_fun <- function(index, nm) {
        grid::pushViewport(grid::viewport(
            x = grid::unit(0, "npc"),
            y = grid::unit(0.5, "npc"),
            width = grid::unit(profile_width_mm, "mm"),
            height = profile_panel_height,
            just = base::c("left", "center"),
            clip = "off"
        ))
        profile_panel_fun(index, nm)
        grid::popViewport()

        if (has_text && nm %in% base::names(text_row_groups)) {
            text_height <- text_heights_mm[[nm]]
            if (is.finite(text_height) && text_height > 0) {
                grid::pushViewport(grid::viewport(
                    x = grid::unit(
                        profile_width_mm + textbox_gap_mm, "mm"
                    ),
                    y = grid::unit(0.5, "npc"),
                    width = grid::unit(textbox_width_mm, "mm"),
                    height = grid::unit(text_height, "mm"),
                    just = base::c("left", "center"),
                    clip = "off"
                ))
                text_panel_fun(index, nm)
                grid::popViewport()
            }
        }
    }

    anno_zoom_fun(
        align_to = profile_row_groups,
        which = "row",
        side = "right",
        panel_fun = panel_fun,
        size = panel_size,
        gap = panel_gap,
        link_width = grid::unit(link_width_mm, "mm"),
        link_gp = link_gp,
        width = grid::unit(
            link_width_mm + content_width_mm, "mm"
        ),
        extend = panel_extend,
        internal_line = TRUE
    )
}

#' Plot modules (horizontal layout)
#'
#' @description Plot WGCNA modules' heatmaps, mean expression profiles,
#' and enrichment terms, aligned horizontally.
#'
#' Different from running WGCNA, the input data should have the
#' replicates merged, instead of having multiple samples per group, time and
#' feature (gene).
#'
#' If certain time points are missing in some groups, NA values are added.
#'
#' @param module A data frame with columns "Feature" and "Module"
#' @param se_obj_merged A SummarizedExperiment object, with one value for each
#' feature at each time point in each group (replicates merged), or a named list
#' of such objects (one per dataset). The colData of each object should contain
#' columns "Sample", "Group", and "Time". The object can be produced by
#' `split_groups()`, `merge_replicates()` and `merge_groups()`.
#' @param assay The assay to use in the SummarizedExperiment object: a numeric
#'   index or character name, e.g. 1 or "orig" for original data, 2 or "norm"
#'   for time 0 normalised data. (Required, no default.)
#' @param scale Whether to scale the data (z-score) across samples for each
#' feature (default is TRUE)
#' @param ylabel Y axis label. A single character string shared by all
#'   datasets. Default: NULL, i.e. "Z-score of log2 abundance" when
#'   scale = TRUE, "Log2 abundance" otherwise.
#' @param profile_width Width of the mean expression profile panels in cm
#' (default: 3)
#' @param profile_link_width Width of the link between heatmap
#' and mean expression profile panels in cm (default: 1)
#' @param enrich_list A named enrichment list for modules,
#'   as produced by `enrichGO_list()$all`. Each named element
#'   should be a data.frame with columns `Cluster`, `Description`,
#'   and the rank column. Accepts any enrichment result with compatible format.
#'   When `enrich_category` is a vector, each category
#'   is shown as a separate column in the annotation textbox.
#'   If NULL (default), no enrichment terms will be displayed.
#' @param enrich_category Category name(s) of enrichment to display:
#'   a single string (e.g. `"BP"`), or a vector for multiple
#'   categories (e.g. `c("BP", "CC", "Features")`).
#'   When `"Features"` is included, `mark_features` are shown as a column
#'   in the textbox. Categories must be present in `enrich_list`, except
#'   for `"Features"`.
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
#'   Set to `NA` for `"Features"`.
#' @param enrich_wrap_width Optional maximum enrichment term width in cm.
#'   The default `NULL` disables wrapping.
#' @param fontsize Base font size (default: 8)
#' @param group_annotation Top annotation for group-split heatmap columns:
#'   `"name"` (default), `"colour"`, `"both"`, or `"none"`.
#'   In `"both"` mode, group names are above the colour strip.
#' @param return_legends Logical; if `TRUE`, return a list containing the
#'   Heatmap/HeatmapList and the native ComplexHeatmap legend objects.
#'   (default: `FALSE`)
#' @param legend_position Legend position: `"right"` (default) or
#'   `"bottom"`.
#' @param profile_scale Y-axis scaling for the mean expression profiles:
#'   `"dataset"` (default) uses one range across modules within each
#'   dataset, `"shared"` uses one range across all modules and datasets,
#'   and `"free"` uses an independent range for each module.
#' @param heatmap_width Width of the heatmap body in cm (default: 8)
#' @param heatmap_height Height of the heatmap body in cm (default: 8).
#'   Set to `NULL` to auto-size.
#' @param mark_features Features to label. Accepts three forms:
#'     - A **character vector**: all marked in black by default on the left.
#'     - A **named list** of character vectors, e.g.
#'         `list("Set 1" = c("gene1"), "Set 2" = c("gene2"))`.
#'         Names become categories with auto-assigned colours (via
#'         `ggsci::pal_jco()`).
#'     - A data frame with `Feature` and `Set` columns, such as the output
#'         of `extract_consensus_hubs()`. `Set` values become categories.
#'     Auto-routed: shown as `anno_mark` (left
#'     side, connecting lines) unless `"Features"` is present in
#'     `enrich_category`, in which case an `anno_multicol_textbox` (right side)
#'     is built instead, using the same colours.
#'     Set to `NULL` (default) to skip marking.
#' @param feature_colours Optional named character vector of colours for
#'   marked feature sets. Names must match feature set names in
#'   `mark_features`. If `mark_features` is a plain character vector, this
#'   parameter is ignored. Omitted sets retain their automatic colours.
#'   Default: `NULL`.
#' @param width Width of the saved image (default is 16 (cm))
#' @param height Height of the saved image (default is 12 (cm))
#' @param res Resolution of the saved image (except pdf and svg formats)
#'   (default: 300 (ppi))
#' @param device Image file format(s) for saving. Can be a character
#'   vector with one or more of `"png"`, `"pdf"`, `"svg"`, `"tiff"`, `"jpeg"`
#'   (default: `"png"`)
#' @param save Directory to save the plot, no saving if is NULL
#' (default is NULL)
#' @param suffix Suffix for the saved image file name (default: "")
#' @param names Optional character vector of dataset names for the plot
#'   subtitles when multi-dataset. Default: NULL, i.e.
#'   the names of `se_obj_merged`.
#' @param selected_modules Optional character or numeric module labels for
#'   which profiles and enrichment panels should be displayed. The default
#'   `NULL` keeps all modules.
#' @import ggplot2
#' @import SummarizedExperiment
#'
#' @returns A combined plot of heatmap, mean expression profile, and enrichment
#' terms for each WGCNA module. If `return_legends = TRUE`, returns a list with
#' elements `heatmap` and `legends` instead.
#' @export
#' @examples
#' data(example_obj)
#' example_obj <- normalise_to_start(example_obj)
#' example_obj_list <- split_groups(example_obj)
#' example_obj_merged_list <- merge_replicates(example_obj_list)
#' example_obj_merged <- merge_groups(example_obj_merged_list)
#'
#' data(example_net)
#' # select two modules for demonstration
#' example_module <- WGCNA_module(example_net) |>
#'     dplyr::filter(Module %in% c("1", "2"))
#'
#' data(example_go)
#'
#' plot_modules_h(example_module |> dplyr::filter(Module != '0'),
#'     example_obj_merged, assay = 2, scale = TRUE,
#'     ylabel = "Z-score of log2 expression",
#'     enrich_list = example_go$all,
#'     enrich_category = "BP",
#'     heatmap_width = 6, heatmap_height = 4)
#'
#' plot_modules_h(example_module |> dplyr::filter(Module != '0'),
#'     list(Exp1 = example_obj_merged, Exp2 = example_obj_merged),
#'     assay = 2, scale = TRUE,
#'     enrich_list = example_go$all,
#'     enrich_category = "BP",
#'     heatmap_width = 6, heatmap_height = 4)
#' @references https://github.com/junjunlab/ClusterGVis
plot_modules_h <- function(
    module,
    se_obj_merged,
    scale = TRUE,
    assay,
    ylabel = NULL,
    profile_width = 3,
    profile_link_width = 1,
    enrich_list = NULL,
    enrich_category = NULL,
    enrich_rank_by = "p.adjust",
    enrich_top_n = 3,
    enrich_threshold = 0.05,
    enrich_wrap_width = NULL,
    fontsize = 8,
    group_annotation = c("name", "colour", "both", "none"),
    return_legends = FALSE,
    profile_scale = c("dataset", "shared", "free"),
    heatmap_width = 8,
    heatmap_height = 8,
    mark_features = NULL,
    width = 16,
    height = 12,
    res = 300,
    suffix = "",
    device = "png",
    save = NULL,
    names = NULL,
    selected_modules = NULL,
    legend_position = c("right", "bottom"),
    feature_colours = NULL
) {
    .check_df(module, "module")
    .check_logical(scale, "scale")
    .check_positive(profile_width, "profile_width")
    .check_positive(profile_link_width, "profile_link_width")
    .check_positive(fontsize, "fontsize")
    .check_character(group_annotation, "group_annotation")
    group_annotation <- match.arg(group_annotation)
    .check_logical(return_legends, "return_legends")
    .check_character(profile_scale, "profile_scale")
    profile_scale <- match.arg(profile_scale)
    legend_position <- match.arg(legend_position)
    fontfamily <- "sans"
    .check_positive(heatmap_width, "heatmap_width")
    if (!is.null(heatmap_height))
        .check_positive(heatmap_height, "heatmap_height")
    .check_positive(width, "width")
    .check_positive(height, "height")
    .check_positive_int(res, "res")
    if (!is.null(selected_modules) &&
        !is.character(selected_modules) &&
        !is.numeric(selected_modules)) {
        stop("'selected_modules' must be character or numeric.")
    }
    if (!is.null(selected_modules) && length(selected_modules) == 0L) {
        stop("'selected_modules' must contain at least one module.")
    }
    if (!is.null(enrich_list) && !is.null(enrich_wrap_width) &&
        (!is.numeric(enrich_wrap_width) ||
            length(enrich_wrap_width) != 1L ||
            !is.finite(enrich_wrap_width) || enrich_wrap_width <= 0)) {
        stop("'enrich_wrap_width' must be one positive number or NULL.")
    }
    if (is.null(ylabel)) {
        ylabel <- if (scale) "Z-score of log2 abundance" else "Log2 abundance"
    }
    .check_character(ylabel, "ylabel")
    .check_character(suffix, "suffix")
    .check_character(device, "device")

    # Transform hub/annotation tables to the named-list format.
    if (is.data.frame(mark_features)) {
        required_mark_cols <- c("Feature", "Set")
        if (!all(required_mark_cols %in% colnames(mark_features))) {
            stop("Data-frame 'mark_features' must contain 'Feature' and ",
                "'Set' columns.")
        }
        if (anyNA(mark_features$Feature) || anyNA(mark_features$Set) ||
            any(as.character(mark_features$Feature) == "") ||
            any(as.character(mark_features$Set) == "")) {
            stop("Data-frame 'mark_features' must not contain missing or ",
                "empty 'Feature' or 'Set' values.")
        }
        mark_features <- split(
            as.character(mark_features$Feature),
            as.character(mark_features$Set)
        )
    }

    # Enrichment parameters validation and recycling
    if (!is.null(enrich_category)) {
        .check_character(enrich_category, "enrich_category")
    }
    if (!is.null(enrich_list)) {
        .check_list(enrich_list, "enrich_list")
    }
    if (!is.null(enrich_category) &&
        (anyNA(enrich_category) || any(enrich_category == ""))) {
        stop("'enrich_category' must not contain missing or empty strings.")
    }
    valid_categories <- c(
        if (is.null(enrich_list)) character(0) else names(enrich_list),
        "Features"
    )
    invalid_categories <- setdiff(enrich_category, valid_categories)
    if (length(invalid_categories) > 0L) {
        stop("'enrich_category' contains categories not present in ",
            "'enrich_list': ", paste(invalid_categories, collapse = ", "),
            ". Use 'Features' for marked feature annotations.")
    }
    has_feature_text <- is.character(enrich_category) &&
        !is.null(mark_features) && length(mark_features) > 0L &&
        "Features" %in% enrich_category
    if (!is.null(enrich_list)) {
        if (anyNA(enrich_category) || any(enrich_category == "")) {
            stop("'enrich_category' must not contain missing or empty strings.")
        }
        if (anyDuplicated(enrich_category)) {
            stop("'enrich_category' values must be unique.")
        }
        n_cats <- length(enrich_category)
        non_hub <- enrich_category != "Features"
        n_enrich <- sum(non_hub)
        enrich_rank_full <- rep(NA_character_, n_cats)
        enrich_top_full <- rep(NA_integer_, n_cats)
        enrich_threshold_full <- rep(NA_real_, n_cats)

        if (n_enrich > 0L) {
            .check_character(enrich_rank_by, "enrich_rank_by")
            rank_input <- enrich_rank_by
            if (length(rank_input) == n_cats) {
                rank_input <- rank_input[non_hub]
            }
            rank_non_hub <- .recycle_enrich_param(
                rank_input, "enrich_rank_by", n_enrich)
            if (anyNA(rank_non_hub) || any(rank_non_hub == "")) {
                stop("'enrich_rank_by' must not contain missing or empty ",
                    "strings for enrichment categories.")
            }
            enrich_rank_full[non_hub] <- rank_non_hub

            top_input <- enrich_top_n
            if (length(top_input) == n_cats) {
                top_input <- top_input[non_hub]
            }
            top_non_hub <- .recycle_enrich_param(
                top_input, "enrich_top_n", n_enrich)
            if (!is.numeric(top_non_hub) ||
                any(!is.finite(top_non_hub))) {
                stop("'enrich_top_n' must contain positive integers.")
            }
            for (v in top_non_hub) {
                .check_positive_int(v, "enrich_top_n")
            }
            enrich_top_full[non_hub] <- top_non_hub

            if (!is.null(enrich_threshold)) {
                threshold_input <- enrich_threshold
                if (length(threshold_input) == n_cats) {
                    threshold_input <- threshold_input[non_hub]
                }
                threshold_non_hub <- .recycle_enrich_param(
                    threshold_input, "enrich_threshold", n_enrich)
                for (v in threshold_non_hub[!is.na(threshold_non_hub)]) {
                    .check_pval(v, "enrich_threshold")
                }
                enrich_threshold_full[non_hub] <- threshold_non_hub
            }
        }
        enrich_rank_by <- enrich_rank_full
        enrich_top_n <- enrich_top_full
        enrich_threshold <- enrich_threshold_full
    } else {
        enrich_category <- if (has_feature_text) "Features" else character(0)
        enrich_rank_by <- if (has_feature_text) NA_character_ else character(0)
        enrich_top_n <- if (has_feature_text) NA_integer_ else integer(0)
        enrich_threshold <- if (has_feature_text) NA_real_ else numeric(0)
    }
    if (!is.null(mark_features)) {
        if (!is.character(mark_features) && !is.list(mark_features)) {
            stop("'mark_features' must be a character vector or named list.")
        }
        if (is.list(mark_features) && length(mark_features) > 0 &&
                is.null(base::names(mark_features))) {
            stop("'mark_features' list must have names.")
        }
    }
    if (!is.null(feature_colours)) {
        colour_names <- base::names(feature_colours)
        if (!is.list(mark_features) ||
            is.null(base::names(mark_features)) ||
            !is.character(feature_colours) ||
            is.null(colour_names)) {
            stop("'feature_colours' requires named feature sets and a ",
                "named character vector.")
        }
        colour_values <- c(colour_names, unname(feature_colours))
        if (anyNA(colour_values) || any(colour_values == "") ||
            anyDuplicated(colour_names)) {
            stop("'feature_colours' must have unique, non-empty names and ",
                "values.")
        }
        tryCatch(
            grDevices::col2rgb(unname(feature_colours)),
            error = function(e) stop(
                "'feature_colours' must contain valid colours."
            )
        )
        unknown_colours <- setdiff(
            colour_names, base::names(mark_features)
        )
        if (length(unknown_colours) > 0L) {
            stop("'feature_colours' contains unknown feature set(s): ",
                paste(unknown_colours, collapse = ", "))
        }
    }
    if (!is.null(save)) .check_character(save, "save")

    se_obj_list <- .as_se_list(se_obj_merged, arg = "se_obj_merged")
    n_datasets <- length(se_obj_list)
    nms <- base::names(se_obj_list)
    for (nm in nms) {
        .check_se_merged(se_obj_list[[nm]], arg = nm)
    }
    if (length(ylabel) != 1L) {
        stop("The 'ylabel' must be a single character string.")
    }
    if (!is.null(names)) {
        .check_character(names, "names")
        if (length(names) != n_datasets) {
            stop("'names' must have one entry per dataset.")
        }
    }
    disp <- if (n_datasets > 1L) {
        if (is.null(names)) nms else names
    } else {
        NULL
    }
    resolved <- vapply(nms, function(nm) {
        .match_assay(assay, se_obj_list[[nm]])
    }, character(1))
    if (!all(resolved == resolved[1])) {
        stop("The 'assay' must be the same across all datasets.")
    }
    assay <- unname(resolved[1])

    module$Module <- .module_labels(module$Module)
    # Exclude grey module (M0)
    module <- module[!module$Module %in% c(
        "M0", "grey", "gray"
    ), , drop = FALSE]
    if (n_datasets > 1L) {
        module <- .warn_unshared_features(module, se_obj_list, assay, disp)
    }
    # Reorder modules by size (largest first)
    mod_sizes <- table(module$Module)
    module$Module <- factor(module$Module,
        levels = base::names(sort(mod_sizes, decreasing = TRUE))
    )
    stopifnot(nrow(module) > 0)
    all_module_levels <- levels(module$Module)
    selected_mode <- !is.null(selected_modules)
    if (selected_mode) {
        selected_modules <- as.character(.module_labels(selected_modules))
        missing_modules <- setdiff(selected_modules, all_module_levels)
        if (length(missing_modules) > 0L) {
            stop("'selected_modules' contains unknown module(s): ",
                paste(missing_modules, collapse = ", "))
        }
        profile_modules <- all_module_levels[
            all_module_levels %in% selected_modules
        ]
    } else {
        profile_modules <- all_module_levels
    }

    data_module_longs <- list()
    groups <- list()
    sp_infos <- list()
    mats <- list()
    mods <- list()
    for (i in seq_len(n_datasets)) {
        data_module_longs[[i]] <-
            .plot_modules_input(module = module,
            se_obj_merged = se_obj_list[[i]], scale = scale, assay = assay)
        groups[[i]] <- data_module_longs[[i]]$Group |> levels()
        sp_infos[[i]] <-
            data_module_longs[[i]] |>
                dplyr::select(Group, Time) |>
                dplyr::distinct()
        data_module_wider <- data_module_longs[[i]] |>
            tidyr::pivot_wider(id_cols = c(Feature, Module),
                names_from = c(Group, Time),
                values_from = Abundance) |>
            dplyr::arrange(Module) |>
            tibble::column_to_rownames("Feature")
        mats[[i]] <- data_module_wider |>
            dplyr::select(-Module) |>
            as.matrix()
        mods[[i]] <- data_module_wider$Module
    }

    selected_profile_height <- NULL
    if (selected_mode && !is.null(heatmap_height)) {
        n_profile_modules <- length(profile_modules)
        available_profile_height <- heatmap_height / n_profile_modules
        if (available_profile_height <= 0) {
            stop("The heatmap is too short for selected profiles.")
        }
        selected_profile_height <- min(
            profile_width, available_profile_height
        )
    }

    profile_values <- lapply(data_module_longs, function(x) {
        x |>
            dplyr::group_by(Module, Group, Time) |>
            dplyr::summarise(
                Abundance = mean(Abundance, na.rm = TRUE),
                .groups = "drop"
            ) |>
            dplyr::pull(Abundance)
    }) |>
        lapply(function(values) values[is.finite(values)])
    get_profile_limits <- function(values) {
        if (length(values) == 0L) return(NULL)
        limits <- range(values)
        if (diff(limits) == 0) {
            padding <- max(abs(limits[1]) * 0.05, 0.5)
            limits <- limits + c(-padding, padding)
        }
        limits
    }
    if (profile_scale == "shared") {
        profile_y_limits <- rep(list(get_profile_limits(
            unlist(profile_values, use.names = FALSE)
        )), n_datasets)
    } else if (profile_scale == "dataset") {
        profile_y_limits <- lapply(profile_values, get_profile_limits)
    } else {
        profile_y_limits <- rep(list(NULL), n_datasets)
    }

    all_vals <- unlist(lapply(mats, function(m) as.vector(m)))
    min_val <- min(all_vals, na.rm = TRUE)
    max_val <- max(all_vals, na.rm = TRUE)
    col_fun <- circlize::colorRamp2(
        c(min_val, 0, max_val),
        c("#3C5488FF", "white", "#E64B35FF")
    )

    # Pre-process mark_features
    mark_state <- list(
        feat_col = character(0), set_col = character(0), lgd = NULL
    )
    if (!is.null(mark_features) && length(mark_features) > 0) {
        if (is.list(mark_features) && !is.null(base::names(mark_features))) {
            cats <- base::names(mark_features)
            n_cats <- length(cats)
            cat_cols <- stats::setNames(
                ggsci::pal_jco()(max(n_cats, 2))[seq_len(n_cats)], cats
            )
            if (!is.null(feature_colours)) {
                cat_cols[base::names(feature_colours)] <-
                    unname(feature_colours)
            }
            all_feats <- unlist(mark_features, use.names = FALSE)
            feat_col <- rep(unname(cat_cols), times = lengths(mark_features))
            base::names(feat_col) <- unlist(mark_features, use.names = FALSE)
            # De-duplicate: use first category if a feature appears in multiple
            feat_col <- feat_col[!duplicated(base::names(feat_col))]
            mark_state$feat_col <- feat_col
            mark_state$set_col <- cat_cols
            mark_state$lgd <- ComplexHeatmap::Legend(
                title = "Feature set", labels = cats,
                legend_gp = grid::gpar(fill = unname(cat_cols)),
                labels_gp = grid::gpar(fontsize = fontsize),
                title_gp = grid::gpar(fontsize = fontsize)
            )
        } else {
            all_feats <- as.character(mark_features)
            mark_state$feat_col <- stats::setNames(
                rep("black", length(all_feats)), all_feats
            )
            mark_state$set_col <- stats::setNames(
                "black", "Features"
            )
        }
    }

    enrich_legend_colours <- new.env(parent = emptyenv())

    # Use a shared mapping so the same group has the same colour across
    # datasets, including when group sets or their ordering differ.
    all_groups <- unique(unlist(groups, use.names = FALSE))
    group_colours <- get_custom_palette(all_groups) |>
        stats::setNames(all_groups)

    group_name_height <- grid::unit(max(4, 0.75 * fontsize), "mm")
    group_colour_height <- grid::unit(max(2, 0.4 * fontsize), "mm")
    make_group_name_anno <- function(group_levels) {
        ComplexHeatmap::anno_block(
            labels = group_levels,
            labels_gp = grid::gpar(fontsize = fontsize,
                fontface = "bold", fontfamily = fontfamily),
            gp = grid::gpar(fill = NA, col = NA),
            height = group_name_height
        )
    }
    make_group_colour_anno <- function(group_levels) {
        ComplexHeatmap::anno_block(
            gp = grid::gpar(
                fill = unname(group_colours[group_levels]), col = NA
            ),
            height = group_colour_height
        )
    }
    make_group_top_anno <- function(group_levels) {
        if (group_annotation == "none") {
            return(NULL)
        }
        if (group_annotation == "name") {
            return(ComplexHeatmap::HeatmapAnnotation(
                Group = make_group_name_anno(group_levels),
                show_annotation_name = FALSE
            ))
        }
        if (group_annotation == "colour") {
            return(ComplexHeatmap::HeatmapAnnotation(
                Group = make_group_colour_anno(group_levels),
                show_annotation_name = FALSE
            ))
        }
        ComplexHeatmap::HeatmapAnnotation(
            Group_name = make_group_name_anno(group_levels),
            Group_colour = make_group_colour_anno(group_levels),
            show_annotation_name = FALSE
        )
    }

    heatmaps <- lapply(seq_len(n_datasets), function(i) {
        # Pretty time breaks for profile x-axis
        time_vals <- as.numeric(levels(data_module_longs[[i]]$Time))
        t_range <- range(time_vals)
        time_breaks <- pretty(t_range, n = 4)
        time_breaks <- time_breaks[time_breaks >= t_range[1] &
            time_breaks <= t_range[2]]
        t_range_exp <- t_range + c(-1, 1) * 0.05 * diff(t_range)
        time_breaks_norm <- (time_breaks - t_range_exp[1]) / diff(t_range_exp)
        last_mod <- utils::tail(profile_modules, 1)
        make_profile_panel_fun <- function(profile_grobs, isolate = FALSE) {
            function(index, nm) {
                if (isolate) {
                    grid::pushViewport(grid::viewport())
                    on.exit(grid::popViewport(), add = TRUE)
                }
                g <- profile_grobs[[nm]]
                g <- grid::grid.grabExpr(grid::grid.draw(g))
                grid::grid.rect(gp = grid::gpar(
                    fill = "transparent", col = "grey50", lwd = 1
                ))
                grid::grid.draw(g)
                if (nm == last_mod) {
                    grid::grid.xaxis(at = time_breaks_norm,
                        label = FALSE,
                        gp = grid::gpar(col = "grey50",
                            fontsize = fontsize,
                            fontfamily = fontfamily))
                    grid::grid.text(as.character(time_breaks),
                        x = grid::unit(time_breaks_norm, "npc"),
                        y = grid::unit(-1.2, "lines"),
                        just = c("center", "top"),
                        gp = grid::gpar(col = "black",
                            fontsize = fontsize,
                            fontfamily = fontfamily))
                    grid::grid.text("Time",
                        x = grid::unit(0.5, "npc"),
                        y = grid::unit(-2.4, "lines"),
                        just = c("center", "top"),
                        gp = grid::gpar(fontsize = fontsize,
                            fontfamily = fontfamily))
                }
            }
        }
        ggplot2_profile_list <- list()
        for (m in levels(data_module_longs[[i]]$Module)) {
            profile_plot <- data_module_longs[[i]] |>
                dplyr::filter(Module == m) |>
                dplyr::mutate(Time = as.numeric(as.character(Time))) |>
                ggplot(aes(x = Time, y = Abundance, group = Feature)) +
                stat_summary(aes(group = Group, color = Group),
                    fun = mean, geom = "line", linewidth = 1, na.rm = TRUE
                ) +
                scale_x_continuous(expand = c(0.05, 0)) + # default
                scale_color_manual(values = group_colours[groups[[i]]]) +
                labs(x = NULL, y = NULL, title = NULL) +
                theme_custom(panel_border = TRUE,
                    base_family = fontfamily) +
                theme(
                    text = element_text(family = fontfamily,
                        size = fontsize),
                    axis.text = element_blank(),
                    axis.title = element_blank(),
                    axis.ticks = element_blank(),
                    panel.border = element_rect(
                        fill = NA, colour = "grey50", linewidth = 1
                    ),
                    legend.position = "none",
                    plot.margin = margin(0, 0, 0, 0, unit = "cm")
                )
            if (!is.null(profile_y_limits[[i]])) {
                profile_plot <- profile_plot +
                    coord_cartesian(ylim = profile_y_limits[[i]])
            }
            ggplot2_profile_list[[m]] <- profile_plot
        }

        profile_align_to <- if (selected_mode) {
            align_to <- lapply(profile_modules, function(mod) {
                which(as.character(mods[[i]]) == mod)
            })
            base::names(align_to) <- profile_modules
            align_to
        } else {
            mods[[i]]
        }
        if (!selected_mode) {
            profile_size <- 1 / length(unique(mods[[i]]))
            profile_gap <- unit(0, "cm")
        }

        selected_profile_spec <- NULL
        text_panel_fun <- NULL
        text_row_groups <- NULL
        text_heights_mm <- numeric(0)
        textbox_width_mm <- 0
        combined_anno <- NULL
        if (selected_mode) {
            profile_grobs <- lapply(profile_modules, function(mod) {
                ggplot2_profile_list[[mod]]
            })
            base::names(profile_grobs) <- profile_modules
            profile_panel_fun <- make_profile_panel_fun(
                profile_grobs = profile_grobs
            )
            selected_profile_spec <- list(
                row_groups = profile_align_to,
                panel_fun = profile_panel_fun
            )
            anno_ggplot2_profile <- NULL
        } else {
            anno_ggplot2_profile <- ComplexHeatmap::anno_zoom(
                align_to = profile_align_to,
                which = "row",
                panel_fun = make_profile_panel_fun(
                    profile_grobs = ggplot2_profile_list,
                    isolate = TRUE
                ),
                size = profile_size,
                gap = profile_gap,
                width = unit(profile_width, "cm"),
                side = "right",
                link_gp = grid::gpar(
                    fill = "grey98", col = "grey50", lwd = 1
                ),
                link_width = unit(profile_link_width, "cm")
            )
        }

        # ---- Multi-category annotation (rightmost dataset only) ----
        if (i == n_datasets &&
            (!is.null(enrich_list) || has_feature_text)) {
            mod_levels <- if (selected_mode) {
                profile_modules
            } else {
                levels(data_module_longs[[i]]$Module)
            }

            # Build combined multi-column textbox from all enrichment categories
            enr_cats <- if (is.null(enrich_list)) {
                character(0)
            } else {
                setdiff(enrich_category, "Features")
            }
            n_enr <- length(enr_cats)

            # First pass: collect per-category per-module term data
            cat_terms <- list() # cate -> list of module data.frames
            for (ci in seq_len(n_enr)) {
                cate <- enr_cats[ci]
                if (is.null(enrich_list[[cate]])) next
                param_i <- match(cate, enrich_category)
                top_n <- enrich_top_n[param_i]
                rank_by <- enrich_rank_by[param_i]
                enr_df <- as.data.frame(enrich_list[[cate]])
                if (!rank_by %in% colnames(enr_df)) {
                    stop("'enrich_rank_by' column '", rank_by,
                        "' not found in enrich_list category '", cate, "'.")
                }
                if (!is.numeric(enr_df[[rank_by]])) {
                    stop("'enrich_rank_by' column '", rank_by,
                        "' in enrich_list category '", cate,
                        "' must be numeric.")
                }

                termanno <- enr_df |>
                    as.data.frame() |>
                    dplyr::mutate(
                        id = as.character(.module_labels(Cluster))
                    ) |>
                    dplyr::filter(id %in% mod_levels) |>
                    dplyr::group_by(id) |>
                    dplyr::slice_min(
                        order_by = !!sym(rank_by),
                        n = top_n, with_ties = FALSE
                    ) |>
                    dplyr::mutate(term = Description) |>
                    dplyr::ungroup()

                # Colour-coding
                if ("col" %in% colnames(termanno)) {
                    termanno$col <- termanno[["col"]]
                } else if (!is.na(enrich_threshold[param_i]) &&
                    rank_by %in% colnames(termanno)) {
                    termanno$col <- ifelse(
                        termanno[[rank_by]] < enrich_threshold[param_i],
                        "black", "grey50")
                } else {
                    termanno$col <- "black"
                }
                termanno$fontsize <- fontsize
                termanno$fontfamily <- fontfamily
                termanno <- termanno |>
                    dplyr::select(id, term, col, fontsize, fontfamily)

                if (nrow(termanno) == 0) next

                enrich_legend_colours[[cate]] <- intersect(
                    c("black", "grey50"), unique(termanno$col)
                )

                term_list <- lapply(mod_levels, function(mod) {
                    tmp <- termanno[termanno$id == mod, ]
                    if (nrow(tmp) == 0) {
                        data.frame(
                            text = character(0), col = character(0),
                            fontsize = numeric(0), fontfamily = character(0),
                            stringsAsFactors = FALSE
                        )
                    } else {
                        data.frame(
                            text = tmp$term, col = tmp$col,
                            fontsize = tmp$fontsize,
                            fontfamily = tmp$fontfamily,
                            stringsAsFactors = FALSE
                        )
                    }
                })
                base::names(term_list) <- mod_levels
                cat_terms[[cate]] <- term_list
            }

            # Feature annotations as a column in the multi-column textbox
            if (has_feature_text) {
                hres <- .build_marked_features_anno(
                    mark_features, module, mark_state,
                    fontsize, mod_levels
                )
                if (!is.null(hres)) {
                    hres <- lapply(hres, function(df) {
                        if (nrow(df) > 0L) {
                            df$fontfamily <- fontfamily
                        }
                        df
                    })
                    cat_terms[["Features"]] <- hres
                } else {
                    message("Features: no marked features found ",
                        "in module data")
                }
            }

            # Drop modules with no terms in any category
            has_terms <- vapply(mod_levels, function(mod) {
                any(vapply(cat_terms, function(ct) {
                    mod %in% base::names(ct) && nrow(ct[[mod]]) > 0
                }, logical(1)))
            }, logical(1))
            if (!selected_mode && any(has_terms)) {
                cat_terms <- lapply(cat_terms, function(ct)
                    ct[mod_levels[has_terms]])
            }

            # Build a single multi-column textbox from all categories
            if (length(cat_terms) > 0) {
                text <- .reformat_cat_terms(cat_terms)
                textbox_args <- list(
                    word_wrap = !is.null(enrich_wrap_width),
                    max_width = if (is.null(enrich_wrap_width)) {
                        NULL
                    } else {
                        unit(enrich_wrap_width, "cm")
                    },
                    line_space = unit(0.15 * fontsize, "pt"),
                    term_space = unit(0.75 * fontsize, "pt"),
                    text_padding = unit(0.65 * fontsize, "pt")
                )
                if (!selected_mode) {
                    align_to <- factor(
                        rep(mod_levels,
                            each = round(nrow(mats[[i]]) / length(mod_levels))
                        ),
                        levels = mod_levels
                    )
                    combined_anno <- do.call(
                        .anno_multicol_textbox,
                        c(
                            list(
                                align_to = align_to,
                                text = text,
                                background_gp = grid::gpar(
                                    fill = "grey98", col = "grey50"
                                ),
                                by = "anno_link",
                                side = "right",
                                show_headers = TRUE,
                                header_gp = grid::gpar(
                                    fontsize = fontsize,
                                    fontface = "bold",
                                    fontfamily = fontfamily
                                )
                            ),
                            if (is.null(enrich_wrap_width)) {
                                list()
                            } else {
                                textbox_args
                            }
                        )
                    )
                } else {
                    text_row_groups <- profile_align_to[names(text)]
                    base::names(text_row_groups) <- names(text)

                    dev_was_open <- length(grDevices::dev.list()) > 0L
                    if (!dev_was_open) {
                        tmp <- tempfile(fileext = ".pdf")
                        grDevices::pdf(tmp)
                    }
                    textbox_layout <- tryCatch(
                        .multicol_textbox_layout(
                            text, textbox_args = textbox_args
                        ),
                        finally = if (!dev_was_open) {
                            grDevices::dev.off()
                            unlink(tmp)
                        }
                    )
                    textbox_width_mm <- textbox_layout$width_mm
                    text_heights_mm <- textbox_layout$heights_mm
                    text_panel_fun <- .multicol_textbox_panel_fun(
                        textbox_layout,
                        background_gp = grid::gpar(
                            fill = "grey98", col = "grey50", lty = 1,
                            lwd = 1
                        ),
                        show_headers = TRUE,
                        header_gp = grid::gpar(
                            fontsize = fontsize, fontface = "bold",
                            fontfamily = fontfamily
                        )
                    )

                }

            }
        }

        if (!is.null(selected_profile_spec)) {
            selected_anno <- .make_selected_anno_zoom(
                profile_row_groups = selected_profile_spec$row_groups,
                profile_width_mm = profile_width * 10,
                link_width_mm = profile_link_width * 10,
                heatmap_height_mm = if (is.null(heatmap_height)) {
                    NULL
                } else {
                    heatmap_height * 10
                },
                profile_height_mm = if (is.null(heatmap_height)) {
                    NULL
                } else {
                    selected_profile_height * 10
                },
                profile_panel_fun = selected_profile_spec$panel_fun,
                text_panel_fun = text_panel_fun,
                text_row_groups = text_row_groups,
                text_heights_mm = text_heights_mm,
                textbox_width_mm = textbox_width_mm
            )
            right_ann <- ComplexHeatmap::rowAnnotation(
                Selected = selected_anno
            )
        } else if (!is.null(combined_anno)) {
            right_ann <- ComplexHeatmap::rowAnnotation(
                profile = anno_ggplot2_profile,
                Enrichment = combined_anno
            )
        } else {
            right_ann <- ComplexHeatmap::rowAnnotation(
                profile = anno_ggplot2_profile
            )
        }

        # ---- anno_mark for key features (leftmost dataset only) ----
        # Skip anno_mark if feature annotations are already shown as textbox
        left_ann <- NULL
        if (i == 1 && length(mark_state$feat_col) > 0 && !has_feature_text) {
            feats_present <- intersect(
                base::names(mark_state$feat_col), rownames(mats[[i]])
            )
            if (length(feats_present) > 0) {
                at_idx <- match(feats_present, rownames(mats[[i]]))
                cols <- mark_state$feat_col[feats_present]
                left_ann <- ComplexHeatmap::rowAnnotation(
                    mark = ComplexHeatmap::anno_mark(
                        at = at_idx,
                        labels = feats_present,
                        side = "left",
                        labels_gp = grid::gpar(
                            col = cols, fontsize = fontsize,
                            fontfamily = fontfamily
                        ),
                        link_gp = grid::gpar(col = "grey50", lwd = 1),
                        padding = unit(2, "mm")
                    )
                )
            }
        }

        ht_args <- list(
            mats[[i]],
            name = disp[i], col = col_fun,
            cluster_rows = FALSE, cluster_columns = FALSE,
            show_row_names = FALSE, show_column_names = FALSE,
            row_split = mods[[i]],
            row_gap = unit(0, "mm"), row_title_rot = 0,
            row_title_gp = grid::gpar(fontsize = fontsize,
                fontface = "bold", fontfamily = fontfamily),
            layer_fun = function(j, i, x, y, w, h, fill, slice_r, slice_c) {
                grid::grid.rect(gp = grid::gpar(
                    fill = "transparent", col = "grey50", lwd = 1
                ))
            }, # draw frames around the modules
            column_split = sp_infos[[i]]$Group,
            column_title = if (n_datasets > 1) disp[i] else NULL,
            column_title_gp =
                grid::gpar(fontsize =
                    if (n_datasets > 1) (fontsize + 2) else fontsize,
                    fontface = if (n_datasets > 1) "bold" else "plain",
                    fontfamily = fontfamily),
            column_title_side = "top",
            column_names_gp = grid::gpar(fontsize = fontsize,
                fontfamily = fontfamily),
            top_annotation = if (group_annotation != "none") {
                make_group_top_anno(groups[[i]])
            } else {
                NULL
            },
            right_annotation = right_ann,
            use_raster = FALSE,
            show_heatmap_legend = i == n_datasets,
            heatmap_legend_param = list(
                title = ylabel,
                labels_gp = grid::gpar(fontsize = fontsize,
                    fontfamily = fontfamily),
                title_gp = grid::gpar(fontsize = fontsize,
                    fontface = "bold", fontfamily = fontfamily)
            ),
            width = grid::unit(heatmap_width, "cm"),
            height = if (is.null(heatmap_height)) NULL
                else grid::unit(heatmap_height, "cm")
        )
        if (!is.null(left_ann)) ht_args$left_annotation <- left_ann
        row_title_col <- if (selected_mode) {
            ifelse(
                levels(data_module_longs[[i]]$Module) %in%
                    profile_modules,
                "black", "grey65"
            )
        } else {
            "black"
        }
        ht_args$row_title_gp <- grid::gpar(
            fontsize = fontsize, fontface = "bold",
            fontfamily = fontfamily, col = row_title_col
        )
        do.call(ComplexHeatmap::Heatmap, ht_args)
    })

    anno_ggplot2_profile_lgd <- ComplexHeatmap::Legend(
        title = "Group", labels = all_groups,
        legend_gp = grid::gpar(
            fill = unname(group_colours)
        ),
        labels_gp = grid::gpar(fontsize = fontsize,
            fontfamily = fontfamily),
        title_gp = grid::gpar(fontsize = fontsize,
            fontface = "bold", fontfamily = fontfamily)
    )

    ht_list <- heatmaps[[1]]
    if (n_datasets > 1) {
        for (i in 2:n_datasets) {
            ht_list <- ht_list + heatmaps[[i]]
        }
    }

    build_legends <- function() {
        lgd_list <- list(anno_ggplot2_profile_lgd)
        n_cats <- length(enrich_category)
        for (ci in seq_len(n_cats)) {
            cate <- enrich_category[ci]
            pt <- enrich_threshold[ci]
            if (cate == "Features") {
                if (!is.null(mark_features) && length(mark_features) > 0) {
                    if (is.list(mark_features)) {
                        hub_cols <- mark_state$set_col[
                            base::names(mark_features)
                        ]
                        hub_lgd <- ComplexHeatmap::Legend(
                            title = "Features",
                            labels = base::names(mark_features),
                            type = "lines",
                            legend_gp = grid::gpar(col = unname(hub_cols),
                                lwd = 2),
                            labels_gp = grid::gpar(fontsize = fontsize,
                                fontfamily = fontfamily),
                            title_gp = grid::gpar(fontsize = fontsize,
                                fontface = "bold", fontfamily = fontfamily)
                        )
                        lgd_list <- c(lgd_list, list(hub_lgd))
                    } else {
                        feature_colour <- mark_state$set_col[["Features"]]
                        hub_lgd <- ComplexHeatmap::Legend(
                            title = "Features",
                            labels = "Features",
                            type = "lines",
                            legend_gp = grid::gpar(
                                col = feature_colour, lwd = 2
                            ),
                            labels_gp = grid::gpar(fontsize = fontsize,
                                fontfamily = fontfamily),
                            title_gp = grid::gpar(fontsize = fontsize,
                                fontface = "bold", fontfamily = fontfamily)
                        )
                        lgd_list <- c(lgd_list, list(hub_lgd))
                    }
                }
            } else if (!is.null(enrich_list) &&
                !is.null(enrich_list[[cate]]) &&
                !is.na(pt) && length(enrich_legend_colours[[cate]]) > 0L) {
                legend_colours <- enrich_legend_colours[[cate]]
                legend_labels <- c(
                    black = sprintf("%s < %.2f", enrich_rank_by[ci], pt),
                    grey50 = sprintf("%s >= %.2f", enrich_rank_by[ci], pt)
                )[legend_colours]
                go_lgd <- ComplexHeatmap::Legend(
                    title = cate,
                    labels = unname(legend_labels),
                    type = "lines",
                    legend_gp = grid::gpar(col = legend_colours, lwd = 2),
                    labels_gp = grid::gpar(fontsize = fontsize,
                        fontfamily = fontfamily),
                    title_gp = grid::gpar(fontsize = fontsize,
                        fontface = "bold", fontfamily = fontfamily)
                )
                lgd_list <- c(lgd_list, list(go_lgd))
            }
        }
        if (!is.null(mark_state$lgd) && !has_feature_text) {
            lgd_list <- c(lgd_list, list(mark_state$lgd))
        }
        lgd_list
    }

    lgd_list <- build_legends()

    draw_ht <- function() {

        pad_bottom <- 10
        ComplexHeatmap::draw(ht_list,
            heatmap_legend_list = lgd_list,
            heatmap_legend_side = legend_position, merge_legends = TRUE,
            # top, right, bottom, left
            padding = unit(c(2, 2, pad_bottom, 2), "mm")
        )
    }

    if (!is.null(save)) {
        if (suffix != "") {
            suffix <- paste0(suffix, "_")
        }

        prefix <- if (n_datasets > 1L) "WGCNA_h_multi_" else "WGCNA_h_"
        file_base <- file.path(save, paste0(prefix, suffix, Sys.Date()))
        for (ext in device) {
            if (ext == "pdf") {
                grDevices::cairo_pdf(
                    file = paste0(file_base, ".", ext),
                    width = width / 2.54,
                    height = height / 2.54,
                    family = fontfamily
                )
            } else if (ext == "svg") {
                grDevices::svg(
                    filename = paste0(file_base, ".", ext),
                    width = width / 2.54,
                    height = height / 2.54,
                    family = fontfamily
                )
            } else if (ext == "png") {
                grDevices::png(paste0(file_base, ".", ext),
                    width = width, height = height,
                    units = "cm", res = res)
            } else if (ext == "tiff") {
                grDevices::tiff(paste0(file_base, ".", ext),
                    width = width, height = height,
                    units = "cm", res = res)
            } else if (ext == "jpeg") {
                grDevices::jpeg(paste0(file_base, ".", ext),
                    width = width, height = height,
                    units = "cm", res = res)
            } else {
                message("Unsupported device: ", ext, ", skipping.")
                next
            }
            draw_ht()
            grDevices::dev.off()
        }
    }

    draw_ht()
    result <- if (return_legends) {
        list(heatmap = ht_list, legends = lgd_list)
    } else {
        ht_list
    }
    invisible(result)
}
