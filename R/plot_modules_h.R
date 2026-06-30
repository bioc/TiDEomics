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
#' in `.anno_multicol_textbox`. Accepts a named list (per-category colours
#' from `ggsci::pal_jco()`) or a character vector (all black). Colours
#' are taken from `mark_state$feat_col`.
#'
#' @param mark_features A named list of character vectors, mapping
#'   category names to feature IDs, or a plain character vector of
#'   feature IDs (all shown in black).
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
        feats <- as.character(mark_features)
    } else {
        # Named list: map each feature to its category
        feats <- stats::setNames(
            rep(names(mark_features), lengths(mark_features)),
            unlist(mark_features, use.names = FALSE)
        )
    }

    marked_in_mod <- intersect(names(feats), module_df$Feature)
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
#' feature at each time point in each group (replicates merged). The colData of
#' the object should contain columns "Sample", "Group", and "Time". The object
#' can be produced by `split_groups()`, `merge_replicates()` and
#' `merge_groups()`.
#' @param assay The assay index in the SummarizedExperiment object to use
#' (default is 2, time 0 normalised data)
#' @param scale Whether to scale the data (z-score) across samples for each
#' feature (default is TRUE)
#' @param ylabel Y axis label prefix (default is "Log2 abundance")
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
#'   categories (e.g. `c("BP", "CC", "Hub features")`).
#'   When `"Hub features"` is included, a textbox is auto-built from
#'   `mark_features` with colors matching `anno_mark`.
#' @param enrich_rank_by Column to rank enrichment terms. A single string
#'   (default: `"p.adjust"`) is recycled for all categories in
#'   `enrich_category`. Supply a character vector of matching length for
#'   per-category control. Also used for colour-coding: terms with values
#'   below `enrich_p_threshold` are shown in black, terms above in grey.
#' @param enrich_top_n Top N enrichment terms per module. A single integer
#'   (default: 3) is recycled for all categories in `enrich_category`.
#'   Supply an integer vector of matching length for per-category control.
#' @param enrich_p_threshold P-value / adjusted p-value threshold for
#'   colour-coding enrichment terms. A single value (default: 0.05) is
#'   recycled for all categories; supply a vector of matching length for
#'   per-category control. Terms below threshold are drawn in black,
#'   terms above in grey70. Set to `NULL` to skip colour-coding (all
#'   terms black). Set to `NA` for categories using a pre-assigned
#'   `col` column (e.g. `"Hub features"`).
#' @param fontsize Base font size (default: 8)
#' @param heatmap_width Width of the heatmap body in cm (default: 8)
#' @param heatmap_height Height of the heatmap body in cm (default: 8)
#' @param mark_features Features to highlight. Accepts two forms:
#'     - A **character vector**: all marked in black on the left.
#'     - A **named list** of character vectors, e.g.
#'         `list("Hub 1" = c("gene1"), "Hub 2" = c("gene2"))`.
#'         Names become categories with auto-assigned colours (via
#'         `ggsci::pal_jco()`).
#'     Auto-routed: shown as `anno_mark` (left
#'     side, connecting lines) unless `"Hub features"` is present in
#'     `enrich_category`, in which case an `anno_textbox` (right side)
#'     is built instead, using the same colours.
#'     Set to `NULL` (default) to skip marking.
#' @param width Width of the saved image (default is 16 (cm))
#' @param height Height of the saved image (default is 12 (cm))
#' @param res Resolution of the saved image (except pdf format)
#' (default: 300 (ppi))
#' @param device Image file format(s) for saving. Can be a character
#'   vector with one or more of `"png"`, `"pdf"`, `"tiff"`, `"jpeg"`
#'   (default: `"png"`)
#' @param save Directory to save plot, NULL (default) for no saving.
#' @param suffix Suffix for the saved image file name (default: "")
#' @import ggplot2
#' @import SummarizedExperiment
#'
#' @returns A combined plot of heatmap, mean expression profile, and enrichment
#' terms for each WGCNA module
#' @export
#' @examples
#' library(org.Mm.eg.db)
#'
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
#' # set cutoff to 1 to show all results for demonstration
#' example_go_list = enrichGO_list(example_module, OrgDb = org.Mm.eg.db,
#'     universe = example_module$Feature,
#'     pvalueCutoff = 1, qvalueCutoff = 1,
#'     category = "BP", simplify = FALSE)
#'
#' plot_modules_h(example_module |> dplyr::filter(Module != '0'),
#'     example_obj_merged, scale = TRUE,
#'     ylabel = "Z-score of log2 expression",
#'     enrich_list = example_go_list$all, enrich_category = "BP",
#'     heatmap_width = 6, heatmap_height = 4)
#' @references https://github.com/junjunlab/ClusterGVis
plot_modules_h <- function(
    module,
    se_obj_merged,
    scale = TRUE,
    assay = 2,
    ylabel = "Log2 abundance",
    profile_width = 3,
    profile_link_width = 1,
    enrich_list = NULL,
    enrich_category = "BP",
    enrich_rank_by = "p.adjust",
    enrich_top_n = 3,
    enrich_p_threshold = 0.05,
    fontsize = 8,
    heatmap_width = 8,
    heatmap_height = 8,
    mark_features = NULL,
    width = 16,
    height = 12,
    res = 300,
    suffix = "",
    device = "png",
    save = NULL
) {
    .check_df(module, "module")
    .check_character(enrich_category, "enrich_category")
    .check_logical(scale, "scale")
    .check_positive(profile_width, "profile_width")
    .check_positive(profile_link_width, "profile_link_width")
    .check_positive_int(enrich_top_n, "enrich_top_n")
    .check_positive(fontsize, "fontsize")
    .check_positive(heatmap_width, "heatmap_width")
    .check_positive(heatmap_height, "heatmap_height")
    .check_positive(width, "width")
    .check_positive(height, "height")
    .check_positive_int(res, "res")
    .check_character(ylabel, "ylabel")
    .check_character(enrich_rank_by, "enrich_rank_by")
    .check_character(suffix, "suffix")
    .check_character(device, "device")
    .check_se_merged(se_obj_merged, "se_obj_merged")
    assay <- .match_assay(assay, se_obj_merged)
    if (!is.null(enrich_list))
        .check_list(enrich_list, "enrich_list")
    if (!is.null(enrich_p_threshold)) {
        valid_vals <- enrich_p_threshold[!is.na(enrich_p_threshold)]
        for (v in valid_vals) .check_pval(v, "enrich_p_threshold")
    }
    if (!is.null(mark_features)) {
        if (!is.character(mark_features) && !is.list(mark_features)) {
            stop("'mark_features' must be a character vector or named list.")
        }
        if (is.list(mark_features) && is.null(names(mark_features))) {
            stop("'mark_features' list must have names.")
        }
    }
    if (!is.null(save)) .check_character(save, "save")

    module$Module <- .module_labels(module$Module)
    # Exclude grey module (M0)
    module <- module[!module$Module %in% c(
        "M0", "grey", "gray"
    ), , drop = FALSE]
    # Reorder modules by size (largest first)
    mod_sizes <- table(module$Module)
    module$Module <- factor(module$Module,
        levels = names(sort(mod_sizes, decreasing = TRUE))
    )
    stopifnot(nrow(module) > 0)
    data_module_long <- .plot_modules_input(module = module,
        se_obj_merged = se_obj_merged, scale = scale, assay = assay)

    n_modules <- length(unique(data_module_long$Module))

    stopifnot(n_modules == length(levels(data_module_long$Module)))

    groups <- data_module_long$Group |> levels()

    sp_info <- data_module_long |>
        dplyr::select(Group, Time) |>
        dplyr::distinct()

    data_module_wider <- data_module_long |>
        tidyr::pivot_wider(id_cols = c(Feature, Module),
            names_from = c(Group, Time), values_from = Abundance) |>
        tibble::column_to_rownames("Feature") |>
        as.data.frame() |>
        dplyr::arrange(Module)
    mat <- data_module_wider |>
        dplyr::select(-Module) |>
        as.matrix()

    min_val <- min(mat, na.rm = TRUE)
    max_val <- max(mat, na.rm = TRUE)
    col_fun <- circlize::colorRamp2(
        c(min_val, 0, max_val),
        c("#3C5488FF", "white", "#E64B35FF")
    )

    # Pre-process mark_features
    mark_state <- list(feat_col = character(0), lgd = NULL)
    if (!is.null(mark_features) && length(mark_features) > 0) {
        if (is.list(mark_features) && !is.null(names(mark_features))) {
            cats <- names(mark_features)
            n_cats <- length(cats)
            cat_cols <- stats::setNames(
                ggsci::pal_jco()(max(n_cats, 2))[seq_len(n_cats)], cats
            )
            all_feats <- unlist(mark_features, use.names = FALSE)
            feat_col <- rep(unname(cat_cols), times = lengths(mark_features))
            names(feat_col) <- unlist(mark_features, use.names = FALSE)
            # De-duplicate: first category wins if a feature appears in multiple
            feat_col <- feat_col[!duplicated(names(feat_col))]
            mark_state$feat_col <- feat_col
            mark_state$lgd <- ComplexHeatmap::Legend(
                title = "Hub type", labels = cats,
                legend_gp = grid::gpar(fill = unname(cat_cols)),
                labels_gp = grid::gpar(fontsize = fontsize),
                title_gp = grid::gpar(fontsize = fontsize)
            )
        } else {
            all_feats <- as.character(mark_features)
            mark_state$feat_col <- stats::setNames(
                rep("black", length(all_feats)), all_feats
            )
        }
    }
    ggplot2_profile_list <- list()
    for (i in levels(data_module_long$Module)) {
        ggplot2_profile_list[[i]] <- ggplot(
            data_module_long |> dplyr::filter(Module == i),
            aes(x = Time, y = Abundance, group = Feature)
        ) +
            stat_summary(aes(group = Group, color = Group),
                fun = mean, geom = "line", linewidth = 1
            ) +
            scale_color_manual(values = get_custom_palette(groups)) +
            labs(x = NULL, y = NULL, title = NULL) +
            theme_custom(panel_border = TRUE) +
            theme(
                text = element_text(size = fontsize),
                axis.text = element_blank(),
                axis.title = element_blank(),
                axis.ticks = element_blank(),
                legend.position = "none", # duplicated legends for each module
                plot.margin = margin(0, 0, 0, 0, unit = "cm")
            )
    }

    anno_ggplot2_profile <- ComplexHeatmap::anno_zoom(
        align_to = data_module_wider$Module,
        which = "row",
        panel_fun = function(index, nm) {
            g <- ggplot2_profile_list[[nm]]
            g <- grid::grid.grabExpr(grid::grid.draw(g))
            grid::pushViewport(grid::viewport())
            grid::grid.rect()
            grid::grid.draw(g)
            grid::popViewport()
        },
        size = 1 / n_modules,
        gap = unit(0, "cm"),
        width = unit(profile_width, "cm"),
        side = "right",
        link_gp = grid::gpar(fill = "grey95", col = "grey50"),
        link_width = unit(profile_link_width, "cm")
    )

    right_ann <- ComplexHeatmap::rowAnnotation(profile = anno_ggplot2_profile)

    # ---- Multi-category annotation ----
    if (!is.null(enrich_list)) {
        n_cats <- length(enrich_category)
        .check_len <- function(val, name, n) {
            if (length(val) > 1 && length(val) != n) {
                message(
                    "  ", name, " has length ", length(val),
                    " but enrich_category has ", n, " categories. ",
                    "First ", n, " values used."
                )
            }
        }
        .check_len(enrich_top_n, "enrich_top_n", n_cats)
        .check_len(enrich_rank_by, "enrich_rank_by", n_cats)
        .check_len(enrich_p_threshold, "enrich_p_threshold", n_cats)
        enrich_top_n <- rep(enrich_top_n, length.out = n_cats)
        enrich_rank_by <- rep(enrich_rank_by, length.out = n_cats)
        if (is.null(enrich_p_threshold)) {
            enrich_p_threshold <- rep(NA_real_, n_cats)
        } else {
            enrich_p_threshold <- rep(enrich_p_threshold,
                length.out = n_cats
            )
        }
        mod_levels <- levels(data_module_long$Module)
        mat_nrow <- nrow(mat)

        # Build combined multi-column textbox from all enrichment categories
        enr_cats <- setdiff(enrich_category, "Hub features")
        n_enr <- length(enr_cats)

        # First pass: collect per-category per-module term data
        cat_terms <- list() # cate -> list of module data.frames
        for (ci in seq_len(n_enr)) {
            cate <- enr_cats[ci]
            if (is.null(enrich_list[[cate]])) next
            top_n <- enrich_top_n[ci]
            rank_by <- enrich_rank_by[ci]
            p_thresh <- enrich_p_threshold[ci]

            termanno <- enrich_list[[cate]] |>
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
            } else if (!is.na(p_thresh) && !is.null(p_thresh) &&
                rank_by %in% colnames(termanno)) {
                termanno$col <- ifelse(
                    termanno[[rank_by]] < p_thresh, "black", "grey70")
            } else {
                termanno$col <- "black"
            }
            termanno$fontsize <- fontsize
            termanno <- termanno |>
                dplyr::select(id, term, col, fontsize)

            if (nrow(termanno) == 0) next

            term_list <- lapply(mod_levels, function(mod) {
                tmp <- termanno[termanno$id == mod, ]
                if (nrow(tmp) == 0) {
                    data.frame(
                        text = character(0), col = character(0),
                        fontsize = numeric(0),
                        stringsAsFactors = FALSE
                    )
                } else {
                    data.frame(
                        text = tmp$term, col = tmp$col,
                        fontsize = tmp$fontsize,
                        stringsAsFactors = FALSE
                    )
                }
            })
            names(term_list) <- mod_levels
            cat_terms[[cate]] <- term_list
        }

        # Hub features as a column in the multi-column textbox
        if ("Hub features" %in% enrich_category &&
            !is.null(mark_features) && length(mark_features) > 0) {
            hres <- .build_marked_features_anno(
                mark_features, module, mark_state,
                fontsize, mod_levels
            )
            if (!is.null(hres)) {
                cat_terms[["Hub features"]] <- hres
            }
        }

        # Drop modules with no terms in any category
        has_terms <- vapply(mod_levels, function(mod) {
            any(vapply(cat_terms, function(ct) {
                mod %in% names(ct) && nrow(ct[[mod]]) > 0
            }, logical(1)))
        }, logical(1))
        if (any(has_terms)) {
            cat_terms <- lapply(cat_terms, function(ct)
                ct[mod_levels[has_terms]])
        }

        # Build a single multi-column textbox from all categories
        if (length(cat_terms) > 0) {
            align_to <- factor(
                rep(mod_levels,
                    each = round(mat_nrow / length(mod_levels))
                ),
                levels = mod_levels
            )
            text <- .reformat_cat_terms(cat_terms)

            combined_anno <- .anno_multicol_textbox(
                align_to = align_to, text = text,
                background_gp = grid::gpar(
                    fill = "grey95", col = "grey50"
                ),
                by = "anno_link", side = "right",
                show_headers = TRUE,
                header_gp = grid::gpar(
                    fontsize = fontsize, fontface = "bold"
                )
            )

            right_ann <- ComplexHeatmap::rowAnnotation(
                profile = anno_ggplot2_profile,
                Enrichment = combined_anno
            )
        }
    }

    # ---- anno_mark for key features ----
    # Skip anno_mark if Hub features is already shown as textbox
    has_hub_text <- "Hub features" %in% enrich_category
    left_ann <- NULL
    if (length(mark_state$feat_col) > 0 && !has_hub_text) {
        feats_present <- intersect(
            names(mark_state$feat_col), rownames(mat)
        )
        if (length(feats_present) > 0) {
            at_idx <- match(feats_present, rownames(mat))
            cols <- mark_state$feat_col[feats_present]
            left_ann <- ComplexHeatmap::rowAnnotation(
                mark = ComplexHeatmap::anno_mark(
                    at = at_idx,
                    labels = feats_present,
                    side = "left",
                    labels_gp = grid::gpar(
                        col = cols,
                        fontsize = fontsize - 1
                    ),
                    link_gp = grid::gpar(col = "grey40", lwd = 0.5),
                    padding = unit(2, "mm")
                )
            )
        }
    }

    ht_args <- list(
        mat,
        name = "heatmap", col = col_fun,
        cluster_rows = FALSE, cluster_columns = FALSE,
        show_row_names = FALSE, show_column_names = FALSE,
        row_split = data_module_wider$Module,
        row_gap = unit(0, "mm"), row_title_rot = 0,
        row_title_gp = grid::gpar(fontsize = fontsize),
        layer_fun = function(j, i, x, y, w, h, fill, slice_r, slice_c) {
            grid::grid.rect(gp = grid::gpar(
                fill = "transparent", col = "grey50"
            ))
        }, # draw frames around the modules
        column_split = sp_info$Group,
        column_title_gp = grid::gpar(fontsize = fontsize),
        right_annotation = right_ann,
        use_raster = FALSE,
        heatmap_legend_param = list(
            title = ylabel,
            title_gp = grid::gpar(fontsize = fontsize, fontface = "bold"),
            labels_gp = grid::gpar(fontsize = fontsize)
        ),
        width = grid::unit(heatmap_width, "cm"),
        height = grid::unit(heatmap_height, "cm")
    )
    if (!is.null(left_ann)) ht_args$left_annotation <- left_ann
    heatmap <- do.call(ComplexHeatmap::Heatmap, ht_args)

    anno_ggplot2_profile_lgd <- ComplexHeatmap::Legend(
        title = "Group", labels = groups,
        legend_gp = grid::gpar(
            fill = get_custom_palette(groups)
        ),
        labels_gp = grid::gpar(fontsize = fontsize),
        title_gp = grid::gpar(fontsize = fontsize)
    )

    draw_ht <- function() {
        lgd_list <- list(anno_ggplot2_profile_lgd)
        n_cats <- length(enrich_category)
        if (is.null(enrich_p_threshold)) {
            p_thresh <- rep(NA_real_, n_cats)
        } else {
            p_thresh <- rep(enrich_p_threshold, length.out = n_cats)
        }
        for (ci in seq_len(n_cats)) {
            cate <- enrich_category[ci]
            pt <- p_thresh[ci]
            if (cate == "Hub features") {
                if (!is.null(mark_features) && length(mark_features) > 0) {
                    if (is.list(mark_features)) {
                        n_hub <- length(mark_features)
                        hub_cols <-
                            ggsci::pal_jco()(max(n_hub, 2))[seq_len(n_hub)]
                        hub_lgd <- ComplexHeatmap::Legend(
                            title = "Hub features",
                            labels = names(mark_features),
                            type = "lines",
                            legend_gp = grid::gpar(col = unname(hub_cols),
                                lwd = 2),
                            labels_gp = grid::gpar(fontsize = fontsize),
                            title_gp = grid::gpar(fontsize = fontsize)
                        )
                        lgd_list <- c(lgd_list, list(hub_lgd))
                    } else {
                        # Character vector: black
                        hub_lgd <- ComplexHeatmap::Legend(
                            title = "Hub features",
                            labels = "Hub features",
                            type = "lines",
                            legend_gp = grid::gpar(col = "black", lwd = 2),
                            labels_gp = grid::gpar(fontsize = fontsize),
                            title_gp = grid::gpar(fontsize = fontsize)
                        )
                        lgd_list <- c(lgd_list, list(hub_lgd))
                    }
                }
            } else if (!is.na(pt) && !is.null(pt)) {
                go_lgd <- ComplexHeatmap::Legend(
                    title = cate,
                    labels = c(
                        sprintf("p.adj < %.2f", pt),
                        sprintf("p.adj >= %.2f", pt)
                    ),
                    type = "lines",
                    legend_gp = grid::gpar(col = c("black", "grey70"), lwd = 2),
                    labels_gp = grid::gpar(fontsize = fontsize),
                    title_gp = grid::gpar(fontsize = fontsize)
                )
                lgd_list <- c(lgd_list, list(go_lgd))
            }
        }
        has_hub_text <- "Hub features" %in% enrich_category
        if (!is.null(mark_state$lgd) && !has_hub_text) {
            lgd_list <- c(lgd_list, list(mark_state$lgd))
        }

        pad_right <- 2
        ComplexHeatmap::draw(heatmap,
            heatmap_legend_list = lgd_list,
            heatmap_legend_side = "right", merge_legends = TRUE,
            padding = unit(c(2, 2, pad_right, 2), "mm")
        )
    }

    if (!is.null(save)) {
        if (suffix != "") {
            suffix <- paste0(suffix, "_")
        }

        for (ext in device) {
            if (ext == "pdf") {
                grDevices::pdf(
                    file = file.path(save, paste0("WGCNA_h_", suffix,
                        Sys.Date(), ".", ext)),
                    width = width / 2.54,
                    height = height / 2.54
                )
            } else if (ext == "png") {
                grDevices::png(file.path(save, paste0("WGCNA_h_",
                    suffix, Sys.Date(), ".", ext)),
                    width = width, height = height,
                    units = "cm", res = res)
            } else if (ext == "tiff") {
                grDevices::tiff(file.path(save, paste0("WGCNA_h_",
                    suffix, Sys.Date(), ".", ext)),
                    width = width, height = height,
                    units = "cm", res = res)
            } else if (ext == "jpeg") {
                grDevices::jpeg(file.path(save, paste0("WGCNA_h_",
                    suffix, Sys.Date(), ".", ext)),
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
}
