#' Plot modules (vertical layout)
#' @description Plot WGCNA modules' mean expression profiles and heatmaps,
#' align vertically.
#'
#' Different from running WGCNA, the input data should have the
#' replicates merged, instead of having multiple samples per group, time and
#' feature (gene).
#'
#' If certain time points are missing in some groups, NA values are added.
#'
#' @param module A data frame with columns "Feature" and "Module"
#' @param se_obj_merged A SummarizedExperiment object, with one value for each
#' feature at each time point in each group (replicates merged), or a named
#' list of such objects (one per dataset). The colData of
#' the object should contain columns "Sample", "Group", and "Time". The object
#' can be produced by `split_groups()`, `merge_replicates()` and
#' `merge_groups()`.
#' @param assay The assay to use in the SummarizedExperiment object: a numeric
#'   index or character name, e.g. 1 or "orig" for original data, 2 or "norm"
#'   for time 0 normalised data. (Required, no default.)
#' @param scale Whether to scale the data (z-score) across samples for each
#' feature (default is TRUE)
#' @param ylabel Y axis label. A single character string shared by all
#'   datasets, or a character vector with one entry per dataset
#'   (multi-input only). Default: NULL, i.e. "Z-score of log2 abundance"
#'   when scale = TRUE, "Log2 abundance" otherwise.
#' @param suffix Suffix for the saved image file name (default is an empty
#' string)
#' @param device Image file format(s) for saving. Can be a character
#'   vector, e.g. `c("png", "pdf")`, to save in multiple formats
#'   (default: `"png"`).
#' @param save Directory to save the plot, no saving if is NULL
#' (default is NULL)
#' @param fontsize Font size for the plot (default is 8)
#' @param width Width of the saved image (default is 12 (cm))
#' @param height Height of the saved image (default is 8 (cm))
#' @param height_ratio Ratio of height of heatmap to the line plot
#' (default is 2)
#' @param res Resolution of the saved image (except pdf format)
#' (default is 300 (ppi))
#' @param names Optional character vector of dataset names for the plot
#'   subtitles when multi-dataset. Default: NULL, i.e.
#'   the names of `se_obj_merged`.
#'
#' @import ggplot2
#' @import patchwork
#'
#' @returns For a single SummarizedExperiment, two plots aligned vertically by
#' groups: the top one is a line plot of module feature mean expression
#' profiles, the bottom one is a heatmap of feature expression across time
#' points. For a named list of datasets, the per-dataset plot pairs are
#' stacked vertically, with each dataset's line plot annotated with its name
#' as a subtitle and an overall title "Module profiles: ...".
#' @export
#' @examples
#' data(example_obj)
#' example_obj <- normalise_to_start(example_obj)
#' example_obj_list <- split_groups(example_obj)
#' example_obj_merged_list <- merge_replicates(example_obj_list)
#' example_obj_merged <- merge_groups(example_obj_merged_list)
#'
#' data(example_net)
#' example_module <- WGCNA_module(example_net)
#'
#' plot_modules_v(example_module |> dplyr::filter(Module != '0'),
#'     example_obj_merged, assay = 2, scale = TRUE,
#'     ylabel = "Z-score of log2 (expression)",
#'     height_ratio = 2,
#'     fontsize = 6)
#'
#' plot_modules_v(example_module |> dplyr::filter(Module != '0'),
#'     list(Exp1 = example_obj_merged, Exp2 = example_obj_merged),
#'     assay = 2, scale = TRUE, height_ratio = 2, fontsize = 6)
plot_modules_v <- function(
    module,
    se_obj_merged,
    scale = TRUE,
    assay,
    ylabel = NULL,
    suffix = "",
    device = "png",
    save = NULL,
    width = 12,
    height = 8,
    height_ratio = 2,
    fontsize = 8,
    res = 300,
    names = NULL
) {
    .check_df(module, "module")
    .check_character(device, "device")
    .check_character(suffix, "suffix")
    .check_logical(scale, "scale")
    if (is.null(ylabel)) {
        ylabel <- if (scale) "Z-score of log2 abundance" else "Log2 abundance"
    }
    .check_character(ylabel, "ylabel")
    .check_positive(width, "width")
    .check_positive(height, "height")
    .check_positive(height_ratio, "height_ratio")
    .check_positive(fontsize, "fontsize")
    .check_positive_int(res, "res")
    if (!is.null(save)) .check_character(save, "save")
    se_obj_list <- .as_se_list(se_obj_merged, arg = "se_obj_merged")
    multi <- length(se_obj_list) > 1L
    nms <- base::names(se_obj_list)
    for (nm in nms) {
        .check_se_merged(se_obj_list[[nm]], arg = nm)
    }
    if (length(ylabel) != 1L &&
            length(ylabel) != length(se_obj_list)) {
        stop("The 'ylabel' must be length 1 or length equal to the ",
            "number of datasets (", length(se_obj_list), ").")
    }
    if (!is.null(names)) {
        .check_character(names, "names")
        if (length(names) != length(se_obj_list)) {
            stop("'names' must have one entry per dataset.")
        }
    }
    disp <- if (multi) {
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

    if (multi) {
        module <- .warn_unshared_features(module, se_obj_list, assay, disp)
    }

    panel_list <- list()
    for (i in seq_along(se_obj_list)) {
        data_module_long <- .plot_modules_input(module = module,
            se_obj_merged = se_obj_list[[i]], scale = scale, assay = assay)

        p <- data_module_long |>
            dplyr::mutate(Time = as.numeric(as.character(Time))) |>
            ggplot(aes(x = Time, y = Abundance, group = Feature)) +
            geom_line(alpha = 0.3, color = "grey") +
            stat_summary(aes(group = Group, color = Group),
                fun = mean, geom = "line",
                linewidth = 1
            ) +
            scale_color_manual(values =
                get_custom_palette(levels(data_module_long$Group))) +
            facet_grid(~Module, scales = "free") +
            ylab(if (length(ylabel) == 1L) ylabel else ylabel[i]) +
            theme_custom(panel_border = TRUE, legend_position = "right",
                base_size = fontsize) +
            theme(
                axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1),
                strip.text = element_text(size = fontsize + 2, face = "plain"),
                strip.background = element_rect(fill = NA, colour = "black",
                    linewidth = 0.5)
            )
        if (multi) {
            p <- p + labs(subtitle = disp[i])
        }

        q <- data_module_long |>
            ggplot(aes(x = Time, y = Feature, fill = Abundance)) +
            geom_tile() +
            scale_y_discrete(
                expand = expansion(add = c(30, 0))
            ) +
            scale_fill_gradientn(
                colours = grDevices::colorRampPalette(c("#3C5488FF",
                    "white", "#E64B35FF"))(100),
                na.value = "grey"
            ) +
            ggh4x::facet_grid2(vars(Group), vars(Module),
                scales = "free", axes = "all",
                remove_labels = "x", independent = "y",
                strip =
                    ggh4x::strip_themed(background_y = ggh4x::elem_list_rect(
                    fill = get_custom_palette(levels(data_module_long$Group))
                ))
            ) +
            theme_custom(legend_position = "right", base_size = fontsize) +
            theme(
                axis.text = element_blank(),
                axis.line = element_blank(),
                axis.ticks.length = unit(0, "cm"),
                strip.text = element_text(size = fontsize + 2, face = "plain"),
                strip.background = element_rect(fill = NA, colour = "black",
                    linewidth = 0.5)
            )

        panel_list[[2 * i - 1]] <- p
        panel_list[[2 * i]] <- q
    }

    heights <- rep(c(1, height_ratio), length(se_obj_list))
    pq <- patchwork::wrap_plots(panel_list, ncol = 1, heights = heights)
    if (multi) {
        pq <- pq + patchwork::plot_annotation(
            title = paste0("Module profiles: ", paste(disp, collapse = ", ")))
    }

    if (!is.null(save)) {
        if (suffix != "") {
            suffix <- paste0(suffix, "_")
        }
        prefix <- if (multi) "WGCNA_v_multi_" else "WGCNA_v_"
        for (ext in device) {
            ggsave(paste0(prefix, suffix, Sys.Date(), ".", ext),
                plot = pq,
                device = ext,
                path = save,
                dpi = res,
                width = width,
                height = height,
                units = "cm",
                limitsize = FALSE
            )
        }
    }

    return(pq)
}
