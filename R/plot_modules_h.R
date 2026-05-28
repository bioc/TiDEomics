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
#' @param ylabel Y axis label prefix (default is "Abundance")
#' @param suffix Suffix for the saved image file name (default is an empty
#' string)
#' @param device Image file format(s) for saving. Can be a character
#'   vector with one or more of `"png"`, `"pdf"`, `"tiff"`, `"jpeg"`
#'   (default: `"png"`)
#' @param save Directory to save the plot, no saving if is NULL
#' (default is NULL)
#' @param profile_width Width of the mean expression profile plot
#' (default is 3 (cm))
#' @param profile_link_width Width of the link in cm between heatmap
#' and mean expression profile plot (default is 1)
#' @param enrich_list A named list of enrichment results for WGCNA modules,
#'   as produced by `enrichGO_list()$all`. Each named element should be a
#'   data.frame with columns `Cluster`, `Description`, and the rank column.
#'   If NULL (default), no enrichment terms will be displayed.
#' @param enrich_category Name of the enrichment category to display
#'   (default is `"BP"`).
#' @param enrich_rank_by Column name to rank enrichment terms by
#'   (default is `"p.adjust"`).
#' @param enrich_top_n Number of top enrichment terms to display per module
#'   (default is 3)
#' @param fontsize Font size (default is 8)
#' @param heatmap_width Width of the heatmap body in cm (default is 8)
#' @param heatmap_height Height of the heatmap body in cm (default is 8)
#' @param width Width of the saved image (default is 16 (cm))
#' @param height Height of the saved image (default is 12 (cm))
#' @param res Resolution of the saved image (except pdf format)
#' (default is 300 (ppi))
#'
#' @import magrittr
#' @import ggplot2
#' @import SummarizedExperiment
#'
#' @returns A combined plot of heatmap, mean expression profile, and enrichment
#' terms for each WGCNA module
#' @export
#' @examples
#' library(magrittr)
#' library(org.Mm.eg.db)
#'
#' data(example)
#' example_obj <- normalise_to_start(example_obj)
#' example_obj_list <- split_groups(example_obj)
#' example_obj_merged_list <- merge_replicates(example_obj_list)
#' example_obj_merged <- merge_groups(example_obj_merged_list)
#'
#' data(example_net)
#' # select two modules for demonstration
#' example_module <- WGCNA_module(example_net) %>%
#'     dplyr::filter(Module %in% c("1", "2"))
#' # set cutoff to 1 to show all results for demonstration
#' example_go_list = enrichGO_list(example_module, OrgDb = org.Mm.eg.db,
#'     universe = example_module$Feature,
#'     pvalueCutoff = 1, qvalueCutoff = 1,
#'     category = "BP", simplify = FALSE)
#' # plot_GO(example_go_list$all, plot_dotplot = TRUE,
#' #     plot_emapplot = FALSE, plot_cnetplot = FALSE)
#'
#' plot_modules_h(example_module %>% dplyr::filter(Module != '0'),
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
    suffix = "",
    device = "png",
    save = NULL,
    profile_width = 3,
    profile_link_width = 1,
    enrich_list = NULL,
    enrich_category = "BP",
    enrich_rank_by = "p.adjust",
    enrich_top_n = 3,
    fontsize = 8,
    heatmap_width = 8,
    heatmap_height = 8,
    width = 16,
    height = 12,
    res = 300
) {
    data_module_long <- .plot_modules_input(module = module,
        se_obj_merged = se_obj_merged, scale = scale, assay = assay)

    n_modules <- length(unique(data_module_long$Module))

    stopifnot(n_modules == length(levels(data_module_long$Module)))

    groups <- data_module_long$Group %>% levels()

    sp_info <- data_module_long %>%
        dplyr::select(Group, Time) %>%
        dplyr::distinct()

    # plot modules with complexheatmap
    data_module_wider <- data_module_long %>%
        tidyr::pivot_wider(id_cols = c(Feature, Module),
            names_from = c(Group, Time), values_from = Abundance) %>%
        tibble::column_to_rownames("Feature") %>%
        as.data.frame() %>%
        dplyr::arrange(Module)
    mat <- data_module_wider %>%
        dplyr::select(-Module) %>%
        as.matrix()

    ggplot2_profile_list <- list()
    for (i in levels(data_module_long$Module)) {
        ggplot2_profile_list[[i]] <- ggplot(
            data_module_long %>% dplyr::filter(Module == i),
            aes(x = Time, y = Abundance, group = Feature)
        ) +
            stat_summary(aes(group = Group, color = Group),
                fun = mean, geom = "line",
                linewidth = 1
            ) +
            scale_color_manual(values = get_custom_palette(groups)) +
            labs(
                x = NULL, y = NULL,
                title = NULL
            ) +
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

    ggplot.panel.arg <- c(1, 0, profile_width,
        "grey95", "grey50", profile_link_width)
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
        gap = unit(as.numeric(ggplot.panel.arg[2]), "cm"),
        width = unit(as.numeric(ggplot.panel.arg[3]), "cm"),
        side = "right",
        link_gp = grid::gpar(
            fill = ggplot.panel.arg[4],
            col = ggplot.panel.arg[5]
        ),
        link_width = unit(ggplot.panel.arg[6], "cm")
    )

    # legend for mean expression profile
    anno_ggplot2_profile_lgd <- ComplexHeatmap::Legend(
        title = "Group",
        labels = groups,
        legend_gp = grid::gpar(fill = get_custom_palette(groups)),
        labels_gp = grid::gpar(fontsize = fontsize),
        title_gp = grid::gpar(fontsize = fontsize)
    )

    min_val <- min(mat, na.rm = TRUE)
    max_val <- max(mat, na.rm = TRUE)

    # GO terms of the specified category
    if (!is.null(enrich_list[[enrich_category]])) {
        termanno <- enrich_list[[enrich_category]] %>%
            as.data.frame() %>%
            dplyr::mutate(id = as.character(Cluster)) %>%
            dplyr::filter(id %in% levels(data_module_long$Module)) %>%
            dplyr::group_by(id) %>%
            dplyr::slice_min(order_by = .data[[enrich_rank_by]],
                n = enrich_top_n, with_ties = FALSE) %>%
            dplyr::mutate(term = Description) %>%

            dplyr::select(id, term)

        termanno <- termanno |>
            dplyr::ungroup() |>
            dplyr::mutate(col = "black", fontsize = fontsize)

        # to list
        term.list <- lapply(seq_len(length(unique(termanno$id))), function(x) {
            tmp <- termanno[which(termanno$id == unique(termanno$id)[x]), ]
            df <- data.frame(
                text = tmp$term,
                col = tmp$col,
                fontsize = tmp$fontsize
            )
            return(df)
        })

        # add names
        termAnno.arg <- c("grey95", "grey50")
        names(term.list) <- unique(termanno$id)
        textbox <- ComplexHeatmap::anno_textbox(
            align_to = factor(
                rep(levels(data_module_wider$Module),
                    each = round(nrow(mat) / n_modules, digits = 0)
                ),
                levels = levels(data_module_wider$Module)
            ),
            text = term.list,
            word_wrap = FALSE,
            add_new_line = TRUE,
            side = "right",
            background_gp = grid::gpar(
                fill = termAnno.arg[1],
                col = termAnno.arg[2]
            ),
            by = "anno_link",
        )

        p <- ComplexHeatmap::Heatmap(mat,
            name = "heatmap",
            col = circlize::colorRamp2(c(min_val, 0, max_val),
                c("#3C5488FF", "white", "#E64B35FF")),
            cluster_rows = FALSE, cluster_columns = FALSE,
            show_row_names = FALSE, show_column_names = FALSE,
            row_split = data_module_wider$Module,
            row_gap = unit(0, "mm"),
            row_title_rot = 0,
            row_title_gp = grid::gpar(fontsize = fontsize),
            layer_fun = function(j, i, x, y, w, h, fill, slice_r, slice_c) {
                grid::grid.rect(gp = grid::gpar(
                    fill = "transparent",
                    col = "grey50"
                ))
            }, # draw frames around the modules
            column_split = sp_info$Group,
            column_title_gp = grid::gpar(fontsize = fontsize),
            right_annotation = ComplexHeatmap::rowAnnotation(
                anno_ggplot2 = anno_ggplot2_profile,
                textbox = textbox
            ),
            width = unit(heatmap_width, "cm"),
            height = unit(heatmap_height, "cm"),
            use_raster = FALSE,
            heatmap_legend_param = list(
                title = ylabel,
                title_gp = grid::gpar(fontsize = fontsize),
                labels_gp = grid::gpar(fontsize = fontsize)
            )
        )
    } else {
        if (!is.null(enrich_list)) { # & is.null(enrich_list[[category]])
            message("GO category '", enrich_category, 
                "' not found in enrich_list, ",
                "skipping enrichment term annotation")
        }

        p <- ComplexHeatmap::Heatmap(mat,
            name = "heatmap",
            col = circlize::colorRamp2(c(min_val, 0, max_val),
                c("#3C5488FF", "white", "#E64B35FF")),
            cluster_rows = FALSE, cluster_columns = FALSE,
            show_row_names = FALSE, show_column_names = FALSE,
            row_split = data_module_wider$Module,
            row_gap = unit(0, "mm"),
            row_title_rot = 0,
            row_title_gp = grid::gpar(fontsize = fontsize),
            layer_fun = function(j, i, x, y, w, h, fill, slice_r, slice_c) {
                grid::grid.rect(gp = grid::gpar(
                    fill = "transparent",
                    col = "grey50"
                ))
            }, # draw frames around the modules
            column_split = sp_info$Group,
            column_title_gp = grid::gpar(fontsize = fontsize),
            right_annotation = ComplexHeatmap::rowAnnotation(anno_ggplot2 =
                anno_ggplot2_profile),
            width = unit(heatmap_width, "cm"),
            height = unit(heatmap_height, "cm"),
            use_raster = FALSE,
            heatmap_legend_param = list(
                title = ylabel,
                title_gp = grid::gpar(fontsize = fontsize),
                labels_gp = grid::gpar(fontsize = fontsize)
            )
        )
    }

    # add legend for mean expression profile
    p <- grid::grid.grabExpr(ComplexHeatmap::draw(p,
        heatmap_legend_list = list(anno_ggplot2_profile_lgd),
        heatmap_legend_side = "right",
        merge_legends = TRUE
    ))
    # cannot use show_annotation_legend here because will produce
    # duplicated legends for each module

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
            on.exit(try(grDevices::dev.off(), silent = TRUE), add = TRUE)
            grid::grid.draw(p)
            grDevices::dev.off()
        }
    }

    grid::grid.newpage()
    grid::grid.draw(p)
}
