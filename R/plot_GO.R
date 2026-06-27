#' Plot GO enrichment
#' @description Visualization of GO enrichment analysis with dotplot,
#' cnetplot, or emapplot in clusterProfiler
#'
#' @param go_list A list of enriched GO results, output from `enrichGO_list()`
#' @param plot_dotplot Whether to plot dotplot (default is TRUE)
#' @param plot_cnetplot Whether to plot cnetplot (default is FALSE)
#' @param plot_emapplot Whether to plot emapplot (default is FALSE)
#' @param showCategory_dotplot showCategory parameter of
#' `clusterProfiler::dotplot()` (default is 5)
#' @param showCategory_cnetplot showCategory parameter of
#' `clusterProfiler::cnetplot()` (default is 5)
#' @param showCategory_emapplot showCategory parameter of
#' `clusterProfiler::emapplot()` (default is 5)
#' @param label Title for the plots (default is "genes")
#' @param fontsize Font size for the plots (default is 8)
#' @param ... Additional parameters to pass to the plotting functions
#' `clusterProfiler::dotplot()`, `clusterProfiler::cnetplot()`, or
#' `clusterProfiler::emapplot()`
#'
#' @import ggplot2
#'
#' @returns A list of plots visualizing the GO enrichment results.
#' @export
#' @examples
#' library(org.Mm.eg.db)
#' data(example_net)
#' # select two modules for demonstration
#' example_module <- WGCNA_module(example_net) |>
#'     dplyr::filter(Module %in% c("1", "2"))
#' # set cutoff to 1 to show all results for demonstration
#' example_go_list = enrichGO_list(example_module, OrgDb = org.Mm.eg.db,
#'     universe = example_module$Feature,
#'     pvalueCutoff = 1, qvalueCutoff = 1,
#'     category = "BP", simplify = FALSE)
#' plot_GO(example_go_list$all, plot_dotplot = TRUE,
#'     plot_emapplot = FALSE, plot_cnetplot = FALSE)
plot_GO <- function(
    go_list,
    plot_dotplot = FALSE,
    plot_cnetplot = FALSE,
    plot_emapplot = FALSE,
    showCategory_dotplot = 5,
    showCategory_cnetplot = 5,
    showCategory_emapplot = 5,
    fontsize = 8,
    label = "features",
    ...
) {
    .check_character(label, "label")
    .check_list(go_list, "go_list", "enrichGO_list")
    .check_logical(plot_dotplot, "plot_dotplot")
    .check_logical(plot_cnetplot, "plot_cnetplot")
    .check_logical(plot_emapplot, "plot_emapplot")
    .check_positive_int(showCategory_dotplot, "showCategory_dotplot")
    .check_positive_int(showCategory_cnetplot, "showCategory_cnetplot")
    .check_positive_int(showCategory_emapplot, "showCategory_emapplot")
    .check_positive(fontsize, "fontsize")
    if (length(intersect(names(go_list), c("BP", "MF", "CC"))) == 0) {
        if ("all" %in% names(go_list) || "simplified" %in% names(go_list)) {
            stop("The input list contains 'all' or 'simplified' sublists. ",
            "Please specify one of them, e.g. go_list$all, ",
            "to visualize the GO enrichment results.")
        } else {
            stop("No valid GO category ('BP', 'MF', 'CC') found in ",
            "the input list. Please check to ensure it contains the expected ",
            "GO enrichment results.")
        }
    }

    if (!plot_dotplot && !plot_cnetplot && !plot_emapplot) {
        message("Please specify at least one plot type to visualize ",
        "the GO enrichment results.")
    }

    p_list <- list()

    # dotplot
    if (plot_dotplot) {
        for (cate in intersect(names(go_list), c("BP", "MF", "CC"))) {
            if (dim(go_list[[cate]])[1] > 0) {
                p_list[[paste0("dotplot_", cate)]] <-
                    clusterProfiler::dotplot(go_list[[cate]],
                        showCategory = showCategory_dotplot,
                        title = paste0("GO ", cate, " in ", label),
                        ...
                    ) +
                    theme_custom(base_size = fontsize, panel_border = TRUE) +
                    theme(
                        axis.text.x = element_text(angle = 45, hjust = 1),
                        panel.grid = element_blank()
                    )
            } else {
                message(sprintf("No significant GO %s terms found.", cate))
            }
        }
    }

    # gene-GO network plot
    if (plot_cnetplot) {
        for (cate in intersect(names(go_list), c("BP", "MF", "CC"))) {
            if (dim(go_list[[cate]])[1] > 0) {
                clus_names <- go_list[[cate]] |>
                    as.data.frame() |>
                    dplyr::pull(Cluster) |>
                    unique()

                p_list[[paste0("cnetplot_", cate)]] <-
                    clusterProfiler::cnetplot(go_list[[cate]],
                        showCategory = showCategory_cnetplot,
                        ...
                    ) +
                    scale_fill_manual(values = get_custom_palette(clus_names)) +
                    theme_custom(base_size = fontsize) +
                    theme(
                        axis.line = element_blank(),
                        axis.title = element_blank(),
                        axis.text = element_blank(),
                        axis.ticks = element_blank()
                    ) +
                    ggtitle(paste0("GO ", cate, " in ", label))
            } else {
                message(sprintf("No significant GO %s terms found.", cate))
            }
        }
    }

    # emapplot
    if (plot_emapplot) {
        for (cate in intersect(names(go_list), c("BP", "MF", "CC"))) {
            if (dim(go_list[[cate]])[1] > 0) {
                clus_names <- go_list[[cate]] |>
                    as.data.frame() |>
                    dplyr::pull(Cluster) |>
                    unique()

                p_list[[paste0("emapplot_", cate)]] <-
                    clusterProfiler::emapplot(
                        enrichplot::pairwise_termsim(go_list[[cate]]),
                        showCategory = showCategory_emapplot,
                        node_label_size = fontsize - 5, # default: 5
                        size_category = 1.5, # default: 1
                        ...
                    ) +
                    scale_fill_manual(values = get_custom_palette(clus_names)) +
                    theme_custom(base_size = fontsize) +
                    theme(
                        axis.line = element_blank(),
                        axis.title = element_blank(),
                        axis.text = element_blank(),
                        axis.ticks = element_blank()
                    ) +
                    ggtitle(paste0("GO ", cate, " in ", label))
            } else {
                message(sprintf("No significant GO %s terms found.", cate))
            }
        }
    }

    return(p_list)
}
