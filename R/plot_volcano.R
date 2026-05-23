#' Volcano plot of DE results
#' 
#' @description This function creates a volcano plot to visualize 
#' the results of differential expression (DE) analysis performed by 
#' `DE_between_group()` or `DE_between_time()`. The plot displays 
#' log2 fold change on the x-axis and -log10 adjusted p-value on the 
#' y-axis, with DE features highlighted based on specified thresholds.
#' If no thresholds are provided, the function will use the thresholds that
#' were applied when running `DE_between_group()` or `DE_between_time()`.
#'
#' @param DE_out Output from `DE_between_group()` or `DE_between_time()`
#' @param group1 (Required for `DE_between_group()` output) 
#' Name of the first group for comparison
#' @param group2 (Required for `DE_between_group()` output) 
#' Name of the second group for comparison
#' @param time (Required for `DE_between_group()` output) 
#' Time point for comparison
#' @param group (Required for `DE_between_time()` output) 
#' Name of the group for comparison
#' @param time1 (Required for `DE_between_time()` output) 
#' First time point for comparison
#' @param time2 (Required for `DE_between_time()` output) 
#' Second time point for comparison
#' @param logFC_thres (Optional) Log2 fold change threshold for 
#' highlighting DE features, default is NULL (using threshold when 
#' `DE_between_group()` or `DE_between_time()` was run)
#' @param adjP_thres (Optional) Adjusted p-value threshold for 
#' highlighting DE features, default is NULL (using threshold when 
#' `DE_between_group()` or `DE_between_time()` was run)
#' @param label (Optional) Whether to label names of DE features on the plot 
#' (default is FALSE)
#' @param fontsize (Optional) Font size for the plot (default is 8)
#' @param ... Additional arguments for ggrepel::geom_text_repel() 
#' when label = TRUE
#'
#' @import ggplot2
#' 
#' @returns A volcano plot of DE results
#' @export
#'
#' @examples
#' data("example")
#' example_obj <- normalise_to_start(example_obj)
#'
#' DE_between_group_out <- DE_between_group(example_obj, assay = 2)
#' plot_volcano(DE_between_group_out, group1 = "untreated", 
#'     group2 = "IFNbeta", time = 24,
#'     logFC_thres = 0.5, adjP_thres = 0.05, label = TRUE)
plot_volcano <- function(DE_out, 
    group1, group2, time, 
    group, time1, time2,
    logFC_thres = NULL, adjP_thres = NULL, 
    label = FALSE, fontsize = 8, ...) {
        
    if (!("all_list" %in% names(DE_out)) | !("de_list" %in% names(DE_out))) {
        stop("DE_out must be the output of DE_between_group() or ",
            "DE_between_time()")
    }

    if (missing(group) | missing(time1) | missing(time2)) {
        if (missing(group1) | missing(group2) | missing(time)) {
            stop("Please provide either 'group', 'time1', and 'time2' for ",
            "output of DE_between_time(), or ",
            "'group1', 'group2', and 'time' for output of DE_between_group().")
        } 
        # DE_between_group
        volcano_tb <- DE_out$all_list[[paste0(group2, "-",
            group1)]][[as.character(time)]]
        if (is.null(volcano_tb)) {
            stop("No DE results found for the specified group and time. ",
                "Please check the input parameters.")
        }

        if (is.null(logFC_thres) & is.null(adjP_thres)) {
            de_tb <- DE_out$de_list[[paste0(group2, "-", group1)]] %>% 
                dplyr::filter(Time == time)
        } else {
            de_tb <- volcano_tb %>% 
                dplyr::filter(Time == time & abs(logFC) >= logFC_thres &
                adj.P.Val <= adjP_thres)
        }

        volcano_tb <- volcano_tb %>%
            dplyr::mutate(Color = ifelse(Feature %in% de_tb$Feature, 
                ifelse(logFC > 0, "Red", "Blue"), "Grey"))

        p <- volcano_tb %>%
            ggplot(aes(x = logFC, y = -log10(adj.P.Val))) +
            geom_point(aes(color = Color)) +
            theme_custom(base_size = fontsize) +
            scale_color_identity() +
            xlab(paste0("Log2 fold change (", group2, "-", group1, 
                ", time ", time, ")")) +
            ylab("-Log10 adjusted p-value") 

    } else {
        # DE_between_time
        volcano_tb <- 
            DE_out$all_list[[group]][[paste0("t", time2, "-t", time1)]]
        if (is.null(volcano_tb)) {
            stop("No DE results found for the specified group and time. ",
                "Please check the input parameters.")
        }
        
        if (is.null(logFC_thres) & is.null(adjP_thres)) {
            de_tb <- 
                DE_out$de_list[[group]][[paste0("t", time2, "-t", time1)]] 
        } else { # filter DE features based on thresholds provided
            de_tb <- volcano_tb %>% 
                dplyr::filter(
                    abs(logFC) >= logFC_thres & adj.P.Val <= adjP_thres)
        }

        volcano_tb <- volcano_tb %>%
            dplyr::mutate(Color = ifelse(Feature %in% de_tb$Feature, 
                ifelse(logFC > 0, "Red", "Blue"), "Grey"))

        p <- volcano_tb %>%
            ggplot(aes(x = logFC, y = -log10(adj.P.Val))) +
            geom_point(aes(color = Color)) +
            theme_custom(base_size = fontsize) +
            scale_color_identity() +
            xlab(paste0("Log2 fold change (", time2, "-", time1, 
                ", group ", group, ")")) +
            ylab("-Log10 adjusted p-value") 
    }
    
    if (!is.null(logFC_thres)) {
        p <- p + geom_vline(xintercept = c(-logFC_thres, logFC_thres), 
            linetype = "dashed", color = "black")
    }
    if (!is.null(adjP_thres)) {
        p <- p + geom_hline(yintercept = -log10(adjP_thres), 
            linetype = "dashed", color = "black")
    }

    if (label) {
        p <- p + ggrepel::geom_text_repel(data = de_tb, 
            aes(label = Feature), size = fontsize * 0.35, max.overlaps = Inf,
            ...)
    }

    print(p)
}
