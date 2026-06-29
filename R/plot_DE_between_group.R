#' DE number between groups
#' @description Plot the number of differentially expressed features between
#' groups over time, using the output of `DE_between_group()`.
#'
#' DE features were pre-filtered with `DE_between_group()`, but by specifying
#' `adjP_thres` and `logFC_thres`, users can re-filter the DE features for
#' plotting.
#'
#' @param DE_between_group_out Output of `DE_between_group()`, a list with
#' `all_list`, `de_list`, `ref_groups`, `all_groups` elements.
#' @param group A character vector specifying which groups to use
#' as the reference (Cond1) for plotting. If NULL, all groups in the output
#' will be plotted. (default is NULL)
#' @param fontsize Font size for the plot (default is 8)
#' @param adjP_thres (Optional) Threshold for adjusted p-value to consider a
#' feature as differentially expressed, for re-filtering
#' the DE features. (default is NULL, no re-filtering)
#' @param logFC_thres (Optional) Threshold for log2 fold change to consider a
#' feature as differentially expressed, for re-filtering
#' the DE features. (default is NULL, no re-filtering)
#'
#' @import ggplot2
#'
#' @returns A series of line plots showing the number of DE features over time
#' for each group comparison.
#' @export
#' @examples
#' data(example_obj)
#' example_obj <- normalise_to_start(example_obj)
#'
#' DE_between_group_out <- DE_between_group(example_obj, assay = 2)
#' plot_DE_between_group(DE_between_group_out, fontsize = 8)
plot_DE_between_group <- function(
    DE_between_group_out,
    group = NULL,
    fontsize = 8,
    adjP_thres = NULL, logFC_thres = NULL
) {
    .check_list(DE_between_group_out, "DE_between_group_out")
    .check_positive(fontsize, "fontsize")
    if (!is.null(group)) .check_character(group, "group")
    if (!is.null(adjP_thres)) .check_pval(adjP_thres, "adjP_thres")
    if (!is.null(logFC_thres)) .check_nonneg(logFC_thres, "logFC_thres")
    if (!is.list(DE_between_group_out) ||
        !all(c("all_list", "de_list") %in% names(DE_between_group_out))) {
        stop("'DE_between_group_out' must be the output of DE_between_group().")
    }

    all_list <- DE_between_group_out$all_list

    if (!is.null(adjP_thres) && !is.null(logFC_thres)) {
        message("Re-filtering DE features with adjP_thres = ", adjP_thres,
            " and logFC_thres = ", logFC_thres, ".")

        # Extract DE features with new thresholds
        de_list <- list()

        for (label_comparison in names(all_list)) {
            de_list[[label_comparison]] <- list()

            for (tb_name in names(all_list[[label_comparison]])) {
                tb <- all_list[[label_comparison]][[tb_name]]
                tb_de <- tb |>
                    dplyr::select(Feature, Comparison, Time, Cond1, Cond2,
                        logFC, adj.P.Val) |>
                    dplyr::filter(adj.P.Val < adjP_thres) |>
                    dplyr::filter(logFC > logFC_thres | logFC < -logFC_thres)
                de_list[[label_comparison]][[tb_name]] <- tb_de
            }

            de_list[[label_comparison]] <-
                de_list[[label_comparison]] |>
                dplyr::bind_rows()
        }
    } else {
        de_list <- DE_between_group_out$de_list
    }

    all_ref_groups <- DE_between_group_out$ref_groups
    all_groups <- DE_between_group_out$all_groups

    if (is.null(group)) {
        group <- all_ref_groups
    } else if (!all(group %in% all_ref_groups)) {
        stop("Specified group(s) not found in DE_between_group output: ",
            paste(setdiff(group, all_ref_groups), collapse = ", "))
    }

    de_num_limma <- de_list |>
        lapply(function(x) {
            x |>
                dplyr::group_by(Time) |>
                dplyr::summarise(N = dplyr::n(), .groups = "drop") |>
                dplyr::mutate(
                    Cond1 = unique(x$Cond1),
                    Cond2 = unique(x$Cond2),
                    Time = as.numeric(as.character(Time))
                )
        }) |>
        dplyr::bind_rows()

    # When certain comparisons are valid but have no DE features,
    # the time point is not included in the plot, so add it with N = 0.
    # Derive expected time points from all_list (all comparisons,
    # not just those with DE features).
    all_comparisons_t <- all_list |>
        lapply(function(x) {
            lapply(x, function(tb) {
                data.frame(
                    Cond1 = unique(tb$Cond1),
                    Cond2 = unique(tb$Cond2),
                    Time = as.numeric(as.character(unique(tb$Time)))
                )
            }) |>
                dplyr::bind_rows()
        }) |>
        dplyr::bind_rows() |>
        dplyr::distinct()

    if (nrow(all_comparisons_t) > nrow(de_num_limma)) {
        missing_comparisons <- dplyr::anti_join(all_comparisons_t,
            de_num_limma,
            by = c("Cond1", "Cond2", "Time"))
        missing_comparisons$N <- 0
        de_num_limma <- dplyr::bind_rows(de_num_limma, missing_comparisons)
    }

    p_list <- list()
    for (i in group) {
        p_list[[i]] <- de_num_limma |>
            dplyr::filter(Cond1 == i) |>
            dplyr::mutate(Group = Cond2) |>
            ggplot(aes(x = Time, y = N, group = Group, color = Group)) +
            geom_point() +
            geom_line() +
            labs(
                title = paste0("DE feature number (vs. ", i, ")"),
                y = "DE feature number"
            ) +
            scale_color_manual(values = get_custom_palette(all_groups)) +
            theme_custom(base_size = fontsize)
    }

    return(p_list)
}
