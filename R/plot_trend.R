#' Plot feature abundance over time
#'
#' @description Plot trend of feature abundances / expression over time by
#' mean and standard deviation (SD) for each group. Accepts either a
#' SummarizedExperiment object (computing mean/SD internally via
#' `calc_mean_sd()`) or a pre-computed table from
#' `calc_mean_sd()`.
#'
#' @param se_obj A SummarizedExperiment object, or a data.frame from
#'   `calc_mean_sd()` with columns: Feature, Time, Group, Mean, SD.
#' @param assay The assay to use in the SummarizedExperiment object: a numeric
#'   index or character name, e.g. 1 or "orig" for original data, 2 or "norm"
#'   for time 0 normalised data. (Default: 1.)
#' @param groups Groups to be plotted, if NULL, all groups will be used
#'   (default is NULL)
#' @param features Features to be plotted, if NULL, an error will be raised
#' @param title Title of the plot. "Feature" if NULL (default).
#' @param ylabel Y axis label of the plot (default is "Log2 abundance")
#' @param errorbar Whether to plot error bars (default is TRUE)
#' @param fontsize Font size for the plot (default is 8)
#'
#' @import ggplot2
#' @returns Plot of feature abundances over time by mean and SD
#' @export
#' @examples
#' data(example_obj)
#' plot_trend(example_obj,
#'     features = sample(rownames(example_obj), 4))
plot_trend <- function(se_obj, assay = 1, groups = NULL, features,
    title = NULL, ylabel = "Log2 abundance", errorbar = TRUE, fontsize = 8) {

    if (is.null(title)) title <- "Feature"
    if (!inherits(se_obj, "SummarizedExperiment") && !is.data.frame(se_obj)) {
        stop("'se_obj' must be a SummarizedExperiment or a data.frame.")
    }
    .check_character(features, "features")
    .check_character(title, "title")
    .check_character(ylabel, "ylabel")
    .check_logical(errorbar, "errorbar")
    .check_positive(fontsize, "fontsize")
    if (!is.null(groups)) .check_character(groups, "groups")

    if (length(features) == 0) {
        stop("'features' must not be empty.")
    }

    # Accept either an SE object or a pre-computed table_mean_sd
    if (inherits(se_obj, "SummarizedExperiment")) {
        assay <- .match_assay(assay, se_obj)
        # Subset to features to be plotted only
        keep <- intersect(rownames(se_obj), features)
        se_obj <- se_obj[keep, , drop = FALSE]
        tbl_list <- calc_mean_sd(se_obj)
        table_mean_sd <- tbl_list[[assay]]
    } else if (is.data.frame(se_obj)) {
        table_mean_sd <- se_obj
    } else {
        stop("se_obj must be a SummarizedExperiment object or a data.frame ",
            "produced by calc_mean_sd().")
    }

    if (is.null(groups)) {
        groups <- unique(table_mean_sd$Group)
        message(sprintf("Group not specified. Plotting all groups: %s",
            paste(groups, collapse = ", ")))
    } else if (!all(groups %in% unique(table_mean_sd$Group))) {
        stop("At least one of the specified groups is not found ",
            "in the input.")
    }

    p <- table_mean_sd |>
        dplyr::filter(Group %in% groups & Feature %in% features) |>
        dplyr::mutate(Feature = factor(Feature, levels = features)) |>
        ggplot(aes(x = Time, y = Mean, group = Group, color = Group)) +
        geom_line(linewidth = 0.7) +
        geom_point(size = 0.3) +
        scale_color_manual(values = get_custom_palette(groups)) +
        theme_custom(
            panel_border = TRUE, legend_position = "bottom",
            base_size = fontsize
        ) +
        ggtitle(paste0(title, " in ", paste(groups, collapse = ", "))) +
        xlab("Time") +
        ylab(ylabel) +
        facet_wrap(~Feature, scales = "free_y")

    if (errorbar) {
        p <- p + geom_errorbar(aes(ymin = Mean - SD, ymax = Mean + SD),
            width = 0.5,
            linewidth = 0.7,
            position = position_dodge(0.05)
        )
    }

    return(p)
}
