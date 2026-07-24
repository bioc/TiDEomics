#' Plot variance decomposition
#'
#' @description Based on `PALMO::variancefeaturePlot()` function
#'
#' @param var_decomp A data frame output from `decomp_variance()`
#' @param rank Rank the plot by selected variable (default is "Group")
#' @param features Features to be plotted (default is NULL)
#' @param top_n Top n features ranked by `rank` to be plotted when `features`
#' is not specified (default is 20)
#' @param show_ylab Whether to show y axis text (default is TRUE)
#' @param fontsize Font size for the plot (default is 8)
#'
#' @import ggplot2
#'
#' @returns A stacked bar plot showing the percentage of variance explained
#' by each model component (Group, Time, Residual, plus Subject and
#' interaction terms when present). Features
#' are ranked by the specified `rank` variable (Group or Time) and only the
#' top n features are plotted if `features` is not specified.
#' @export
#' @examples
#' data(example_obj)
#' example_obj <- normalise_to_start(example_obj)
#'
#' var_decomp <- decomp_variance(example_obj, assay = 1)
#' plot_variance(var_decomp, rank = "Time", top_n = 20)
plot_variance <- function(var_decomp, rank = "Group",
    features = NULL,
    top_n = 20, show_ylab = TRUE, fontsize = 8) {
    .check_df(var_decomp, "var_decomp")
    .check_character(rank, "rank")
    .check_logical(show_ylab, "show_ylab")
    .check_positive_int(top_n, "top_n")
    .check_positive(fontsize, "fontsize")

    if (!rank %in% colnames(var_decomp)) {
        avail <- setdiff(colnames(var_decomp),
            c("Feature", "mean", "median", "sd", "max"))
        stop(sprintf("Rank column '%s' not in decomp output. Available: %s",
            rank, paste(avail, collapse = ", ")))
    }
    rank <- sym(rank)

    if (is.null(features)) {
        if (top_n > nrow(var_decomp)) {
            top_n <- nrow(var_decomp)
            message("Specified top_n is larger than the number of features. ",
                "Plotting all ", top_n, " features.")
        } else {
            message("Features not specified. Plotting top ", top_n,
                " features ", "ranked by ", as.character(rank), ".")
        }
        features <- var_decomp |>
            dplyr::arrange(dplyr::desc(!!rank), Feature) |>
            dplyr::slice(seq(1, top_n)) |>
            dplyr::pull(Feature)
    } else if (!all(features %in% var_decomp$Feature)) {
        stop("At least one of the specified features is not found in the ",
        "input.")
    }

    # Dynamic column detection: all variance components (not metadata columns)
    meta_cols <- c("Feature", "mean", "median", "sd", "max")
    comp_cols <- setdiff(colnames(var_decomp), meta_cols)
    var_decomp_longer <- var_decomp |>
        dplyr::select(dplyr::all_of(c("Feature", comp_cols))) |>
        dplyr::arrange(dplyr::desc(!!rank)) |>
        tidyr::pivot_longer(cols = -Feature, names_to = "Composition",
            values_to = "Percentage") |>
        dplyr::mutate(Composition = factor(Composition,
            levels = unique(c(as.character(rank), comp_cols)) |> rev()
        )) |>
        dplyr::mutate(Feature = factor(Feature, levels = rev(unique(Feature))))

    p <- var_decomp_longer |>
        dplyr::filter(Feature %in% features) |>
        ggplot(aes(x = Feature, y = Percentage, fill = Composition)) +
        geom_bar(stat = "identity", position = "stack") +
        labs(x = "Features", y = "% Variance explained") +
        scale_fill_manual(values = c(
            "Group" = "#FED439FF",
            "Time" = "#709AE1FF",
            "Subject" = "#00A087FF",
            "Group:Time" = "#E69F00",
            "Subject:Time" = "#E69F00",
            "Residual" = "grey"
        )) +
        theme_custom(base_size = fontsize) +
        theme(
            axis.text.x = element_text(angle = 90, hjust = 0.5, vjust = 1),
            legend.position = "right"
        ) +
        coord_flip()

    if (!show_ylab) {
        p <- p + theme(axis.text.y = element_blank())
    }

    return(p)
}
