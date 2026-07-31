#' Plot coefficient of variation (CV)
#'
#' @description Plot the distribution of coefficient of variation (CV) for
#' each feature with replicates. The CV is calculated as the standard
#' deviation divided by the mean of the abundance values.
#'
#' Note: CV is only meaningful for positive-valued data.
#' @param se_obj A SummarizedExperiment object, produced by `create_input()`
#' function, containing the abundance data and associated sample information.
#' @param fontsize Font size for the plot (default is 8)
#' (default is 8).
#' @import SummarizedExperiment
#' @import ggplot2
#' @returns A plot showing the distribution of coefficient of variation (CV)
#' for each feature, grouped by Time and Group. If there are no replicates, a
#' message will be printed indicating that CV cannot be calculated.
#' @export
#' @examples
#' data(example_obj)
#' plot_cv(example_obj)
plot_cv <- function(se_obj, fontsize = 8) {
    .check_se(se_obj)
    .check_positive(fontsize, "fontsize")
    if (unique(se_obj$Replicate) |> length() < 2) {
        message("CV cannot be calculated without replicates.")
        return(NULL)
    }

    if (any(assays(se_obj)[[1]] <= 0, na.rm = TRUE)) {
        message(
            "Data contains non-positive values. CV is only meaningful for ",
            "positive data."
        )
    }

    calc_cv <- function(x, na.rm = TRUE) stats::sd(x, na.rm = na.rm) /
        mean(x, na.rm = na.rm)

    cv_tb <- assays(se_obj)[[1]] |>
        as.data.frame() |>
        tibble::rownames_to_column("Feature") |>
        tidyr::pivot_longer(cols = -Feature) |>
        merge(colData(se_obj), by.x = "name", by.y = "Sample") |>
        as.data.frame() |>
        dplyr::mutate(Time = as.factor(Time)) |>
        dplyr::group_by(Time, Group, Feature) |>
        dplyr::summarise(CV = calc_cv(value), .groups = "keep")

    p <- cv_tb |>
        dplyr::filter(!is.na(CV) & !is.infinite(CV)) |>
        ggplot(aes(x = Time, y = CV)) +
        geom_violin() +
        geom_boxplot(
            position = "dodge", color = "black",
            outliers = FALSE, width = 0.3
        ) +
        facet_grid(Group ~ Time, scales = "free_x") +
        ylab("Coefficient of variation") +
        theme_custom(base_size = fontsize)
    return(p)
}
