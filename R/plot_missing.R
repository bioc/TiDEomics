#' Plot missing rate
#'
#' @description Plot the ratio of missing values for each sample
#' @param se_obj A SummarizedExperiment object, produced by `create_input()`
#' function, containing the abundance data and associated sample information.
#' @param fontsize (Optional) An integer specifying the font size for the plot
#' (default is 8).
#' @param signif (Optional) A logical value indicating whether to perform
#' significance testing between groups and add significance annotations
#' to the plot (default is FALSE).
#' @param ... Additional arguments to be passed to
#' `ggpubr::stat_compare_means()` when `signif` is TRUE,
#' for customizing the significance annotations.
#' @import SummarizedExperiment
#' @import ggplot2
#' @returns A plot showing the missing value ratio for each sample, with a
#' dashed line indicating the global missing value ratio across all samples.
#' And a boxplot comparing the missing value ratio between groups, with optional
#' significance annotations.
#' @export
#' @examples
#' # simulate data with random missing values
#' na_data <- matrix(rnorm(1500), nrow = 100, ncol = 150)
#' na_data[sample(length(na_data), size = 2000)] <- NA
#' na_data <- data.frame(Feature = paste0("Feature", 1:100), na_data)
#' colnames(na_data)[-1] <- paste0("Sample", 1:150)
#'
#' na_obj <- create_input(na_data,
#'     data.frame(Sample = paste0("Sample", 1:150),
#'     Time = rep(rep(1:10, each = 5), 3),
#'     Group = rep(c("A", "B", "C"), each = 50),
#'     Replicate = rep(1:5, 30)))
#' plot_missing(na_obj, signif = TRUE)
plot_missing <- function(se_obj, fontsize = 8, signif = FALSE, ...) {
    .check_se(se_obj)
    .check_logical(signif, "signif")
    .check_positive(fontsize, "fontsize")
    global_missing_rate <- sum(is.na(assays(se_obj)[[1]])) /
        dim(assays(se_obj)[[1]])[1] / dim(assays(se_obj)[[1]])[2]

    missing_tb <- assays(se_obj)[[1]] |>
        as.data.frame() |>
        dplyr::summarise(dplyr::across(dplyr::everything(),
            ~ sum(is.na(.)) / length(.))) |>
        tidyr::pivot_longer(cols = dplyr::everything(),
            values_to = "Missing") |>
        merge(colData(se_obj), by.x = "name", by.y = "Sample") |>
        as.data.frame() |>
        dplyr::mutate(Time = as.factor(Time))

    overview <- missing_tb |>
        ggplot(aes(x = Time, y = Missing, fill = Replicate)) +
        geom_bar(stat = "identity", position = "dodge", color = "black") +
        geom_hline(yintercept = global_missing_rate, linetype = "dashed") +
        facet_grid(Group ~ Time, scales = "free_x") +
        ggtitle(paste0(
            "NA ratio (out of ", dim(assays(se_obj)[[1]])[1],
            " features) (global: ",
            global_missing_rate |> round(digits = 2), ")"
        )) +
        theme_custom(base_size = fontsize) +
        theme(panel.grid.major.y = element_line(color = "grey")) +
        ylab("Missing rate")

    p <- missing_tb |>
        ggplot(aes(x = Group, y = Missing)) +
        geom_boxplot(size = 0.2, outlier.shape = NA) +
        geom_point(aes(color = Time), size = 0.7,
        position = position_jitter(width = 0.3)) +
        theme_custom(base_size = fontsize)

    if (signif) {
        p <- p +
            ggpubr::stat_compare_means(
                comparisons = utils::combn(unique(se_obj$Group) |>
                as.vector(), m = 2) |>
                as.data.frame() |> as.list(), ...)
    }
    res <- list(overview = overview, comparison = p)
    return(res)
}
