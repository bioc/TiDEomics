#' Plot number of identified features
#'
#' @description Plot the number of identified features (i.e., features 
#' with non-missing values) for each sample, with a dashed line indicating 
#' the average number across all samples.
#' @param se_obj A SummarizedExperiment object, produced by `create_input()` 
#' function, containing the abundance data and associated sample information.
#' @param fontsize (Optional) An integer specifying the font size for the plot 
#' (default is 8).
#' @param signif (Optional) A logical value indicating whether to perform 
#' significance testing between groups and add significance annotations to 
#' the plot (default is FALSE).
#' @param ... Additional arguments to be passed to `ggsignif::geom_signif()` 
#' function when `signif` is TRUE, for customizing the significance annotations.
#' @import SummarizedExperiment
#' @import ggplot2
#' @returns A plot showing the ID number for each sample, with a dashed line 
#' indicating the average number across all samples.
#' @export
#' @examples
#' # simulate data with random missing values
#' na_data <- matrix(rnorm(1000), nrow = 100, ncol = 100)
#' na_data[sample(length(na_data), size = 1000)] <- NA
#' na_data <- data.frame(Feature = paste0("Feature", 1:100), na_data)
#' colnames(na_data)[-1] <- paste0("Sample", 1:100)
#'
#' na_obj <- create_input(na_data,
#'     data.frame(Sample = paste0("Sample", 1:100),
#'     Time = rep(rep(1:10, each = 5), 2),
#'     Group = rep(c("A", "B"), each = 50),
#'     Replicate = rep(1:5, 20)))
#' plot_ID(na_obj)
#' plot_missing(na_obj)
plot_ID <- function(se_obj, fontsize = 8, signif = FALSE, ...) {
    .check_se(se_obj)
    .check_logical(signif, "signif")
    .check_positive(fontsize, "fontsize")
    global_id <- sum(!is.na(assays(se_obj)[[1]])) /
        dim(assays(se_obj)[[1]])[2]

    id_tb <- assays(se_obj)[[1]] |>
        dplyr::summarise(dplyr::across(dplyr::everything(), 
            ~ sum(!is.na(.)))) |>
        tidyr::pivot_longer(cols = dplyr::everything(), values_to = "ID") |>
        merge(colData(se_obj), by.x = "name", by.y = "Sample") |>
        as.data.frame() |>
        dplyr::mutate(Time = as.factor(Time))

    overview <- id_tb |>
        ggplot(aes(x = Time, y = ID, fill = Replicate)) +
        geom_bar(stat = "identity", position = "dodge", color = "black") +
        geom_hline(yintercept = global_id, linetype = "dashed") +
        facet_grid(Group ~ Time, scales = "free_x") +
        ggtitle(paste0(
            "# ID (total: ", dim(assays(se_obj)[[1]])[1],
            " features) (average: ",
            global_id |> round(digits = 2), ")"
        )) +
        theme_custom(base_size = fontsize) +
        theme(panel.grid.major.y = element_line(color = "grey")) +
        ylab("# ID")

    p <- id_tb |>
        ggplot(aes(x = Group, y = ID)) +
        geom_boxplot(size = 0.2, outlier.shape = NA) +
        geom_point(aes(color = Time), size = 0.7, 
            position = position_jitter(width = 0.3)) +
        theme_custom(base_size = fontsize)

    if (signif) {
        p <- p +
            ggsignif::geom_signif(
                comparisons = utils::combn(unique(se_obj$Group) |> 
                as.vector(), m = 2) |>
                as.data.frame() |> as.list(), ...)
    }
    res <- list(overview = overview, comparison = p)
    print(res)
    return(invisible(res))
}
