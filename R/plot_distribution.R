#' Abundance distribution plot
#'
#' @description This function generates density plots to visualise the 
#' distribution of abundance values across samples. The user can choose 
#' to facet the plot by "Group", "Time", or "Sample" to compare 
#' distributions across different conditions or samples. If no faceting 
#' variable is specified, a single density plot will be generated for all 
#' samples combined.
#'
#' @param se_obj A SummarizedExperiment object, produced by `create_input()` 
#' function, containing the abundance data and associated sample information.
#' @param facet_by (Optional) A character string specifying the variable to 
#' facet the plot by. It can be one of "Group", "Time", or "Sample". If NULL, 
#' no faceting will be applied and a single density plot will be generated for 
#' all samples (default is NULL).
#' @param fontsize (Optional) An integer specifying the font size for the plot 
#' (default is 8).
#'
#' @import SummarizedExperiment
#' @import magrittr
#' @import ggplot2
#'
#' @returns A plot showing the density distribution of abundance values, 
#' optionally faceted by the specified variable.
#' @export
#' @examples
#' data("example")
#' plot_distribution(example_obj)
#' plot_distribution(example_obj, facet_by = "Group")
plot_distribution <- function(se_obj, facet_by = NULL, fontsize = 8) {
    if (is.null(facet_by)) {
        (assays(se_obj)[[1]] %>%
            tidyr::pivot_longer(cols = dplyr::everything()) %>%
            ggplot(aes(x = value, group = name)) +
            geom_density() +
            ggtitle("Abundance distribution") +
            labs(x = "Log2 Abundance", y = "Density") +
            theme_custom(base_size = fontsize)) %>% print()
    } else {
        if (!facet_by %in% c("Group", "Time", "Sample")) {
            stop("Please specify 'facet_by' as either 'Group', 'Time', ", 
            "or 'Sample'.")
        }

        if (facet_by == "Group") {
            (assays(se_obj)[[1]] %>%
                tidyr::pivot_longer(cols = dplyr::everything()) %>%
                merge(., colData(se_obj), by.x = "name", by.y = "Sample") %>%
                as.data.frame() %>%
                dplyr::mutate(Time = as.factor(Time)) %>%
                ggplot(aes(x = value, y = Time, fill = Replicate)) +
                ggridges::geom_density_ridges(color = "black", alpha = 0.3) +
                facet_wrap(~Group, nrow = 1) +
                ggtitle("Abundance distribution") +
                xlab("Log2 Abundance") +
                theme_custom(base_size = fontsize)) %>% print()
        } else if (facet_by == "Time") {
            (assays(se_obj)[[1]] %>%
                tidyr::pivot_longer(cols = dplyr::everything()) %>%
                merge(., colData(se_obj), by.x = "name", by.y = "Sample") %>%
                as.data.frame() %>%
                dplyr::mutate(Time = as.factor(Time)) %>%
                ggplot(aes(x = value, y = Group, fill = Replicate)) +
                ggridges::geom_density_ridges(color = "black", alpha = 0.3) +
                facet_wrap(~Time, nrow = 1) +
                ggtitle("Abundance distribution") +
                xlab("Log2 Abundance") +
                theme_custom(base_size = fontsize)) %>% print()
        } else if (facet_by == "Sample") {
            (assays(se_obj)[[1]] %>%
                tidyr::pivot_longer(cols = dplyr::everything()) %>%
                ggplot(aes(x = value)) +
                geom_density() +
                facet_wrap(~name) +
                ggtitle("Abundance distribution") +
                labs(x = "Log2 Abundance", y = "Density") +
                theme_custom(panel_border = TRUE, base_size = fontsize)) %>% 
                print()
        }
    }
}
