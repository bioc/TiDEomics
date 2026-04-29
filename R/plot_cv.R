#' Plot coefficient of variation (CV)
#'
#' @description Plot the distribution of coefficient of variation (CV) for 
#' each feature with replicates. The CV is calculated as the standard 
#' deviation divided by the mean of the abundance values.
#' @param se_obj A SummarizedExperiment object, produced by `create_input()` 
#' function, containing the abundance data and associated sample information.
#' @param fontsize (Optional) An integer specifying the font size for the plot 
#' (default is 8).
#' @import SummarizedExperiment
#' @import magrittr
#' @import ggplot2
#' @importFrom dplyr summarise mutate
#' @returns A plot showing the distribution of coefficient of variation (CV) 
#' for each feature, grouped by Time and Group. If there are no replicates, a 
#' message will be printed indicating that CV cannot be calculated.
#' @export
#' @examples
#' data("example")
#' plot_cv(example_obj)
plot_cv <- function(se_obj, fontsize = 8) {
    if (unique(se_obj$Replicate) %>% length() < 2) {
        message("CV cannot be calculated without replicates.")
        return(NULL)
    }

    calc_cv <- function(x, na.rm = TRUE) sd(x, na.rm = na.rm) / 
        mean(x, na.rm = na.rm)

    cv_tb <- assays(se_obj)[[1]] %>%
        mutate(Feature = row.names(.)) %>%
        tidyr::pivot_longer(cols = -Feature) %>%
        merge(., colData(se_obj), by.x = "name", by.y = "Sample") %>%
        as.data.frame() %>%
        mutate(Time = as.factor(Time)) %>%
        group_by(Time, Group, Feature) %>%
        summarise(CV = calc_cv(value), .groups = "keep")

    (cv_tb %>%
        filter(!is.na(CV)) %>%
        ggplot(aes(x = Time, y = CV)) +
        geom_violin() +
        geom_boxplot(
            position = "dodge", color = "black",
            outliers = FALSE, width = 0.3
        ) +
        facet_grid(Group ~ Time, scales = "free_x") +
        ylab("Coeffecient of variation") +
        theme_custom(base_size = fontsize)) %>% print()
}
