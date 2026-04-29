#' Plot variance decomposition
#'
#' @description Based on `PALMO::variancefeaturePlot()` function
#'
#' @param var_decomp A data frame of variance decomposition result from
#' `variance_decomp()` function
#' @param rank Rank the plot by "Group" or "Time"
#' @param features Features to be plotted (default is NULL)
#' @param top_n Top n features ranked by `rank` to be plotted when `features`
#' is not specified (default is 20)
#' @param ylab Whether to show y axis text (default is TRUE)
#' @param fontsize Font size for the plot (default is 8)
#'
#' @importFrom dplyr arrange desc rename select mutate filter pull
#' @import ggplot2
#' @import magrittr
#' @importFrom dplyr slice
#'
#' @returns A plot of variance decomposition for each feature, showing the
#' percentage of variance explained by Group, Time, and Residual. The features
#' are ranked by the specified `rank` variable (Group or Time) and only the
#' top n features are plotted if `features` is not specified.
#' @export
#' @examples
#' data("example")
#' example_obj <- normalise_to_start(example_obj)
#'
#' var_decomp <- decomp_variance(example_obj, assay = 1)
#' plot_variance(var_decomp, rank = "Time", top_n = 20)
plot_variance <- function(var_decomp, rank = c("Group", "Time"),
    features = NULL,
    top_n = 20, ylab = TRUE, fontsize = 8) {

    rank <- sym(rank) # or Group
    other <- setdiff(c("Group", "Time"), as.character(rank))

    if (is.null(features)) {
        if (top_n > nrow(var_decomp)) {
            top_n <- nrow(var_decomp)
            message("Specified top_n is larger than the number of features. ",
                "Plotting all ", top_n, " features.")
        } else {
            message("Features not specified. Plotting top ", top_n,
                " features ", "ranked by ", as.character(rank), ".")
        }
        features <- var_decomp %>%
            arrange(desc(!!rank)) %>%
            slice(seq(1, top_n)) %>%
            pull(Feature)
    } else if (!all(features %in% var_decomp$Feature)) {
        stop("At least one of the specified features is not found in the ",
        "input.")
    }

    var_decomp_longer <- var_decomp %>%
        dplyr::select(Feature, Group, Time, Residual) %>%
        arrange(desc(!!rank)) %>%
        tidyr::pivot_longer(cols = -Feature, names_to = "Composition",
            values_to = "Percentage") %>%
        mutate(Composition = factor(Composition,
            levels = c(as.character(rank), other, "Residual") %>% rev()
        )) %>%
        mutate(Feature = factor(Feature, levels = rev(unique(Feature))))

    p <- var_decomp_longer %>%
        filter(Feature %in% features) %>%
        ggplot(aes(x = Feature, y = Percentage, fill = Composition)) +
        geom_bar(stat = "identity", position = "stack") +
        labs(x = "Features", y = "% Variance explained") +
        scale_fill_manual(values = c(
            "Group" = "#FED439FF", # pal_simpsons()(2)[1],
            "Time" = "#709AE1FF", # pal_simpsons()(2)[2],
            "Residual" = "grey"
        )) +
        theme_custom(base_size = fontsize) +
        theme(
            axis.text.x = element_text(angle = 90, hjust = 0.5, vjust = 1),
            legend.position = "right"
        ) +
        coord_flip()

    if (!ylab) {
        p <- p + theme(axis.text.y = element_blank())
    }

    print(p)
}
