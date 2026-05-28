utils::globalVariables(c(
    "Feature", "Group", "Time", "Replicate", "Batch", "Sample", "Subject",
    "Row.names", "Comparison", "Pattern", "Module",
    "Cond1", "Cond2", "logFC", "adj.P.Val",
    "One", "Three", "PTID", "Residual",
    "Composition", "Percentage", "Day", "Count",
    "Abundance", "Cluster", "id", "p.adjust",
    "Description", "term", "PC1", "PC2",
    "x_start", "x_end", "y_start", "y_end",
    "Mean", "SD", "V1", "V2", "ID", "N", "Sample_new",
    "Exp_ratio", "SYMBOL", "CV", "value", "Missing",
    "Adjusted.P.value", "name", ".", "Breakpoint", "Freq", "Color",
    "Database", "MeanConn", "Power", "SignedR2",
    "Combined.Score", "Genes", "Term"
))

#' Custom ggplot2 theme
#'
#' A minimal theme used by all TiDEomics plotting functions, built on
#' `theme_minimal()`. Exported so users can apply it to their own plots
#' for consistent styling.
#'
#' @param base_size Base font size (default: 8).
#' @param panel_border Logical, whether to draw a border around the panel
#'   (default: FALSE).
#' @param legend_position Legend position (default: "right").
#'
#' @return A ggplot2 theme object.
#' @import ggplot2

#' @export
#'
#' @examples
#' library(ggplot2)
#' ggplot(mtcars, aes(wt, mpg)) +
#'     geom_point() +
#'     theme_custom(base_size = 10)
theme_custom <- function(
    base_size = 8, panel_border = FALSE,
    legend_position = "right"
) {
    half_line <- base_size / 2
    if (panel_border) {
        panel.border <- element_rect(color = "black", fill = NA,
            linewidth = 0.5)
        axis.line <- element_blank()
    } else {
        panel.border <- element_blank()
        axis.line <- element_line(color = "black", linewidth = 0.5)
    }
    theme_minimal(base_size = base_size) +
        theme(
            text = element_text(size = base_size, colour = "black"),
            plot.title = element_text(size = base_size + 1, face = "bold",
                hjust = 0, margin = margin(b = half_line)),
            plot.subtitle = element_text(size = base_size,
                margin = margin(b = half_line)),
            axis.title = element_text(size = base_size + 1),
            axis.text = element_text(size = base_size),
            axis.ticks = element_line(linewidth = 0.5),
            axis.line = axis.line,
            panel.background = element_rect(fill = NA, colour = NA),
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            panel.border = panel.border,
            legend.background = element_blank(),
            legend.key = element_blank(),
            legend.title = element_text(size = base_size),
            legend.text = element_text(size = base_size),
            legend.position = legend_position,
            strip.text = element_text(size = base_size, face = "bold",
                lineheight = 1.5),
            plot.margin = margin(6, 6, 6, 6)
        )
}
