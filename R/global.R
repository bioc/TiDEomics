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

#' Resolve assay argument
#'
#' Accepts assay as a name or index, validates it against the object,
#' and returns the resolved value. Open to any assay present in the object.
#'
#' @param assay A numeric index or character name.
#' @param se_obj A SummarizedExperiment object.
#' @return The resolved assay identifier (name if available, otherwise index).
#' @keywords internal
.match_assay <- function(assay, se_obj) {
    if (length(assay) != 1 || !(is.numeric(assay) || is.character(assay))) {
        stop("'assay' must be a single numeric index or character name.")
    }
    # Resolve numeric indices to names when names are available
    assay_names <- names(SummarizedExperiment::assays(se_obj))
    if (!is.null(assay_names) && length(assay_names) > 0 && is.numeric(assay)) {
        if (!assay %in% seq_along(assay_names)) {
            stop("Assay index ", assay, " out of range. Available: ",
                paste(seq_along(assay_names), assay_names,
                    sep = " = ", collapse = ", "))
        }
        assay <- assay_names[assay]
    }
    # Validate character names
    if (is.character(assay)) {
        if (is.null(assay_names)) {
            stop("Assays are not named. Use a numeric index instead.")
        }
        if (!assay %in% assay_names) {
            stop("Assay '", assay, "' not found. Available: ",
                paste(assay_names, collapse = ", "))
        }
    }
    assay
}

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
    .check_character(legend_position, "legend_position")
    legend_position <- match.arg(legend_position,
        c("right", "left", "bottom", "top", "none"))
    .check_logical(panel_border, "panel_border")
    .check_positive(base_size, "base_size")
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
