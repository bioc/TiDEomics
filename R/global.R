utils::globalVariables(c(
    "Feature", "Group", "Time", "Replicate", "Batch", "Sample",
    "Row.names", "Comparison", "Pattern", "Module",
    "Cond1", "Cond2", "logFC", "adj.P.Val",
    "One", "Three", "Feature", "PTID", "Residual",
    "Composition", "Percentage", "Day", "Count",
    "Sample", "Abundance", "Cluster", "id", "p.adjust",
    "Description", "term", "PC1", "PC2",
    "x_start", "x_end", "y_start", "y_end",
    "Mean", "SD", "V1", "V2", "ID", "N", "Sample_new",
    "Exp_ratio", "SYMBOL", "CV", "value", "Missing",
    "Adjusted.P.value", "name", ".", "Breakpoint", "Freq", "Color"
))

theme_custom <- function(
    base_size = 8, panel_border = FALSE,
    legend_position = "right", x_text_angle = 0
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
            # axis.text.x = element_text(angle = x_text_angle,
            #     hjust = ifelse(x_text_angle > 5, 1, NULL)),
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
