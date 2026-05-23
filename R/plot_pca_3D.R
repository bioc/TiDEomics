#' Plot PCA in 3D
#'
#' @description Plot principal component analysis (PCA) results in 3D. This 
#' function takes the output PCA object of the `plot_pca` function and 
#' visualizes the samples in a 3D space defined by the specified principal 
#' components.
#'
#' The samples are coloured by Group and sized by Time.
#'
#' @param pca A PCA object returned by the `plot_pca()` function.
#' @param pcs A numeric vector specifying which three principal components 
#' to plot (default is 1:3)
#'
#' @import magrittr
#'
#' @returns An interactive 3D PCA plot showing the distribution of samples 
#' in the space defined by the specified principal components. The samples 
#' are coloured by Group and sized by Time.
#' @export
#' @examples
#' data("example")
#' PC = plot_pca(example_obj, morepc = seq(1, 3))
#' plot_pca_3D(PC, pcs = seq(1, 3))
plot_pca_3D <- function(pca, pcs = seq(1, 3)) {
    if (!requireNamespace("plotly", quietly = TRUE))
        stop("Package 'plotly' is required for 3D PCA. ", 
        "Install with: install.packages('plotly')")
    if (length(pcs) != 3) {
        stop("Please specify three principal components to plot.")
    }
    if (!all(pcs %in% seq(1, length(pca$variance)))) {
        stop("Invalid PCs specified. Please choose from 1 to ", 
            length(pca$variance))
    }

    pc <- as.data.frame(pca$rotated)
    pc$Sample <- rownames(pc)
    pc <- merge(pc, pca$metadata, by = "Sample")

    xlab <- sprintf("PC%d: %.2f%%", pcs[1], 
        pca$variance[[paste0("PC", pcs[1])]])
    ylab <- sprintf("PC%d: %.2f%%", pcs[2], 
        pca$variance[[paste0("PC", pcs[2])]])
    zlab <- sprintf("PC%d: %.2f%%", pcs[3], 
        pca$variance[[paste0("PC", pcs[3])]])

    # plot 3 PCs
    p <- plotly::plot_ly(pc,
        x = pc[[paste0("PC", pcs[1])]],
        y = pc[[paste0("PC", pcs[2])]],
        z = pc[[paste0("PC", pcs[3])]],
        color = ~Group, size = ~Time, text = ~Sample,
        type = "scatter3d", mode = "markers",
        colors = get_custom_palette(levels(pca$metadata$Group))
    ) %>%
        plotly::layout(
            scene = list(
                xaxis = list(title = xlab),
                yaxis = list(title = ylab),
                zaxis = list(title = zlab)
            )
        )
    return(p)
}
