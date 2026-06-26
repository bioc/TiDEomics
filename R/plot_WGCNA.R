#' Plot WGCNA results
#'
#' @description Plot WGCNA module eigengenes and module-group/time correlations
#' based on output of `run_WGCNA()`
#'
#' @param net WGCNA network object output by `run_WGCNA()`
#' @param fontsize Font size for plots (default is 8)
#'
#'
#' @returns Plots of WGCNA module dendrogram, module eigengenes, pairwise
#' scatterplots of eigengenes, clustering of module eigengenes, and
#' module-trait correlation heatmap
#' @export
#' @examples
#' data("example")
#' example_obj <- normalise_to_start(example_obj)
#'
#' # wgcna_input <- prepare_WGCNA(example_obj, assay = 2, powers = seq(1, 30),
#' #     networkType = "signed", RsquaredCut = 0.8)
#' # wgcna_input$fitIndices
#' # picked_power <- wgcna_input$powerEstimate
#' # example_net <- run_WGCNA(wgcna_input,
#' #    power = picked_power,
#' #    minModuleSize = 10, # only 100 genes in the example data
#' #    numericLabels = TRUE)
#' data("example_net")
#' plot_WGCNA(example_net, fontsize = 8)
#' @references https://github.com/edo98811/WGCNA_official_documentation/blob/main/FemaleLiver-03-relateModsToExt.R
plot_WGCNA <- function(net, fontsize = 8) {
    .check_list(net, "net", "run_WGCNA")
    if (!"colors" %in% names(net)) stop("'net' must contain a 'colors' element.")
    .check_positive(fontsize, "fontsize")
    if ("numericLabels" %in% names(net$parameters)) {
        if (net$parameters$numericLabels == TRUE) {
            moduleColors <- WGCNA::labels2colors(net$colors)
        } else {
            moduleColors <- net$colors
        }
    } else {
        moduleColors <- net$colors
    }
    WGCNA::plotDendroAndColors(net$dendrograms[[1]],
        moduleColors[net$blockGenes[[1]]],
        main = "Feature dendrogram and module colors",
        dendroLabels = FALSE, hang = 0.03,
        addGuide = TRUE, guideHang = 0.05
    )

    ## plot eigengenes
    # Module eigengene is defined as the first principal component of
    # the expression matrix of the corresponding module.
    MEs0 <- WGCNA::moduleEigengenes(net$input_data, net$colors)$eigengenes
    MEs <- WGCNA::orderMEs(MEs0)

    ann_row <- net$sample_info |>
        as.data.frame()
    rownames(ann_row) <- NULL
    ann_row <- ann_row |>
        tibble::column_to_rownames("Sample") |>
        dplyr::select(Time, Group) |>
        dplyr::mutate(Time = as.numeric(Time))
    stopifnot(identical(row.names(MEs), row.names(ann_row)))

    col_time_func <- circlize::colorRamp2(
        seq(min(ann_row$Time), max(ann_row$Time), length.out = 3),
        c("#dadaeb", "#9e9ac8", "#54278f")
    )

    ann_colors <- list(
        Group = get_custom_palette(levels(net$sample_info$Group)),
        Time = col_time_func
    )

    ComplexHeatmap::Heatmap(as.matrix(MEs),
        column_title = "WGCNA module eigengenes",
        col = grDevices::colorRampPalette(c("#3C5488FF", "white",
            "#E64B35FF"))(100),
        show_row_names = FALSE,
        show_column_names = TRUE,
        column_title_gp = grid::gpar(fontsize = fontsize + 2,
            fontface = "bold"),
        column_names_gp = grid::gpar(fontsize = fontsize),
        left_annotation = ComplexHeatmap::rowAnnotation(
            df = ann_row,
            col = ann_colors,
            gp = grid::gpar(fontsize = fontsize),
            annotation_name_gp = grid::gpar(fontsize = fontsize,
                fontface = "bold"),
            annotation_legend_param = list(
                title_gp = grid::gpar(fontsize = fontsize, fontface = "bold"),
                labels_gp = grid::gpar(fontsize = fontsize)
            )
        ),
        heatmap_legend_param = list(
            title = "r",
            title_gp = grid::gpar(fontsize = fontsize, fontface = "bold"),
            labels_gp = grid::gpar(fontsize = fontsize)
        )
    ) |> print()

    # Pairwise scatterplots of eigengenes
    WGCNA::plotMEpairs(MEs,
        cex.labels = 1,
        gap = 1 / 10,
        clusterMEs = TRUE
    )

    dissimME <- (1 - t(stats::cor(MEs, method = "p", use = "p"))) / 2
    hclustdatME <- stats::hclust(stats::as.dist(dissimME), method = "average")
    graphics::par(mfrow = c(1, 1))
    graphics::plot(hclustdatME, 
        main = "Clustering based on the module eigengenes")

    ## module-trait correlation

    # test correlation of modules to groups
    datTraits <- ann_row |>
        dplyr::select(Group, Time) |>
        dplyr::mutate(Time = as.numeric(Time)) |>
        WGCNA::binarizeCategoricalColumns(
            convertColumns = c("Group"),
            dropFirstLevelVsAll = FALSE,
            includePairwise = TRUE
        )

    # Compute correlation and p-values
    moduleTraitCor <- stats::cor(MEs, datTraits, use = "p")
    moduleTraitPvalue <- WGCNA::corPvalueStudent(moduleTraitCor,
        nSamples = nrow(net$input_data)) |>
        cut(
            breaks = c(-Inf, 0.001, 0.01, 0.05, Inf),
            labels = c("***", "**", "*", "")
        ) # label with significance stars

    textMatrix <- moduleTraitPvalue
    dim(textMatrix) <- dim(moduleTraitCor)

    # Display correlation heatmap
    graphics::par(mfrow = c(1, 1))
    WGCNA::labeledHeatmap(
        Matrix = moduleTraitCor,
        xLabels = names(datTraits),
        yLabels = names(MEs),
        ySymbols = names(MEs),
        setStdMargins = FALSE,
        cex.text = 0.7,
        cex.lab = 0.7,
        zlim = c(-1, 1),
        textMatrix = textMatrix,
        colors = grDevices::colorRampPalette(c("#3C5488FF", "white",
            "#E64B35FF"))(100),
        main = "Module-Group Correlation"
    )
}
