#' Prepare data and choose power for WGCNA
#'
#' @description Prepare input data for WGCNA (format, QC),
#' then choose the appropriate soft thresholding power for
#' network construction, by analysing scale free topology
#' with different soft thresholding powers.
#'
#' @param se_obj A SummarizedExperiment object. Data
#' normalised to time point 0 can be in the second assay slot, created by
#' `normalise_to_start()`. Features can be pre-filtered, e.g. by residual
#' variance calculated by `decomp_variance()`, to remove noisy features
#' before running WGCNA.
#' @param assay Which assay slot of the SummarizedExperiment object to use
#' for WGCNA input (default is 2, which is where the time 0 normalised data is
#' stored by `normalise_to_start()`)
#' @param networkType (Optional) Parameter of `WGCNA::pickSoftThreshold()`
#' (default is "signed")
#' @param RsquaredCut (Optional) Parameter of `WGCNA::pickSoftThreshold()`
#' (default is 0.8)
#' @param MeanConnectivity (Optional) Line of mean connectivity
#' (default is 100)
#' @param powers (Optional) Parameter of `WGCNA::pickSoftThreshold()`
#' (default is `c(seq(1, 10, by = 1), seq(12, 20, by = 2))`)
#' @param ... Additional parameters to be passed to `WGCNA::pickSoftThreshold()`
#'
#' @importFrom graphics par plot text abline
#' @importFrom grDevices dev.off
#'
#' @returns A list containing results of the scale-free topology
#' fit indices with different powers, suggested power, network type and prepared
#' input data used for reuse in `run_WGCNA()`
#' @export
#' @examples
#' data("example")
#' example_obj <- normalise_to_start(example_obj)
#'
#' wgcna_input <- prepare_WGCNA(example_obj, assay = 2, powers = seq(1, 30),
#'     networkType = "signed", RsquaredCut = 0.8)
#' wgcna_input$fitIndices
#' picked_power <- wgcna_input$powerEstimate
#' # example_net <- run_WGCNA(wgcna_input,
#' #    power = picked_power,
#' #    minModuleSize = 10, # only 100 genes in the example data
#' #    numericLabels = TRUE)
#' # plot_WGCNA(example_net, fontsize = 8)
#' @references https://github.com/edo98811/WGCNA_official_documentation/
prepare_WGCNA <- function(
    se_obj, assay = 2,
    networkType = "signed",
    RsquaredCut = 0.8,
    MeanConnectivity = 100,
    powers = NULL, ...
) {
    WGCNA::allowWGCNAThreads()
    data_wgcna <- WGCNA_input(se_obj, assay = assay)

    if (is.null(powers)) {
        powers <- c(seq(1, 10, by = 1), seq(12, 20, by = 2))
    } else if (!all(powers %% 1 == 0) | !all(powers > 0)) {
        stop("Powers must be a vector of positive integers.")
    }

    cor <- WGCNA::cor
    sft <- WGCNA::pickSoftThreshold(data_wgcna,
        powerVector = powers,
        verbose = 5,
        networkType = networkType,
        RsquaredCut = RsquaredCut,
        ...
    )
    cor <- stats::cor

    # for reuse in `run_WGCNA()`
    sft$data <- data_wgcna
    sft$sample_info <- colData(se_obj)
    sft$networkType <- networkType

    par(mfrow = c(2, 1))
    cex1 <- 0.9

    # Scale-free topology fit index as a function of the soft-thresholding power
    plot(sft$fitIndices[, 1], -sign(sft$fitIndices[, 3]) * sft$fitIndices[, 2],
        xlab = "Soft Threshold (power)",
        ylab = "Scale Free Topology Model Fit, signed R^2", type = "n",
        main = paste("Scale independence")
    )
    text(sft$fitIndices[, 1], -sign(sft$fitIndices[, 3]) * sft$fitIndices[, 2],
        labels = sft$fitIndices[, 1], cex = cex1, col = "black"
    )
    # this line corresponds to using an R^2 cut-off of h
    abline(h = RsquaredCut, col = "red")

    # Mean connectivity as a function of the soft-thresholding power
    plot(sft$fitIndices[, 1], sft$fitIndices[, 5],
        xlab = "Soft Threshold (power)",
        ylab = "Mean Connectivity", type = "n",
        main = paste("Mean connectivity")
    )
    text(sft$fitIndices[, 1], sft$fitIndices[, 5],
        labels = sft$fitIndices[, 1], cex = cex1, col = "black"
    )
    abline(h = MeanConnectivity, col = "red")

    return(sft)
}
