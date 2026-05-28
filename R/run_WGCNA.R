#' Weighted gene co-expression network analysis
#'
#' @description Run WGCNA based on output of `prepare_WGCNA()`, to
#' identify co-expression modules of features with `WGCNA::blockwiseModules()`
#'
#' @param wgcna_input A list output by `prepare_WGCNA()`, containing the
#' prepared input data and networkType parameter used in power selection.
#' @param power Soft-thresholding power to be used in
#' `WGCNA::blockwiseModules()`, selected automatically or manually based
#' on the output of `prepare_WGCNA()`
#' @param numericLabels Whether to use numeric labels for modules in the output
#' (default is TRUE)
#' @param ... Additional parameters to be passed to `WGCNA::blockwiseModules()`
#'
#' @returns The built network and parameters of `WGCNA::blockwiseModules()`,
#' and the input data and sample information for use in `plot_WGCNA()`.
#' @export
#' @examples
#' data("example")
#' example_obj <- normalise_to_start(example_obj)
#'
#' wgcna_input <- prepare_WGCNA(example_obj, assay = 2, powers = seq(1, 30),
#'     networkType = "signed", RsquaredCut = 0.8)
#' wgcna_input$fitIndices
#' picked_power <- wgcna_input$powerEstimate
#' example_net <- run_WGCNA(wgcna_input,
#'     power = picked_power,
#'     minModuleSize = 10, # only 100 genes in the example data
#'     numericLabels = TRUE)
#' # plot_WGCNA(example_net, fontsize = 8)
#' # use_data(example_net)
#' @references https://github.com/edo98811/WGCNA_official_documentation/blob/main/FemaleLiver-03-relateModsToExt.R
run_WGCNA <- function(wgcna_input, power, numericLabels = TRUE, ...) {
    WGCNA::allowWGCNAThreads()

    .cor_orig <- cor
    cor <- WGCNA::cor
    net <- WGCNA::blockwiseModules(wgcna_input$data,
        power = power,
        networkType = wgcna_input$networkType,
        numericLabels = numericLabels,
        ...
    )
    cor <- .cor_orig

    net$input_data <- wgcna_input$data
    net$sample_info <- wgcna_input$sample_info # colData(se_obj)

    net$parameters <- list(
        power = power,
        networkType = wgcna_input$networkType,
        numericLabels = numericLabels,
        ...
    )

    return(net)
}
