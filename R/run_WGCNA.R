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
#' data(example_obj)
#'
#' wgcna_input <- prepare_WGCNA(example_obj, assay = 2, powers = seq(1, 30),
#'     networkType = "signed", RsquaredCut = 0.8)
#' wgcna_input$plot
#' picked_power <- wgcna_input$custom_powerEstimate
#' example_net <- run_WGCNA(wgcna_input,
#'     power = picked_power,
#'     minModuleSize = 10, # only 100 genes in the example data
#'     numericLabels = TRUE)
#' @references https://github.com/edo98811/WGCNA_official_documentation/blob/main/FemaleLiver-03-relateModsToExt.R
run_WGCNA <- function(wgcna_input, power, numericLabels = TRUE, ...) {
    .check_list(wgcna_input, "wgcna_input", "prepare_WGCNA")
    .check_logical(numericLabels, "numericLabels")
    .check_positive_int(power, "power")

    dots <- list(...)
    network_type <- wgcna_input$networkType
    if ("networkType" %in% names(dots)) {
        if (!identical(dots$networkType, network_type)) {
            stop("'networkType' does not match the networkType stored in ",
                "'wgcna_input'.")
        }
        dots$networkType <- NULL
    }

    # Check minModuleSize
    minModuleSize <- dots$minModuleSize
    if (is.null(minModuleSize)) {
        # Default is min(20, ncol(datExpr)/2)
        minModuleSize <- eval(formals(WGCNA::blockwiseModules)$minModuleSize,
            envir = list(datExpr = wgcna_input$data))
    }
    .check_positive_int(minModuleSize, "minModuleSize")

    if (ncol(wgcna_input$data) < minModuleSize) {
        stop(
            "Too few features (", ncol(wgcna_input$data),
            ") for WGCNA analysis (need >= ", minModuleSize, ")."
        )
    }

    WGCNA::allowWGCNAThreads()

    .cor_orig <- cor
    cor <- WGCNA::cor
    on.exit(cor <- .cor_orig)
    net <- do.call(WGCNA::blockwiseModules, c(
        list(
            datExpr = wgcna_input$data,
            power = power,
            networkType = network_type,
            numericLabels = numericLabels
        ),
        dots
    ))

    net$input_data <- wgcna_input$data
    net$sample_info <- wgcna_input$sample_info # colData(se_obj)

    net$parameters <- c(
        list(
            power = power,
            networkType = network_type,
            numericLabels = numericLabels
        ),
        dots
    )

    return(net)
}
