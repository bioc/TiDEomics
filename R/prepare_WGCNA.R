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
#' @param assay The assay to use in the SummarizedExperiment object: a numeric
#'   index or character name, e.g. 1 or "orig" for original data, 2 or "norm"
#'   for time 0 normalised data. (Default: 2.)
#' @param networkType Parameter of `WGCNA::pickSoftThreshold()`
#' (default is "signed")
#' @param RsquaredCut Parameter of `WGCNA::pickSoftThreshold()`
#' (default is 0.8)
#' @param MeanConnectivity Line of mean connectivity
#' (default is 100)
#' @param powers Parameter of `WGCNA::pickSoftThreshold()`
#' (default: NULL, auto-assigned as `seq(1, 20)`)
#' @param fontsize Font size for the plot (default is 8)
#' @param ... Additional parameters to be passed to `WGCNA::pickSoftThreshold()`
#'
#' @returns A list containing results of the scale-free topology
#' fit indices with different powers, suggested power, network type and prepared
#' input data used for reuse in `run_WGCNA()`
#' @export
#' @examples
#' data(example_obj)
#' example_obj <- normalise_to_start(example_obj)
#'
#' wgcna_input <- prepare_WGCNA(example_obj, assay = 2, powers = seq(1, 30),
#'     networkType = "signed", RsquaredCut = 0.8)
#' wgcna_input$fitIndices
#' picked_power <- wgcna_input$powerEstimate
#' @references https://github.com/edo98811/WGCNA_official_documentation/
prepare_WGCNA <- function(
    se_obj, assay = 2,
    networkType = "signed",
    RsquaredCut = 0.8,
    MeanConnectivity = 100,
    powers = NULL,
    fontsize = 8, ...
) {
    .check_se(se_obj)
    .check_character(networkType, "networkType")
    networkType <- match.arg(networkType,
        c("unsigned", "signed", "signed hybrid"))
    .check_pval(RsquaredCut, "RsquaredCut")
    .check_nonneg(MeanConnectivity, "MeanConnectivity")
    .check_positive(fontsize, "fontsize")
    if (!is.null(powers)) .check_positive(powers, "powers")

    assay <- .match_assay(assay, se_obj)
    WGCNA::allowWGCNAThreads()
    data_wgcna <- .WGCNA_input(se_obj, assay = assay)

    if (is.null(powers)) {
        powers <- seq(1, 20)
    } else if (length(powers) == 0 || !all(powers %% 1 == 0) ||
            !all(powers > 0)) {
        stop("Powers must be a non-empty vector of positive integers.")
    }

    .cor_orig <- cor
    cor <- WGCNA::cor
    on.exit(cor <- .cor_orig)
    sft <- WGCNA::pickSoftThreshold(data_wgcna,
        powerVector = powers,
        verbose = 5,
        networkType = networkType,
        RsquaredCut = RsquaredCut,
        ...
    )

    # for reuse in `run_WGCNA()`
    sft$data <- data_wgcna
    sft$sample_info <- colData(se_obj)
    sft$networkType <- networkType

    # Diagnostic plots: scale-free topology fit + mean connectivity
    fit_tb <- data.frame(
        Power       = sft$fitIndices[, 1],
        SignedR2    = -sign(sft$fitIndices[, 3]) * sft$fitIndices[, 2],
        MeanConn    = sft$fitIndices[, 5]
    )

    p1 <- ggplot(fit_tb, aes(x = Power, y = SignedR2)) +
        geom_text(aes(label = Power), size = fontsize * 0.35) +
        geom_hline(yintercept = RsquaredCut, colour = "red", linewidth = 0.5) +
        labs(x = "Soft Threshold (power)",
            y = expression(Signed~R^2),
            title = "Scale independence") +
        theme_custom(base_size = fontsize)

    p2 <- ggplot(fit_tb, aes(x = Power, y = MeanConn)) +
        geom_text(aes(label = Power), size = fontsize * 0.35) +
        geom_hline(yintercept = MeanConnectivity, colour = "red",
                linewidth = 0.5) +
        labs(x = "Soft Threshold (power)",
            y = "Mean Connectivity",
            title = "Mean connectivity") +
        theme_custom(base_size = fontsize)

    sft$plot <- p1 / p2
    return(sft)
}
