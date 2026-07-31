#' Prepare WGCNA input data
#'
#' @description Prepare input data for WGCNA, including transposing the data
#' to have samples in rows and features in columns, and removing bad samples
#' and features. Features can be pre-filtered, e.g. by residual variance
#' calculated by `decomp_variance()`, to remove noisy features before
#' preparing the data for WGCNA.
#'
#' @param se_obj A SummarizedExperiment object. Data
#' normalised to time point 0 can be in the second assay slot, created by
#' `normalise_to_start()`.
#' @param assay The assay to use in the SummarizedExperiment object: a numeric
#'   index or character name, e.g. 1 or "orig" for original data, 2 or "norm"
#'   for time 0 normalised data. (Default: 2.)
#'
#' @import SummarizedExperiment
#'
#' @returns A data frame with samples in rows and features in columns, filtered
#' to remove bad samples and features, ready for use in `prepare_WGCNA()`.
#' @keywords internal
.WGCNA_input <- function(se_obj, assay) {
    assay <- .match_assay(assay, se_obj)

    # WGCNA requires the rows as samples, and columns as features e.g. genes.
    data_wgcna <- assays(se_obj)[[assay]] |>
        t() |>
        as.data.frame()

    # remove bad features and samples
    data_wgcna <- data_wgcna[WGCNA::goodSamples(data_wgcna),
        WGCNA::goodGenes(data_wgcna)]

    stopifnot(sum(WGCNA::goodSamples(data_wgcna) == FALSE) == 0)
    stopifnot(sum(WGCNA::goodGenes(data_wgcna) == FALSE) == 0)

    return(data_wgcna)
}
