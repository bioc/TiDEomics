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
#' @param assay Which assay slot of the SummarizedExperiment object to use 
#' for WGCNA input.
#'
#' @import SummarizedExperiment
#' @import magrittr
#' @importFrom dplyr select mutate
#'
#' @returns A data frame with samples in rows and features in columns, filtered 
#' to remove bad samples and features, ready for use in `prepare_WGCNA()`.
#' @keywords internal
.WGCNA_input <- function(se_obj, assay) {
    if (assay > length(assays(se_obj))) {
        stop("The input object does not have assay 2: time 0 normalised ", 
        "data. Please run `normalise_to_start()` to create the time 0 ", 
        "normalised data in assay 2.")
    }

    # WGCNA requires the rows as samples, and columns as features e.g. genes.
    data_wgcna <- assays(se_obj)[[assay]] %>%
        t() %>% #
        as.data.frame()

    goodgenes <- 
        colnames(data_wgcna)[which(WGCNA::goodGenes(data_wgcna) == TRUE)]

    # remove bad features and samples
    data_wgcna <- data_wgcna[WGCNA::goodSamples(data_wgcna), 
        WGCNA::goodGenes(data_wgcna)]

    stopifnot(sum(WGCNA::goodSamples(data_wgcna) == FALSE) == 0)
    stopifnot(sum(WGCNA::goodGenes(data_wgcna) == FALSE) == 0)

    # goodSamples(data_wgcna) %>% summary 
    # goodGenes(data_wgcna) %>% summary 

    return(data_wgcna)
}
