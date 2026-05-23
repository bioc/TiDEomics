#' Calculate mean and SD
#' @description Calculate mean and standard deviation for each group 
#' and time point
#'
#' @param se_obj A SummarizedExperiment object with assays containing 
#' expression data. The first assay should be the original data, and 
#' the second (if present) should be normalized data.
#'
#' @import SummarizedExperiment
#' @import magrittr
#' @returns A list containing two data frames: one for original values 
#' and one for normalized values (if available). Each data frame includes 
#' columns for Mean, SD, Group, Time, and Feature.
#' @export
#' @examples
#' data("example")
#' table_mean_sd_list <- calc_mean_sd(example_obj)
#' table_mean_sd <- table_mean_sd_list$norm0
#' plot_trend(table_mean_sd, 
#'     features = sample(unique(table_mean_sd$Feature), 4))
calc_mean_sd <- function(se_obj) {
    # function to calculate mean and sd for each group and time point
    mean_sd <- function(se_obj, i, j, assay) {
        # assay 1 is original, 2 is normalised to time 0
        # table = expression matrix, i = group, j = time
        sp <- se_obj[, se_obj$Group == i & se_obj$Time == j]

        table_mean <- assays(sp)[[assay]] %>%
            apply(1, function(x) mean(x, na.rm = TRUE)) %>%
            as.data.frame()
        colnames(table_mean) <- "Mean"
        table_mean$SD <- assays(sp)[[assay]] %>%
            apply(1, function(x) stats::sd(x, na.rm = TRUE))
        table_mean$Group <- i
        table_mean$Time <- j
        table_mean$Feature <- row.names(table_mean)
        return(table_mean)
    }

    # norm and non-norm values
    table_mean_sd_orig <- data.frame(
        Mean = NULL, SD = NULL,
        Group = NULL, Time = NULL, Feature = NULL
    )
    table_mean_sd_norm <- data.frame(
        Mean = NULL, SD = NULL,
        Group = NULL, Time = NULL, Feature = NULL
    )
    for (i in unique(colData(se_obj)$Group)) {
        for (j in sort(unique(colData(se_obj)$Time))) {
            table_orig <- mean_sd(se_obj, i, j, assay = 1) # not normalised
            table_mean_sd_orig <- rbind(table_mean_sd_orig, table_orig)

            if (length(assays(se_obj)) == 2) {
                table_norm <- mean_sd(se_obj, i, j, assay = 2) # t0 normalised
                table_mean_sd_norm <- rbind(table_mean_sd_norm, table_norm)
            }
        }
    }

    # replace NaN with NA
    table_mean_sd_orig[is.nan(table_mean_sd_orig$Mean), "Mean"] <- NA
    table_mean_sd_orig[is.nan(table_mean_sd_orig$SD), "SD"] <- NA
    table_mean_sd_orig$Group <- factor(table_mean_sd_orig$Group,
        levels = levels(colData(se_obj)$Group))

    if (length(assays(se_obj)) == 2) {
        table_mean_sd_norm[is.nan(table_mean_sd_norm$Mean), "Mean"] <- NA
        table_mean_sd_norm[is.nan(table_mean_sd_norm$SD), "SD"] <- NA
        table_mean_sd_norm$Group <- factor(table_mean_sd_norm$Group,
            levels = levels(colData(se_obj)$Group))
    }

    return(list(orig = table_mean_sd_orig, norm0 = table_mean_sd_norm))
}
