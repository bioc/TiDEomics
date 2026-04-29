#' Normalise to time 0
#' @description Normalise to starting time point, to make the mean of 
#' starting time point replicates 0
#'
#' @param se_obj A SummarizedExperiment object created by `create_input()`
#'
#' @import SummarizedExperiment
#' @import magrittr
#'
#' @returns A SummarizedExperiment object with normalised data in the second 
#' assay slot.
#' @export
#' @examples
#' data("example")
#' example_obj <- normalise_to_start(example_obj)
#' example_obj_list <- split_groups(example_obj)
#' example_obj_merged_list <- merge_replicates(example_obj_list)
#' example_obj_merged <- merge_groups(example_obj_merged_list)
normalise_to_start <- function(se_obj) {
    d_list_0norm <- list()

    for (i in unique(colData(se_obj)$Group)) {
        input <- se_obj[, se_obj$Group == i]

        time_series <- sort(unique(colData(se_obj)$Time))

        mean_by_time <- vapply(time_series, function(t) {
            cols_t <- input[, input$Time == t]
            if (dim(cols_t)[2] == 0) {
                rep(NA_real_, nrow(cols_t))
            } else {
                rowMeans(assays(cols_t)[[1]], na.rm = TRUE)
            }
        }, FUN.VALUE = numeric(nrow(input)))
        # This function checks that all values of FUN are 
        # compatible with the FUN.VALUE,
        # in that they must have the same length and type.

        colnames(mean_by_time) <- as.character(time_series)

        first_non_na <- apply(mean_by_time, 1, function(row) {
            pos <- which(!is.na(row))
            if (length(pos) == 0) NA_real_ else row[pos[1]]
        })

        d_list_0norm[[as.character(i)]] <- sweep(assays(input)[[1]], 
            1, first_non_na, FUN = "-")
    }
    d_0norm <- do.call(cbind, d_list_0norm)

    se_obj@assays@data[[2]] <- d_0norm

    return(se_obj)
}
