#' Impute missing values
#' @description Impute missing values for each group of samples in a list of 
#' SummarizedExperiment objects, with the minimum value in the group.
#'
#' The input samples can contain replicates, or merged replicates (mean of 
#' replicates).
#'
#' @param se_obj_list A list of SummarizedExperiment objects, created by 
#' `split_groups()` or `merge_replicates()`, each corresponds to one group 
#' of samples.
#'
#' @import SummarizedExperiment
#' @import magrittr
#' @importFrom dplyr bind_rows
#'
#' @returns A list of SummarizedExperiment objects with the missing values 
#' imputed for each group. Each object in the list corresponds to one group 
#' of samples.
#' @export
#' @examples
#' data("example")
#' example_obj <- normalise_to_start(example_obj)
#' example_obj_list <- split_groups(example_obj)
#' example_obj_merged_list <- merge_replicates(example_obj_list)
#' example_obj_merged_imp_list <- impute_groups(example_obj_merged_list)
impute_groups <- function(se_obj_list) {
    se_obj_imp_list <- list()

    for (i in names(se_obj_list)) {
        input <- se_obj_list[[i]]

        # calculate mean of replicates for one group
        assay_list <- list()
        for (assay in seq(1, length(assays(input)))) {
            M_na <- assays(input)[[assay]]
            M_imp <- M_na
            M_imp[is.na(M_imp)] <- min(M_imp, na.rm = TRUE)

            assay_list[[assay]] <- M_imp
        }

        se_obj_imp <- input
        assays(se_obj_imp) <- assay_list

        se_obj_imp_list[[i]] <- se_obj_imp
    }

    return(se_obj_imp_list)
}
