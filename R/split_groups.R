#' Split groups
#' @description Split SummarizedExperiment object by groups
#'
#' @param se_obj A SummarizedExperiment object created by `create_input()`
#'
#' @import SummarizedExperiment
#'
#' @returns A list of SummarizedExperiment objects for each group in the input 
#' object. Each object in the list corresponds to one group of samples.
#' @export
#' @examples
#' data("example")
#' example_obj <- normalise_to_start(example_obj)
#' example_obj_list <- split_groups(example_obj)
#' example_obj_merged_list <- merge_replicates(example_obj_list)
#' example_obj_merged <- merge_groups(example_obj_merged_list)
split_groups <- function(se_obj) {
    se_obj_list <- list()

    for (i in unique(se_obj$Group)) {
        input <- se_obj[, se_obj$Group == i]
        se_obj_list[[i]] <- input
    }

    return(se_obj_list)
}
