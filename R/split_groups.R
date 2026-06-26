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
split_groups <- function(se_obj) {
    .check_se(se_obj)
    if (!"Group" %in% colnames(colData(se_obj))) {
        stop("'se_obj' must have a 'Group' column in colData.")
    }
    se_obj_list <- list()

    for (i in unique(se_obj$Group)) {
        input <- se_obj[, se_obj$Group == i]
        se_obj_list[[i]] <- input
    }

    return(se_obj_list)
}
