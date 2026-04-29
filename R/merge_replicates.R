#' Merge replicates
#' @description Calculate mean of replicates for each feature at each time 
#' point for each group
#'
#' @param se_obj_list A list of SummarizedExperiment objects created by 
#' `split_groups()`, each corresponds to one group of samples.
#'
#' @import SummarizedExperiment
#' @import magrittr
#'
#' @returns A list of SummarizedExperiment objects containing the mean of 
#' replicates for each feature at each time point for each group. Each object 
#' in the list corresponds to one group of samples.
#' @export
#' @examples
#' data("example")
#' example_obj <- normalise_to_start(example_obj)
#' example_obj_list <- split_groups(example_obj)
#' example_obj_merged_list <- merge_replicates(example_obj_list)
#' example_obj_merged <- merge_groups(example_obj_merged_list)
merge_replicates <- function(se_obj_list) {
    se_obj_merged_list <- list()

    for (i in names(se_obj_list)) {
        input <- se_obj_list[[i]]

        # calculate mean of replicates for one group
        assay_list <- list()
        for (assay in seq(1, length(assays(input)))) {
            d_mean <- data.frame(row.names = row.names(input))
            for (j in sort(unique(colData(input)$Time))) {
                sp <- input[, input$Time == j]

                j_mean <- assays(sp)[[assay]] %>%
                    apply(1, function(x) mean(x, na.rm = TRUE))
                j_mean[is.nan(j_mean)] <- NA

                d_mean[[as.character(j)]] <- j_mean
            }

            assay_list[[assay]] <- d_mean
        }

        # coldata for merged replicates
        coldata_group <- data.frame(
            row.names = colnames(d_mean),
            Sample = colnames(d_mean),
            Group = i,
            Time = colnames(d_mean) %>% as.numeric()
        )

        se_obj_merged <- SummarizedExperiment(
            assays = assay_list,
            colData = coldata_group
        )

        se_obj_merged_list[[i]] <- se_obj_merged
    }

    return(se_obj_merged_list)
}
