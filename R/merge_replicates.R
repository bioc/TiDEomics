#' Merge replicates
#' @description Calculate mean of replicates for each feature at each time
#' point for each group.
#'
#' When a Subject column is present, merging is done in two stages:
#' 1. Average replicates within each Group-Subject-Time combination.
#' 2. Average across subjects within each Group-Time combination.
#' This gives equal weight to each subject regardless of replicate count.
#'
#' Without a Subject column, all samples at the same Group-Time are averaged
#' together in a single step.
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
    has_subject <- "Subject" %in% colnames(colData(se_obj_list[[1]]))

    if (has_subject) {
        message("Subject column detected. Merging within-subject first, ", 
        "then across subjects. ")
    }

    for (i in names(se_obj_list)) {
        input <- se_obj_list[[i]]

        assay_list <- list()
        for (assay in seq(1, length(assays(input)))) {
            d_mean <- data.frame(row.names = row.names(input))

            for (j in sort(unique(colData(input)$Time))) {
                sp <- input[, input$Time == j]

                if (has_subject) {
                    # Stage 1: average within each Subject per Group, Time
                    subjects_j <- unique(colData(sp)$Subject)
                    subj_means <- vapply(subjects_j, function(s) {
                        sp_s <- sp[, sp$Subject == s]
                        rowMeans(assays(sp_s)[[assay]], na.rm = TRUE)
                    }, FUN.VALUE = numeric(nrow(sp)))
                    # Stage 2: average across subjects
                    j_mean <- rowMeans(subj_means, na.rm = TRUE)
                } else {
                    # Average all samples at per Group, Time
                    j_mean <- assays(sp)[[assay]] %>%
                        apply(1, function(x) mean(x, na.rm = TRUE))
                }

                j_mean[is.nan(j_mean)] <- NA
                d_mean[[as.character(j)]] <- j_mean
            }

            assay_list[[assay]] <- d_mean
        }

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
