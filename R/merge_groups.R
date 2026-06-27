#' Merge groups into one object
#' @description Merge SummarizedExperiment object of different groups into one
#'
#' @param se_obj_list A list of SummarizedExperiment objects, such as output of 
#' `split_groups()` and `merge_replicates()`. Names of the list components are 
#' the group names.
#'
#' @import SummarizedExperiment
#'
#' @returns A SummarizedExperiment object containing all samples from the input 
#' list. The 'Group' column in the colData will indicate the group of each 
#' sample. New sample names will be prefixed with the group name.
#' @export
#' @examples
#' data(example_obj)
#' example_obj <- normalise_to_start(example_obj)
#' example_obj_list <- split_groups(example_obj)
#' example_obj_merged_list <- merge_replicates(example_obj_list)
#' example_obj_merged <- merge_groups(example_obj_merged_list)
merge_groups <- function(se_obj_list) {
    .check_se_list(se_obj_list, "se_obj_list")
    .check_se_list_merged(se_obj_list, "se_obj_list")
    assay_list <- list()
    coldata_list <- list()

    for (i in names(se_obj_list)) {
        se_obj <- se_obj_list[[i]]

        # assay data
        for (j in seq_along(assays(se_obj))) {
            df <- assays(se_obj)[[j]]
            colnames(df) <- paste0(i, "_", colnames(df))

            if (j > length(assay_list)) {
                assay_list[[j]] <- df
            } else {
                assay_list[[j]] <- dplyr::full_join(
                    assay_list[[j]] |> 
                    as.data.frame() |>
                    tibble::rownames_to_column(".rowname"),
                    df |> 
                    as.data.frame() |>
                    tibble::rownames_to_column(".rowname"),
                    by = ".rowname"
                ) |>
                    tibble::column_to_rownames(".rowname") |>
                    as.matrix()
            }
        }

        # colData
        cd <- colData(se_obj) |> as.data.frame()
        cd$Group <- i # unnecessary
        cd$Sample <- paste0(i, "_", cd$Sample)
        coldata_list[[i]] <- cd
    }

    cd_all <- do.call(rbind, coldata_list)
    rownames(cd_all) <- cd_all$Sample
    # keep order
    cd_all$Group <- factor(cd_all$Group, levels = names(se_obj_list)) 

    # keep format consistent
    if (!("Replicate" %in% names(cd_all))) {
        cd_all$Replicate <- 1
    }
    if (!("Batch" %in% names(cd_all))) {
        cd_all$Batch <- 1
    }

    se_obj_merged <- SummarizedExperiment(
        assays = assay_list,
        colData = cd_all
    )

    return(se_obj_merged)
}
