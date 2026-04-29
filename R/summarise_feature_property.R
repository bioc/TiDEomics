#' Summarise feature properties
#'
#' @description This function takes a list of merged `SummarizedExperiment` 
#' objects, which are the output of the `calc_feature_property()` function, 
#' and summarizes the feature properties across all groups. It extracts the 
#' row data from each `SummarizedExperiment` object, combines them into a 
#' single data frame, and returns this summary data frame. Each row in the 
#' resulting data frame represents a feature, and the columns contain the 
#' properties of that feature for each group, including the feature name, 
#' group name, proportion of NA values for that feature in that group, 
#' randomness p-value, and maximum fold change over time.
#'
#' @param se_obj_merged_list A list of merged `SummarizedExperiment` objects, 
#' output of `calc_feature_property()` function
#' @import SummarizedExperiment
#' @import magrittr
#'
#' @returns A data frame summarizing the feature properties across all groups, 
#' with each row representing a feature and columns containing the properties 
#' including the feature name, group name, and the proportion of NA values for 
#' that feature in that group, randomness p-value, and maximum fold change 
#' over time.
#' @export
#' @examples
#' data("example")
#' example_obj <- normalise_to_start(example_obj)
#' example_obj_list <- split_groups(example_obj)
#' example_obj_merged_list <- merge_replicates(example_obj_list)
#'
#' example_obj_merged_list <- 
#'     calc_feature_property(example_obj_merged_list, threshold = 0)
#' property_random_fc <- summarise_feature_property(example_obj_merged_list)
summarise_feature_property <- function(se_obj_merged_list) {
    if (rowData(se_obj_merged_list[[1]]) %>% colnames() %>% length() < 5) {
        stop("The input object does not contain the expected properties. ", 
        "Please make sure to input the output of `calc_feature_property` ", 
        "function.")
    }

    property_list <- list()
    for (i in names(se_obj_merged_list)) {
        property_list[[i]] <- rowData(se_obj_merged_list[[i]]) %>% 
            as.data.frame()
    }
    property_random_fc <- do.call(rbind, property_list)
    row.names(property_random_fc) <- NULL

    return(property_random_fc)
}
