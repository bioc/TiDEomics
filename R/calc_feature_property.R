#' Calculate feature property
#' @description Calculate randomness p-value, maximum fold change along time
#' course, and ratio of expressed time points for each feature in each group,
#' based on the mean of replicates at each time point.
#'
#' @param se_obj_merged_list A list of SummarizedExperiment objects created
#' by `merge_replicates()`, each containing the mean of replicates for each
#' feature at each time point for one group.
#' @param threshold A numeric value to determine whether a feature is
#' "expressed" in samples. (default is NULL, meaning any non-NA value
#' is considered expressed).
#'
#' @import SummarizedExperiment
#' @import magrittr
#'
#' @returns A list of SummarizedExperiment objects, each corresponding to
#' one group, with updated rowData containing the following columns: Feature,
#' Group, Exp_threshold, T_exp (number of time points with values >
#' Exp_threshold), T_total
#' (total number of time points), P_trend (p-value from Bartels' test for
#' randomness against a trend), Max_FC (maximum fold change across time
#' points), Max_FC_time (difference between the time points at which maximum
#' and minimum expression occur; positive if max occurs after min, negative
#' otherwise), Exp_ratio (T_exp / T_total).
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
calc_feature_property <- function(se_obj_merged_list, threshold = NULL) {
    # Validate that input has been through merge_replicates
    for (i in names(se_obj_merged_list)) {
        n_col <- ncol(assays(se_obj_merged_list[[i]])[[1]])
        n_time <- length(unique(colData(se_obj_merged_list[[i]])$Time))
        if (n_col != n_time) {
            stop("Input has multiple values per time point, replicates ",
                "have not been merged. ",
                "Run merge_replicates() before calc_feature_property().")
        }
    }

    for (i in names(se_obj_merged_list)) {
        d_mean <- assays(se_obj_merged_list[[i]])[[1]]

        # threshold to consider a feature as "expressed"
        # If NULL, any non-NA value is considered expressed.
        filter_func <- function(row) {
            if (is.null(threshold)) {
                return(sum(!is.na(row)))
            } else {
                return(sum(row > threshold))
            }
        }

        property_tb <- data.frame(
            row.names = rownames(d_mean),
            Feature = rownames(d_mean),
            Group = i,
            Exp_threshold = ifelse(is.null(threshold), NA, threshold),
            # non-missing count
            T_exp = apply(d_mean, 1, function(row) filter_func(row)),
            # total time points
            T_total = ncol(d_mean),
            # features with at least 1 value
            One = apply(d_mean, 1, function(row) (filter_func(row) >= 1)),
            # features with at least 3 values
            Three = apply(d_mean, 1, function(row) (filter_func(row) >= 3))
        )

        d_mean_1 <- d_mean[property_tb$One, ]
        d_mean_n <- d_mean[property_tb$Three, ]

        # bartel's test for randomness, against a trend, not oscillation
        # need at least 3 values
        random_pv <- apply(
            d_mean_n, 1,
            function(x) {
                randtests::bartels.rank.test(stats::na.omit(x),
                    alternative = "left.sided")$p.value
            }
        ) %>%
            as.data.frame() %>%
            set_colnames("P_trend")
        # alternative "two.sided": the null hypothesis of randomness is
        # tested against nonrandomness. "left.sided": the null hypothesis
        # of randomness is tested against a trend. "right.sided": the null
        # hypothesis of randomness is tested against a systematic oscillation.

        # max fold change and the time points of max and min
        max_fc <- apply(d_mean_1, 1, function(x) {
            max(stats::na.omit(x)) - min(stats::na.omit(x))
        }) %>%
            as.data.frame() %>%
            set_colnames("Max_FC")

        time_vals <- as.numeric(colnames(d_mean_1))
        max_fc_day <- apply(d_mean_1, 1, function(x) {
            time_vals[which.max(x)] - time_vals[which.min(x)]
        }) %>%
            as.data.frame() %>%
            set_colnames("Max_FC_time")

        property_tb <- merge(property_tb, random_pv,
            by = "row.names", all.x = TRUE) %>%
            tibble::column_to_rownames("Row.names") %>%
            merge(max_fc, by = "row.names", all.x = TRUE) %>%
            tibble::column_to_rownames("Row.names") %>%
            merge(max_fc_day, by = "row.names", all.x = TRUE) %>%
            tibble::column_to_rownames("Row.names") %>%
            dplyr::select(-One, -Three)

        property_tb$Exp_ratio <- property_tb$T_exp / property_tb$T_total

        rowData(se_obj_merged_list[[i]]) <- property_tb
    }

    return(se_obj_merged_list)
}
