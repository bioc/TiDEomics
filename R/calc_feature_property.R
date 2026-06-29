#' Calculate feature property
#' @description Calculate per-feature temporal properties for each group
#' from the mean of replicates at each time point.
#'
#' @param se_obj_merged_list A list of SummarizedExperiment objects created
#' by `merge_replicates()`, each containing the mean of replicates for each
#' feature at each time point for one group.
#' @param threshold A numeric value to determine whether a feature is
#' "expressed" in samples. (default is NULL, meaning any non-NA value
#' is considered expressed).
#'
#' @import SummarizedExperiment
#'
#' @returns A list of SummarizedExperiment objects, each corresponding to
#' one group, with updated rowData containing the following columns: Feature,
#' Group, Exp_threshold, T_exp (number of time points with values >
#' Exp_threshold), T_total
#' (total number of time points), P_trend (p-value from Bartels' test for
#' randomness against a trend), Max_FC (maximum fold change across time
#' points), Max_FC_time (difference between the time points at which maximum
#' and minimum expression occur; positive if max occurs after min, negative
#' otherwise), Exp_ratio (T_exp / T_total), Rho_time (Spearman correlation
#' with time; positive = up trend), Peak_ratio (number of local maxima
#' / T_total; 0 = no peaks detected),
#' Log2_CV (log2(1 + CV) of expression across time points),
#' AUC (area under the time-0-normalised curve from assay 2; NA if assay 2
#' is not available).
#' @export
#' @examples
#' data(example_obj)
#' example_obj <- normalise_to_start(example_obj)
#' example_obj_list <- split_groups(example_obj)
#' example_obj_merged_list <- merge_replicates(example_obj_list)
#'
#' example_obj_merged_list <-
#'     calc_feature_property(example_obj_merged_list, threshold = 0)
calc_feature_property <- function(se_obj_merged_list, threshold = NULL) {
    .check_se_list(se_obj_merged_list, "se_obj_merged_list")
    .check_se_list_merged(se_obj_merged_list, "se_obj_merged_list")
    if (!is.null(threshold)) .check_numeric(threshold, "threshold")
    .prop_to_df <- function(data, fun, colname) {
        apply(data, 1, fun) |>
            as.data.frame() |>
            stats::setNames(colname)
    }

    for (i in names(se_obj_merged_list)) {
        n_col <- ncol(assays(se_obj_merged_list[[i]])[[1]])
        n_time <- length(unique(colData(se_obj_merged_list[[i]])$Time))
        if (n_col != n_time) {
            stop("Input has multiple values per time point, replicates ",
                "have not been merged. ",
                "Run merge_replicates() before calc_feature_property().")
        }
    }

    filter_func <- function(row) {
        if (is.null(threshold)) {
            return(sum(!is.na(row)))
        } else {
            return(sum(row > threshold))
        }
    }

    for (i in names(se_obj_merged_list)) {
        d_mean <- assays(se_obj_merged_list[[i]])[[1]]

        property_tb <- data.frame(
            row.names = rownames(d_mean),
            Feature = rownames(d_mean),
            Group = i,
            Exp_threshold = ifelse(is.null(threshold), NA, threshold),
            T_exp = apply(d_mean, 1, function(row) filter_func(row)),
            T_total = ncol(d_mean),
            One = apply(d_mean, 1, function(row) (filter_func(row) >= 1)),
            Three = apply(d_mean, 1, function(row) (filter_func(row) >= 3))
        )

        d_mean_1 <- d_mean[property_tb$One, ]
        d_mean_n <- d_mean[property_tb$Three, ]
        time_vals <- as.numeric(colnames(d_mean_1))

        random_pv <- apply(
            d_mean_n, 1,
            function(x) {
                randtests::bartels.rank.test(stats::na.omit(x),
                    alternative = "left.sided")$p.value
            }
        ) |>
            as.data.frame() |>
            stats::setNames("P_trend")

        max_fc <- .prop_to_df(d_mean_1, function(x) {
            max(stats::na.omit(x)) - min(stats::na.omit(x))
        }, "Max_FC")

        max_fc_day <- .prop_to_df(d_mean_1, function(x) {
            time_vals[which.max(x)] - time_vals[which.min(x)]
        }, "Max_FC_time")

        rho_time <- apply(d_mean_n, 1, function(x) {
            idx <- !is.na(x)
            if (sum(idx) < 3) return(NA_real_)
            stats::cor(time_vals[idx], x[idx], method = "spearman")
        }) |>
            as.data.frame() |>
            stats::setNames("Rho_time")

        peak_ratio <- apply(d_mean_n, 1, function(x) {
            idx <- !is.na(x)
            if (sum(idx) < 3) return(NA_real_)
            v <- x[idx]
            n <- length(v)
            peaks <- 0
            for (j in 2:(n - 1)) {
                if (v[j] > v[j - 1] && v[j] > v[j + 1]) peaks <- peaks + 1
            }
            peaks / n
        }) |>
            as.data.frame() |>
            stats::setNames("Peak_ratio")

        temporal_cv <- .prop_to_df(d_mean_n, function(x) {
            m <- mean(x, na.rm = TRUE)
            s <- stats::sd(x, na.rm = TRUE)
            if (m == 0 || is.na(m)) return(NA_real_)
            log2(1 + s / abs(m))
        }, "Log2_CV")

        has_assay2 <- length(assays(se_obj_merged_list[[i]])) >= 2
        if (has_assay2) {
            d_norm <- assays(se_obj_merged_list[[i]])[[2]]
            d_norm_1 <- d_norm[property_tb$One, , drop = FALSE]
            auc_val <- apply(d_norm_1, 1, function(x) {
                idx <- !is.na(x)
                if (sum(idx) < 2) return(NA_real_)
                v <- x[idx]
                t_v <- time_vals[idx]
                sum(diff(t_v) * (v[-1] + v[-length(v)]) / 2)
            }) |>
                as.data.frame() |>
                stats::setNames("AUC")
        } else {
            message("Assay 2 (time-0-normalised) not available. ",
                "AUC will be NA. Run normalise_to_start() first.")
            auc_val <- data.frame(
                row.names = rownames(d_mean_1),
                AUC = rep(NA_real_, nrow(d_mean_1)))
        }

        prop_cols <- list(random_pv, max_fc, max_fc_day, rho_time,
            peak_ratio, temporal_cv, auc_val)
        for (col_df in prop_cols) {
            property_tb <- merge(property_tb, col_df,
                by = "row.names", all.x = TRUE) |>
                tibble::column_to_rownames("Row.names")
        }
        property_tb <- property_tb |>
            dplyr::select(-One, -Three)

        property_tb$Exp_ratio <- property_tb$T_exp / property_tb$T_total

        rowData(se_obj_merged_list[[i]]) <- property_tb
    }

    return(se_obj_merged_list)
}
