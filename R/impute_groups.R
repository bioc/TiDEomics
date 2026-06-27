#' Impute missing values
#' @description Impute missing values for each group of samples in a list of
#' SummarizedExperiment objects. By default, replaces NA with the minimum
#' non-NA value per group / subject. 
#' A custom function can be supplied for other strategies.
#'
#' The input samples can contain replicates, or merged replicates (mean of
#' replicates).
#'
#' @param se_obj_list A list of SummarizedExperiment objects, created by
#' `split_groups()` or `merge_replicates()`, each corresponds to one group
#' of samples.
#' @param fun Function applied to the non-NA values to generate the
#'   replacement value. Default: `min`. Use `function(x) min(x) / 2`
#'   for half-minimum, `median` for median imputation, etc.
#' @param impute_by `"group"` (default): compute replacement value from
#'   all non-NA values in the group. `"subject"`: compute per-subject
#'   replacement value (requires Subject column). If a subject is all
#'   NA, fall back to group-level replacement.
#'
#' @import SummarizedExperiment
#'
#' @returns A list of SummarizedExperiment objects with missing values
#' imputed. Each object corresponds to one group of samples.
#' @export
#' @examples
#' data(example_obj)
#' example_obj <- normalise_to_start(example_obj)
#' example_obj_list <- split_groups(example_obj)
#' example_obj_merged_list <- merge_replicates(example_obj_list)
#' example_obj_merged_imp_list <- impute_groups(example_obj_merged_list)
impute_groups <- function(se_obj_list, fun = min,
    impute_by = c("group", "subject")) {
    .check_se_list(se_obj_list, "se_obj_list")
    .check_se_list_merged(se_obj_list, "se_obj_list")
    if (!is.function(fun)) stop("'fun' must be a function, e.g. min or median.")
    impute_by <- match.arg(impute_by)
    se_obj_imp_list <- list()

    has_subject <- "Subject" %in% colnames(colData(se_obj_list[[1]]))

    if (impute_by == "subject" && !has_subject) {
        warning("impute_by = 'subject' but no Subject column found. ",
            "Falling back to group-level imputation.")
        impute_by <- "group"
    }

    if (impute_by == "subject") {
        message("Imputing per-subject: each subject's missing values ", 
            "replaced by subject-level statistic.")
    }

    for (i in names(se_obj_list)) {
        input <- se_obj_list[[i]]

        # Report missing value status
        n_na <- sum(is.na(assays(input)[[1]]))
        n_total <- length(as.matrix(assays(input)[[1]]))
        if (n_na > 0) {
            message(sprintf("Group %s: %d missing values (%.1f%%).",
                i, n_na, 100 * n_na / n_total))
        } else {
            message(sprintf("Group %s: no missing values.", i))
        }

        assay_list <- list()
        for (assay in seq_along(assays(input))) {
            M_na <- assays(input)[[assay]]
            M_imp <- M_na

            if (impute_by == "subject") {
                subjects <- unique(colData(input)$Subject)
                for (s in subjects) {
                    s_cols <- which(colData(input)$Subject == s)
                    s_vals <- M_na[, s_cols, drop = FALSE]
                    s_non_na <- s_vals[!is.na(s_vals)]
                    s_imp_val <- if (length(s_non_na) > 0) {
                        fun(s_non_na)
                    } else {
                        fun(M_imp[!is.na(M_imp)])
                    }
                    na_mask <- is.na(M_imp[, s_cols, drop = FALSE])
                    M_imp[, s_cols][na_mask] <- s_imp_val
                }
            } else {
                M_imp[is.na(M_imp)] <- fun(M_imp[!is.na(M_imp)])
            }

            assay_list[[assay]] <- M_imp
        }

        se_obj_imp <- input
        assays(se_obj_imp) <- assay_list
        se_obj_imp_list[[i]] <- se_obj_imp
    }

    return(se_obj_imp_list)
}
