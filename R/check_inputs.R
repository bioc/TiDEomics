# Input validation helpers
# These are internal functions used across the package to validate
# function arguments with descriptive error messages.

#' Check that x is a SummarizedExperiment
#' @param x Object to check.
#' @param arg Name of the argument (for error message).
#' @keywords internal
.check_se <- function(x, arg = "se_obj") {
    if (!methods::is(x, "SummarizedExperiment")) {
        stop("'", arg, "' must be a SummarizedExperiment object.")
    }
}

#' Check that x is a list of SummarizedExperiment objects
#' @param x Object to check.
#' @param arg Name of the argument (for error message).
#' @keywords internal
.check_se_list <- function(x, arg = "se_obj_list") {
    if (!is.list(x) || length(x) == 0) {
        stop("'", arg, "' must be a non-empty list.")
    }
    is_se <- vapply(x, methods::is, logical(1), "SummarizedExperiment")
    if (!all(is_se)) {
        bad <- which(!is_se)
        stop("'", arg, "[[", bad[1],
            "]]' is not a SummarizedExperiment object.")
    }
}

#' Check that x is a data.frame
#' @param x Object to check.
#' @param arg Name of the argument (for error message).
#' @keywords internal
.check_df <- function(x, arg = "data") {
    if (!is.data.frame(x) && !is.matrix(x)) {
        stop("'", arg, "' must be a data.frame or matrix.")
    }
    if (nrow(x) == 0L) {
        stop("'", arg, "' must have at least 1 row.")
    }
}

#' Check that x is a character vector
#' @param x Object to check.
#' @param arg Name of the argument (for error message).
#' @keywords internal
.check_character <- function(x, arg = deparse(substitute(x))) {
    if (!is.character(x)) {
        stop("'", arg, "' must be a character vector.")
    }
    if (length(x) == 0L) {
        stop("'", arg, "' must be non-empty.")
    }
}

#' Check that x is numeric
#' @param x Object to check.
#' @param arg Name of the argument (for error message).
#' @keywords internal
.check_numeric <- function(x, arg = deparse(substitute(x))) {
    if (!is.numeric(x)) {
        stop("'", arg, "' must be numeric.")
    }
}

#' Check that x is a logical flag
#' @param x Object to check.
#' @param arg Name of the argument (for error message).
#' @keywords internal
.check_logical <- function(x, arg = deparse(substitute(x))) {
    if (!is.logical(x) || length(x) != 1 || is.na(x)) {
        stop("'", arg, "' must be TRUE or FALSE.")
    }
}

#' Check that x is a single positive integer
#' @param x Object to check.
#' @param arg Name of the argument (for error message).
#' @keywords internal
.check_positive_int <- function(x, arg = deparse(substitute(x))) {
    if (length(x) != 1 || !is.numeric(x) || is.na(x) ||
            x < 1 || x != floor(x)) {
        stop("'", arg, "' must be a positive integer (>= 1).")
    }
}

#' Check that x is a single number between 0 and 1
#' @param x Object to check.
#' @param arg Name of the argument (for error message).
#' @keywords internal
.check_pval <- function(x, arg = deparse(substitute(x))) {
    if (length(x) != 1 || !is.numeric(x) || is.na(x) || x < 0 || x > 1) {
        stop("'", arg, "' must be a number between 0 and 1.")
    }
}

#' Check that x is a single non-negative number
#' @param x Object to check.
#' @param arg Name of the argument (for error message).
#' @keywords internal
.check_nonneg <- function(x, arg = deparse(substitute(x))) {
    if (length(x) != 1 || !is.numeric(x) || is.na(x) || x < 0) {
        stop("'", arg, "' must be a non-negative number.")
    }
}

#' Check that x is a single positive number
#' @param x Object to check.
#' @param arg Name of the argument (for error message).
#' @keywords internal
.check_positive <- function(x, arg = deparse(substitute(x))) {
    if (!is.numeric(x) || length(x) == 0 || anyNA(x) || any(x <= 0)) {
        stop("'", arg, "' must be positive.")
    }
}

#' Check that x is a list
#' @param x Object to check.
#' @param arg Name of the argument (for error message).
#' @param expected_from Optional: name of the function whose output is expected,
#'   added to the error message.
#' @keywords internal
.check_list <- function(x, arg = deparse(substitute(x)),
                        expected_from = NULL) {
    if (!is.list(x)) {
        msg <- paste0("'", arg, "' must be a list")
        if (!is.null(expected_from)) {
            msg <- paste0(msg, ", the output of ", expected_from, "().")
        } else {
            msg <- paste0(msg, ".")
        }
        stop(msg)
    }
}

#' Check that an SE list has merged data (one value per Group x Time)
#'
#' After `merge_replicates()`, each SE should have exactly one sample
#' per combination of Group and Time.
#'
#' @param se_obj_list A list of SummarizedExperiment objects.
#' @param arg Name of the argument (for error message).
#' @keywords internal
.check_se_list_merged <- function(se_obj_list,
    arg = "se_obj_list") {
    for (nm in names(se_obj_list)) {
        .check_se_merged(se_obj_list[[nm]], paste0(arg, "[[", nm, "]]"))
    }
}

#' Check that an SE has merged data (one value per Group x Time)
#'
#' After `merge_replicates()`, each SE should have exactly one sample
#' per combination of Group and Time.
#'
#' @param se_obj A SummarizedExperiment object.
#' @param arg Name of the argument (for error message).
#' @keywords internal
.check_se_merged <- function(se_obj, arg = "se_obj") {
    cd <- SummarizedExperiment::colData(se_obj)
    n_combos <- length(unique(paste(cd$Group, cd$Time)))
    if (ncol(se_obj) != n_combos) {
        stop("'", arg, "' must contain merged data (one value per ",
            "Group x Time combination). Run merge_replicates() on the ",
            "output of split_groups() first.")
    }
}

#' Check that an SE list has no missing values (imputed)
#'
#' After `impute_groups()`, no assay should contain NA values.
#'
#' @param se_obj_list A list of SummarizedExperiment objects.
#' @param arg Name of the argument (for error message).
#' @keywords internal
.check_se_list_no_na <- function(se_obj_list,
    arg = "se_obj_list") {
    for (nm in names(se_obj_list)) {
        if (anyNA(SummarizedExperiment::assay(se_obj_list[[nm]]))) {
            stop("'", arg, "' contains missing values. ",
                "Run impute_groups() first.")
        }
    }
}

#' Check that an SE has calc_feature_property results in rowData
#'
#' @param se_obj A SummarizedExperiment object.
#' @param arg Name of the argument (for error message).
#' @keywords internal
.check_se_has_properties <- function(se_obj,
    arg = "se_obj") {
    required <- c("Exp_ratio", "P_trend", "Max_FC")
    missing_cols <- setdiff(required,
                            colnames(SummarizedExperiment::rowData(se_obj)))
    if (length(missing_cols) > 0) {
        stop("'", arg, "' missing rowData columns: ",
            paste(missing_cols, collapse = ", "),
            ". Run calc_feature_property() first.")
    }
}

#' Check that an SE list has calc_feature_property results
#'
#' @param se_obj_list A list of SummarizedExperiment objects.
#' @param arg Name of the argument (for error message).
#' @keywords internal
.check_se_list_has_properties <- function(se_obj_list,
    arg = "se_obj_list") {
    for (nm in names(se_obj_list)) {
        .check_se_has_properties(se_obj_list[[nm]], arg)
    }
}

#' Check that a data.frame is from summarise_Trendy
#'
#' @param x A data.frame.
#' @param arg Name of the argument (for error message).
#' @keywords internal
.check_df_trendy_summary <- function(x, arg = deparse(substitute(x))) {
    required <- c("Feature", "Group", "Pattern")
    missing_cols <- setdiff(required, colnames(x))
    if (length(missing_cols) > 0) {
        stop("'", arg, "' must be the output of summarise_Trendy(). ",
            "Missing columns: ", paste(missing_cols, collapse = ", "))
    }
}

#' Check that a data.frame is from summarise_feature_property
#'
#' @param x A data.frame.
#' @param arg Name of the argument (for error message).
#' @keywords internal
.check_df_feature_property <- function(x,
    arg = deparse(substitute(x))) {
    required <- c("Feature", "Group", "Exp_ratio", "P_trend", "Max_FC")
    missing_cols <- setdiff(required, colnames(x))
    if (length(missing_cols) > 0) {
        stop("'", arg, "' must be the output of ",
            "summarise_feature_property(). ",
            "Missing columns: ", paste(missing_cols, collapse = ", "))
    }
}

#' Warn if assay data looks like raw counts rather than log-transformed values
#'
#' limma expects log-transformed input. Raw counts (integers, large max)
#' produce unreliable results. This check warns if the data has
#' characteristics of un-logged counts: all non-negative, >90% of values
#' are integer-like, and max value > 100.
#'
#' @param mat A numeric matrix of expression values.
#' @param assay_name Name of the assay (for the warning message).
#' @keywords internal
.check_limma_input <- function(mat, assay_name) {
    vals <- as.numeric(as.matrix(mat))
    vals <- vals[!is.na(vals)]

    if (length(vals) == 0) return(invisible(NULL))

    all_nonneg <- all(vals >= 0)
    if (!all_nonneg) return(invisible(NULL))

    int_like <- mean(abs(vals - round(vals)) < 1e-8)
    large_max <- max(vals) > 100

    if (int_like > 0.9 && large_max) {
        warning(
            "Assay '", assay_name, "' looks like raw counts (non-negative, ",
            "mostly integers, max = ", round(max(vals)), "). ",
            "limma expects log-transformed values. ",
            "Use log2(x + 1) or similar transformation before ",
            "creating the input with create_input(). ",
            "Set trend = TRUE for RNA-seq count-derived data.",
            call. = FALSE
        )
    }
    invisible(NULL)
}
