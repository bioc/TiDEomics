#' Normalise to starting time point
#' @description Normalise to starting time point, to make the mean of
#' starting time point samples 0.
#'
#' Two modes are available:
#' - `by_subject = FALSE` (default): computes the group-level mean at the
#'   first non-NA time point per feature and subtracts it. All samples within
#'   a group share the same baseline.
#' - `by_subject = TRUE`: computes each subject's value at the first non-NA
#'   time point per feature and subtracts that subject-specific baseline.
#'   Removes between-subject baseline differences. Use when each subject has
#'   their own initial condition. Requires a Subject column in colData.
#'
#' @param se_obj A SummarizedExperiment object created by `create_input()`
#' @param by_subject If FALSE (default), use group-level baseline. If TRUE,
#'   use subject-level baseline (requires Subject column in colData).
#'
#' @import SummarizedExperiment
#'
#' @returns A SummarizedExperiment object with normalised data in the second
#' assay slot.
#' @export
#' @examples
#' data(example_obj)
#' example_obj <- normalise_to_start(example_obj)
normalise_to_start <- function(se_obj, by_subject = FALSE) {
    .check_se(se_obj)
    .check_logical(by_subject, "by_subject")
    if (by_subject && !"Subject" %in% colnames(colData(se_obj))) {
        stop("by_subject = TRUE requires a 'Subject' column in colData. ",
            "Set by_subject = FALSE for group-level normalisation, or ",
            "use create_input(..., subject_col = '<column_name>') to ",
            "add a Subject column.")
    }
    d_list_0norm <- list()

    if (by_subject) {
        # ---- Subject-level normalisation ----
        subjects <- unique(colData(se_obj)$Subject)
        message(
            "Normalising to subject-level baseline at each feature's first ",
            "non-NA time point for ", length(subjects), " subjects."
        )

        for (subj in subjects) {
            input <- se_obj[, se_obj$Subject == subj]
            grp <- unique(colData(input)$Group)

            input_grp <- input[, input$Group == grp]
            time_series <- sort(unique(colData(input_grp)$Time))

            # Per-time-point means for this subject
            mean_by_time <- vapply(time_series, function(t) {
                cols_t <- input_grp[, input_grp$Time == t]
                if (dim(cols_t)[2] == 0) {
                    rep(NA_real_, nrow(cols_t))
                } else {
                    rowMeans(assays(cols_t)[[1]], na.rm = TRUE)
                }
            }, FUN.VALUE = numeric(nrow(input_grp)))
            colnames(mean_by_time) <- as.character(time_series)

            first_non_na <- apply(mean_by_time, 1, function(row) {
                pos <- which(!is.na(row))
                if (length(pos) == 0) NA_real_ else row[pos[1]]
            })

            d_list_0norm[[paste0(subj)]] <- sweep(
                assays(input_grp)[[1]], 1, first_non_na, FUN = "-"
            )
        }
    } else {
        # ---- Group-level normalisation (original behaviour) ----
        message(
            "Normalising to group baseline at each feature's ",
            "first non-NA time point."
        )

        for (i in unique(colData(se_obj)$Group)) {
            input <- se_obj[, se_obj$Group == i]

            time_series <- sort(unique(colData(input)$Time))

            mean_by_time <- vapply(time_series, function(t) {
                cols_t <- input[, input$Time == t]
                if (dim(cols_t)[2] == 0) {
                    rep(NA_real_, nrow(cols_t))
                } else {
                    rowMeans(assays(cols_t)[[1]], na.rm = TRUE)
                }
            }, FUN.VALUE = numeric(nrow(input)))
            colnames(mean_by_time) <- as.character(time_series)

            first_non_na <- apply(mean_by_time, 1, function(row) {
                pos <- which(!is.na(row))
                if (length(pos) == 0) NA_real_ else row[pos[1]]
            })

            n_later <- sum(
                apply(mean_by_time, 1, function(row) {
                    pos <- which(!is.na(row))
                    length(pos) > 0 && pos[1] > 1
                }),
                na.rm = TRUE
            )
            if (n_later > 0) {
                message(
                    "Group ", i, ": ", n_later,
                    " feature(s) has NA in the first time point, and used",
                    " first non-NA at a later time point to normalise."
                )
            }

            d_list_0norm[[as.character(i)]] <- sweep(assays(input)[[1]],
                1, first_non_na, FUN = "-")
        }
    }

    d_0norm <- do.call(cbind, unname(d_list_0norm))
    assay(se_obj, "norm") <- as.matrix(d_0norm)

    return(se_obj)
}
