#' Create object
#'
#' @description Check format of input data and sample annotation, create a
#' SummarizedExperiment object
#'
#' @param data A data frame with rows as features (e.g., genes, proteins)
#' and columns as samples. The first column should contain feature identifiers
#' (e.g., gene symbols) and named as 'Feature'. The data should have been log
#' transformed (and normalised if needed). Missing values are allowed in some
#' of the downstream analyses.
#' @param sample_ann A data frame containing sample annotations with required
#' columns: 'Sample', 'Group', and 'Time'. 'Replicate' and 'Batch' are optional
#' columns and will be auto-generated if not provided. 'Time',
#' 'Replicate' and 'Batch' should be numeric.
#' @param subject_col Optional: name of a column in `sample_ann` identifying
#' biological subjects measured repeatedly across time points (e.g.,
#' `"PatientID"`). If provided, the column is renamed to 'Subject' and used
#' for repeated-measures analyses downstream. If NULL (default), all samples
#' are treated as independent, e.g. for cell culture experiments.
#' experiments.
#'
#' @import SummarizedExperiment
#' @import magrittr
#'
#' @returns A SummarizedExperiment object containing the input data and sample
#' annotations.
#' @export
#' @examples
#' data <- data.frame(Feature = c("Gene1", "Gene2"),
#'                    Sample1 = c(1, 2), Sample2 = c(3, 4))
#' sample_ann <- data.frame(
#'     Sample = c("Sample1", "Sample2"),
#'     Group = c("A", "A"),
#'     Time = c(0, 1),
#'     Replicate = c(1, 1),
#'     Batch = c(1, 1)
#' )
#' se_obj <- create_input(data, sample_ann)
create_input <- function(data, sample_ann, subject_col = NULL) {
    # check data format
    if (!is.data.frame(data)) {
        stop("Input data must be a data frame.")
    }
    if (colnames(data)[1] != "Feature") {
        stop("The first column of data must be named 'Feature' ",
        "and contain feature identifiers.")
    }

    # check sample annotation
    required_cols <- c("Sample", "Group", "Time")
    if (!all(required_cols %in% names(sample_ann))) {
        stop(
            "Sample annotation is missing one or more required columns: ",
            paste(required_cols, collapse = ", ")
        )
    }

    # ---- Subject column handling ----
    if (!is.null(subject_col)) {
        if (!subject_col %in% colnames(sample_ann)) {
            stop("subject_col = '", subject_col, 
                "' not found in sample_ann. Available columns: ", 
                paste(colnames(sample_ann), collapse = ", "))
        }
        colnames(sample_ann)[colnames(sample_ann) == subject_col] <- "Subject"

        # Report Subject statistics
        n_subjects <- dplyr::n_distinct(sample_ann$Subject)
        subj_times <- sample_ann %>%
            dplyr::group_by(Subject, Group) %>%
            dplyr::summarise(n_tp = dplyr::n_distinct(Time), .groups = "drop")
        n_groups_per_subj <- sample_ann %>%
            dplyr::group_by(Subject) %>%
            dplyr::summarise(n_grp = dplyr::n_distinct(Group), .groups = "drop")

        message(sprintf(
            "Subject column '%s' recognized:", subject_col
        ))
        message(sprintf(
            "  %d unique subjects", n_subjects
        ))
        message(sprintf(
            "  Median %d time points per subject (range: %d-%d)",
            stats::median(subj_times$n_tp),
            min(subj_times$n_tp), max(subj_times$n_tp)
        ))

        cross_subj <- n_groups_per_subj$Subject[n_groups_per_subj$n_grp > 1]
        if (length(cross_subj) > 0) {
            stop(
            length(cross_subj), " subject(s) appear in multiple groups: ",
            paste(cross_subj, collapse = ", "),
            " Downstream functions assume subjects are unique to one group."
            )
        }

        if (all(subj_times$n_tp == 1)) {
            warning(
                "Subject column '", subject_col, "' found but every value ",
                "appears at only one time point, ", 
                "no repeated measures detected. ",
                "Treating samples as independent."
            )
            sample_ann$Subject <- NULL
        }
    } else {
        message(
            "No Subject column specified. Samples treated as independent. ",
            "For repeated-measures designs, set subject_col to the column ",
            "identifying biological subjects."
        )
    }

    # ---- Replicate handling ----
    if (!("Replicate" %in% names(sample_ann))) {
        message("No 'Replicate' column provided. Auto-generating replicate IDs.")
        if ("Subject" %in% names(sample_ann)) {
            sample_ann <- sample_ann %>%
                dplyr::group_by(Subject, Group, Time) %>%
                dplyr::arrange(Sample, .by_group = TRUE) %>%
                dplyr::mutate(Replicate = dplyr::row_number()) %>%
                dplyr::ungroup() %>%
                as.data.frame()
            comb_n <- table(paste(sample_ann$Subject, 
                sample_ann$Group, sample_ann$Time))
            message("  Replicate IDs assigned within each Group, Subject, ",
                "and Time. ", sum(comb_n > 1), 
                " combinations with >1 replicate."
            )
        } else {
            sample_ann <- sample_ann %>%
                dplyr::group_by(Group, Time) %>%
                dplyr::arrange(Sample, .by_group = TRUE) %>%
                dplyr::mutate(Replicate = dplyr::row_number()) %>%
                dplyr::ungroup() %>%
                as.data.frame()
            comb_n <- table(paste(sample_ann$Group, sample_ann$Time))
            message(
                "  Replicate IDs assigned within each Group and Time. ",
                sum(comb_n > 1), " combinations with >1 replicate."
            )
        }
    }

    if (!("Batch" %in% names(sample_ann))) {
        sample_ann$Batch <- 1
        sample_ann$Batch <- factor(sample_ann$Batch,
            levels = sort(as.numeric(unique(sample_ann$Batch)))
        )
    }
    if (!is.numeric(sample_ann$Time)) {
        sample_ann$Time <- as.numeric(as.character(sample_ann$Time))
        if (any(is.na(sample_ann$Time))) {
            stop("The 'Time' column must be numeric.")
        }
    }

    # Validate replicate uniqueness
    if ("Subject" %in% names(sample_ann)) {
        rep_check <- sample_ann %>%
            dplyr::group_by(Subject, Group, Time) %>%
            dplyr::summarise(
                rep_unique = (length(unique(Replicate)) == length(Replicate)),
                .groups = "drop"
            )
        if (any(rep_check$rep_unique != TRUE)) {
            stop("Duplicate 'Replicate' IDs found for the same ",
            "Group, Subject, and Time combination.")
        }
    } else {
        rep_check <- sample_ann %>%
            dplyr::group_by(Group, Time) %>%
            dplyr::summarise(
                rep_unique = (length(unique(Replicate)) == length(Replicate)),
                .groups = "drop"
            )
        if (any(rep_check$rep_unique != TRUE)) {
            stop("Duplicate 'Replicate' IDs found for the same ",
            "Group and Time combination.")
        }
    }

    # check if all samples have annotation
    if (!all(sample_ann$Sample %in% colnames(data)[-1])) {
        stop("Some samples in sample annotation are not present in the data.")
    }
    if (!all(colnames(data)[-1] %in% sample_ann$Sample)) {
        stop("Some samples in the data are not present in the sample ",
        "annotation.")
    }

    if (is.null(levels(sample_ann$Group))) {
        sample_ann$Group <- factor(sample_ann$Group,
            levels = sort(unique(sample_ann$Group)))
        message("Converting 'Group' column to factor. Default order ",
        "is alphabetical.")
    }
    if (is.null(levels(sample_ann$Replicate))) {
        sample_ann$Replicate <- factor(sample_ann$Replicate,
            levels = sort(as.numeric(unique(sample_ann$Replicate)))
        )
        message("Converting 'Replicate' column to factor. Default order ",
        "is numerical.")
    }
    if (is.null(levels(sample_ann$Batch))) {
        sample_ann$Batch <- factor(sample_ann$Batch,
            levels = sort(as.numeric(unique(sample_ann$Batch)))
        )
        message("Converting 'Batch' column to factor. Default order ",
        "is numerical.")
    }

    # Arrange in order: Group, Subject (if present), Time, Replicate, Batch
    if ("Subject" %in% names(sample_ann)) {
        sample_ann <- sample_ann %>%
            dplyr::arrange(Group, Subject, Time, Replicate, Batch) %>%
            droplevels()
    } else {
        sample_ann <- sample_ann %>%
            dplyr::arrange(Group, Time, Replicate, Batch) %>%
            droplevels()
    }
    sample_order <- sample_ann$Sample

    row.names(data) <- data$Feature
    data <- data[, c("Feature", sample_order)]

    stopifnot(identical(colnames(data)[-1], sample_ann$Sample))

    se_obj <- SummarizedExperiment(
        assays = data[, -1],
        colData = sample_ann
    )

    return(se_obj)
}
