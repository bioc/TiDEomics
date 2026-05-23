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
#' columns and will be set to 1 for all samples if not provided. 'Time',
#' 'Replicate' and 'Batch' should be numeric.
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
create_input <- function(data, sample_ann) {
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
    if (!("Replicate" %in% names(sample_ann))) {
        sample_ann$Replicate <- 1
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
    # no duplicate replicates ID for each group and time point
    rep_check <- sample_ann %>%
        dplyr::group_by(Group, Time) %>%
        dplyr::summarise(rep_unique =
            (length(unique(Replicate)) == length(Replicate))) %>%
        dplyr::ungroup()
    if (any(rep_check$rep_unique != TRUE)) {
        stop("Duplicate 'Replicate' IDs found for the same 'Group' ",
        "and 'Time' combination.")
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

    sample_ann <- sample_ann %>%
        dplyr::arrange(Group, Time, Replicate, Batch) %>%
        droplevels()
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
