#' Flatten nested differential expression results
#'
#' @description
#' Flattens the nested DE result structure returned by
#' [DE_between_group()] or [DE_between_time()] into a named list of
#' data frames, one entry per contrast-time combination.
#' This is useful for exporting results to other tools such as
#' `DeeDeeExperiment`.
#'
#' @param de_list A list of DE results. Accepts individual output of
#' [DE_between_group()] or [DE_between_time()] with an `all_list`
#' element; or a named list (e.g. `list(between_group = ...,
#' between_time = ...)`), each containing an `all_list` element.
#'
#' @returns A named list of data frames with columns renamed for
#'   `DeeDeeExperiment` compatibility (`log2FoldChange`, `pvalue`, `padj`).
#'   Returns an empty list if `de_list` is empty.
#'
#' @export
#'
#' @examples
#' data(tutorial_data)
#' data(tutorial_sample_info)
#' tide <- prepare_tide(tutorial_data, tutorial_sample_info,
#'     keep = "threshold", residual_threshold = 100)
#' tide$DE <- DE_between_group(tide$se, assay = "norm",
#'     filter = 1, trend = TRUE)
#' de_flat <- flatten_DE(tide$DE)
#' str(de_flat, max.level = 1)
flatten_DE <- function(de_list) {
    .check_list(de_list, "de_list")
    if (length(de_list) == 0) return(list())

    .rename_de_cols <- function(df) {
        colnames(df)[colnames(df) == "logFC"]     <- "log2FoldChange"
        colnames(df)[colnames(df) == "P.Value"]   <- "pvalue"
        colnames(df)[colnames(df) == "adj.P.Val"] <- "padj"
        if ("Feature" %in% colnames(df))
            rownames(df) <- df$Feature
        df
    }

    .extract_contrast <- function(contrast_data, base_name) {
        if (is.data.frame(contrast_data))
            return(stats::setNames(
                list(.rename_de_cols(contrast_data)), base_name))
        if (is.list(contrast_data)) {
            out <- list()
            for (time in names(contrast_data)) {
                df <- contrast_data[[time]]
                if (is.data.frame(df)) {
                    t_prefix <- if (grepl("^t", time)) "_" else "_T"
                    out[[paste0(base_name, t_prefix, time)]] <-
                        .rename_de_cols(df)
                }
            }
            return(out)
        }
        list()
    }

    # Case 1: DE_between_group / DE_between_time direct output
    if (!is.null(de_list$all_list) && is.list(de_list$all_list)) {
        out <- list()
        for (contrast in names(de_list$all_list))
            out <- c(out, .extract_contrast(
                de_list$all_list[[contrast]], contrast))
        return(out)
    }

    # Case 2: Nested format, prepare_tide element style
    out <- list()
    for (nm in names(de_list)) {
        chunk <- de_list[[nm]]
        if (!is.list(chunk)) next
        if (!is.null(chunk$all_list) && is.list(chunk$all_list)) {
            for (contrast in names(chunk$all_list))
                out <- c(out, .extract_contrast(chunk$all_list[[contrast]],
                    paste0(nm, "_", contrast)))
        }
    }
    out
}
