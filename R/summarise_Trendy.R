#' Summarise Trendy results
#' @description Summarise Trendy results into data frame
#'
#' @param res_list A list of Trendy analysis results, output of `run_Trendy()`
#' @param ... Additional arguments to be passed to the `Trendy::topTrendy()`
#'
#' @returns A data frame of summary results of Trendy, including breakpoints,
#' segment slopes and p-values for each fitted feature in each group. A
#' "Pattern" column is added to show the combination of segment trends (up,
#' down, stable) for each feature.
#' @export
#' @examples
#' data("example_res_list")
#' trendy_summary <- summarise_Trendy(example_res_list)
summarise_Trendy <- function(res_list, ...) {
    .check_list(res_list, "res_list", "run_Trendy")
    res_summary <- list()
    for (i in names(res_list)) {
        if (is.null(res_list[[i]])) {
            message("Skipping group '", i,
                "': Trendy analysis was not performed (insufficient data).")
            next
        }
        trendy.summary <- .summarise_Trendy_one_group(res_list[[i]], ...)
        trendy.summary$Group <- i
        trendy.summary <- trendy.summary |>
            dplyr::select(Group, dplyr::everything())
        res_summary[[i]] <- trendy.summary
    }
    if (length(res_summary) == 0) {
        message("No groups with valid Trendy results.")
        return(NULL)
    }
    res_summary_df <- do.call(rbind, res_summary) |> as.data.frame()

    # convert segment trends from 1, -1, 0 to up, down, stable
    df <- res_summary_df

    seg_cols <- grep("^Segment[0-9]+\\.Trend$", colnames(df), value = TRUE)
    df[seg_cols] <- df[seg_cols] |> lapply(function(col) {
        c("1" = "up", "0" = "stable", "-1" = "down")[as.character(col)]
    })

    # all combinations of segment trends
    patterns <- apply(df[seg_cols], 1, function(row) {
        prs <- which(!is.na(row))
        paste(row[prs], collapse = "_")
    })
    df_patterns <- df |>
        dplyr::mutate(Pattern = patterns)

    return(df_patterns)
}

#' Summarise Trendy results (one group)
#' @description Summarise Trendy results into data frame (one group)
#'
#' @param res Result of the Trendy analysis, output of `run_Trendy()`
#' @param ... Additional arguments to be passed to the `Trendy::topTrendy`
#'
#'
#' @returns A data frame of summary results of Trendy, including breakpoints,
#' segment slopes and p-values for each fitted feature.
#' @keywords internal
.summarise_Trendy_one_group <- function(res, ...) {
    res.top <- Trendy::topTrendy(res, ...)
    trendy.summary <- Trendy::formatResults(res.top)

    # bug in colname when maxK = 1
    if ("topTrendyData.Breakpoints.featureNames..." %in%
        colnames(trendy.summary)) {
        colnames(trendy.summary)[
            colnames(trendy.summary) ==
                "topTrendyData.Breakpoints.featureNames..."
        ] <- "Breakpoint"
    }

    return(trendy.summary)
}
