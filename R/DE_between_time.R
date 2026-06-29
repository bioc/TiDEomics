#' DE between time points
#' @description Differential expression analysis between time points within
#' each group by limma
#'
#' @param se_obj A SummarizedExperiment object created by `create_input`
#' @param group (Optional) A character vector specifying which groups to
#' analyse. If NULL, all groups in the 'Group' column will be used. (default
#' is NULL)
#' @param filter (Optional) Minimum number of replicates required in both
#' conditions for a feature to be tested. If NULL, the minimum number of
#' replicates across all groups and time points will be used. (default is NULL)
#' @param assay Assay to use: `"orig"` for original data, `"norm"` for
#' normalised-to-start data. Numeric indices (1, 2) are also accepted.
#' No default, must be specified explicitly. The selected assay should
#' contain log-transformed, normalised values (e.g. log2-CPM for RNA-seq,
#' log2-intensity for proteomics). A warning is issued if the data appears
#' to be un-logged raw counts.
#' @param adjP_thres (Optional) Threshold for adjusted p-value to consider a
#' feature as differentially expressed (default is 0.05)
#' @param logFC_thres (Optional) Threshold for log2 fold change to consider a
#' feature as differentially expressed (default is 1)
#' @param trend (Optional) Logical, passed to `limma::eBayes()`.
#'   Set to `TRUE` for RNA-seq count-derived data to model the mean-variance
#'   trend. Leave as `FALSE` (default) for microarray, proteomics,
#'   metabolomics, or other log-intensity data where the mean-variance
#'   relationship is typically flat.
#'
#' @import SummarizedExperiment
#'
#' @returns A list with: `all_list` (nested list of DE results per group
#'   and time comparison, all features); `de_list` (nested list of
#'   significant features only); `fit_list` (nested list of limma
#'   `MArrayLM` fit objects, for use with `limma::plotSA()`);
#'   `time_series` (vector of all available time points). The output
#'   can be passed to `plot_DE_between_time()` for visualisation.
#' @export
#' @examples
#' data(example_obj)
#' example_obj <- normalise_to_start(example_obj)
#'
#' DE_between_time_out <- DE_between_time(example_obj, assay = 1)
DE_between_time <- function(se_obj, group = NULL, filter = NULL,
    assay, adjP_thres = 0.05, logFC_thres = 1,
    trend = FALSE) {
    .check_se(se_obj)
    .check_pval(adjP_thres, "adjP_thres")
    .check_nonneg(logFC_thres, "logFC_thres")
    .check_logical(trend, "trend")
    if (!is.null(filter)) {
        .check_positive_int(filter, "filter")
    }

    if (is.null(group)) {
        group <- unique(se_obj$Group)
    } else if (!all(group %in% unique(se_obj$Group))) {
        stop("At least one of the specified groups is not found in the ",
        "'Group' column of the input object.")
    }

    if (is.null(filter)) {
        # minimum number of replicates in any group at any time point
        filter <- min(table(se_obj$Group, se_obj$Time))
        message("Non-NA replicate number filter not specified. ",
        "Using minimum number of replicates across all groups and ",
        "time points: ", filter)
    }

    assay <- .match_assay(assay, se_obj)
    .check_limma_input(assay(se_obj, assay), assay)

    # compare all time point by all time point
    outlist_limma <- list()
    fitlist_limma <- list()

    time_series_all <- NULL

    for (i in group) {
        input <- se_obj[, se_obj$Group == i]
        time_series <- sort(unique(colData(input)$Time)) |> as.character()
        time_series_all <- unique(c(time_series_all, time_series))

        outlist_limma[[i]] <- list()
        fitlist_limma[[i]] <- list()

        for (cond1 in time_series) {
            for (cond2 in time_series) {
                if (as.numeric(cond1) < as.numeric(cond2)) {
                    # comparisons between conditions
                    label_comparison <- paste0("t", cond2, "-t", cond1)
                    d_cond1 <- input[, input$Time == as.numeric(cond1)]
                    d_cond2 <- input[, input$Time == as.numeric(cond2)]

                    # require at least 'filter' values in both conditions
                    cond1_non_na <-
                        rowSums(!is.na(assays(d_cond1)[[assay]]))
                    cond2_non_na <-
                        rowSums(!is.na(assays(d_cond2)[[assay]]))

                    if (filter > dim(assays(d_cond1)[[assay]])[2] ||
                        filter > dim(assays(d_cond2)[[assay]])[2]) {
                        stop("Filter value is larger than the number of ",
                        "replicates in one or both conditions.")
                    }
                    rows_selected <- (cond1_non_na >= filter &
                        cond2_non_na >= filter)

                    n <- length(rows_selected)
                    n_keep <- sum(rows_selected)
                    if (n == 0 || n_keep == 0) {
                        message("Comparing group ", i, " time ", cond2, " to ",
                            cond1, ": no features pass the filter. Skipping.")
                        next
                    }
                    pct <- round(n_keep / n, 3) * 100
                    message("Comparing group ", i, " time ", cond2, " to ",
                        cond1, ": keeping ", n_keep, " of ", n, " features (",
                        sprintf("%.1f", pct), "%)")

                    d_cond1_filter <- d_cond1[rows_selected, ]
                    d_cond2_filter <- d_cond2[rows_selected, ]
                    tb_compare <- cbind(d_cond1_filter, d_cond2_filter)

                    # limma
                    design <- stats::model.matrix(~ 0 +
                        as.character(tb_compare$Time))
                    # numbers are non-valid names
                    colnames(design) <- colnames(design) |>
                        gsub(pattern = "as.character(tb_compare$Time)",
                        fixed = TRUE, replacement = "t")
                    rownames(design) <- colnames(tb_compare)

                    # Paired analysis if Subject column present
                    has_subject <- "Subject" %in% colnames(colData(tb_compare))
                    if (has_subject) {
                        corfit <- limma::duplicateCorrelation(
                            assays(tb_compare)[[assay]],
                            design = design,
                            block = colData(tb_compare)$Subject
                        )
                        fit1 <- limma::lmFit(assays(tb_compare)[[assay]],
                            design = design,
                            block = colData(tb_compare)$Subject,
                            correlation = corfit$consensus.correlation,
                            maxit = 2000
                        )
                    } else {
                        fit1 <- limma::lmFit(assays(tb_compare)[[assay]],
                            design = design,
                            maxit = 2000
                        )
                    }
                    contrast <- limma::makeContrasts(
                        contrasts = paste0("t", cond2, "-t", cond1),
                        levels = colnames(design)
                    )
                    fit2 <- limma::contrasts.fit(fit1, contrasts = contrast)
                    fit3 <- limma::eBayes(fit2, trend = trend)
                    fitlist_limma[[i]][[label_comparison]] <- fit3
                    modtest <- limma::topTable(fit3, number = Inf,
                        sort.by = "none")
                    modtest <- modtest[, !colnames(modtest) %in%
                        c("AveExpr", "t", "B")]
                    limma_result <- merge(assays(tb_compare)[[assay]],
                        modtest,
                        by = "row.names", all = TRUE
                    )
                    rownames(limma_result) <- limma_result$Row.names
                    limma_result <- limma_result |> dplyr::select(-Row.names)

                    gene_df <- data.frame(
                        row.names = rownames(tb_compare),
                        Feature = rownames(tb_compare)
                    )

                    # setup columns
                    gene_df$Comparison <- label_comparison
                    gene_df$Group <- i
                    gene_df$Cond1 <- factor(cond1, levels = time_series)
                    gene_df$Cond2 <- factor(cond2, levels = time_series)

                    # combine with gene df and sort
                    out_limma <- merge(gene_df, limma_result,
                        by = "row.names", all = TRUE) |>
                        dplyr::select(-Row.names)

                    outlist_limma[[i]][[label_comparison]] <- out_limma
                }
            }
        }
    }

    # list of DE tables
    de_list_limma <- list()

    for (i in group) {
        de_list_limma[[i]] <- list()

        for (tb_name in names(outlist_limma[[i]])) { # example: t1-t0
            tb <- outlist_limma[[i]][[tb_name]]
            tb_de <- tb |>
                dplyr::select(Feature, Comparison, Group, Cond1, Cond2,
                    logFC, adj.P.Val) |>
                dplyr::filter(adj.P.Val < adjP_thres) |>
                dplyr::filter(logFC > logFC_thres | logFC < -logFC_thres)
            de_list_limma[[i]][[tb_name]] <- tb_de
        }

    }

    return(list(all_list = outlist_limma, de_list = de_list_limma,
                fit_list = fitlist_limma, time_series = time_series_all))
}
