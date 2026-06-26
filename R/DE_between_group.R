#' DE between groups
#' @description Differential expression analysis between groups at each time
#' point by limma. The function compares each pair of groups at each time
#' point, and returns a nested list of DE analysis results for each pair of
#' groups and each time point including all features, as well as a list of
#' filtered DE results for each pair of groups based on the specified
#' thresholds including only significant features. Use
#' `plot_DE_between_group()` to visualise the number of DE features between
#' groups over time.
#'
#' @details
#' Only time points present in both groups are compared. No multiple-testing
#' correction is applied across the pairwise group-by-time comparisons; each
#' comparison is independent. Apply your own correction (e.g.,
#' `stats::p.adjust()`) across combined results if needed.
#'
#' @param se_obj A SummarizedExperiment object
#' @param group (Optional) A character vector specifying which group to be
#' compared to. If NULL, all groups in the 'Group' column will be compared to.
#' (default is NULL)
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
#'   pair and time point, all features); `de_list` (significant features
#'   only); `fit_list` (nested list of limma `MArrayLM` fit objects,
#'   for use with `limma::plotSA()`); `ref_groups` (group to be compared to);
#'   `all_groups` (all groups). The output can be passed to
#'   `plot_DE_between_group()` for visualisation.
#' @export
#' @examples
#' data("example")
#' example_obj <- normalise_to_start(example_obj)
#'
#' DE_between_group_out <- DE_between_group(example_obj, assay = 2)
DE_between_group <- function(
    se_obj, group = NULL, filter = NULL, assay,
    adjP_thres = 0.05, logFC_thres = 1, trend = FALSE
) {
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
        stop("At least one of the specified groups is not found in ",
        "the 'Group' column of the input object.")
    }

    if (is.null(filter)) {
        # minimum number of replicates in any group at any time point
        filter <- min(table(se_obj$Group, se_obj$Time))
        message("Non-NA replicate number filter not specified. ",
        "Using minimum number of replicates across all groups ",
        "and time points: ", filter)
    }

    assay <- .match_assay(assay, se_obj)
    .check_limma_input(assay(se_obj, assay), assay)

    # compare all time point by all time point
    outlist_limma <- list()
    fitlist_limma <- list()

    groups <- unique(as.character(se_obj$Group))

    for (i in group) { # all groups if unselected

        for (j in groups) {
            if (i == j) next

            time_series_1 <-
                sort(unique(colData(se_obj[, se_obj$Group == i])$Time)) |>
                as.character()
            time_series_2 <-
                sort(unique(colData(se_obj[, se_obj$Group == j])$Time)) |>
                as.character()
            time_series <- intersect(time_series_1, time_series_2)
            # only compare time points present in both groups

            cond1 <- i
            cond2 <- j

            dropped1 <- setdiff(time_series_1, time_series_2)
            dropped2 <- setdiff(time_series_2, time_series_1)
            if (length(dropped1) > 0 || length(dropped2) > 0) {
                message("Comparing ", cond2, " vs ", cond1,
                    ": time points only in ", cond1, ": ",
                    if (length(dropped1) > 0)
                        paste(dropped1, collapse = ", ") else "none",
                    "; only in ", cond2, ": ",
                    if (length(dropped2) > 0)
                        paste(dropped2, collapse = ", ") else "none")
            }
            label_comparison <- paste0(cond2, "-", cond1)
            outlist_limma[[label_comparison]] <- list()

            for (t in time_series) {
                outlist_limma[[label_comparison]][[t]] <- list()

                input <- se_obj[, se_obj$Time == t]

                # comparisons between conditions
                d_cond1 <- input[, input$Group == cond1]
                d_cond2 <- input[, input$Group == cond2]

                # require at least 'filter' values in both conditions
                cond1_non_na <- rowSums(!is.na(assays(d_cond1)[[assay]]))
                cond2_non_na <- rowSums(!is.na(assays(d_cond2)[[assay]]))

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
                    message("Comparing group ", cond2, " to ", cond1,
                        " at Time ", t,
                        ": no features pass the filter. Skipping.")
                    next
                }
                pct <- round(n_keep / n, 3) * 100
                message("Comparing group ", cond2, " to ", cond1,
                    " at Time ", t, ": keeping ", n_keep, " of ", n,
                    " features (", sprintf("%.1f", pct), "%)")

                d_cond1_filter <- d_cond1[rows_selected, ]
                d_cond2_filter <- d_cond2[rows_selected, ]
                tb_compare <- cbind(d_cond1_filter, d_cond2_filter)

                # limma
                design <- stats::model.matrix(~ 0 +
                    as.character(tb_compare$Group))
                # numbers are non-valid names
                colnames(design) <- colnames(design) |>
                    gsub(pattern = "as.character(tb_compare$Group)",
                        fixed = TRUE, replacement = "Group")
                rownames(design) <- colnames(tb_compare)
                fit1 <- limma::lmFit(assays(tb_compare)[[assay]],
                    design = design,
                    maxit = 2000
                )
                contrast <- limma::makeContrasts(
                    contrasts = paste0("Group", cond2, "-Group", cond1),
                    levels = colnames(design)
                )
                fit2 <- limma::contrasts.fit(fit1, contrasts = contrast)
                fit3 <- limma::eBayes(fit2, trend = trend)
                fitlist_limma[[label_comparison]][[t]] <- fit3
                modtest <- limma::topTable(fit3, number = Inf, sort.by = "none")
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
                gene_df$Time <- t
                gene_df$Cond1 <- factor(cond1, levels = groups)
                gene_df$Cond2 <- factor(cond2, levels = groups)

                # combine with gene df and sort
                out_limma <- merge(gene_df, limma_result,
                    by = "row.names", all = TRUE) |>
                    dplyr::select(-Row.names)

                outlist_limma[[label_comparison]][[t]] <- out_limma
            }
        }
    }

    # list of DE tables

    # default criteria: p.adj < 0.05, |log2FC| > 1
    de_list_limma <- list()

    for (i in group) {

        for (j in groups) {
            if (i == j) next

            label_comparison <- paste0(j, "-", i)
            de_list_limma[[label_comparison]] <- list()

            for (tb_name in names(outlist_limma[[label_comparison]])) {
                tb <- outlist_limma[[label_comparison]][[tb_name]]
                tb_de <- tb |>
                    dplyr::select(Feature, Comparison, Time, Cond1, Cond2,
                        logFC, adj.P.Val) |>
                    dplyr::filter(adj.P.Val < adjP_thres) |>
                    dplyr::filter(logFC > logFC_thres | logFC < -logFC_thres)
                de_list_limma[[label_comparison]][[tb_name]] <- tb_de
            }

            de_list_limma[[label_comparison]] <-
                de_list_limma[[label_comparison]] |>
                dplyr::bind_rows()
        }
    }

    return(list(all_list = outlist_limma, de_list = de_list_limma,
                fit_list = fitlist_limma, ref_groups = group, all_groups = groups))
}
