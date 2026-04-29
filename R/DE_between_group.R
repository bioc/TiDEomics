#' DE between groups
#' @description Differential expression analysis between groups at each time
#' point by limma. The function compares each pair of groups at each time
#' point, and returns a nested list of DE analysis results for each pair of
#' groups and each time point including all features, as well as a list of
#' filtered DE results for each pair of groups based on the specified
#' thresholds including only significant features. The function can also plot
#' the number of DE features between groups over time.
#'
#' @param se_obj A SummarizedExperiment object
#' @param group (Optional) A character vector specifying which group to be
#' compared to. If NULL, all groups in the 'Group' column will be compared to.
#' (default is NULL)
#' @param filter (Optional) Minimum number of replicates required in both
#' conditions for a feature to be tested. If NULL, the minimum number of
#' replicates across all groups and time points will be used. (default is NULL)
#' @param assay Assay index to use, where 1 is the original data and 2 is
#' normalised to time 0 (if available) (default is 1)
#' @param adjP_thres (Optional) Threshold for adjusted p-value to consider a
#' feature as differentially expressed (default is 0.05)
#' @param logFC_thres (Optional) Threshold for log2 fold change to consider a
#' feature as differentially expressed (default is 1)
#' @param plot (Optional) Whether to plot the number of DE features between
#' groups over time (default is TRUE)
#' @param fontsize (Optional) Font size for the plot (default is 8)
#'
#' @import SummarizedExperiment
#' @import magrittr
#' @import ggplot2
#' @importFrom dplyr filter select bind_rows anti_join group_by summarise
#' @importFrom dplyr distinct
#'
#' @returns A list containing two elements: 'all_list' is a nested list of DE
#' results for each pair of groups and each time point including all features;
#' 'de_list' is a list of filtered DE results for each pair of groups based on
#' the specified thresholds including only significant features.
#' @export
#' @examples
#' data("example")
#' example_obj <- normalise_to_start(example_obj)
#'
#' DE_between_group_out <- DE_between_group(example_obj, assay = 2)
DE_between_group <- function(
    se_obj, group = NULL, filter = NULL, assay = c(1, 2),
    adjP_thres = 0.05, logFC_thres = 1, plot = TRUE, fontsize = 8
) {
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

    if (length(assay) > 1 | is.null(assay)) {
        assay <- 1
        message("Using assay 1 (original data) for DE analysis.")
    } else if (!(assay %in% seq_along(assays(se_obj)))) {
        assay_avail <- paste(seq_along(assays(se_obj)), collapse = ", ")
        stop(sprintf("Invalid assay index specified. Please choose from: %s",
            assay_avail))
    }

    # compare all time point by all time point
    outlist_limma <- list()

    groups <- levels(se_obj$Group)

    for (i in group) { # all groups if unselected

        for (j in groups) {
            if (i == j) next

            time_series_1 <-
                sort(unique(colData(se_obj[, se_obj$Group == i])$Time)) %>%
                as.character()
            time_series_2 <-
                sort(unique(colData(se_obj[, se_obj$Group == j])$Time)) %>%
                as.character()
            time_series <- intersect(time_series_1, time_series_2)
            # only compare time points present in both groups

            cond1 <- i
            cond2 <- j
            label_comparison <- paste0(cond2, "-", cond1)
            outlist_limma[[label_comparison]] <- list()

            for (t in time_series) {
                outlist_limma[[label_comparison]][[t]] <- list()

                input <- se_obj[, se_obj$Time == t]

                # comparisons between conditions
                d_cond1 <- input[, input$Group == cond1]
                d_cond2 <- input[, input$Group == cond2]

                # require at least 'filter' values in both conditions
                cond1_non_na <- rowSums(!is.na(assays(d_cond1)[[assay]])) %>%
                    as.data.frame()
                cond2_non_na <- rowSums(!is.na(assays(d_cond2)[[assay]])) %>%
                    as.data.frame()

                if (filter > dim(assays(d_cond1)[[assay]])[2] |
                    filter > dim(assays(d_cond2)[[assay]])[2]) {
                    stop("Filter value is larger than the number of ",
                    "replicates in one or both conditions.")
                }
                rows_selected <- (cond1_non_na$. >= filter &
                    cond2_non_na$. >= filter)

                n <- length(rows_selected)
                n_keep <- sum(rows_selected)
                pct <- round(n_keep / n, 3) * 100
                message("Comparing group ", cond2, " to ", cond1, " at Time ",
                    t, ": keeping ", n_keep, " of ", n, " features (",
                    sprintf("%.1f", pct), "%)")

                d_cond1_filter <- d_cond1[rows_selected, ]
                d_cond2_filter <- d_cond2[rows_selected, ]
                tb_compare <- cbind(d_cond1_filter, d_cond2_filter)

                # limma
                design <- stats::model.matrix(~ 0 +
                    as.character(tb_compare$Group))
                # numbers are non-valid names
                colnames(design) <- colnames(design) %>%
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
                fit3 <- limma::eBayes(fit2)
                modtest <- limma::topTable(fit3, number = Inf, sort.by = "none")
                limma_result <- merge(assays(tb_compare)[[assay]],
                    as.data.frame(modtest[, -c(2, 3, 6)]),
                    by = "row.names", all = TRUE
                )
                rownames(limma_result) <- limma_result$Row.names
                limma_result <- limma_result %>% dplyr::select(-Row.names)

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
                    by = "row.names", all = TRUE) %>%
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
                tb_de <- tb %>%
                    dplyr::select(Feature, Comparison, Time, Cond1, Cond2,
                        logFC, adj.P.Val) %>%
                    filter(adj.P.Val < adjP_thres) %>%
                    filter(logFC > logFC_thres | logFC < -logFC_thres)
                de_list_limma[[label_comparison]][[tb_name]] <- tb_de
            }

            de_list_limma[[label_comparison]] <-
                de_list_limma[[label_comparison]] %>%
                bind_rows()
        }
    }

    if (plot) {
        de_num_limma <- de_list_limma %>%
            lapply(function(x) {
                x %>%
                    group_by(Time) %>%
                    summarise(N = n()) %>%
                    mutate(
                        Cond1 = unique(x$Cond1),
                        Cond2 = unique(x$Cond2),
                        Time = as.numeric(as.character(Time))
                    )
            }) %>%
            do.call(rbind, .) %>%
            as.data.frame()

        # when certain comparisons are valid but have no DE features
        # the time point is not included in the plot, so add it with n = 0
        all_comparisons_t <- outlist_limma %>%
            lapply(function(x) { # list of comparisons
                lapply(x, function(tb) { # list of time points
                    data.frame(
                        Cond1 = unique(tb$Cond1),
                        Cond2 = unique(tb$Cond2),
                        Time = as.numeric(as.character(unique(tb$Time)))
                    )
                }) %>%
                    bind_rows()
            }) %>%
            bind_rows() %>%
            distinct()
        if (dim(all_comparisons_t)[1] > dim(de_num_limma)[1]) {
            missing_comparisons <- anti_join(all_comparisons_t, de_num_limma,
                by = c("Cond1", "Cond2", "Time")) # rows in x and not in y
            missing_comparisons$N <- 0
            de_num_limma <- bind_rows(de_num_limma, missing_comparisons)
        }

        for (i in group) {
            (de_num_limma %>%
                filter(Cond1 == i) %>%
                mutate(Group = Cond2) %>%
                ggplot(aes(x = Time, y = N, group = Group, color = Group)) +
                geom_point() +
                geom_line() +
                labs(
                    title = paste0("DE feature number (vs. ", i, ")"),
                    y = "DE feature number"
                ) +
                scale_color_manual(values = get_custom_palette(groups)) +
                theme_custom(base_size = fontsize)
            ) %>% print()
        }
    }

    return(list(all_list = outlist_limma, de_list = de_list_limma))
}
