#' GO enrichment with ranked gene list
#'
#' @description Gene ontology enrichment analysis of a ranked gene list
#' with clusterProfiler. Genes can be ranked by variance decomposition results.
#'
#' @param rank_table A data frame with at least two columns: 'Feature' for gene
#' names and one or more columns of variables for ranking the genes,
#' e.g. output of `decomp_variance()` containing variance decomposition results
#' @param gene_rank_by Variable in `rank_table` to rank the genes by,
#' e.g. "Time", "Group" in the output of `decomp_variance()`
#' @param OrgDb Organism database, e.g. org.Hs.eg.db, org.Mm.eg.db
#' @param keyType (Optional) Available options are
#' `AnnotationDbi::keytypes(OrgDb)` (default is "SYMBOL")
#' @param go_rank_by (Optional) Variable in the GO enrichment result to
#' rank the GO terms by (default is "p.adjust", other options include
#' "pvalue", "qvalue", "NES", "setSize", "enrichmentScore", etc.)
#' @param category (Optional) GO category to analyze (default is all three of
#' BP, MF, CC)
#' @param pvalueCutoff  (Optional) Parameter of `clusterProfiler::gseGO()`
#' (default is 0.05)
#' @param pAdjustMethod (Optional) Parameter of `clusterProfiler::gseGO()`
#' (default is "BH")
#' @param ... additional arguments passed to `clusterProfiler::gseGO()`
#'
#' @returns A list of `gseaResult` objects containing the GSEA results
#' @export
#'
#' @examples
#' library(org.Mm.eg.db)
#' data(example_obj)
#' example_obj <- normalise_to_start(example_obj)
#'
#' var_decomp <- decomp_variance(example_obj, assay = 1)
#' example_go_rank <- enrichGO_rank(var_decomp, gene_rank_by = "Time",
#'     OrgDb = org.Mm.eg.db, keyType = "SYMBOL", category = "BP")
enrichGO_rank <- function(
    rank_table,
    gene_rank_by,
    OrgDb,
    keyType = "SYMBOL",
    go_rank_by = "p.adjust",
    category = NULL,
    pvalueCutoff = 0.05,
    pAdjustMethod = "BH",
    ...
) {
    .check_df(rank_table, "rank_table")
    .check_character(gene_rank_by, "gene_rank_by")
    .check_character(keyType, "keyType")
    .check_character(go_rank_by, "go_rank_by")
    pAdjustMethod <- match.arg(pAdjustMethod,
        c("holm", "hochberg", "hommel", "bonferroni",
            "BH", "BY", "fdr", "none"))
    .check_pval(pvalueCutoff, "pvalueCutoff")
    if (!inherits(OrgDb, "OrgDb")) {
        stop("'OrgDb' must be an OrgDb object, e.g. org.Hs.eg.db.")
    }

    if (is.null(category)) {
        category <- c("BP", "MF", "CC")
        message("GO category not specified. Using all three: BP, MF, CC.")
    } else {
        category <- match.arg(category, c("BP", "MF", "CC"),
            several.ok = TRUE)
    }

    if (is.null(gene_rank_by) || length(gene_rank_by) == 0 ||
        length(gene_rank_by) > 1) {
        stop("Please specify one variable for ranking the genes.")
    }

    if (!(gene_rank_by %in% colnames(rank_table))) {
        stop("Gene ranking variable not found in the input: ",
            paste(setdiff(gene_rank_by, colnames(rank_table)), collapse = ", "))
    }

    # gene set enrichment with ranked gene list
    rank_list <- rank_table |>
        dplyr::arrange(dplyr::desc(.data[[gene_rank_by]])) |>
        dplyr::pull(.data[[gene_rank_by]], name = Feature)

    if (any(is.na(rank_list)) || any(is.nan(rank_list)) ||
            any(is.infinite(rank_list))) {
        message("NA / NaN / Inf values found in the gene ranking variable. ",
            "Those genes will be removed.")
        rank_list <- rank_list[!is.na(rank_list)]
        rank_list <- rank_list[!is.infinite(rank_list)]
    }

    gse_list <- list()
    for (ont in category) {
        gse_rank <- clusterProfiler::gseGO(
            geneList = rank_list,
            OrgDb = OrgDb,
            keyType = keyType,
            ont = ont,
            pvalueCutoff = pvalueCutoff,
            pAdjustMethod = pAdjustMethod,
            ...
        )
        if (is.null(gse_rank) || is.null(gse_rank@result) ||
            nrow(gse_rank@result) == 0) {
            message("No significant GO terms found for category: ", ont)
            next
        }

        message("Removing NA ID gene sets for ", ont, ".")

        gse_rank@result <- gse_rank@result |>
            dplyr::filter(!is.na(ID))

        if (go_rank_by %in% colnames(gse_rank@result)) {
            gse_rank@result <- gse_rank@result |>
                dplyr::arrange(.data[[go_rank_by]])
        } else {
            message("GO term ranking variable not found in the result: ",
                paste(setdiff(go_rank_by, colnames(gse_rank@result)),
                collapse = ", "), ". Returning results in original order.")
        }

        gse_list[[ont]] <- gse_rank
    }

    if (length(gse_list) == 0) {
        message("No significant GO terms found in any category.")
        return(NULL)
    }

    return(gse_list)
}
