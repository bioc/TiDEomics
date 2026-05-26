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
#' @returns A `gseaResult`object containing the GO enrichment results
#' @export
#'
#' @examples
#' library(org.Mm.eg.db)
#' data("example")
#' example_obj <- normalise_to_start(example_obj)
#'
#' var_decomp <- decomp_variance(example_obj, assay = 1)
#' example_go_rank <- enrichGO_rank(var_decomp, gene_rank_by = "Time",
#'     OrgDb = org.Mm.eg.db, keyType = "SYMBOL", category = "BP")
#' enrichplot::gseaplot2(example_go_rank, geneSetID = 1:2)
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
    if (is.null(category)) {
        category <- c("BP", "MF", "CC")
        message("GO category not specified. Using all three: BP, MF, CC.")
    } else if (!all(category %in% c("BP", "MF", "CC"))) {
        stop("Invalid GO category. Please choose from 'BP', 'MF', 'CC'.")
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
    rank_list <- rank_table %>%
        dplyr::arrange(dplyr::desc(.data[[gene_rank_by]])) %>%
        dplyr::pull(.data[[gene_rank_by]], name = Feature)

    if (any(is.na(rank_list)) || any(is.nan(rank_list)) ||
            any(is.infinite(rank_list))) {
        message("NA / NaN / Inf values found in the gene ranking variable. ",
            "Those genes will be removed.")
        rank_list <- rank_list[!is.na(rank_list)]
        rank_list <- rank_list[!is.infinite(rank_list)]
    }

    gse_rank <- clusterProfiler::gseGO(
        geneList = rank_list,
        OrgDb = OrgDb,
        keyType = keyType,
        ont = category,
        pvalueCutoff = pvalueCutoff,
        pAdjustMethod = pAdjustMethod,
        ...
    )

    message("Removing NA ID gene sets.")

    gse_rank@result <- gse_rank@result %>%
        dplyr::filter(!is.na(ID))

    if (is.null(gse_rank) || is.null(gse_rank@result) ||
        dim(gse_rank@result)[1] == 0) {
        message("No significant GO terms found.")
        return(NULL)
    }

    if (go_rank_by %in% colnames(gse_rank@result)) {
        gse_rank@result <- gse_rank@result %>%
            dplyr::arrange(.data[[go_rank_by]])
    } else {
        message("GO term ranking variable not found in the result: ",
            paste(setdiff(go_rank_by, colnames(gse_rank@result)),
            collapse = ", "), ". Returning results in original order.")
    }

    return(gse_rank)
}
