#' GO enrichment with gene sets
#' @description Gene ontology enrichment analysis of multiple gene sets
#'   with clusterProfiler, allowing for using different or same background
#'   genes for each gene set
#'
#' @param gene_list A named list of gene vectors, or a data.frame
#'   with `Feature` and `Module` columns from `WGCNA_module()`.
#' @param keyType (Optional) Available options are
#'   `AnnotationDbi::keytypes(OrgDb)` (default is "SYMBOL")
#' @param OrgDb Organism database, e.g. org.Hs.eg.db, org.Mm.eg.db
#' @param universe Background genes for all input gene sets,
#' used if `universe_list` is not provided
#' @param universe_list Background genes for each input gene set,
#' a list of gene vectors with the same names as `gene_list`
#' @param pAdjustMethod (Optional) Parameter of `clusterProfiler::enrichGO()`
#'   (default is "BH")
#' @param pvalueCutoff (Optional) Parameter of `clusterProfiler::enrichGO()`
#'   (default is 0.05)
#' @param qvalueCutoff (Optional) Parameter of `clusterProfiler::enrichGO()`
#'   (default is 0.05)
#' @param category (Optional) GO category to analyze (default is all three of
#'   BP, MF, CC)
#' @param simplify (Optional) Whether to simplify the GO terms by removing
#'   redundant terms with `clusterProfiler::simplify()` function.
#'   (default is FALSE)
#' @param simplify_cutoff (Optional) Parameter of `clusterProfiler::simplify()`,
#'   cutoff for similarity when simplifying GO terms (default is 0.7)
#' @param simplify_by (Optional) Parameter of `clusterProfiler::simplify()`,
#'   method to choose representative term when simplifying GO terms
#'   (default is "p.adjust")
#' @param simplify_select_fun (Optional) Parameter of
#'   `clusterProfiler::simplify()`, function to select representative term when
#'   simplifying GO terms (default is `min`)
#' @param simplify_measure (Optional) Parameter of
#'   `clusterProfiler::simplify()`,
#'   method to calculate similarity when simplifying GO terms
#'   (default is "Wang")
#' @param ... additional arguments passed to `clusterProfiler::enrichGO()`
#'
#' @returns A nested list of GO enrichment results with sublists:
#' 'all' including all terms and 'simplified' including simplified terms
#' (if `simplify = TRUE`), both containing further sublists
#' for each GO category (BP, MF, CC).
#' @export
#' @examples
#' library(magrittr)
#' library(org.Mm.eg.db)
#' library(clusterProfiler)
#' data(example_net)
#' # select two modules for demonstration
#' example_module <- WGCNA_module(example_net) %>%
#'     dplyr::filter(Module %in% c("1", "2"))
#' # set cutoff to 1 to show all results for demonstration
#' example_go_list = enrichGO_list(example_module, OrgDb = org.Mm.eg.db,
#'     universe = example_module$Feature,
#'     pvalueCutoff = 1, qvalueCutoff = 1,
#'     category = "BP", simplify = FALSE)
#' # plot_GO(example_go_list$all, plot_dotplot = TRUE,
#' #     plot_emapplot = FALSE, plot_cnetplot = FALSE)
enrichGO_list <- function(gene_list, keyType = "SYMBOL",
    OrgDb,
    universe = NULL,
    universe_list = NULL,
    pAdjustMethod = "BH",
    pvalueCutoff = 0.05,
    qvalueCutoff = 0.05,
    category = NULL,
    simplify = FALSE,
    simplify_cutoff = 0.7,
    simplify_by = "p.adjust",
    simplify_select_fun = min,
    simplify_measure = "Wang",
    ...
) {
    if (is.null(category)) {
        category <- c("BP", "MF", "CC")
        message("GO category not specified. Using all three: BP, MF, CC.")
    } else if (!all(category %in% c("BP", "MF", "CC"))) {
        stop("Invalid GO category. Please choose from 'BP', 'MF', 'CC'.")
    }

    gene_list <- .prepare_gene_list(gene_list)

    if (!is.null(universe) && !is.null(universe_list)) {
        stop("Please provide only one of universe or universe_list.")
    }

    if (is.null(universe_list)) {
        if (is.null(universe)) {
            message("Background genes not specified, ",
                "using default (all genes in the database). ")
            universe_list <- lapply(gene_list, function(x) NULL)
        } else {
            universe_list <- lapply(gene_list, function(x) universe)
        }
    }

    if (is.null(names(universe_list))) {
        stop("Input universe_list must be a named list.")
    }
    if (!all(names(gene_list) %in% names(universe_list))) {
        stop("Names of gene_list and universe_list must match.")
    }

    go_list <- list()
    go_list_simplify <- list()
    for (cate in category) {
        message("Performing GO enrichment for category: ", cate)

        for (clus in names(gene_list)) {
            if (length(gene_list[[clus]]) == 0) {
                message("Gene list ", clus, " is empty. Skipping.")
                next
            } else {
                message("Processing gene list: ", clus)
            }

            go_list[[cate]][[clus]] <- clusterProfiler::enrichGO(
                gene = gene_list[[clus]],
                OrgDb = OrgDb,
                ont = cate,
                keyType = keyType,
                universe = universe_list[[clus]],
                pAdjustMethod = pAdjustMethod,
                pvalueCutoff = pvalueCutoff,
                qvalueCutoff = qvalueCutoff,
                ...
            )
        }

        go_list[[cate]] <- go_list[[cate]][
            !vapply(go_list[[cate]], is.null, logical(1))
        ]

        if (simplify && length(go_list[[cate]]) > 0) {
            message("Simplifying GO terms for category: ", cate)

            n_before <- length(go_list[[cate]])

            for (clus in names(go_list[[cate]])) {
                go_list_simplify[[cate]][[clus]] <-
                    clusterProfiler::simplify(go_list[[cate]][[clus]],
                        cutoff = simplify_cutoff,
                        by = simplify_by,
                        select_fun = simplify_select_fun,
                        measure = simplify_measure
                    )
            }

            go_list_simplify[[cate]] <- go_list_simplify[[cate]][
                !vapply(go_list_simplify[[cate]], is.null, logical(1))
            ]

            n_after <- length(go_list_simplify[[cate]])
            if (n_after < n_before) {
                message(n_before - n_after,
                    " input gene set(s) dropped after simplification ", 
                    "(no terms retained).")
            }
        }
    }

    message("Merging GO enrichment results across gene lists ",
        "for each category.")

    go_list_merged <- list()
    for (cate in category) {
        if (length(go_list[[cate]]) == 0) {
            go_list_merged[[cate]] <- NULL
        } else {
            go_list_merged[[cate]] <-
                clusterProfiler::merge_result(go_list[[cate]])
        }
    }

    if (simplify) {
        go_list_merged_simplify <- list()
        for (cate in category) {
            if (length(go_list_simplify[[cate]]) == 0) {
                go_list_merged_simplify[[cate]] <- NULL
            } else {
                go_list_merged_simplify[[cate]] <-
                    clusterProfiler::merge_result(go_list_simplify[[cate]])
            }
        }
        return(list("all" = go_list_merged,
            "simplified" = go_list_merged_simplify))
    } else {
        return(list("all" = go_list_merged))
    }
}
