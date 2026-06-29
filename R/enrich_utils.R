# Utility for enrichment functions

#' Prepare input gene list for enrichment functions
#' @description Convert a WGCNA_module() output data.frame to a named gene
#' list, or pass through an existing named list unchanged.
#' @param x A data.frame with 'Feature' and 'Module' columns
#' (as returned by WGCNA_module()), or a named list of gene vectors.
#' @return A named list of gene vectors, where names correspond to module
#' names.
#' @keywords internal
.prepare_gene_list <- function(x) {
    if (is.data.frame(x)) {
        if (!all(c("Feature", "Module") %in% colnames(x))) {
            stop("Input data.frame must have 'Feature' and 'Module' columns ",
                "(as returned by WGCNA_module()).")
        }
        x$Module <- droplevels(as.factor(x$Module))
        x <- split(x$Feature, x$Module)
        # drop empty groups from unused factor levels
        x <- x[lengths(x) > 0]
    }
    if (is.null(names(x))) {
        stop("Input gene_list must be a named list or a data.frame ",
            "from WGCNA_module().")
    }
    return(x)
}

#' Hypergeometric enrichment test
#'
#' @description One-tailed hypergeometric (Fisher's exact) test of query
#'   genes against named gene sets.
#'
#' @param query Character vector of query genes.
#' @param gene_sets Named list of gene vectors (gene set name -> genes,
#'   already filtered to universe).
#' @param universe Character vector of all background genes.
#' @param min_overlap Minimum number of overlapping genes to report
#'   (default: 1).
#'
#' @return A data.frame with columns `Term`, `N`, `m`, `k`, `q`,
#'   `pvalue`, or NULL if no gene set passes filters.
#' @keywords internal
.hypergeometric_test <- function(query, gene_sets, universe,
    min_overlap = 1L) {
    query <- intersect(query, universe)
    k <- length(query)
    N <- length(universe)
    if (k == 0) return(NULL)

    res <- lapply(names(gene_sets), function(gs_name) {
        gs <- intersect(gene_sets[[gs_name]], universe)
        m <- length(gs)
        q <- length(intersect(query, gs))
        if (q < min_overlap) return(NULL)
        data.frame(
            Term    = gs_name,
            N       = N,
            m       = m,
            k       = k,
            q       = q,
            pvalue  = stats::phyper(q - 1, m, N - m, k,
                                    lower.tail = FALSE),
            stringsAsFactors = FALSE
        )
    })
    res <- res[!vapply(res, is.null, logical(1))]
    if (length(res) == 0) return(NULL)
    do.call(rbind, res)
}
