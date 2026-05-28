#' Gene-set enrichment via enrichR
#'
#' @description Enrichment analysis of multiple gene sets against databases
#'   from enrichR (e.g. DSigDB for drug signatures,
#'   ChEA for TF targets, KEGG for pathways, DrugBank for drug targets).
#'
#' @details
#' **Multiple testing:** P-values are adjusted (Benjamini--Hochberg)
#' independently for each module within each database. No correction is
#' applied across databases or across modules.
#'
#' Requires an internet connection.
#'
#' @param gene_list A named list of gene vectors, or a data.frame
#'   with `Feature` and `Module` columns from `WGCNA_module()`.
#'   Gene identifiers must match the convention of the chosen site
#'   (e.g. human gene symbols for Enrichr, fly symbols for FlyEnrichr).
#' @param databases Character vector of enrichR database names to query.
#'   Available databases can be checked with `enrichR::listEnrichrDbs()`.
#'   (default: `c("DSigDB", "DrugMatrix")`)
#' @param site enrichR site to query. See `enrichR::listEnrichrSites()`.
#'   (default: `"Enrichr"`)
#' @param universe Background genes for all input gene sets,
#'   used if `universe_list` is not provided.
#' @param universe_list Background genes for each input gene set,
#'   a list of gene vectors with the same names as `gene_list`.
#' @param pvalueCutoff Adjusted p-value cutoff for filtering enriched
#'   terms (default: 0.05)
#' @param include_overlap Parameter passed to `enrichR::enrichr()`. 
#' If `TRUE`, databases are downloaded during each query to 
#' output 'Overlap' when analysing with a background. (default: `FALSE`)
#'
#' @returns A named list of data.frames, one per database. Each data.frame
#'   has columns `Cluster`, `Description`, `p.adjust` 
#'   (Adjusted.P.value from enrichR output),
#'   `Odds.Ratio`, `Combined.Score`, `Genes`, and any additional
#'   columns returned by the enrichR API. Compatible with
#'   `plot_modules_h(enrich_list = result, enrich_category = "DSigDB")`.
#' @export
#'
#' @examples
#' data(example_net)
#' library(dplyr)
#' example_module <- WGCNA_module(example_net) %>%
#'     dplyr::filter(Module %in% c("1", "2"))
#' # Use high pvalueCutoff for demonstration
#' # enrichr_out <- enrichR_list(example_module, 
#' #     databases = c("KEGG_2019_Mouse"),
#' #     universe = example_module$Feature, pvalueCutoff = 0.5)
enrichR_list <- function(
    gene_list,
    databases = c("DSigDB", "DrugMatrix"),
    site = "Enrichr",
    universe = NULL,
    universe_list = NULL,
    pvalueCutoff = 0.05,
    include_overlap = FALSE
) {
    if (!requireNamespace("enrichR", quietly = TRUE)) {
        stop("Package 'enrichR' is required. Install with: ",
            "install.packages('enrichR')")
    }
    loadNamespace("enrichR")
    getNamespace("enrichR")$.onAttach(NULL, "enrichR")

    site <- match.arg(site, c("Enrichr", "FlyEnrichr", "WormEnrichr",
        "YeastEnrichr", "FishEnrichr", "OxEnrichr"))

    enrichR::setEnrichrSite(site)

    gene_list <- .prepare_gene_list(gene_list)

    dbs_available <- tryCatch(
        enrichR::listEnrichrDbs(),
        error = function(e) {
            message("enrichR server unreachable: ", e$message)
            return(NULL)
        }
    )

    missing <- setdiff(databases, dbs_available$libraryName)
    if (length(missing) > 0) {
        warning("Database(s) not available: ",
            paste(missing, collapse = ", "))
        databases <- intersect(databases, dbs_available$libraryName)
    }
    if (length(databases) == 0) {
        stop("None of the requested databases are available.")
    }

    if (!is.null(universe) & !is.null(universe_list)) {
        stop("Please provide only one of universe or universe_list.")
    }

    if (is.null(universe_list)) {
        if (is.null(universe)) {
            message("Background genes not specified, ",
                "using default (all genes in the database).")
            universe_list <- lapply(gene_list, function(x) NULL)
        } else {
            universe_list <- lapply(gene_list, function(x) universe)
        }
    }

    if (!is.null(names(universe_list)) &&
        !all(names(gene_list) %in% names(universe_list))) {
        stop("Names of gene_list and universe_list must match.")
    }

    # Query all databases per gene list to minimise API calls
    db_tables <- list()
    for (i in names(gene_list)) {
        if (length(gene_list[[i]]) == 0) next
        enr <- enrichR::enrichr(gene_list[[i]],
            databases = databases,
            background = universe_list[[i]],
            include_overlap = include_overlap
        )
        for (db in databases) {
            if (is.null(enr[[db]]) || nrow(enr[[db]]) == 0) next
            db_tables[[db]] <- rbind(
                db_tables[[db]],
                enr[[db]] %>%
                    dplyr::filter(Adjusted.P.value < pvalueCutoff) %>%
                    dplyr::mutate(Cluster = i)
            )
        }
    }

    for (db in names(db_tables)) {
        db_tables[[db]] <- db_tables[[db]] %>%
            dplyr::rename(Description = Term,
                        p.adjust = Adjusted.P.value) %>%
            dplyr::select(Cluster, Description, p.adjust,
                Odds.Ratio, Combined.Score, Genes,
                dplyr::everything()) %>%
            dplyr::arrange(Cluster, p.adjust)
        rownames(db_tables[[db]]) <- NULL
    }

    if (length(db_tables) == 0) {
        message("No enriched terms found.")
        return(list())
    }

    return(db_tables)
}
