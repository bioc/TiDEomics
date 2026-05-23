#' Enrich for drug targets
#'
#' @description Enrichment analysis for drug targets of multiple gene sets
#' using enrichR, allowing for using different or same background genes for
#' each gene set
#'
#' @param gene_list A list of gene vectors, need to be gene symbols
#' @param drug_dbs A character vector of drug-related databases in enrichR
#' to use for enrichment analysis. Available databases can be checked with
#' `enrichR::listEnrichrDbs()`. (default is c("DSigDB", "DrugMatrix"))
#' @param universe Background genes for all input gene sets,
#' used if `universe_list` is not provided
#' @param universe_list Background genes for each input gene set,
#' a list of gene vectors with the same names as `gene_list`,
#' need to be gene symbols.
#' @param pvalueCutoff (Optional) Adjusted p-value cutoff for filtering
#' enriched terms (default is 0.05)
#'
#' @returns A data frame of enriched drug targets for each gene set in the
#' input list.
#' @export
#' @examples
#' library(magrittr)
#' data(example_net)
#' example_module <- data.frame(Module = as.factor(example_net$colors)) %>%
#'     tibble::rownames_to_column("Feature") %>% dplyr::arrange(Module)
#' example_module_list <- example_module %>% dplyr::filter(Module != 0) %>%
#'     split(as.character(.$Module)) %>%
#'     lapply(`[[`, "Feature")
#' # drugs_tb = enrich_drug_list(example_module_list, drug_dbs = c("DSigDB"),
#' #     universe = example_module$Feature, pvalueCutoff = 0.1)
enrich_drug_list <- function(
    gene_list,
    drug_dbs = c("DSigDB", "DrugMatrix"),
    universe = NULL,
    universe_list = NULL,
    pvalueCutoff = 0.05
) {
    options(enrichR.base.address = "https://maayanlab.cloud/Enrichr/")

    dbs_available <- tryCatch(
        enrichR::listEnrichrDbs(),
        error = function(e) {
            warning("enrichR server unreachable: ", e$message)
            return(NULL)
        }
    )
    if (is.null(dbs_available)) {
        return(data.frame(Module = character(), Database = character(),
            Term = character(), Adjusted.P.value = numeric(),
            Genes = character(), stringsAsFactors = FALSE))
    }

    if (!all(drug_dbs %in% dbs_available$libraryName)) {
        message("The specified drug database(s) ",
        paste(setdiff(drug_dbs, dbs_available$libraryName), collapse = ", "),
        " are not available ",
        "in enrichR. Please check the available databases with ",
        "`enrichR::listEnrichrDbs()`")
        return(NULL)
    }

    if (!is.null(universe) & !is.null(universe_list)) {
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

    drugs_tb <- list()
    for (i in names(gene_list)) {
        drugs_tb[[i]] <- list()
        enr <- enrichR::enrichr(gene_list[[i]],
            databases = drug_dbs,
            background = universe_list[[i]]
        )
        for (db in drug_dbs) {
            drugs_tb[[i]][[db]] <- enr[[db]] %>%
                dplyr::filter(Adjusted.P.value < pvalueCutoff) %>%
                dplyr::mutate(
                    Module = i,
                    Database = db
                )
        }
        drugs_tb[[i]] <- do.call(rbind, drugs_tb[[i]])
    }
    drugs_tb_all <- do.call(rbind, drugs_tb)

    drugs_tb_all <- drugs_tb_all %>%
        dplyr::select("Module", "Database", "Term", "Adjusted.P.value",
            "Genes", dplyr::everything()) %>%
        dplyr::arrange("Module", "Database", "Adjusted.P.value")

    return(drugs_tb_all)
}
