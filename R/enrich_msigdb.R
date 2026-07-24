#' Gene set enrichment via MSigDB
#'
#' @description Enrichment analysis against MSigDB gene sets via
#'   hypergeometric test. Supports MSigDB categories, subcategories,
#'   and specified gene sets. Requires the `msigdbr` package.
#'
#' Use cases:
#' - Drug perturbations: `category = "C2", subcategory = "CGP"`
#'   (chemical and genetic perturbations)
#' - TF targets: `category = "C3", subcategory = "TFT:GTRD"`
#' - Hallmarks: `category = "H", subcategory = NULL` (no subcategory)
#' - GO: `category = "C5", subcategory = "BP"` (GO biological process)
#'
#' Multiple testing: P-values are adjusted (Benjamini-Hochberg)
#' per module across the gene-set tests within that module.
#'
#' @param gene_list A named list of gene vectors, or a data.frame
#'   with `Feature` and `Module` columns from `WGCNA_module()`.
#' @param universe Background genes (required).
#' @param category MSigDB category (default: `NULL`).
#'   See `msigdbr::msigdbr_collections(db_species = "HS"/"MM")`
#'   for available options.
#'   Not required when `gene_sets` is provided; in that case all
#'   collections are searched to locate the named gene sets.
#' @param subcategory MSigDB subcategory (default: `NULL`).
#'   Set to NULL to include all subcategories of `category`.
#' @param gene_sets Optional character vector of specific MSigDB gene set
#'   names to include (e.g., `c("HALLMARK_APOPTOSIS",
#'   "HALLMARK_P53_PATHWAY")`). When provided, `msigdbr` is queried across
#'   all collections to find these gene sets, so `category` and
#'   `subcategory` are not required. If NULL (default), all gene sets in
#'   the category/subcategory are used.
#' @param species Species of input genes (default: `"Homo sapiens"`).
#'   See `msigdbr::msigdbr_species()` for available options.
#' @param db_species Species of MSigDB genesets (default: `"HS"`).
#'   Use "MM" for mouse-native gene sets. If species and db_species mismatch,
#'   `db_species` MSigDB will be used with ortholog mapping to `species`.
#' @param name Name for the output list element. If NULL, auto-derived
#'   from `category` (e.g., `"C2"` names the output element `"C2"`). If
#'   `category` is also NULL, it uses `"msigdb"`.
#' @param pvalueCutoff Adjusted p-value cutoff (default: 0.05).
#' @param pAdjustMethod P-value adjustment method (default: `"BH"`).
#' @param minGSSize Minimum gene set size (default: 10).
#' @param maxGSSize Maximum gene set size (default: 500).
#' @returns A named list with one element: a data.frame with columns
#'   `Cluster`, `ID`, `Description`, `GeneRatio`, `BgRatio`, `pvalue`,
#'   `p.adjust`, `Count`, `geneID`.
#' @export
#'
#' @examples
#' if (requireNamespace("msigdbr", quietly = TRUE)) {
#'     data(example_net)
#'     example_module <- WGCNA_module(example_net) |>
#'         dplyr::filter(Module %in% c("1", "2"))
#'     hallmark_msigdb <- enrich_msigdb(example_module, category = "MH",
#'         species = "Mus musculus", db_species = "MM",
#'         minGSSize = 1,
#'         pvalueCutoff = 0.9,
#'         universe = WGCNA_module(example_net, exclude_grey = FALSE)$Feature)
#' }
enrich_msigdb <- function(
    gene_list,
    universe,
    category = NULL,
    subcategory = NULL,
    gene_sets = NULL,
    species = "Homo sapiens",
    db_species = "HS",
    name = NULL,
    pvalueCutoff = 0.05,
    pAdjustMethod = "BH",
    minGSSize = 10,
    maxGSSize = 500
) {
    .check_pval(pvalueCutoff, "pvalueCutoff")
    .check_character(pAdjustMethod, "pAdjustMethod")
    .check_character(species, "species")
    .check_character(db_species, "db_species")
    .check_positive_int(minGSSize, "minGSSize")
    .check_positive_int(maxGSSize, "maxGSSize")
    .check_character(universe, "universe")
    if (!is.null(category)) .check_character(category, "category")
    if (!is.null(subcategory)) .check_character(subcategory, "subcategory")
    if (!is.null(gene_sets)) .check_character(gene_sets, "gene_sets")
    if (!is.null(name)) .check_character(name, "name")
    db_species <- match.arg(db_species, c("HS", "MM"))

    if (is.null(category) && is.null(gene_sets)) {
        stop("Either 'category' or 'gene_sets' must be specified.")
    }

    if (!requireNamespace("msigdbr", quietly = TRUE)) {
        stop("Package 'msigdbr' is required. ",
            "Install with: BiocManager::install('msigdbr')")
    }

    # --- Load MSigDB gene sets ---
    # When gene_sets is provided, query all collections to find them
    # (category/subcategory are optional filters)
    if (!is.null(gene_sets)) {
        msig <- msigdbr::msigdbr(species = species, db_species = db_species)
    } else {
        msig <- msigdbr::msigdbr(species = species, db_species = db_species,
            collection = category, subcollection = subcategory)
    }

    if (nrow(msig) == 0) {
        stop("No gene sets found for specified collection, ",
            "species '", species, "'.")
    }

    # --- Filter to specific gene sets if requested ---
    if (!is.null(gene_sets)) {
        found <- intersect(gene_sets, msig$gs_name)
        missing <- setdiff(gene_sets, msig$gs_name)
        if (length(missing) > 0) {
            pst_missing <- paste(missing, collapse = ", ")
            warning("Gene set(s) not found: ", pst_missing)
        }
        if (length(found) == 0) {
            stop("None of the requested gene sets were found ",
                "in species '", species, "'.")
        }
        msig <- msig[msig$gs_name %in% found, ]
        message("Using ", length(found), " gene set(s).")
    }

    term2gene <- msig[, c("gs_name", "gene_symbol")]
    colnames(term2gene) <- c("Term", "Gene")

    if (is.null(name)) {
        name <- if (is.null(category)) "msigdb" else category
    }

    gene_list <- .prepare_gene_list(gene_list)
    universe <- unique(universe)

    # --- Build gene sets from term2gene, filtered to universe ---
    term2gene <- term2gene[term2gene$Gene %in% universe, , drop = FALSE]
    if (nrow(term2gene) == 0) {
        warning("No MSigDB genes found in universe for ", name, ".")
        return(stats::setNames(list(NULL), name))
    }
    gs_list <- split(term2gene$Gene, term2gene$Term)
    gs_sizes <- lengths(gs_list)
    keep <- gs_sizes >= minGSSize & gs_sizes <= maxGSSize
    gs_list <- gs_list[keep]
    if (length(gs_list) == 0) {
        warning("No gene sets with size in [", minGSSize, ", ",
            maxGSSize, "] for ", name, ".")
        return(stats::setNames(list(NULL), name))
    }

    # --- Per-module hypergeometric test ---
    all_results <- list()
    for (clus in names(gene_list)) {
        if (length(gene_list[[clus]]) == 0) {
            message("Gene list ", clus, " is empty. Skipping.")
            next
        }
        message("Processing gene list: ", clus)

        hg <- .hypergeometric_test(
            query     = gene_list[[clus]],
            gene_sets = gs_list,
            universe  = universe
        )
        if (is.null(hg)) next

        hg$Cluster     <- clus
        hg$ID          <- hg$Term
        hg$Description <- hg$Term
        hg$GeneRatio   <- paste0(hg$q, "/", hg$k)
        hg$BgRatio     <- paste0(hg$m, "/", hg$N)
        hg$Count       <- hg$q
        hg$geneID      <- vapply(hg$Term, function(t) {
            paste(intersect(gene_list[[clus]],
                            intersect(gs_list[[t]], universe)),
                    collapse = ";")
        }, character(1))
        # Per-module multiple-testing correction
        hg$p.adjust <- stats::p.adjust(hg$pvalue,
            method = pAdjustMethod)

        if (nrow(hg) > 0) {
            all_results[[clus]] <- hg
        }
    }

    if (length(all_results) == 0) {
        warning("No enriched terms found for ", name, ".")
        return(stats::setNames(list(NULL), name))
    }

    combined <- do.call(rbind, all_results)
    rownames(combined) <- NULL

    # Filter by cutoff (p.adjust already computed per-module above)
    combined <- combined[!is.na(combined$p.adjust) &
        combined$p.adjust <= pvalueCutoff, , drop = FALSE]

    if (nrow(combined) == 0) {
        warning("No enriched terms pass cutoffs for ", name, ".")
        return(stats::setNames(list(NULL), name))
    }

    combined <- combined[order(combined$Cluster, combined$p.adjust),
        c("Cluster", "ID", "Description", "GeneRatio",
        "BgRatio", "pvalue", "p.adjust",
        "Count", "geneID"), drop = FALSE]

    stats::setNames(list(combined), name)
}
