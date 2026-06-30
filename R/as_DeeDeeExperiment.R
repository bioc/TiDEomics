#' Convert a TiDEomics result to DeeDeeExperiment
#'
#' @description Wraps the output of `prepare_tide()` into a
#'   `DeeDeeExperiment` object.
#'
#' @param tide A named list returned by `prepare_tide()`.
#' @returns A `DeeDeeExperiment` object if the `DeeDeeExperiment`
#'   package is available, otherwise an error.
#' @export
#' @examples
#' data(tutorial_data)
#' data(tutorial_sample_info)
#' tide <- prepare_tide(tutorial_data, tutorial_sample_info,
#'     keep = "threshold", residual_threshold = 100)
#' if (requireNamespace("DeeDeeExperiment", quietly = TRUE)) {
#'     dde <- as_DeeDeeExperiment(tide)
#' }
as_DeeDeeExperiment <- function(tide) {
    .check_list(tide, "tide")
    required <- c("se", "DE", "enrichment", "variance", "feature_property",
                "filter_summary", "WGCNA")
    missing <- setdiff(required, names(tide))
    if (length(missing) > 0) {
        stop("'tide' must be the output of prepare_tide(). Missing elements: ",
            paste(missing, collapse = ", "))
    }
    if (!requireNamespace("DeeDeeExperiment", quietly = TRUE)) {
        stop(
            "Package 'DeeDeeExperiment' is required. ",
            "Install with: BiocManager::install('DeeDeeExperiment')"
        )
    }
    if (!requireNamespace("SingleCellExperiment", quietly = TRUE)) {
        stop("Package 'SingleCellExperiment' is required.")
    }
    if (!"package:SingleCellExperiment" %in% search())
        attachNamespace(getNamespace("SingleCellExperiment"))

    # Build SingleCellExperiment from SummarizedExperiment
    sce <- SingleCellExperiment::SingleCellExperiment(
        assays    = SummarizedExperiment::assays(tide$se),
        colData   = SummarizedExperiment::colData(tide$se),
        rowData   = SummarizedExperiment::rowData(tide$se),
        metadata  = tide$se@metadata
    )

    # ---- DE results: flatten nested structure ----
    # TiDEomics fits one limma model per time point (1-coefficient MArrayLM).
    # DeeDeeExperiment requires >=2 coefficients, so we pass data.frames
    # with columns renamed to the DeeDeeExperiment convention.
    .flatten_de <- function(de_list) {
        if (length(de_list) == 0) return(list())
        # Helper: extract data.frames from a contrast
        # (could be df or list of df per time)
        .extract_contrast <- function(contrast_data, base_name) {
            if (is.data.frame(contrast_data))
                return(stats::setNames(
                    list(.rename_de_cols(contrast_data)), base_name))
            if (is.list(contrast_data)) {
                out <- list()
                for (time in names(contrast_data)) {
                    df <- contrast_data[[time]]
                    if (is.data.frame(df))
                        # Avoid double "T" when time already starts with "t"
                        t_prefix <- if (grepl("^t", time)) "_" else "_T"
                        out[[paste0(base_name, t_prefix, time)]] <-
                            .rename_de_cols(df)
                }
                return(out)
            }
            list()
        }
        # Case 1: Flat DE_between_group() output - all_list at top level
        if (!is.null(de_list$all_list) && is.list(de_list$all_list)) {
            out <- list()
            for (contrast in names(de_list$all_list))
                out <- c(out, .extract_contrast(
                    de_list$all_list[[contrast]], contrast))
            return(out)
        }
        # Case 2: Nested structure - x$DE$source$all_list$contrast$time
        out <- list()
        for (nm in names(de_list)) {
            chunk <- de_list[[nm]]
            if (!is.list(chunk)) next
            if (!is.null(chunk$all_list) && is.list(chunk$all_list)) {
                for (contrast in names(chunk$all_list))
                    out <- c(out, .extract_contrast(chunk$all_list[[contrast]],
                        paste0(nm, "_", contrast)))
            }
        }
        out
    }
    .rename_de_cols <- function(df) {
        colnames(df)[colnames(df) == "logFC"]     <- "log2FoldChange"
        colnames(df)[colnames(df) == "P.Value"]   <- "pvalue"
        colnames(df)[colnames(df) == "adj.P.Val"] <- "padj"
        # Set rownames to feature IDs so DeeDeeExperiment can match to sce
        if ("Feature" %in% colnames(df))
            rownames(df) <- df$Feature
        df
    }
    de_flat <- .flatten_de(tide$DE)
    # DeeDeeExperiment checks !is.null(de_results) before processing;
    # DDE .check_de_results fails on empty list: names(list()) is NULL
    if (length(de_flat) == 0) de_flat <- NULL

    # ---- Enrichment results: unwrap and flatten ----
    # enrichGO_list() returns:
    #   list(all = list(category=df, ...),
    #   unmerged_all = list(category=list(mod=enrichResult)))
    # enrichR_list() / enrich_msigdb() return a flat named list directly.
    .flatten_enrich <- function(enrich_list) {
        if (length(enrich_list) == 0) return(list())
        # Already flat (named list of data.frames / enrichResult objects)
        if (all(vapply(enrich_list, function(x) is.data.frame(x) ||
            inherits(x, "enrichResult"), logical(1))))
            return(enrich_list)
        # enrichGO_list(): prefer unmerged enrichResult objects for DDE
        if (!is.null(enrich_list$unmerged_all) &&
                is.list(enrich_list$unmerged_all)) {
            out <- list()
            for (cate in names(enrich_list$unmerged_all)) {
                for (mod in names(enrich_list$unmerged_all[[cate]])) {
                    nm <- paste0(cate, "_", mod)
                    out[[nm]] <- enrich_list$unmerged_all[[cate]][[mod]]
                }
            }
            return(out)
        }
        # fallback: unwrap enrichGO_list()$all (merge_result data.frames)
        if (!is.null(enrich_list$all) && is.list(enrich_list$all))
            return(.flatten_enrich(enrich_list$all))
        # Nested: x$enrichment$source = list(category=df, ...)
        out <- list()
        for (nm in names(enrich_list)) {
            chunk <- enrich_list[[nm]]
            if (is.data.frame(chunk) || inherits(chunk, "enrichResult")) {
                out[[nm]] <- chunk
            } else if (is.list(chunk)) {
                sub <- .flatten_enrich(chunk)
                for (sub_nm in names(sub)) {
                    key <- if (nm == "all" || is.null(nm) || nm == "")
                        sub_nm else paste0(nm, "_", sub_nm)
                    out[[key]] <- sub[[sub_nm]]
                }
            }
        }
        out
    }
    enrich <- .flatten_enrich(tide$enrichment)
    # Map TiDEomics enrichR column names back to what DeeDeeExperiment expects.
    # enrichR_list() renames Term -> Description,
    # Adjusted.P.value -> p.adjust
    # for plot_modules_h() compatibility; DDE needs the original names.
    .is_enrichr_format <- function(df) {
        is.data.frame(df) &&
            "Description" %in% colnames(df) &&
            "p.adjust" %in% colnames(df) &&
            !"Term" %in% colnames(df) &&
            ("Overlap" %in% colnames(df) || "Genes" %in% colnames(df))
    }

    is_enrichr <- lapply(enrich, .is_enrichr_format)
    for (nm in names(enrich)) {
        if (is_enrichr[[nm]]) {
            enrich[[nm]]$Term <- enrich[[nm]]$Description
            enrich[[nm]]$Adjusted.P.value <- enrich[[nm]]$p.adjust
        }
    }

    # DeeDeeExperiment auto-detects the parser from name prefixes
    # (clusterProfiler_, enrichr_, topGO_, fgsea_, etc.).
    # Add the appropriate prefix if none is present.
    known_prefixes <- c("clusterProfiler_", "enrichr_", "topGO_", "fgsea_",
                        "gsea_", "gProfiler_", "DAVID_", "GeneTonic_")
    has_prefix <- vapply(names(enrich), function(nm) {
        any(vapply(known_prefixes, function(p) startsWith(nm, p), logical(1)))
    }, logical(1))
    if (any(!has_prefix)) {
        # Assign prefix from content:
        # enrichR data.frames -> enrichr_, others -> clusterProfiler_
        for (i in which(!has_prefix)) {
            nm <- names(enrich)[i]
            prefix <- if (is_enrichr[[nm]]) "enrichr_" else "clusterProfiler_"
            names(enrich)[i] <- paste0(prefix, nm)
        }
    }
    # Drop NULL entries and ensure all elements are named
    enrich <- enrich[!vapply(enrich, is.null, logical(1))]
    if (length(enrich) > 0 && is.null(names(enrich)))
        names(enrich) <- rep("", length(enrich))
    enrich <- enrich[names(enrich) != ""]
    # DeeDeeExperiment requires a named list; empty -> NULL
    if (length(enrich) == 0) enrich <- NULL

    dde <- DeeDeeExperiment::DeeDeeExperiment(
        sce            = sce,
        de_results     = de_flat,
        enrich_results = enrich
    )

    dde@metadata$variance <- tide$variance
    dde@metadata$feature_property <- tide$feature_property
    dde@metadata$filter_summary <- tide$filter_summary
    dde@metadata$WGCNA <- tide$WGCNA

    dde
}
