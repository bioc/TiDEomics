#' Flatten nested enrichment results
#'
#' @description
#' Flattens the nested enrichment result structure returned by
#' [enrichGO_list()], [enrichR_list()], or [enrich_msigdb()] into a
#' flat named list of data frames or `enrichResult` objects.
#' This is useful for exporting results to other tools such as
#' `DeeDeeExperiment`.
#'
#' @param enrich_list Output of [enrichGO_list()], [enrichR_list()], or
#'   [enrich_msigdb()]; or a named list combining several such outputs
#'   (e.g. `list(GO = go_res, MSigDB = msigdb_res)`).
#'
#' @returns A flat named list of data frames or `enrichResult` objects
#'   with `DeeDeeExperiment`-compatible names and columns
#'   (e.g. "clusterProfiler_BP_1", "enrichr_DSigDB").
#'   Returns an empty list if `enrich_list` is empty.
#'
#' @export
#'
#' @examples
#' data(example_go)
#' enrich_flat <- flatten_enrich(example_go)
#' str(enrich_flat, max.level = 1)
flatten_enrich <- function(enrich_list) {
    .check_list(enrich_list, "enrich_list")

    .flatten <- function(x) {
        if (length(x) == 0) return(list())

        # Already flat (named list of data.frames / enrichResult objects)
        if (all(vapply(x, function(e) is.data.frame(e) ||
            inherits(e, "enrichResult"), logical(1))))
            return(x)

        # enrichGO_list(): prefer unmerged enrichResult objects
        if (!is.null(x$unmerged_all) && is.list(x$unmerged_all)) {
            out <- list()
            for (cate in names(x$unmerged_all)) {
                for (mod in names(x$unmerged_all[[cate]])) {
                    nm <- paste0(cate, "_", mod)
                    out[[nm]] <- x$unmerged_all[[cate]][[mod]]
                }
            }
            return(out)
        }

        # Fallback: unwrap enrichGO_list()$all (merge_result data.frames)
        if (!is.null(x$all) && is.list(x$all))
            return(.flatten(x$all))

        # Nested: recursively flatten
        out <- list()
        for (nm in names(x)) {
            chunk <- x[[nm]]
            if (is.data.frame(chunk) || inherits(chunk, "enrichResult")) {
                out[[nm]] <- chunk
            } else if (is.list(chunk)) {
                sub <- .flatten(chunk)
                for (sub_nm in names(sub)) {
                    key <- if (nm == "all" || is.null(nm) || nm == "")
                        sub_nm else paste0(nm, "_", sub_nm)
                    out[[key]] <- sub[[sub_nm]]
                }
            }
        }
        out
    }

    enrich <- .flatten(enrich_list)
    if (length(enrich) == 0) return(list())

    # Change TiDEomics enrichR column names back
    # enrichR_list() renames Term -> Description,
    # Adjusted.P.value -> p.adjust for plot_modules_h() compatibility.
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

    # Drop NULL and unnamed entries
    enrich <- enrich[!vapply(enrich, is.null, logical(1))]
    if (length(enrich) > 0 && is.null(names(enrich)))
        names(enrich) <- rep("", length(enrich))
    enrich <- enrich[names(enrich) != ""]
    if (length(enrich) == 0) enrich <- NULL
    enrich
}
