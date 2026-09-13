#' Extract hub features from consensus WGCNA network
#'
#' Identifies hub features within each consensus module using module
#' membership (kME). Consensus hubs are ranked using `WGCNA::consensusKME()`
#' for signed and signed-hybrid networks, and absolute kME for unsigned
#' networks.
#'
#' @param consensus An output object from `run_consensus_WGCNA()`.
#' @param top_n Number of hub features to return per module, for each
#'   dataset and consensus. Default: 1.
#' @param exclude_grey Logical; if `TRUE`, exclude the grey/unassigned
#'   module. Default: `TRUE`.
#'
#' @section Requirements:
#' The `consensus` object must be the output of `run_consensus_WGCNA()`,
#' which augments the raw `WGCNA::blockwiseConsensusModules()` output
#' with `multiExpr` and `parameters`.
#'
#' @returns A data frame with columns `Feature`, `Module`, `Set`, `Rank`,
#'   `Score`, and `ScoreType`. `Set` identifies the consensus hub or a
#'   dataset-specific hub set, for example `"Consensus hub"` or
#'   `"RNA hub"`. Duplicate features across sets are allowed.
#'
#' @examples
#' data(example_cons)
#' hubs <- extract_consensus_hubs(example_cons, top_n = 1)
#' utils::head(hubs)
#'
#' @export
extract_consensus_hubs <- function(
    consensus,
    top_n = 1,
    exclude_grey = TRUE
) {
    .check_list(consensus, "consensus", "run_consensus_WGCNA")
    .check_positive_int(top_n, "top_n")
    .check_logical(exclude_grey, "exclude_grey")

    empty_output <- function() {
        data.frame(
            Feature = character(0), Module = character(0),
            Set = character(0), Rank = integer(0),
            Score = numeric(0), ScoreType = character(0),
            stringsAsFactors = FALSE
        )
    }

    scores <- .consensus_kme_scores(
        consensus, exclude_grey = exclude_grey
    )
    if (length(scores$modules) == 0) {
        message("No modules to extract hubs from.")
        return(empty_output())
    }

    colors <- scores$colors
    feature_names <- scores$features
    modules <- scores$modules
    dataset_names <- names(scores$datasets)
    hub_set_names <- paste0(dataset_names, " hub")
    score_sets <- c(
        list("Consensus hub" = scores$consensus),
        stats::setNames(scores$datasets, hub_set_names)
    )
    score_types <- c(
        "Consensus hub" = scores$consensus_type,
        stats::setNames(rep(scores$dataset_type, length(dataset_names)),
            hub_set_names)
    )
    rows <- lapply(names(score_sets), function(set_name) {
        scores <- score_sets[[set_name]]
        module_rows <- lapply(modules, function(mod) {
            features <- feature_names[colors == mod]
            module_scores <- scores[features]
            keep <- is.finite(module_scores)
            if (!any(keep)) return(NULL)

            features <- features[keep]
            module_scores <- module_scores[keep]
            selected <- utils::head(
                order(module_scores, decreasing = TRUE), top_n
            )

            data.frame(
                Feature = features[selected],
                Module = as.character(mod),
                Set = set_name,
                Rank = seq_along(selected),
                Score = unname(module_scores[selected]),
                ScoreType = unname(score_types[set_name]),
                stringsAsFactors = FALSE
            )
        })
        module_rows[!vapply(module_rows, is.null, logical(1))]
    })
    rows <- do.call(c, rows)

    if (length(rows) == 0L) empty_output() else do.call(rbind, rows)
}

#' Calculate consensus and dataset-specific kME scores
#'
#' Runs `WGCNA::consensusKME()` and extracts assigned-module scores for
#' each dataset and for the consensus network.
#'
#' @param consensus Output from `run_consensus_WGCNA()`.
#' @param exclude_grey Logical; whether to exclude the grey module.
#' @return A list containing module assignments, score vectors, and labels.
#' @keywords internal
.consensus_kme_scores <- function(consensus, exclude_grey = TRUE) {
    .check_list(consensus, "consensus", "run_consensus_WGCNA")
    .check_logical(exclude_grey, "exclude_grey")

    network <- consensus$network
    multiExpr <- consensus$multiExpr
    parameters <- consensus$parameters

    if (!is.list(network)) {
        stop("'consensus$network' must be a list returned by ",
            "WGCNA::blockwiseConsensusModules().")
    }
    if (!is.list(multiExpr) || length(multiExpr) == 0L ||
        is.null(network$multiMEs) ||
        !is.list(network$multiMEs) ||
        length(multiExpr) != length(network$multiMEs)) {
        stop("'consensus$multiExpr' and 'consensus$network$multiMEs' ",
            "must be non-empty lists of the same length.")
    }
    if (!is.list(parameters) || is.null(parameters$networkType) ||
            is.null(parameters$consensusQuantile)) {
        stop("'consensus$parameters' must contain 'networkType' and ",
            "'consensusQuantile'.")
    }

    colors <- network$colors
    modules <- unique(colors)

    if (exclude_grey) {
        modules <- setdiff(modules, c("grey", "0", 0))
    }

    if (length(modules) == 0L) {
        return(list(
            colors = colors,
            modules = modules
        ))
    }

    feature_names <- names(colors)
    dataset_names <- names(multiExpr)

    consensus_kme <- WGCNA::consensusKME(
        multiExpr = multiExpr,
        moduleLabels = colors,
        multiEigengenes = network$multiMEs,
        consensusQuantile = parameters$consensusQuantile,
        signed = parameters$networkType != "unsigned",
        setNames = dataset_names,
        excludeGrey = exclude_grey,
        useRankPvalue = FALSE
    )

    # Dataset-specific assigned-module kME
    kME.list <- lapply(dataset_names, function(dataset) {
        assigned.kME <- rep(NA_real_, length(feature_names))
        names(assigned.kME) <- feature_names

        for (mod in modules) {
            features <- feature_names[colors == mod]
            kme.col <- paste0("kME", mod, ".", dataset)

            if (!kme.col %in% colnames(consensus_kme)) {
                stop("Cannot find ", kme.col,
                    " in consensusKME output.")
            }

            feature_idx <- match(features, consensus_kme$ID)
            assigned.kME[features] <-
                consensus_kme[[kme.col]][feature_idx]
        }

        assigned.kME
    })
    names(kME.list) <- dataset_names

    if (parameters$networkType == "unsigned") {
        # consensusKME() retains signed correlations when signed = FALSE.
        # For unsigned networks, use absolute kME to include anti-correlated
        # features
        score_sets <- lapply(kME.list, abs)
        consensus_scores <- apply(
            do.call(cbind, score_sets),
            1,
            function(scores) {
                scores <- scores[is.finite(scores)]
                if (length(scores) == 0L) return(NA_real_)
                stats::quantile(
                    scores,
                    probs = parameters$consensusQuantile,
                    names = FALSE
                )
            }
        )
        consensus_type <- "quantile(abs(kME))"
        dataset_type <- "abs(kME)"
    } else {
        score_sets <- kME.list
        consensus_scores <- rep(NA_real_, length(feature_names))
        names(consensus_scores) <- feature_names
        for (mod in modules) {
            kme.col <- paste0("consensus.kME", mod)
            if (!kme.col %in% colnames(consensus_kme)) {
                stop("Cannot find ", kme.col,
                    " in consensusKME output.")
            }
            features <- feature_names[colors == mod]
            feature_idx <- match(features, consensus_kme$ID)
            consensus_scores[features] <-
                consensus_kme[[kme.col]][feature_idx]
        }
        consensus_type <- "consensus.kME"
        dataset_type <- "kME"
    }

    list(
        colors = colors,
        features = feature_names,
        modules = modules,
        consensus = consensus_scores,
        datasets = score_sets,
        consensus_type = consensus_type,
        dataset_type = dataset_type
    )
}
