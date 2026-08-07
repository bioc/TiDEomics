#' Extract hub features from consensus WGCNA network
#'
#' Identifies hub features within each consensus module using module
#' membership (kME). Consensus hubs are ranked using consensus kME
#' calculated across datasets according to the consensus quantile used
#' during network construction. Dataset-specific hubs are ranked using
#' dataset-specific assigned-module kME and exclude features selected
#' as consensus hubs from the same module.
#'
#' For unsigned networks, absolute kME values are used to reproduce the
#' behaviour of `WGCNA::blockwiseConsensusModules()`. For signed and
#' signed hybrid networks, signed kME values are retained.
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
#' @returns A named list.
#'     - `Consensus`: consensus hub features.
#'     - One element per dataset containing dataset-specific hub
#'       features.
#' Each element is a character vector containing the top `top_n`
#' features from each module concatenated together.
#'
#' @examples
#' data(example_cons)
#' hubs <- extract_consensus_hubs(example_cons, top_n = 1)
#' names(hubs)
#' hubs$Consensus
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

    if (length(modules) == 0) {
        message("No modules to extract hubs from.")
        empty <- stats::setNames(rep(list(character(0)),
            length(multiExpr)), names(multiExpr))
        return(c(list(Consensus = character(0)), empty))
    }

    feature_names <- names(colors)

    # Dataset-specific assigned-module kME
    kME.list <- lapply(seq_along(multiExpr), function(set_idx) {
        if (!is.list(multiExpr[[set_idx]]) ||
                !is.data.frame(multiExpr[[set_idx]]$data) &&
                !is.matrix(multiExpr[[set_idx]]$data)) {
            stop("'consensus$multiExpr[[", set_idx,"]]$data' must be a ",
                "data.frame or matrix.")
        }
        datExpr <- multiExpr[[set_idx]]$data
        datExpr <- datExpr[, feature_names, drop = FALSE]

        kME.full <- WGCNA::signedKME(
            datExpr = datExpr,
            datME = network$multiMEs[[set_idx]]$data,
            outputColumnName = "kME"
        )

        # If the network is unsigned, take the absolute value of kME
        if (parameters$networkType == "unsigned") {
            kME.full <- abs(kME.full)
        }

        assigned.kME <- rep(NA_real_, length(feature_names))
        names(assigned.kME) <- feature_names

        for (mod in modules) {
            features <- feature_names[colors == mod]

            kme.col <- paste0("kME", mod)

            if (!kme.col %in% colnames(kME.full)) {
                stop("Cannot find ", kme.col,
                    " in signedKME output for dataset ", set_idx, ".")
            }

            assigned.kME[features] <-
                kME.full[features, kme.col, drop = TRUE]
        }

        assigned.kME
    })

    # Consensus kME
    consensus.kME <- apply(
        do.call(cbind, kME.list),
        1,
        stats::quantile,
        probs = parameters$consensusQuantile,
        na.rm = TRUE,
        names = FALSE
    )

    names(consensus.kME) <- feature_names

    out <- list()

    # Consensus hubs
    consensus.features <- unlist(
        lapply(modules, function(mod) {
            features <- feature_names[colors == mod]

            scores <- consensus.kME[features]

            features[order(scores, decreasing = TRUE)][
                seq_len(min(top_n, length(features)))
            ]
        }),
        use.names = FALSE
    )

    out[["Consensus"]] <- consensus.features

    # Dataset-specific hubs (excluding consensus hubs)
    dataset.names <- names(multiExpr)

    for (set_idx in seq_along(kME.list)) {
        dataset.hubs <- unlist(
            lapply(modules, function(mod) {
                features <- feature_names[colors == mod]

                scores <- kME.list[[set_idx]][features]

                ranked <- features[order(scores, decreasing = TRUE)]

                # Exclude consensus hubs for the module
                consensus.mod.features <- consensus.features[
                    colors[consensus.features] == mod
                ]

                ranked <- setdiff(ranked, consensus.mod.features)

                utils::head(ranked, top_n)
            }),
            use.names = FALSE
        )

        out[[dataset.names[set_idx]]] <- dataset.hubs
    }

    out
}
