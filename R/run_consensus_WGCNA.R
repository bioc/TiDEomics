#' Consensus WGCNA across multiple datasets
#'
#' Run blockwise consensus WGCNA on more than one dataset
#' sharing a common feature ID space.
#'
#' @param nets A named list of WGCNA input objects (output from
#'   `prepare_WGCNA`).
#' @param powers A named numeric vector of soft-thresholding powers to be used
#'   in `WGCNA::blockwiseConsensusModules()`. Values may be selected
#'   based on the output of prepare_WGCNA().
#'   Names must match names of `nets`.
#' @param consensusQuantile Quantile at which consensus is to be defined.
#'   (default: 0). See `WGCNA::blockwiseConsensusModules()` for details.
#' @param numericLabels Whether to use numeric labels for modules in the output
#'   (default is TRUE).
#' @param plot Logical. If TRUE (default), plot consensus dendrogram
#'   with module colours for each feature block via
#'   `WGCNA::plotDendroAndColors()`.
#' @param ... Additional parameters to be passed to
#'   `WGCNA::blockwiseConsensusModules`.
#'   The network type is inherited from `nets`, and a
#'   `networkType` passed via `...` is accepted only if it matches
#'   the value stored in `nets`.
#'
#' @returns A list with elements: `network` (the consensus result),
#'   `multiExpr` (the input multi-expression list),
#'   `common_features` (features shared across all datasets),
#'   `names` (dataset names), and `parameters`.
#' @export
#'
#' @examples
#' data(example_obj)
#' data(example_obj_2)
#' inp1 <- prepare_WGCNA(example_obj,
#'     powers = seq(1, 20), assay = 2, networkType = "signed")
#' inp2 <- prepare_WGCNA(example_obj_2,
#'     powers = seq(1, 20), assay = 2, networkType = "signed")
#' cons <- run_consensus_WGCNA(
#'     nets = list(Exp1 = inp1, Exp2 = inp2),
#'     powers = c(Exp1 = 8, Exp2 = 8),
#'     minModuleSize = 10
#' )
run_consensus_WGCNA <- function(
    nets, powers,
    consensusQuantile = 0,
    numericLabels = TRUE,
    plot = TRUE,
    ...
) {
    .check_list(nets, "nets")
    .check_logical(plot, "plot")
    .check_logical(numericLabels, "numericLabels")
    if (length(nets) < 2) stop("'nets' must contain at least 2 networks.")
    if (is.null(names(nets)) || anyNA(names(nets)) ||
            any(names(nets) == "") || anyDuplicated(names(nets))) {
        stop("'nets' must have unique, non-empty names for all elements.")
    }
    bad <- which(!vapply(nets, is.list, logical(1)))
    if (length(bad) > 0) {
        stop("'nets[[", bad[1],
            "]]' must be a list, the output of prepare_WGCNA().")
    }
    if (!is.numeric(powers) || length(powers) == 0L ||
            anyNA(powers) || any(!is.finite(powers)) || any(powers <= 0)) {
        stop("'powers' must contain finite positive numeric values.")
    }
    if (is.null(names(powers)) || anyNA(names(powers)) ||
            any(names(powers) == "") || anyDuplicated(names(powers)) ||
            !identical(sort(names(powers)), sort(names(nets)))) {
        stop("'powers' must have one uniquely named value for each ",
            "dataset in 'nets'.")
    }
    powers <- powers[names(nets)]

    dataset_names <- names(nets)
    n_datasets <- length(nets)
    dots <- list(...)

    required_fields <- c("data", "networkType")
    missing_fields <- lapply(nets, function(x) {
        setdiff(required_fields, names(x))
    })
    if (any(lengths(missing_fields) > 0L)) {
        bad <- which(lengths(missing_fields) > 0L)[1]
        stop("'nets[[", bad, "]]' is missing: ",
            paste(missing_fields[[bad]], collapse = ", "), ".")
    }

    input_network_types <- vapply(nets, function(x) {
        if (!is.character(x$networkType) || length(x$networkType) != 1L) {
            return(NA_character_)
        }
        x$networkType
    }, character(1))
    if (anyNA(input_network_types) || any(input_network_types == "") ||
        length(unique(input_network_types)) != 1L) {
        stop("All elements of 'nets' must use the same non-empty ",
            "'networkType'.")
    }
    if ("networkType" %in% names(dots)) {
        if (!identical(dots$networkType, input_network_types[[1]])) {
            stop("'networkType' does not match the networkType stored in ",
                "'nets'.")
        }
        dots$networkType <- NULL
    }
    networkType <- input_network_types[[1]]

    invalid_data <- vapply(nets, function(x) {
        !is.data.frame(x$data) && !is.matrix(x$data) ||
            nrow(x$data) == 0L
    }, logical(1))
    if (any(invalid_data)) {
        stop("'nets[[", which(invalid_data)[1], "]]$data' must be a ",
            "non-empty data.frame or matrix.")
    }
    feature_lists <- lapply(nets, function(x) colnames(x$data))
    if (any(vapply(feature_lists, is.null, logical(1))) ||
        any(vapply(feature_lists, function(x) {
            anyNA(x) || any(x == "") || anyDuplicated(x)
        }, logical(1)))) {
        stop("Each 'nets' element must contain 'data' with unique, ",
            "non-empty feature names.")
    }

    # WGCNA::blockwiseConsensusModules accepts a vector of powers
    if (length(unique(powers)) > 1) {
        message(
            "Powers differ across datasets (",
            paste(paste0(dataset_names, "=", powers),
                collapse = ", "),
            "). Passing per-dataset powers."
        )
    }

    # Intersection of features across all datasets
    common_features <- Reduce(intersect, feature_lists)

    feature_sets_differ <- any(vapply(
        feature_lists,
        function(features) !setequal(features, common_features),
        logical(1)
    ))
    if (feature_sets_differ) {
        warning(
            "Feature sets differ across datasets. Consensus WGCNA will use ",
            "the shared feature set. The power values may therefore not ",
            "match this feature space; for best comparability, rerun ",
            "prepare_WGCNA() after restricting each dataset to the shared ",
            "features.",
            call. = FALSE
        )
    }

    message(
        "Common features across ", n_datasets, " datasets: ",
        length(common_features)
    )
    for (nm in dataset_names) {
        message("  ", nm, ": ", length(feature_lists[[nm]]),
            " features")
    }

    # Check minModuleSize
    minModuleSize <- dots$minModuleSize
    if (is.null(minModuleSize)) {
        minModuleSize <- formals(WGCNA::blockwiseConsensusModules)$minModuleSize
    }
    .check_positive_int(minModuleSize, "minModuleSize")

    if (length(common_features) < minModuleSize) {
        stop(
            "Too few common features (", length(common_features),
            ") for consensus analysis (need >= ", minModuleSize, ")."
        )
    }

    # Build multiExpr
    multiExpr <- lapply(dataset_names, function(nm) {
        df <- as.data.frame(
            nets[[nm]]$data[, common_features, drop = FALSE]
        )
        list(data = df)
    })
    names(multiExpr) <- dataset_names

    message("Running blockwise consensus WGCNA...")
    WGCNA::allowWGCNAThreads()

    cons_net <- do.call(WGCNA::blockwiseConsensusModules, c(
        list(
            multiExpr          = multiExpr,
            power              = powers,
            consensusQuantile  = consensusQuantile,
            networkType        = networkType,
            numericLabels      = numericLabels,
            verbose            = 3
        ),
        dots
    ))

    result <- list(
        network = cons_net,
        multiExpr = multiExpr,
        common_features = common_features,
        names = dataset_names,
        parameters = c(
            list(
                powers            = powers,
                consensusQuantile = consensusQuantile,
                networkType       = networkType,
                numericLabels     = numericLabels
            ),
            dots
        )
    )
    message("Consensus WGCNA complete.")

    if (plot) {
        module_colors <- if (numericLabels) {
            WGCNA::labels2colors(result$network$colors)
        } else {
            result$network$colors
        }
        n_mods <- length(unique(module_colors))
        n_blocks <- length(result$network$dendrograms)
        for (block in seq_len(n_blocks)) {
            WGCNA::plotDendroAndColors(
                result$network$dendrograms[[block]],
                module_colors[result$network$blockGenes[[block]]],
                main = paste0(
                    "Consensus WGCNA block ", block, "/", n_blocks, ": ",
                    paste(result$names, collapse = " + "),
                    "\n(", n_mods, " modules)"
                ),
                dendroLabels = FALSE,
                hang = 0.03,
                addGuide = TRUE,
                guideHang = 0.05
            )
        }
    }

    result
}
