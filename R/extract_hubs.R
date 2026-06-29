#' Extract hub features from WGCNA modules
#'
#' Identifies hub features within each module using module membership
#' (kME).
#' For unsigned networks, absolute kME values are used to reproduce the
#' behaviour of `WGCNA::blockwiseModules()`. For signed and signed
#' hybrid networks, signed kME values are retained.
#'
#' @param net An output object from `run_WGCNA()`
#' @param top_n Number of hub features to return per module. Default: 3.
#' @param exclude_grey Logical; if `TRUE`, exclude the grey/unassigned
#' module. Default: `TRUE`.
#' @return A character vector containing the top `top_n` features
#' from each module concatenated together.
#' @examples
#' data(example_net)
#' hubs <- extract_hubs(example_net, top_n = 3)
#' @export
extract_hubs <- function(
    net,
    top_n = 3,
    exclude_grey = TRUE
) {
    .check_list(net, "net", "run_WGCNA")
    .check_positive_int(top_n, "top_n")
    .check_logical(exclude_grey, "exclude_grey")
    if (is.null(net$input_data)) {
        stop("net must contain 'input_data' (use run_WGCNA())")
    }
    if (is.null(net$parameters$networkType)) {
        stop(
            "net must contain 'parameters$networkType' ",
            "(use run_WGCNA())"
        )
    }

    colors <- net$colors
    modules <- unique(colors)

    if (exclude_grey) {
        modules <- setdiff(modules, c("grey", "0", 0))
    }

    feature_names <- names(colors)

    # kME
    datExpr <- net$input_data
    datExpr <- datExpr[, feature_names, drop = FALSE]

    kME.full <- WGCNA::signedKME(
        datExpr = datExpr,
        datME = net$MEs,
        outputColumnName = "kME"
    )

    # If the network is unsigned, take the absolute value of kME
    if (net$parameters$networkType == "unsigned") {
        kME.full <- abs(kME.full)
    }

    # Assigned-module kME
    assigned.kME <- rep(NA_real_, length(feature_names))
    names(assigned.kME) <- feature_names

    for (mod in modules) {
        features <- feature_names[colors == mod]

        kme.col <- paste0("kME", mod)

        if (!kme.col %in% colnames(kME.full)) {
            message("Cannot find ", kme.col, " in signedKME output; skipping.")
            next
        }

        assigned.kME[features] <-
            kME.full[features, kme.col, drop = TRUE]
    }

    # Hubs
    hubs <- unlist(
        lapply(modules, function(mod) {
            features <- feature_names[colors == mod]

            scores <- assigned.kME[features]

            ranked <- features[order(scores, decreasing = TRUE)]

            utils::head(ranked, top_n)
        }),
        use.names = FALSE
    )

    return(hubs)
}
