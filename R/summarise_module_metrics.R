#' Summarise WGCNA module metrics
#'
#' Computes per-module metrics for a WGCNA network: module size,
#' proportion of features assigned, and module membership (kME).
#' Unassigned module is excluded. For consensus WGCNA output, metrics
#' are returned for the consensus and for each individual dataset.
#'
#' `MeanKME` is mean module membership
#' (`kME = cor(feature, module eigengene)`). `MeanKME2` is the mean
#' squared module membership (`mean(kME^2)`), equivalent to the mean
#' per-feature proportion of variance explained by the module eigengene.
#' `Q25KME` and `Q75KME` are the 25th and 75th kME percentiles.
#' For consensus output, it is calculated from consensus
#' or individual-dataset kME scores, respectively.
#'
#' @param net A WGCNA network object from `run_WGCNA()` or
#'   `run_consensus_WGCNA()`.
#' @return For an individual WGCNA network, a data frame with columns:
#'   `Module`, `Size`, `Proportion`, `MeanKME`, `MeanKME2`,
#'   `MedianKME`, `SDKME`, `Q25KME`, `Q75KME`, `MinKME`, and `MaxKME`.
#'   For consensus WGCNA, a named list of these data frames, with
#'   `Consensus` followed by one element for each individual dataset.
#' @export
#'
#' @examples
#' data(example_net)
#' summarise_module_metrics(example_net)
#' data(example_cons)
#' summarise_module_metrics(example_cons)
summarise_module_metrics <- function(net) {
    .check_list(net, "net")

    is_individual <- all(c("colors", "input_data") %in% names(net))
    is_consensus <- all(c(
        "network", "multiExpr", "common_features", "names", "parameters"
    ) %in% names(net))

    if ((!is_individual && !is_consensus) ||
            (is_individual && is_consensus)) {
        stop("'net' must be output from run_WGCNA() or ",
            "run_consensus_WGCNA().")
    }

    .summarise_module_scores <- function(
        colors, kme, expr_scaled = NULL, me = NULL, mod_levels
    ) {
        mod_labels <- .module_labels(colors)

        result <- do.call(rbind, lapply(mod_levels, function(m) {
            in_mod <- mod_labels == m
            n_feat <- sum(in_mod)
            mod_num <- gsub("^M", "", m)
            kme_col <- paste0("kME", mod_num)
            kme_vals <- if (is.data.frame(kme)) { # individual WGCNA
                if (kme_col %in% colnames(kme)) {
                    kme[in_mod, kme_col]
                } else rep(NA_real_, n_feat)
            } else kme[in_mod] # consensus WGCNA (consensus/dataset-specific)
            me_col <- paste0("ME", mod_num)
            # MeanKME2 = mean(kME^2): average per-feature variance explained
            kme2 <- if (!is.null(expr_scaled) && !is.null(me) &&
                    me_col %in% colnames(me) && n_feat > 1) {
                # individual WGCNA
                mod_cor <- stats::cor(
                    expr_scaled[, in_mod, drop = FALSE],
                    me[, me_col], use = "pairwise.complete.obs")
                mean(mod_cor^2, na.rm = TRUE)
            } else if (is.null(expr_scaled)) { # consensus WGCNA
                mean(kme_vals^2, na.rm = TRUE)
            } else NA_real_

            kme_quantiles <- stats::quantile(kme_vals,
                probs = c(0.25, 0.75), na.rm = TRUE, names = FALSE)

            data.frame(Module = m, Size = n_feat,
                Proportion = n_feat / length(colors),
                MeanKME   = mean(kme_vals, na.rm = TRUE),
                MeanKME2  = kme2,
                MedianKME = stats::median(kme_vals, na.rm = TRUE),
                SDKME     = stats::sd(kme_vals, na.rm = TRUE),
                Q25KME    = kme_quantiles[1],
                Q75KME    = kme_quantiles[2],
                MinKME    = if (n_feat > 0)
                    min(kme_vals, na.rm = TRUE) else NA_real_,
                MaxKME    = if (n_feat > 0)
                    max(kme_vals, na.rm = TRUE) else NA_real_,
                stringsAsFactors = FALSE)
        }))

        rownames(result) <- NULL
        result
    }

    if (is_individual) {
        colors <- net$colors
        mod_labels <- .module_labels(colors)
        mod_levels <- setdiff(levels(mod_labels), c("M0", "grey", "gray"))

        if (length(mod_levels) == 0) stop("No non-grey modules to evaluate.")

        me <- WGCNA::moduleEigengenes(
            as.data.frame(net$input_data), colors,
                excludeGrey = TRUE)$eigengenes

        # Use absolute or signed kME depending on network type
        is_unsigned <- !is.null(net$parameters$networkType) &&
            net$parameters$networkType == "unsigned"
        if (is_unsigned) {
            kme <- abs(WGCNA::signedKME(
                as.data.frame(net$input_data), me, outputColumnName = "kME"))
        } else {
            kme <- WGCNA::signedKME(
                as.data.frame(net$input_data), me, outputColumnName = "kME")
        }

        expr_scaled <- scale(as.data.frame(net$input_data))
        return(.summarise_module_scores(
            colors, kme, expr_scaled, me, mod_levels
        ))
    }

    scores <- .consensus_kme_scores(net, exclude_grey = TRUE)
    mod_labels <- .module_labels(scores$colors)
    mod_levels <- setdiff(levels(mod_labels), c("M0", "grey", "gray"))

    if (length(mod_levels) == 0) stop("No non-grey modules to evaluate.")

    output <- list(
        Consensus = .summarise_module_scores(
            scores$colors, scores$consensus, mod_levels = mod_levels
        )
    )
    output <- c(output, lapply(scores$datasets, function(dataset_scores) {
        .summarise_module_scores(
            scores$colors, dataset_scores, mod_levels = mod_levels
        )
    }))
    names(output)[-1] <- names(scores$datasets)
    output
}
