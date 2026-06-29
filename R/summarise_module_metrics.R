#' Summarise WGCNA module metrics
#'
#' Computes per-module metrics for a WGCNA network: module size,
#' proportion of features assigned, and module membership (kME).
#' Unassigned module is excluded.
#'
#' `MeanKME` and `MeanKME2` are the mean and mean squared module
#' membership (kME = cor(feature, module eigengene)).
#' `MeanKME2` = mean(kME^2) is the average proportion
#' of per-feature variance explained by the module eigengene.
#'
#' @param net A WGCNA network object from `run_WGCNA()`.
#' @return A data frame with columns: `Module`, `Size`,
#'   `Proportion`, `MeanKME`, `MeanKME2`, `MedianKME`, `SDKME`,
#'   `MinKME`, `MaxKME`.
#' @export
#'
#' @examples
#' data(example_net)
#' summarise_module_metrics(example_net)
summarise_module_metrics <- function(net) {
    .check_list(net, "net", "run_WGCNA")
    if (!"colors" %in% names(net))
        stop("'net' must contain a 'colors' element.")
    if (!"input_data" %in% names(net))
        stop("'net' must contain 'input_data' (use run_WGCNA()).")
    colors <- net$colors
    mod_labels <- .module_labels(colors)
    mod_levels <- setdiff(levels(mod_labels), c("M0", "grey", "gray"))

    if (length(mod_levels) == 0) stop("No non-grey modules to evaluate.")

    me <- WGCNA::moduleEigengenes(
        as.data.frame(net$input_data), colors, excludeGrey = TRUE)$eigengenes

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

    result <- do.call(rbind, lapply(mod_levels, function(m) {
        in_mod <- mod_labels == m
        n_feat <- sum(in_mod)
        mod_num <- gsub("^M", "", m)
        kme_col <- paste0("kME", mod_num)
        kme_vals <- if (kme_col %in% colnames(kme)) {
            kme[in_mod, kme_col]
        } else rep(NA_real_, n_feat)
        me_col <- paste0("ME", mod_num)
        # MeanKME2 = mean(kME^2): average per-feature variance explained
        kme2 <- if (me_col %in% colnames(me) && n_feat > 1) {
            mod_cor <- stats::cor(
                expr_scaled[, in_mod, drop = FALSE],
                me[, me_col], use = "pairwise.complete.obs")
            mean(mod_cor^2, na.rm = TRUE)
        } else NA_real_

        data.frame(Module = m, Size = n_feat,
            Proportion = n_feat / length(colors),
            MeanKME   = mean(kme_vals, na.rm = TRUE),
            MeanKME2  = kme2,
            MedianKME = stats::median(kme_vals, na.rm = TRUE),
            SDKME     = stats::sd(kme_vals, na.rm = TRUE),
            MinKME    = if (n_feat > 0)
                min(kme_vals, na.rm = TRUE) else NA_real_,
            MaxKME    = if (n_feat > 0)
                max(kme_vals, na.rm = TRUE) else NA_real_,
            stringsAsFactors = FALSE)
    }))

    rownames(result) <- NULL
    result
}
