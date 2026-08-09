#' Prepare for plotting WGCNA modules
#' @description Prepare the input data for plotting WGCNA modules
#'
#' Here, different from running WGCNA, the input data should have the
#' replicates merged, instead of having multiple samples per group, time and
#' feature (gene).
#'
#' If certain time points are missing in some groups, NA values are added.
#'
#' @param module A data frame with columns "Feature" and "Module"
#' @param se_obj_merged A SummarizedExperiment object, with one value for each
#' feature at each time point in each group (replicates merged). The colData of
#' the object should contain columns "Sample", "Group", and "Time". The object
#' can be produced by `split_groups()`, `merge_replicates()` and
#' `merge_groups()`.
#' @param scale Whether to scale the data (z-score) across samples for each
#' feature
#' @param assay The assay to use in the SummarizedExperiment object: a numeric
#'   index or character name, e.g. 1 or "orig" for original data, 2 or "norm"
#'   for time 0 normalised data. (Required, no default.)
#'
#' @import SummarizedExperiment
#'
#' @returns A long-format data frame suitable for ggplot2, with columns
#' "Feature", "Module", "Abundance", "Group", and "Time"
#' @keywords internal
.plot_modules_input <- function(module, se_obj_merged, assay, scale) {
    # sample annotation
    sp_info <- colData(se_obj_merged) |>
        as.data.frame() |>
        dplyr::select(c(Sample, Group, Time))

    data_wgcna_merged <- assays(se_obj_merged)[[assay]] |>
        as.data.frame()

    if (scale) {
        data_wgcna_merged <- data_wgcna_merged |>
            t() |>
            scale(center = TRUE, scale = TRUE) |>
            t() |>
            as.data.frame()
    }

    data_wgcna_merged <- data_wgcna_merged |>
        tibble::rownames_to_column("Feature")

    if (!any(module$Feature %in% data_wgcna_merged$Feature)) {
        stop("None of the ", length(unique(module$Feature)),
            " features in the module table are present in 'se_obj_merged'. ",
            "Check that se_obj_merged is derived from the same data as the ",
            "WGCNA input.")
    }

    missing_features <- setdiff(module$Feature, data_wgcna_merged$Feature)
    if (length(missing_features) > 0) {
        msg <- paste0(
            length(missing_features),
            " feature(s) in the module are missing from the assay and ",
            "will be excluded: ",
            paste(utils::head(missing_features, 5), collapse = ", "),
            if (length(missing_features) > 5) " ..."
        )
        warning(msg, call. = FALSE)
    }

    data_wgcna_merged <- data_wgcna_merged |>
        dplyr::filter(Feature %in% module$Feature)

    data_module_long <- data_wgcna_merged |>
        tidyr::pivot_longer(
            c(-Feature),
            names_to = "Sample",
            values_to = "Abundance"
        ) |>
        merge(sp_info, by = "Sample") |>
        dplyr::arrange(Group, Time) |>
        dplyr::select(-Sample) |>
        dplyr::mutate(Time = {
            t <- as.character(Time)
            factor(t, levels = unique(sort(as.numeric(t))))
        }) |>
        dplyr::mutate(Feature = as.character(Feature))
        # geom_tile errors with 'AsIs' character

    # fill in missing time points with NA
    data_module_long_filled <- data_module_long |>
        dplyr::group_by(Feature) |>
        tidyr::complete(Group, Time, fill = list(Abundance = NA)) |>
        dplyr::ungroup()

    data_module_long_filled <- merge(module, data_module_long_filled,
        by.x = "Feature", by.y = "Feature") |>
        dplyr::arrange(Group, Time, Module, Feature) |>
        droplevels()

    return(data_module_long_filled)
}

#' Restrict module features to those shared by all datasets
#'
#' Restrict the module table to the features present in every
#' SummarizedExperiment in `se_obj_list`. When
#' non-shared features are present, it warns with the count per dataset
#' and errors when no feature is shared.
#'
#' @param module A data frame with columns "Feature" and "Module"
#' @param se_obj_list A list of SummarizedExperiment objects
#' @param assay The assay to use in the SummarizedExperiment objects: a numeric
#'   index or character name
#' @param names Optional character vector of dataset names used in warning,
#'   default: NULL
#' @returns The subset of `module` with features present in all datasets
#' @keywords internal
.warn_unshared_features <- function(module, se_obj_list, assay, names = NULL) {
    feature_sets <- lapply(seq_along(se_obj_list), function(i) {
        rownames(SummarizedExperiment::assays(se_obj_list[[i]])[[assay]])
    })
    common_feats <- Reduce(intersect, c(list(module$Feature), feature_sets))
    if (length(common_feats) == 0) {
        stop("No common features across the provided SE objects and ",
            "module table. Check that the SE objects are derived from ",
            "the same data as the module table.")
    }
    ds_names <- if (is.null(names)) {
        paste0("Dataset ", seq_along(se_obj_list))
    } else {
        names
    }
    missing_per_dataset <- vapply(seq_along(se_obj_list), function(i) {
        length(setdiff(module$Feature, feature_sets[[i]]))
    }, integer(1))
    dropped_features <- setdiff(module$Feature, common_feats)
    if (length(dropped_features) > 0) {
        idx_missing <- which(missing_per_dataset > 0)
        counts <- vapply(missing_per_dataset[idx_missing], function(n) {
            if (n == 1) "1 feature" else paste0(n, " features")
        }, character(1))
        parts <- paste0(counts, " (", ds_names[idx_missing], ")")
        if (length(parts) > 1) {
            parts <- paste(
                paste(parts[-length(parts)], collapse = ", "),
                utils::tail(parts, 1), sep = " and "
            )
        }
        verb <- if (sum(missing_per_dataset) == 1L) "is" else "are"
        warning(paste0(parts, " ", verb,
            " not shared by input se_obj_merged, and therefore excluded ",
            "from plotting."), call. = FALSE)
    }
    module[module$Feature %in% common_feats, , drop = FALSE]
}
