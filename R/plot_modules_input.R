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
