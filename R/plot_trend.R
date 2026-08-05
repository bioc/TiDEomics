#' Plot feature abundance over time
#'
#' @description Plot trend of feature abundances / expression over time by
#' mean and standard deviation (SD) for each group. Accepts either a
#' SummarizedExperiment object (computing mean/SD internally via
#' `calc_mean_sd()`), a pre-computed table from `calc_mean_sd()`, or a
#' named list of SummarizedExperiment objects (multiple datasets,
#' e.g. multiomics).
#'
#' @param se_obj A SummarizedExperiment object, a data.frame from
#'   `calc_mean_sd()` with columns: Feature, Time, Group, Mean, SD, or a
#'   named list of SummarizedExperiment objects (multiple datasets,
#'   e.g. multiomics).
#' @param assay The assay to use in the SummarizedExperiment object: a numeric
#'   index or character name, e.g. 1 or "orig" for original data, 2 or "norm"
#'   for time 0 normalised data. (Default: 1.)
#' @param groups Groups to be plotted, if NULL, all groups will be used
#'   (default is NULL)
#' @param features Features to be plotted, if NULL, an error will be raised
#' @param title Title of the plot. "Feature" if NULL (default) for single
#'   input, "Trends across datasets" if NULL for a named list of
#'   SummarizedExperiment objects.
#' @param ylabel Y axis label of the plot (default is "Log2 abundance").
#' @param errorbar Whether to plot error bars (default is TRUE)
#' @param fontsize Font size for the plot (default is 8)
#' @param multi_layout Facet layout for multiple-datasets input:
#'   `"wrap"` (default): wraps features, datasets nested per feature;
#'   `"vertical"`: datasets as rows, features as columns; `"horizontal"`:
#'   features as rows, datasets as columns; `"overlay"`: one panel per
#'   feature, datasets by linetype.
#'   Ignored for single input.
#' @param ncol Number of feature columns per row for `"wrap"` and `"overlay"`
#'   layouts (multiple datasets) or in `facet_wrap()` (single input). If NULL
#'   (default), auto-computed.
#' @param independent_y If TRUE (default), each panel gets its own y-axis
#'   scale. If FALSE, y-axes are shared: within each facet row for
#'   "vertical" and "horizontal" layouts, fully for "wrap" and "overlay".
#'
#' @import ggplot2
#' @returns Plot of feature abundances over time by mean and SD
#' @export
#' @examples
#' data(example_obj)
#' plot_trend(example_obj,
#'     features = sample(rownames(example_obj), 4))
#'
#' # Also accepts a named list of SummarizedExperiment objects
#' se_list <- list(RNA = example_obj, Protein = example_obj)
#' plot_trend(se_list,
#'     features = sample(rownames(example_obj), 2))
plot_trend <- function(se_obj, assay = 1, groups = NULL, features,
    title = NULL, ylabel = "Log2 abundance", errorbar = TRUE, fontsize = 8,
    multi_layout = c("wrap", "vertical", "horizontal", "overlay"),
    ncol = NULL, independent_y = TRUE) {

    .check_character(features, "features")
    if (!is.null(title)) .check_character(title, "title")
    .check_character(ylabel, "ylabel")
    if (length(ylabel) != 1L) {
        stop("'ylabel' must be a single character string.")
    }
    .check_logical(errorbar, "errorbar")
    .check_positive(fontsize, "fontsize")
    if (!is.null(groups)) .check_character(groups, "groups")
    .check_logical(independent_y, "independent_y")
    if (!is.null(ncol)) .check_positive_int(ncol, "ncol")

    if (is.list(se_obj) && !inherits(se_obj, "SummarizedExperiment") &&
        !is.data.frame(se_obj) && length(se_obj) > 0) {
        multi <- TRUE
        .check_se_list(se_obj)
        if (is.null(title)) title <- "Trends across datasets"
        group_vals <- unique(as.character(unlist(
            lapply(se_obj, function(se) se$Group))))
    } else {
        multi <- FALSE
        if (is.null(title)) title <- "Feature"

        # Accept either an SE object or a pre-computed table_mean_sd
        if (inherits(se_obj, "SummarizedExperiment")) {
            assay <- .match_assay(assay, se_obj)
            # Subset to features to be plotted only
            keep <- intersect(rownames(se_obj), features)
            if (length(keep) == 0) {
                stop("None of the requested features found in the data.")
            }
            se_obj <- se_obj[keep, , drop = FALSE]
            tbl_list <- calc_mean_sd(se_obj)
            table_mean_sd <- tbl_list[[assay]]
        } else if (is.data.frame(se_obj)) {
            required <- c("Feature", "Time", "Group", "Mean", "SD")
            missing <- setdiff(required, colnames(se_obj))
            if (length(missing) > 0) {
                stop("se_obj must be a data.frame produced by calc_mean_sd() ",
                    "with columns: ", paste(required, collapse = ", "),
                    ". Missing: ", paste(missing, collapse = ", "))
            }
            table_mean_sd <- se_obj
        } else {
            stop("se_obj must be one or a list of SummarizedExperiment ",
                "object(s) or a data.frame produced by calc_mean_sd().")
        }

        group_vals <- unique(table_mean_sd$Group)
    }

    if (is.null(groups)) {
        groups <- group_vals
        message(sprintf("Group not specified. Plotting all groups: %s",
            paste(groups, collapse = ", ")))
    } else if (!all(groups %in% group_vals)) {
        stop("At least one of the specified groups is not found ",
            "in the input.")
    }

    if (multi) {
        return(.plot_trend_multi(se_obj_list = se_obj, assay = assay,
            groups = groups, features = features, title = title,
            ylabel = ylabel, errorbar = errorbar, fontsize = fontsize,
            multi_layout = multi_layout, ncol = ncol,
            independent_y = independent_y))
    }

    table_mean_sd$Time <- as.numeric(as.character(table_mean_sd$Time))

    p <- table_mean_sd |>
        dplyr::filter(Group %in% groups & Feature %in% features) |>
        dplyr::mutate(Feature = factor(Feature, levels = features)) |>
        ggplot(aes(x = Time, y = Mean, group = Group, color = Group)) +
        geom_line(linewidth = 0.7) +
        geom_point(size = 0.8) +
        scale_color_manual(values = get_custom_palette(groups)) +
        theme_custom(
            panel_border = TRUE, legend_position = "bottom",
            base_size = fontsize
        ) +
        ggtitle(paste0(title, " in ", paste(groups, collapse = ", "))) +
        xlab("Time") +
        ylab(ylabel) +
        facet_wrap(~Feature,
            scales = if (independent_y) "free_y" else "fixed",
            ncol = ncol)

    if (errorbar) {
        p <- p + geom_errorbar(aes(ymin = Mean - SD, ymax = Mean + SD),
            width = 0.5,
            linewidth = 0.5,
            position = position_dodge(0.05)
        )
    }

    return(p)
}


#' Plot feature abundance over time across datasets
#'
#' @description Internal helper for `plot_trend()`: builds a combined
#' mean +/- SD table across datasets and plots with a faceted grid.
#'
#' @inheritParams plot_trend
#' @param se_obj_list A named list of SummarizedExperiment objects.
#'
#' @import ggplot2
#' @returns A ggplot object.
#' @keywords internal
.plot_trend_multi <- function(se_obj_list, assay = 1, groups = NULL,
    features, title = NULL, ylabel = "Log2 abundance", errorbar = TRUE,
    fontsize = 8, multi_layout = c("wrap", "vertical", "horizontal",
        "overlay"), ncol = NULL, independent_y = TRUE) {

    multi_layout <- match.arg(multi_layout)
    independent <- if (isTRUE(independent_y)) "y" else FALSE

    # Resolve the shared assay across all datasets: one assay applies to all
    resolved <- vapply(names(se_obj_list), function(nm) {
        .match_assay(assay, se_obj_list[[nm]])
    }, character(1))
    if (length(unique(resolved)) > 1) {
        stop("'assay' resolves to different names across datasets: ",
            paste(paste0(names(resolved), "=", resolved), collapse = ", "),
            ". Use a character name shared by all datasets.")
    }
    assay <- unname(resolved[1])

    nms <- names(se_obj_list)
    n_datasets <- length(se_obj_list)

    # Subset each dataset to the requested features, dropping empty ones
    # (avoids calc_mean_sd() on empty subsets)
    keeps <- lapply(nms, function(nm) {
        intersect(rownames(se_obj_list[[nm]]), features)
    })
    if (!any(lengths(keeps) > 0)) {
        stop("None of the requested features found in the data.")
    }

    # Build combined mean+/-SD table
    df_list <- lapply(seq_len(n_datasets), function(i) {
        keep <- keeps[[i]]
        if (length(keep) == 0) {
            return(NULL)
        }
        se <- se_obj_list[[i]][keep, , drop = FALSE]
        tbl <- calc_mean_sd(se)[[assay]]
        tbl$Dataset <- nms[i]
        tbl
    })
    df <- do.call(rbind, df_list)

    df$Feature <- factor(df$Feature, levels = features)
    df$Dataset <- factor(df$Dataset, levels = nms)
    df$Time <- as.numeric(as.character(df$Time))

    # groups resolved by plot_trend()
    df <- df[df$Group %in% groups, , drop = FALSE]

    if (is.null(ncol)) {
        ncol <- min(length(features), ceiling(sqrt(length(features) * 2)))
    }

    # Strip fill colours: custom for <= 5 datasets, otherwise default
    strip_fills <- if (n_datasets <= 5) {
        c(
            "#D6EAF8", "#FADBD8", "#D5F5E3",
            "#FCF3CF", "#E8DAEF"
        )[seq_len(n_datasets)]
    } else {
        scales::pal_hue()(n_datasets)
    }

    # ---- Build facets ----
    if (multi_layout == "vertical") {
        p <- ggplot(df, aes(
            x = Time, y = Mean,
            color = Group, group = Group
        )) +
            ggh4x::facet_grid2(Dataset ~ Feature,
                scales = "free_y", independent = independent,
                strip = ggh4x::strip_themed(
                    background_y = ggh4x::elem_list_rect(
                        fill = strip_fills
                    )
                )
            )
    } else if (multi_layout == "horizontal") {
        p <- ggplot(df, aes(
            x = Time, y = Mean,
            color = Group, group = Group
        )) +
            ggh4x::facet_grid2(Feature ~ Dataset,
                scales = "free_y", independent = independent,
                strip = ggh4x::strip_themed(
                    background_x = ggh4x::elem_list_rect(
                        fill = strip_fills
                    )
                )
            )
    } else if (multi_layout == "wrap") {
        p <- ggplot(df, aes(
            x = Time, y = Mean,
            color = Group, group = Group
        )) +
            ggh4x::facet_wrap2(~ Feature + Dataset,
                ncol = ncol * n_datasets,
                scales = if (independent_y) "free_y" else "fixed",
                strip = ggh4x::strip_nested()
            )
    } else {
        # overlay: one panel per feature, datasets via linetype
        p <- ggplot(df, aes(
            x = Time, y = Mean,
            color = Group, group = interaction(Group, Dataset),
            linetype = Dataset
        )) +
            ggh4x::facet_wrap2(~Feature, ncol = ncol,
                scales = if (independent_y) "free_y" else "fixed") +
            labs(linetype = NULL)
    }

    p <- p +
        geom_line(linewidth = 0.7) +
        geom_point(size = 0.8) +
        scale_color_manual(values = get_custom_palette(groups)) +
        theme_custom(
            panel_border = TRUE,
            legend_position = "bottom",
            base_size = fontsize
        ) +
        labs(x = "Time", y = ylabel, title = paste0(title, " in ",
            paste(groups, collapse = ", ")))

    if (errorbar) {
        p <- p + geom_errorbar(
            aes(ymin = Mean - SD, ymax = Mean + SD),
            width = 0.5, linewidth = 0.5,
            position = position_dodge(0.05)
        )
    }

    p
}
