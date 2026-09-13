#' DE over time course
#'
#' @description Differential expression analysis of temporal changes
#'   within and between groups by limma, using time as a factor or a numeric
#'   variable for linear or spline models.
#'
#' @param se_obj A SummarizedExperiment object created by `create_input()`.
#' @param group A character vector specifying which groups to analyse. If NULL,
#'   all groups in the `Group` column are used.
#' @param reference_group A single reference group for a `joint_group` model.
#'   If NULL, the first group in `group` or `colData(se_obj)$Group` is used.
#' @param assay The assay to use in the SummarizedExperiment object: a numeric
#'   index or character name, e.g. 1 or "orig" for original data, 2 or "norm"
#'   for time 0 normalised data. The selected assay should contain
#'   log-transformed, normalised values. A warning is issued if the data
#'   appears to be un-logged raw counts.
#' @param group_mode `"per_group"` or `"joint_group"`.
#'   - `"per_group"` tests each group separately.
#'   - `"joint_group"` fits one Group-by-Time model across the selected groups.
#' @param time_mode `"factor"`, `"numeric"`, or `"spline"`.
#'   - `"factor"` treats time as a categorical variable. The earliest available
#'     time point is the reference, and each later time point is compared with
#'     it. All temporal terms are tested jointly with an omnibus F-test. With
#'     one temporal term, this is equivalent to a two-sided t-test.
#'   - `"numeric"` treats time as continuous and tests one linear trend.
#'   - `"spline"` models time with a natural cubic spline basis, allowing
#'     a smooth non-linear trajectory. The spline coefficients are tested
#'     jointly with an omnibus F-test.
#' @param spline_df An integer of degrees of freedom when
#'   `time_mode = "spline"`, controlling the flexibility of the fitted
#'   trajectory. Larger values allow more complex trajectories but may
#'   increase overfitting. It must be at least 2 and smaller than the
#'   number of time points used by the model. Default: 3.
#' @param trend Logical, passed to `limma::eBayes()`. Use `TRUE` only when the
#'   mean-variance relationship is expected to vary with abundance. The
#'   default `FALSE` is usually appropriate for log-transformed, normalised
#'   data unless diagnostics support modelling a mean-variance trend.
#' @param pAdjustMethod P-value adjustment method passed to
#'   `stats::p.adjust()` for the returned result tables. Must be one of
#'   `stats::p.adjust.methods`. Default: `"BH"`.
#' @details If `colData(se_obj)` contains `Subject`, repeated measurements are
#'   modelled with `limma::duplicateCorrelation()`.
#'   Missingness patterns are defined by whether each Group-Time combination
#'   contains at least one finite value.
#'   In a joint model, `reference_time_effect` tests whether the reference
#'   group changes over time; `group_time_interaction` tests whether the other
#'   groups have different temporal patterns; `global_time_effect` tests all
#'   time and Group-by-Time terms jointly, asking whether any selected group
#'   has a temporal effect; `group_effect` tests selected group differences
#'   at the earliest (reference) time point used by each missingness pattern;
#'   `global_group_effect` tests all Group and Group-by-Time terms jointly,
#'   asking whether group trajectories differ anywhere in the fitted time
#'   course. The joint model is fit on the shared time points
#'   and shared features of all selected groups.
#'   Excluded or non-estimable features are retained in the results with `NA`.
#'   Multiple testing correction is applied within each group for
#'   `per_group` results and separately within each `joint_group` result table,
#'   after combining all missingness patterns.
#'   `ReferenceTime` is the earliest time point available in each fitting
#'   pattern and is used for `group_effect`
#'   in all time modes and for temporal comparisons in factor mode.
#'   It is `NA` for numeric and spline temporal tests.
#'
#' @returns A list with `results`, `fits`, `design`, and `metadata`. In a
#'   `per_group` model, `results` contains `time_effect`. In a
#'   `joint_group` model, `results` contains `reference_time_effect`,
#'   `global_time_effect`, `group_effect`, `global_group_effect`, and
#'   `group_time_interaction`. Result tables report
#'   an omnibus `F` statistic, `P.Value`, and `adj.P.Val` for the selected
#'   terms. When one coefficient is tested, `F = t^2` and the corresponding
#'   limma coefficient estimate (`logFC`) is retained. `ReferenceTime`
#'   identifies the time used for `group_effect`. Features with missing values
#'   are fitted using their available observations when the model remains
#'   estimable. Features with insufficient observations or non-estimable
#'   designs return in the results with `NA` for model coefficients, test
#'   statistics, and p-values.
#' @export
#' @examples
#' data(example_obj)
#'
#' per_group_result <- DE_timecourse(
#'     example_obj, assay = 2, group_mode = "per_group"
#' )
#' joint_result <- DE_timecourse(
#'     example_obj, assay = 2,
#'     group = c("IFNbeta", "IFNgamma"),
#'     group_mode = "joint_group", reference_group = "IFNgamma"
#' )
#' utils::head(joint_result$results$global_time_effect)
DE_timecourse <- function(
    se_obj,
    group = NULL,
    reference_group = NULL,
    assay,
    group_mode = c("per_group", "joint_group"),
    time_mode = c("factor", "numeric", "spline"),
    spline_df = 3,
    trend = FALSE,
    pAdjustMethod = "BH"
) {
    .check_se(se_obj)
    .check_logical(trend, "trend")
    .check_character(pAdjustMethod, "pAdjustMethod")
    if (length(pAdjustMethod) != 1L || is.na(pAdjustMethod) ||
            !nzchar(pAdjustMethod) ||
            !pAdjustMethod %in% stats::p.adjust.methods) {
        stop(
            "'pAdjustMethod' must be one of: ",
            paste(stats::p.adjust.methods, collapse = ", "), "."
        )
    }

    group_mode <- match.arg(group_mode)
    time_mode <- match.arg(time_mode)
    assay <- .match_assay(assay, se_obj)
    assay_data <- assay(se_obj, assay)
    .check_limma_input(assay_data, assay)
    if (any(!is.finite(assay_data) & !is.na(assay_data))) {
        stop("The selected assay must contain finite values or NA.")
    }

    if (time_mode == "spline") {
        if (!is.numeric(spline_df) || length(spline_df) != 1L ||
        !is.finite(spline_df) || spline_df != as.integer(spline_df) ||
        spline_df < 2) {
            stop("'spline_df' must be one integer greater than or equal to 2.")
        }
        spline_df <- as.integer(spline_df)
    }

    sample_data <- as.data.frame(SummarizedExperiment::colData(se_obj))
    if (!all(c("Group", "Time") %in% colnames(sample_data))) {
        stop("'se_obj' must contain 'Group' and 'Time' in colData.")
    }
    sample_data$Group <- as.character(sample_data$Group)
    sample_data$TimeNumeric <- as.numeric(as.character(sample_data$Time))
    if (anyNA(sample_data$TimeNumeric) ||
            any(!is.finite(sample_data$TimeNumeric))) {
        stop("'Time' must be numeric and finite.")
    }

    all_groups <- unique(sample_data$Group)
    if (anyNA(sample_data$Group) || any(!nzchar(sample_data$Group))) {
        stop("'Group' must contain non-missing, non-empty values.")
    }
    if (is.null(group)) {
        group <- all_groups
    } else {
        .check_character(group, "group")
        if (anyDuplicated(group) || !all(group %in% all_groups)) {
            stop("'group' must contain unique groups present in 'se_obj'.")
        }
    }
    if (!is.null(reference_group)) {
        if (group_mode != "joint_group") {
            stop("'reference_group' is only used with 'joint_group'.")
        }
        .check_character(reference_group, "reference_group")
        if (length(reference_group) != 1L || is.na(reference_group) ||
                !reference_group %in% group) {
            stop("'reference_group' must be one selected group.")
        }
        group <- c(reference_group, setdiff(group, reference_group))
    }
    if (group_mode == "joint_group" && length(group) < 2L) {
        stop("'joint_group' requires at least two groups.")
    }

    time_levels <- sort(unique(sample_data$TimeNumeric))
    if (length(time_levels) < 2L) {
        stop("At least two time points are required.")
    }
    if (time_mode == "spline" && spline_df >= length(time_levels)) {
        stop("'spline_df' must be smaller than the number of time points.")
    }

    sample_data$TimeFactor <- factor(
        sample_data$TimeNumeric, levels = time_levels
    )
    sample_data$Group <- factor(sample_data$Group, levels = group)
    expr <- SummarizedExperiment::assay(se_obj, assay)
    expr <- expr[, sample_data$Group %in% group, drop = FALSE]
    sample_data <- sample_data[sample_data$Group %in% group, , drop = FALSE]

    joint_time_levels <- time_levels
    if (group_mode == "joint_group") {
        group_times <- lapply(group, function(current_group) {
            sort(unique(sample_data$TimeNumeric[
                sample_data$Group == current_group
            ]))
        })
        joint_time_levels <- Reduce(intersect, group_times)
        joint_time_levels <- sort(joint_time_levels)
        if (length(joint_time_levels) < 2L) {
            stop(
                "The joint Group-Time model requires at least two time ",
                "points shared by all selected groups."
            )
        }
        if (time_mode == "spline" &&
                spline_df >= length(joint_time_levels)) {
            stop(
                "'spline_df' must be smaller than the number of shared ",
                "time points in the joint model."
            )
        }
    }

    # Common spline spec for different missingness patterns
    make_spline_spec <- function(time_points) {
        time_points <- sort(unique(time_points))
        basis <- splines::ns(time_points, df = spline_df)
        list(
            knots = attr(basis, "knots"),
            boundary_knots = attr(basis, "Boundary.knots")
        )
    }
    spline_spec <- if (time_mode == "spline" &&
            group_mode == "joint_group") {
        make_spline_spec(joint_time_levels)
    } else {
        NULL
    }

    # Features remain in the output even when they cannot be fitted. Invalid
    # features are excluded only from the limma fit and receive NA results.
    build_design <- function(data, spline_spec = NULL) {
        if (time_mode == "factor") {
            if (group_mode == "per_group") {
                return(stats::model.matrix(
                    ~TimeFactor, data = data,
                    contrasts.arg = list(TimeFactor = "contr.treatment")
                ))
            }
            return(stats::model.matrix(
                ~Group * TimeFactor, data = data,
                contrasts.arg = list(
                    Group = "contr.treatment",
                    TimeFactor = "contr.treatment"
                )
            ))
        }
        if (time_mode == "numeric") {
            if (group_mode == "per_group") {
                return(stats::model.matrix(~TimeNumeric, data = data))
            }
            return(stats::model.matrix(
                ~Group * TimeNumeric, data = data,
                contrasts.arg = list(Group = "contr.treatment")
            ))
        }

        if (is.null(spline_spec)) {
            stop("A common spline specification is required.")
        }
        knots <- spline_spec$knots
        boundary_knots <- spline_spec$boundary_knots
        if (group_mode == "per_group") {
            return(stats::model.matrix(
                ~splines::ns(TimeNumeric, knots = knots,
                    Boundary.knots = boundary_knots), data = data
            ))
        }
        stats::model.matrix(
            ~Group * splines::ns(TimeNumeric, knots = knots,
                Boundary.knots = boundary_knots), data = data,
            contrasts.arg = list(Group = "contr.treatment")
        )
    }

    build_design_safely <- function(data, spline_spec, pattern_key) {
        tryCatch(
            build_design(data, spline_spec),
            error = function(err) {
                warning(
                    "Could not build the design for missingness pattern '",
                    pattern_key, "'. Affected features receive NA results: ",
                    conditionMessage(err),
                    call. = FALSE
                )
                NULL
            }
        )
    }

    fit_limma <- function(expression, design, data) {
        if (!"Subject" %in% colnames(data)) {
            return(limma::lmFit(expression, design, maxit = 2000))
        }
        subject <- as.character(data$Subject)
        if (anyNA(subject) || any(!nzchar(subject))) {
            stop("'Subject' must not contain missing or empty values.")
        }
        if (length(unique(subject)) < 2L) {
            stop("'Subject' must contain at least two unique subjects.")
        }
        corfit <- limma::duplicateCorrelation(
            expression, design = design, block = subject
        )
        limma::lmFit(
            expression, design, block = subject,
            correlation = corfit$consensus.correlation, maxit = 2000
        )
    }

    coefficient_columns <- function(design, interaction = FALSE,
            group_effect = FALSE) {
        columns <- colnames(design)
        is_interaction <- grepl("^Group", columns) & grepl(":", columns)
        if (group_effect) {
            return(columns[grepl("^Group", columns) & !is_interaction])
        }
        if (time_mode == "factor") {
            is_time <- grepl("TimeFactor", columns)
        } else if (time_mode == "numeric") {
            is_time <- grepl("TimeNumeric", columns)
        } else {
            is_time <- grepl("TimeNumeric", columns)
        }
        if (interaction) {
            return(columns[is_interaction & is_time])
        }
        columns[is_time & !is_interaction]
    }

    make_group_contrasts <- function(data, spline_spec, reference_time) {
        contrast_data <- data.frame(
            Group = factor(group, levels = levels(data$Group)),
            TimeNumeric = rep(reference_time, length(group)),
            TimeFactor = factor(
                rep(reference_time, length(group)),
                levels = levels(data$TimeFactor)
            )
        )
        profile_design <- build_design(contrast_data, spline_spec)
        profile_difference <- sweep(
            profile_design[-1, , drop = FALSE], 2L,
            profile_design[1, ], FUN = "-"
        )
        contrast_matrix <- t(profile_difference)
        colnames(contrast_matrix) <- group_contrast_names()
        contrast_matrix
    }

    group_contrast_names <- function() {
        paste0("Group", make.names(group[-1L]), "_at_reference_time")
    }

    result_reference_time <- function(time_points) {
        if (time_mode == "factor" && length(time_points) > 0L) {
            return(min(time_points))
        }
        NA_real_
    }

    make_result <- function(fit, columns, test, groups, time_points,
            reference_time = NULL) {
        table <- limma::topTable(
            fit, coef = columns, number = Inf, sort.by = "none",
            adjust.method = "none"
        )
        # limma returns a t-statistic for one coefficient. Convert it to the
        # equivalent one-df F-statistic so every result table uses F.
        if (!"F" %in% colnames(table) && "t" %in% colnames(table)) {
            table$F <- table$t^2
        }
        table$t <- NULL
        table <- tibble::rownames_to_column(table, "Feature")
        table <- cbind(
            data.frame(
                Feature = table$Feature,
                Group = groups,
                Test = test,
                stringsAsFactors = FALSE
            ),
            table[, setdiff(colnames(table), "Feature"), drop = FALSE]
        )
        table$TimePointsTested <- paste(time_points, collapse = ",")
        table$ReferenceTime <- if (is.null(reference_time)) {
            result_reference_time(time_points)
        } else {
            reference_time
        }
        table
    }

    bind_result_rows <- function(result_tables) {
        result_tables <- Filter(function(table) nrow(table) > 0L,
            result_tables)
        if (length(result_tables) == 0L) {
            return(data.frame())
        }
        # result_columns <- unique(unlist(lapply(result_tables, colnames)))
        # result_tables <- lapply(result_tables, function(table) {
        #     missing_columns <- setdiff(result_columns, colnames(table))
        #     for (column in missing_columns) {
        #         table[[column]] <- NA
        #     }
        #     table[, result_columns, drop = FALSE]
        # })
        # do.call(rbind, result_tables)
        result <- as.data.frame(dplyr::bind_rows(result_tables))
        rownames(result) <- NULL
        result
    }

    bind_results <- function(result_tables, adjust_within_group = FALSE) {
        result <- bind_result_rows(result_tables)
        if (nrow(result) == 0L) return(result)
        if ("P.Value" %in% colnames(result)) {
            if (!adjust_within_group) {
                result$adj.P.Val <- stats::p.adjust(
                    result$P.Value, method = pAdjustMethod
                )
            } else {
                result$adj.P.Val <- NA_real_
                for (current_group in unique(result$Group)) {
                    group_rows <- !is.na(result$Group) &
                        result$Group == current_group
                    result$adj.P.Val[group_rows] <- stats::p.adjust(
                        result$P.Value[group_rows], method = pAdjustMethod
                    )
                }
            }
        }
        metadata_columns <- intersect(
            c("Feature", "Group", "Test", "TimePointsTested", "ReferenceTime"),
            colnames(result)
        )
        summary_columns <- intersect(
            c("AveExpr", "F", "P.Value", "adj.P.Val", "B"),
            colnames(result)
        )
        other_columns <- setdiff(
            colnames(result), c(metadata_columns, summary_columns)
        )
        time_columns <- other_columns[grepl("TimeFactor", other_columns)]
        if (length(time_columns) > 0L) {
            time_order <- match(
                as.numeric(sub(".*TimeFactor", "", time_columns)),
                time_levels
            )
            time_columns <- time_columns[order(time_order, na.last = TRUE)]
        }
        other_columns <- c(
            time_columns,
            setdiff(other_columns, time_columns)
        )
        result <- result[, c(
            metadata_columns, other_columns, summary_columns
        ), drop = FALSE]
        result
    }

    complete_result_schema <- function(result, coefficients) {
        summary_columns <- c("AveExpr", "F", "P.Value", "adj.P.Val", "B")
        # limma uses names e.g. `GroupIFNbeta.TimeNumeric` for multiple
        # coefficients, or as `logFC` for single tested coefficient
        schema_coefficients <- if (length(coefficients) > 1L) {
            make.names(coefficients, unique = TRUE)
        } else if (length(coefficients) == 1L) {
            "logFC"
        } else {
            character()
        }
        expected_columns <- c(
            "Feature", "Group", "Test", "TimePointsTested", "ReferenceTime",
            schema_coefficients, summary_columns
        )
        for (column in setdiff(expected_columns, colnames(result))) {
            if (column %in% c(
                "Feature", "Group", "Test", "TimePointsTested"
            )) {
                result[[column]] <- rep(NA_character_, nrow(result))
            } else {
                result[[column]] <- rep(NA_real_, nrow(result))
            }
        }
        extra_columns <- setdiff(colnames(result), expected_columns)
        result[, c(expected_columns, extra_columns), drop = FALSE]
    }

    report_replication <- function(data) {
        comb_counts <- table(data$Group, data$TimeNumeric)
        positive_counts <- as.integer(comb_counts[comb_counts > 0L])
        if (length(positive_counts) == 0L) return(invisible(NULL))
        if (all(positive_counts == 1L)) {
            message(
                "No replicate samples detected: every observed Group-Time ",
                "combination contains one sample."
            )
        } else if (any(positive_counts == 1L)) {
            message(
                "Some Group-Time combinations contain only one sample; ",
                "replicate-level variance cannot be estimated for them."
            )
        }
        invisible(NULL)
    }

    empty_result <- function(features, groups, test, time_points,
            reference_time = NULL) {
        n_features <- length(features)
        result_reference <- if (is.null(reference_time)) {
            result_reference_time(time_points)
        } else {
            reference_time
        }
        data.frame(
            Feature = features,
            Group = rep(groups, length.out = n_features),
            Test = rep(test, length.out = n_features),
            TimePointsTested = rep(
                paste(time_points, collapse = ","),
                length.out = n_features
            ),
            ReferenceTime = rep(
                result_reference,
                length.out = n_features
            ),
            stringsAsFactors = FALSE
        )
    }

    fit_missing_patterns <- function(expression, data, groups,
            spline_spec = NULL) {
        comb_key <- if (group_mode == "joint_group") {
            interaction(data$Group, data$TimeFactor, drop = TRUE, sep = "::")
        } else {
            as.character(data$TimeFactor)
        }
        # row indices of samples for each combination
        comb_indices <- split(seq_len(nrow(data)), comb_key)

        # per feature, 1 or 0 for each combination merged as a binary string
        pattern_keys <- vapply(
            seq_len(nrow(expression)),
            function(feature_index) {
                present <- vapply(comb_indices, function(indices) {
                    any(is.finite(expression[feature_index, indices]))
                }, logical(1))
                paste0(as.integer(present), collapse = "")
            },
            character(1)
        )
        test <- if (group_mode == "per_group") {
            switch(
                time_mode,
                factor = "time_effect",
                numeric = "linear_time_effect",
                spline = "spline_time_effect"
            )
        } else {
            NULL
        }
        time_results <- list()
        interaction_results <- list()
        global_results <- list()
        group_results <- list()
        global_group_results <- list()
        fits <- list()
        designs <- list()

        for (pattern_key in unique(pattern_keys)) {
            rows <- which(pattern_keys == pattern_key)
            bits <- strsplit(pattern_key, "", fixed = TRUE)[[1]]
            available_combs <- names(comb_indices)[bits == "1"]
            sample_indices <- unlist(
                comb_indices[available_combs], use.names = FALSE
            )
            data_fit <- droplevels(data[sample_indices, , drop = FALSE])
            expression_fit <- expression[rows, sample_indices, drop = FALSE]
            groups_supported <- group_mode != "joint_group" || all(
                group %in% unique(as.character(data_fit$Group))
            )
            observed_time_count <- length(unique(data_fit$TimeNumeric))
            spline_supported <- time_mode != "spline" ||
                observed_time_count > spline_df

            design <- if (length(sample_indices) == 0L) {
                NULL
            } else if (!groups_supported) {
                NULL
            } else if (spline_supported) {
                build_design_safely(data_fit, spline_spec, pattern_key)
            } else {
                NULL
            }
            design_key <- pattern_key
            fits[[design_key]] <- NULL
            designs[[design_key]] <- design
            valid_design <- !is.null(design)
            if (valid_design) {
                # reject designs with no residual degrees of freedom or
                # non-estimable model terms before fitting any feature
                residual_df <- nrow(design) - qr(design)$rank
                non_estimable <- limma::nonEstimable(design)
                valid_design <- residual_df > 0L &&
                    length(non_estimable) == 0L
            }

            time_points <- sort(unique(data_fit$TimeNumeric))
            reference_time <- if (length(time_points) > 0L) {
                min(time_points)
            } else {
                NA_real_
            }

            if (!valid_design) {
                valid_indices <- integer(0)
            } else {
                # Recheck each feature after removing its own missing
                # samples. A feature can be invalid even if the common
                # pattern design is valid.
                valid_rows <- vapply(seq_len(nrow(expression_fit)),
                    function(feature_index) {
                        observed <- is.finite(
                            expression_fit[feature_index, ]
                        )
                        observed_times <- unique(
                            data_fit$TimeNumeric[observed]
                        )
                        if (length(observed_times) < 2L) {
                            return(FALSE)
                        }
                        feature_design <- design[observed, , drop = FALSE]
                        feature_rank <- qr(feature_design)$rank
                        feature_rank == ncol(feature_design) &&
                            sum(observed) - feature_rank > 0L
                    }, logical(1))
                valid_indices <- which(valid_rows)
                invalid_indices <- which(!valid_rows)
            }

            # A failed common design and a pattern with no estimable features
            # both receive NA result rows.
            if (length(valid_indices) == 0L) {
                if (group_mode == "per_group") {
                    time_results[[design_key]] <- empty_result(
                        rownames(expression_fit), groups, test, time_points
                    )
                } else {
                    time_results[[design_key]] <- empty_result(
                        rownames(expression_fit), group[1],
                        "reference_time_effect", time_points
                    )
                    interaction_results[[design_key]] <- empty_result(
                        rownames(expression_fit), NA_character_,
                        "group_time_interaction", time_points
                    )
                    global_results[[design_key]] <- empty_result(
                        rownames(expression_fit), NA_character_,
                        "global_time_effect", time_points
                    )
                    group_results[[design_key]] <- empty_result(
                        rownames(expression_fit), NA_character_,
                        "group_effect", time_points, reference_time
                    )
                    global_group_results[[design_key]] <- empty_result(
                        rownames(expression_fit), NA_character_,
                        "global_group_effect", time_points
                    )
                }
                next
            }

            expression_valid <- expression_fit[valid_indices, , drop = FALSE]
            # Fit valid features once, then add explicit NA rows for features
            # that cannot be estimated from this pattern.
            fit_base <- fit_limma(expression_valid, design, data_fit)
            fit <- limma::eBayes(fit_base, trend = trend)
            rownames(fit$coefficients) <- rownames(expression_valid)
            fits[[design_key]] <- fit
            time_columns <- coefficient_columns(design)

            if (group_mode == "per_group") {
                if (length(time_columns) == 0L) {
                    valid_time_result <- data.frame()
                } else {
                    valid_time_result <- make_result(
                        fit, time_columns, test, groups, time_points
                    )
                }
                invalid_time_result <- empty_result(
                    rownames(expression_fit)[invalid_indices], groups, test,
                    time_points
                )
                time_results[[design_key]] <- bind_result_rows(list(
                    valid_time_result, invalid_time_result
                ))
            } else {
                group_effect_columns <- coefficient_columns(
                    design, group_effect = TRUE
                )
                interaction_columns <- coefficient_columns(
                    design, interaction = TRUE
                )
                group_contrasts <- make_group_contrasts(
                    data_fit, spline_spec, reference_time
                )
                group_contrast_fit <- limma::contrasts.fit(
                    fit_base, group_contrasts
                )
                group_contrast_fit <- limma::eBayes(
                    group_contrast_fit, trend = trend
                )
                rownames(group_contrast_fit$coefficients) <-
                    rownames(expression_valid)
                valid_group_result <- make_result(
                    group_contrast_fit, colnames(group_contrasts),
                    "group_effect", NA_character_, time_points,
                    reference_time
                )
                invalid_group_result <- empty_result(
                    rownames(expression_fit)[invalid_indices],
                    NA_character_, "group_effect", time_points,
                    reference_time
                )
                group_results[[design_key]] <- bind_result_rows(list(
                    valid_group_result, invalid_group_result
                ))
                if (length(time_columns) == 0L) {
                    valid_time_result <- data.frame()
                } else {
                    valid_time_result <- make_result(
                        fit, time_columns, "reference_time_effect",
                        group[1], time_points
                    )
                }
                invalid_time_result <- empty_result(
                    rownames(expression_fit)[invalid_indices],
                    group[1], "reference_time_effect", time_points
                )
                time_results[[design_key]] <- bind_result_rows(list(
                    valid_time_result, invalid_time_result
                ))
                if (length(interaction_columns) == 0L) {
                    valid_interaction_result <- data.frame()
                } else {
                    valid_interaction_result <- make_result(
                        fit, interaction_columns, "group_time_interaction",
                        NA_character_, time_points
                    )
                }
                invalid_interaction_result <- empty_result(
                    rownames(expression_fit)[invalid_indices],
                    NA_character_, "group_time_interaction", time_points
                )
                interaction_results[[design_key]] <- bind_result_rows(list(
                    valid_interaction_result, invalid_interaction_result
                ))
                global_columns <- c(time_columns, interaction_columns)
                valid_global_result <- if (length(global_columns) == 0L) {
                    data.frame()
                } else {
                    make_result(
                        fit, global_columns, "global_time_effect",
                        NA_character_, time_points
                    )
                }
                invalid_global_result <- empty_result(
                    rownames(expression_fit)[invalid_indices],
                    NA_character_, "global_time_effect", time_points
                )
                global_results[[design_key]] <- bind_result_rows(list(
                    valid_global_result, invalid_global_result
                ))
                global_group_columns <- c(
                    group_effect_columns, interaction_columns
                )
                valid_global_group_result <- make_result(
                    fit, global_group_columns, "global_group_effect",
                    NA_character_, time_points
                )
                invalid_global_group_result <- empty_result(
                    rownames(expression_fit)[invalid_indices],
                    NA_character_, "global_group_effect", time_points
                )
                global_group_results[[design_key]] <- bind_result_rows(list(
                    valid_global_group_result, invalid_global_group_result
                ))
            }
        }

        fit_output <- if (length(fits) == 1L) fits[[1]] else fits
        design_output <- if (length(designs) == 1L) {
            designs[[1]]
        } else {
            designs
        }
        list(
            time_results = time_results,
            interaction_results = interaction_results,
            global_results = global_results,
            group_results = group_results,
            global_group_results = global_group_results,
            fits = fit_output,
            designs = design_output
        )
    }

    result_list <- list()
    fit_list <- list()
    design_list <- list()
    time_levels_by_group <-
        stats::setNames(vector("list", length(group)), group)

    if (group_mode == "per_group") {
        for (current_group in group) {
            indices <- sample_data$Group == current_group
            data_group <- droplevels(sample_data[indices, , drop = FALSE])
            expr_group <- expr[, indices, drop = FALSE]
            time_levels_by_group[[current_group]] <- sort(
                unique(data_group$TimeNumeric)
            )
            report_replication(data_group)
            if (length(unique(data_group$TimeNumeric)) < 2L) {
                stop(
                    "Group '", current_group,
                    "' requires at least two time points."
                )
            }
            if (time_mode == "spline" &&
                    length(unique(data_group$TimeNumeric)) <= spline_df) {
                stop(
                    "Group '", current_group,
                    "' has too few time points for spline_df = ", spline_df,
                    "."
                )
            }
            group_spline_spec <- if (time_mode == "spline") {
                make_spline_spec(time_levels_by_group[[current_group]])
            } else {
                NULL
            }
            design <- build_design(data_group, group_spline_spec)
            residual_df <- nrow(design) - qr(design)$rank
            if (residual_df <= 0L) {
                stop("The time-course design for group '", current_group,
                    "' has no residual degrees of freedom.")
            }
            non_estimable <- limma::nonEstimable(design)
            if (length(non_estimable) > 0L) {
                stop(
                    "The time-course design for group '", current_group,
                    "' contains non-estimable terms: ",
                    paste(non_estimable, collapse = ", "), "."
                )
            }
            pattern_fit <- fit_missing_patterns(
                expr_group, data_group, current_group, group_spline_spec
            )
            group_result <- bind_result_rows(pattern_fit$time_results)
            result_list[[current_group]] <- complete_result_schema(
                group_result, coefficient_columns(design)
            )
            fit_list[[current_group]] <- pattern_fit$fits
            design_list[[current_group]] <- pattern_fit$designs
        }
        results <- list(
            time_effect = bind_results(result_list, adjust_within_group = TRUE)
        )
    } else {
        for (current_group in group) {
            time_levels_by_group[[current_group]] <- sort(unique(
                sample_data$TimeNumeric[sample_data$Group == current_group]
            ))
        }
        joint_sample <- sample_data$TimeNumeric %in% joint_time_levels
        joint_data <- droplevels(
            sample_data[joint_sample, , drop = FALSE]
        )
        joint_expr <- expr[, joint_sample, drop = FALSE]
        report_replication(joint_data)
        design <- build_design(joint_data, spline_spec)
        residual_df <- nrow(design) - qr(design)$rank
        if (residual_df <= 0L) {
            stop(
                "The joint Group-Time design has no residual degrees of ",
                "freedom."
            )
        }
        non_estimable <- limma::nonEstimable(design)
        if (length(non_estimable) > 0L) {
            stop(
                "The joint Group-Time design contains non-estimable terms: ",
                paste(non_estimable, collapse = ", "), "."
            )
        }
        pattern_fit <- fit_missing_patterns(
            joint_expr, joint_data, group, spline_spec
        )
        group_result <- complete_result_schema(
            bind_results(pattern_fit$group_results),
            group_contrast_names()
        )
        time_result <- complete_result_schema(
            bind_results(pattern_fit$time_results),
            coefficient_columns(design)
        )
        interaction_result <- complete_result_schema(
            bind_results(pattern_fit$interaction_results),
            coefficient_columns(design, interaction = TRUE)
        )
        global_result <- complete_result_schema(
            bind_results(pattern_fit$global_results),
            c(
                coefficient_columns(design),
                coefficient_columns(design, interaction = TRUE)
            )
        )
        global_group_result <- complete_result_schema(
            bind_results(pattern_fit$global_group_results),
            c(
                coefficient_columns(design, group_effect = TRUE),
                coefficient_columns(design, interaction = TRUE)
            )
        )
        results <- list(
            reference_time_effect = time_result,
            global_time_effect = global_result,
            group_effect = group_result,
            global_group_effect = global_group_result,
            group_time_interaction = interaction_result
        )
        fit_list <- list(fit = pattern_fit$fits)
        design_list <- list(design = pattern_fit$designs)
    }

    list(
        results = results,
        fits = fit_list,
        design = design_list,
        metadata = list(
            assay = assay,
            group_mode = group_mode,
            time_mode = time_mode,
            spline_df = if (time_mode == "spline") spline_df else NULL,
            groups = group,
            time_levels = time_levels,
            time_levels_by_group = time_levels_by_group,
            reference_group = if (group_mode == "joint_group") {
                group[1]
            } else {
                NULL
            },
            trend = trend,
            pAdjustMethod = pAdjustMethod,
            features_tested = rownames(expr)
        )
    )
}
