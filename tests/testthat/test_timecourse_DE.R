data("example_obj")

expect_result_statistics <- function(table, require_finite = FALSE) {
    expect_true(all(c("F", "P.Value", "adj.P.Val") %in% names(table)))
    if (require_finite) {
        expect_true(all(is.finite(table$F)))
        expect_true(all(is.finite(table$P.Value)))
    }
}

test_that("DE_timecourse returns global per-group time tests", {
    result <- DE_timecourse(
        example_obj, assay = 2,
        group_mode = "per_group", time_mode = "factor"
    )

    expect_true(is.list(result))
    expect_true(all(c("results", "fits", "design", "metadata") %in%
        names(result)))
    expect_result_statistics(result$results$time_effect)
    expect_true(all(c("Feature", "Group", "Test", "TimePointsTested",
        "ReferenceTime") %in% names(result$results$time_effect)))
    expect_equal(nrow(result$results$time_effect), 400L)
    expect_equal(result$metadata$group_mode, "per_group")
    expect_equal(result$metadata$pAdjustMethod, "BH")
})

test_that("DE_timecourse supports numeric and spline time models", {
    numeric_result <- DE_timecourse(
        example_obj, assay = 2, group = "IFNbeta",
        group_mode = "per_group", time_mode = "numeric"
    )
    spline_result <- DE_timecourse(
        example_obj, assay = 2, group = "IFNbeta",
        group_mode = "per_group", time_mode = "spline", spline_df = 2
    )

    time_results <- list(
        numeric = numeric_result$results$time_effect,
        spline = spline_result$results$time_effect
    )
    for (time_result in time_results) {
        expect_result_statistics(time_result, require_finite = TRUE)
        expect_false("t" %in% names(time_result))
    }
    expect_false("TimeNumeric" %in% names(time_results$numeric))
    expect_true("logFC_TimeSlope" %in% names(time_results$numeric))
    expect_false("logFC" %in% names(time_results$numeric))
    spline_coefficients <- c("SplineBasis_1", "SplineBasis_2")
    expect_true(all(spline_coefficients %in% names(time_results$spline)))
    expect_false(any(grepl(
        "splines", names(time_results$spline), fixed = TRUE
    )))
    expect_equal(length(spline_coefficients), 2L)
    expect_true(all(colSums(
        !is.na(time_results$spline[spline_coefficients])
    ) > 0L))
    expect_false("fitted_range" %in% names(time_results$spline))
    spline_designs <- spline_result$design[["IFNbeta"]]
    if (is.matrix(spline_designs)) spline_designs <- list(spline_designs)
    spline_fits <- spline_result$fits[["IFNbeta"]]
    if (inherits(spline_fits, "MArrayLM")) spline_fits <- list(spline_fits)
    spline_design_names <- unlist(lapply(
        spline_designs,
        function(design) if (is.null(design)) character() else {
            colnames(design)
        }
    ), use.names = FALSE)
    spline_fit_names <- unlist(lapply(
        spline_fits,
        function(fit) if (is.null(fit)) character() else {
            colnames(fit$coefficients)
        }
    ), use.names = FALSE)
    expect_true(any(grepl(
        "splines::ns", spline_design_names, fixed = TRUE
    )))
    expect_true(any(grepl(
        "splines::ns", spline_fit_names, fixed = TRUE
    )))

    numeric_fit <- numeric_result$fits[["IFNbeta"]]
    numeric_table <- limma::topTable(
        numeric_fit, coef = "TimeNumeric", number = Inf,
        sort.by = "none", adjust.method = "none"
    )
    result_rows <- match(
        rownames(numeric_table), numeric_result$results$time_effect$Feature
    )
    expect_equal(
        numeric_result$results$time_effect$F[result_rows],
        numeric_table$t^2
    )
    expect_equal(
        numeric_result$results$time_effect$P.Value[result_rows],
        numeric_table$P.Value
    )
})

test_that("DE_timecourse validates finite inputs and shared spline grids", {
    se_time <- example_obj
    SummarizedExperiment::colData(se_time)$Time[1] <- Inf
    expect_error(
        DE_timecourse(se_time, assay = 2),
        "finite"
    )

    se_assay <- example_obj
    values <- SummarizedExperiment::assay(se_assay, 2)
    values[1, 1] <- Inf
    SummarizedExperiment::assay(se_assay, 2) <- values
    expect_error(
        DE_timecourse(se_assay, assay = 2),
        "finite values or NA"
    )

    sample_data <- as.data.frame(
        SummarizedExperiment::colData(example_obj)
    )
    all_times <- sort(unique(sample_data$Time))
    shared_times <- all_times[seq_len(2L)]
    keep <- sample_data$Group == "IFNbeta" |
        (sample_data$Group == "IFNgamma" &
            sample_data$Time %in% shared_times)
    se_incomplete <- example_obj[, keep]
    expect_error(
        DE_timecourse(
            se_incomplete,
            assay = 2,
            group = c("IFNbeta", "IFNgamma"),
            group_mode = "joint_group",
            time_mode = "spline",
            spline_df = 2
        ),
        "shared"
    )
})

test_that("DE_timecourse uses unique time grids for splines", {
    set.seed(1)
    time <- c(rep(0, 10), 1, 2, 3, 4)
    values <- matrix(
        stats::rnorm(nrow(example_obj) * length(time)),
        nrow = nrow(example_obj),
        dimnames = list(
            rownames(example_obj),
            paste0("sample_", seq_along(time))
        )
    )
    se <- SummarizedExperiment::SummarizedExperiment(
        assays = list(norm = values),
        colData = data.frame(
            Group = rep("A", length(time)),
            Time = time
        )
    )

    result <- suppressMessages(DE_timecourse(
        se, assay = "norm", group = "A",
        group_mode = "per_group", time_mode = "spline", spline_df = 3
    ))

    expect_result_statistics(
        result$results$time_effect, require_finite = TRUE
    )
})

test_that("DE_timecourse validates p-value adjustment methods", {
    expect_error(
        DE_timecourse(
            example_obj, assay = 2, pAdjustMethod = "not_a_method"
        ),
        "'pAdjustMethod' must be one of"
    )
})

test_that("DE_timecourse retains features missing at one time point", {
    se <- example_obj[, colData(example_obj)$Group == "IFNbeta"]
    feature <- rownames(se)[2]
    for (time_mode in c("factor", "numeric", "spline")) {
        values <- SummarizedExperiment::assay(se, 2)
        values[feature, colData(se)$Time == 2] <- NA_real_
        SummarizedExperiment::assay(se, 2) <- values

        result <- DE_timecourse(
            se, assay = 2, group = "IFNbeta",
            group_mode = "per_group", time_mode = time_mode,
            spline_df = 2
        )
        feature_result <- result$results$time_effect[
            result$results$time_effect$Feature == feature, , drop = FALSE
        ]

        expect_equal(nrow(feature_result), 1L)
        expect_true(nzchar(feature_result$TimePointsTested))
        if (time_mode == "factor") {
            expect_equal(feature_result$ReferenceTime, 0)
        } else {
            expect_true(is.na(feature_result$ReferenceTime))
        }
        expect_true(all(is.finite(feature_result$P.Value)))
    }
})

test_that("DE_timecourse reports pattern-specific factor reference times", {
    se <- example_obj[, colData(example_obj)$Group == "IFNbeta"]
    feature <- rownames(se)[2]
    values <- SummarizedExperiment::assay(se, 2)
    values[feature, colData(se)$Time == 0] <- NA_real_
    values[feature, colData(se)$Time %in% c(6, 8, 24)] <- NA_real_
    SummarizedExperiment::assay(se, 2) <- values

    result <- DE_timecourse(
        se, assay = 2, group = "IFNbeta",
        group_mode = "per_group", time_mode = "factor"
    )
    feature_result <- result$results$time_effect[
        result$results$time_effect$Feature == feature, , drop = FALSE
    ]

    expect_equal(nrow(feature_result), 1L)
    expect_equal(feature_result$TimePointsTested, "2,4")
    expect_equal(feature_result$ReferenceTime, 2)
    expect_true(is.finite(feature_result$P.Value))
})

test_that("DE_timecourse retains unestimable features with NA results", {
    se <- example_obj[, colData(example_obj)$Group == "IFNbeta"]
    values <- SummarizedExperiment::assay(se, 2)
    feature <- rownames(se)[3]
    values[feature, ] <- NA_real_
    SummarizedExperiment::assay(se, 2) <- values

    result <- DE_timecourse(
        se, assay = 2, group = "IFNbeta",
        group_mode = "per_group", time_mode = "factor"
    )
    feature_result <- result$results$time_effect[
        result$results$time_effect$Feature == feature, , drop = FALSE
    ]

    expect_equal(nrow(feature_result), 1L)
    expect_true(is.na(feature_result$P.Value))
    expect_true(is.na(feature_result$adj.P.Val))
    expect_true(all(c("F", "P.Value", "adj.P.Val") %in%
        names(result$results$time_effect)))
})

test_that("DE_timecourse retains joint features absent from a group", {
    se <- example_obj[
        , colData(example_obj)$Group %in% c("IFNbeta", "IFNgamma")
    ]
    feature <- rownames(se)[4]
    values <- SummarizedExperiment::assay(se, 2)
    values[feature, colData(se)$Group == "IFNgamma"] <- NA_real_
    SummarizedExperiment::assay(se, 2) <- values

    result <- DE_timecourse(
        se, assay = 2,
        group = c("IFNbeta", "IFNgamma"),
        group_mode = "joint_group", reference_group = "IFNgamma"
    )

    for (table in result$results) {
        feature_result <- table[table$Feature == feature, , drop = FALSE]
        expect_equal(nrow(feature_result), 1L)
        expect_true(is.na(feature_result$P.Value))
    }
    reference_result <- result$results$reference_time_effect
    expect_equal(
        reference_result$Group[reference_result$Feature == feature],
        "IFNgamma"
    )
})

test_that("DE_timecourse fixes treatment contrasts explicitly", {
    old_contrasts <- getOption("contrasts")
    on.exit(options(contrasts = old_contrasts), add = TRUE)
    options(contrasts = c("contr.sum", "contr.poly"))

    result <- DE_timecourse(
        example_obj, assay = 2,
        group = c("IFNbeta", "IFNgamma"),
        group_mode = "joint_group", reference_group = "IFNgamma"
    )

    expect_true(all(is.finite(
        result$results$reference_time_effect$P.Value
    )))
    expect_true(any(grepl("GroupIFNbeta", colnames(result$design$design))))
})

test_that("DE_timecourse reports missing sample replication", {
    se <- example_obj[, colData(example_obj)$Group == "IFNbeta"]
    keep <- !duplicated(colData(se)$Time)
    se <- se[, keep]

    expect_message(
        expect_error(
            DE_timecourse(
                se, assay = 2, group = "IFNbeta",
                group_mode = "per_group", time_mode = "factor"
            ),
            "no residual degrees"
        ),
        "No replicate samples detected"
    )
})

test_that("DE_timecourse reports the default joint reference group", {
    result <- NULL
    expect_message(
        result <- DE_timecourse(
            example_obj, assay = 2,
            group = c("IFNbeta", "IFNgamma"),
            group_mode = "joint_group", time_mode = "numeric"
        ),
        "Using 'IFNbeta' as the reference group"
    )
    expect_equal(result$metadata$reference_group, "IFNbeta")
})

test_that("DE_timecourse returns estimable joint Group-Time tests", {
    result <- DE_timecourse(
        example_obj, assay = 2, group = c("IFNbeta", "IFNgamma"),
        group_mode = "joint_group", time_mode = "factor",
        reference_group = "IFNgamma"
    )

    expect_true(all(c("group_effect", "reference_time_effect",
        "global_time_effect", "group_time_interaction",
        "global_group_effect") %in%
        names(result$results)))
    expect_false("time_effect" %in% names(result$results))
    expect_equal(result$metadata$reference_group, "IFNgamma")
    expect_true(all(result$results$group_effect$ReferenceTime == 0))
    expect_result_statistics(result$results$group_effect)
    expect_result_statistics(result$results$global_group_effect)
    expect_true("logFC_GroupIFNbeta_at_reference_time" %in%
        names(result$results$group_effect))
    expect_true(all(c(
        "GroupIFNbeta", "delta_logFC_GroupIFNbeta_TimeFactor2"
    ) %in%
        names(result$results$global_group_effect)))
    expect_true("logFC_TimeFactor2" %in%
        names(result$results$reference_time_effect))
    expect_true("delta_logFC_GroupIFNbeta_TimeFactor2" %in%
        names(result$results$group_time_interaction))
    expect_true("GroupIFNbeta:TimeFactor2" %in%
        colnames(result$design$design))
    expect_true("GroupIFNbeta:TimeFactor2" %in%
        colnames(result$fits$fit$coefficients))
    expect_equal(nrow(result$results$group_time_interaction), 100L)
    expect_true(any(grepl("TimeFactor", names(
        result$results$group_time_interaction
    ))))
    expect_false(any(grepl(
        ":", names(result$results$group_time_interaction), fixed = TRUE
    )))
    incomplete_result <- DE_timecourse(
        example_obj, assay = 2,
        group_mode = "joint_group", time_mode = "factor"
    )
    expect_true(is.list(incomplete_result))
    expect_true(all(nzchar(
        incomplete_result$results$group_time_interaction$TimePointsTested
    )))
})

test_that("DE_timecourse names single factor contrasts by effect", {
    keep <- colData(example_obj)$Group %in% c("IFNbeta", "IFNgamma") &
        colData(example_obj)$Time %in% c(0, 2)
    result <- DE_timecourse(
        example_obj[, keep], assay = 2,
        group = c("IFNbeta", "IFNgamma"),
        group_mode = "joint_group", time_mode = "factor",
        reference_group = "IFNgamma"
    )

    expect_true("logFC_TimeFactor2" %in%
        names(result$results$reference_time_effect))
    expect_true("logFC_GroupIFNbeta_at_reference_time" %in%
        names(result$results$group_effect))
    expect_true("delta_logFC_GroupIFNbeta_TimeFactor2" %in%
        names(result$results$group_time_interaction))
    expect_true("GroupIFNbeta:TimeFactor2" %in%
        colnames(result$design$design))
    expect_true("GroupIFNbeta:TimeFactor2" %in%
        colnames(result$fits$fit$coefficients))
})

test_that("DE_timecourse tests joint group effects across time modes", {
    for (current_time_mode in c("factor", "numeric", "spline")) {
        result <- DE_timecourse(
            example_obj, assay = 2,
            group = c("IFNbeta", "IFNgamma"),
            group_mode = "joint_group",
            time_mode = current_time_mode,
            spline_df = 2,
            reference_group = "IFNgamma"
        )
        group_result <- result$results$group_effect
        global_group_result <- result$results$global_group_effect

        expect_equal(nrow(group_result), nrow(example_obj))
        expect_equal(unique(group_result$Test), "group_effect")
        expect_true(all(group_result$ReferenceTime == 0))
        expect_result_statistics(group_result, require_finite = TRUE)
        expect_result_statistics(
            global_group_result, require_finite = TRUE
        )
        if (current_time_mode == "spline") {
            spline_names <- c("SplineBasis_1", "SplineBasis_2")
            interaction_names <- c(
                "delta_SplineBasis_GroupIFNbeta_1",
                "delta_SplineBasis_GroupIFNbeta_2"
            )
            expect_true(all(spline_names %in% names(
                result$results$reference_time_effect
            )))
            expect_true(all(interaction_names %in% names(
                result$results$group_time_interaction
            )))
            expect_true(all(c(spline_names, interaction_names) %in% names(
                result$results$global_time_effect
            )))
            expect_true(all(c("GroupIFNbeta", interaction_names) %in%
                names(global_group_result)))
        }
    }

    numeric_result <- DE_timecourse(
        example_obj, assay = 2,
        group = c("IFNbeta", "IFNgamma"),
        group_mode = "joint_group", time_mode = "numeric",
        reference_group = "IFNgamma"
    )
    numeric_fit <- numeric_result$fits$fit
    group_table <- limma::topTable(
        numeric_fit, coef = "GroupIFNbeta", number = Inf,
        sort.by = "none", adjust.method = "none"
    )
    result_rows <- match(
        rownames(group_table),
        numeric_result$results$group_effect$Feature
    )
    group_effect_estimates <- numeric_result$results$group_effect[[
        "logFC_GroupIFNbeta_at_reference_time"
    ]]
    expect_equal(
        group_effect_estimates[result_rows],
        group_table$logFC
    )
    expect_equal(
        numeric_result$results$group_effect$F[result_rows],
        group_table$t^2
    )
    expect_true("logFC_TimeSlope" %in%
        names(numeric_result$results$reference_time_effect))
    expect_true("logFC_TimeSlope" %in%
        names(numeric_result$results$global_time_effect))
    expect_true("delta_logFC_GroupIFNbeta_TimeSlope" %in%
        names(numeric_result$results$group_time_interaction))
})

test_that("DE_timecourse uses omnibus group tests for multiple groups", {
    result <- DE_timecourse(
        example_obj, assay = 2,
        group = c("IFNbeta", "IFNgamma", "LPS"),
        group_mode = "joint_group", time_mode = "numeric",
        reference_group = "IFNgamma"
    )

    expect_true(all(c(
        "logFC_GroupIFNbeta_at_reference_time",
        "logFC_GroupLPS_at_reference_time"
    ) %in% names(result$results$group_effect)))
    expect_true(all(c(
        "GroupIFNbeta", "GroupLPS",
        "delta_logFC_GroupIFNbeta_TimeSlope",
        "delta_logFC_GroupLPS_TimeSlope"
    ) %in% names(result$results$global_group_effect)))
    expect_result_statistics(
        result$results$group_effect, require_finite = TRUE
    )
    expect_result_statistics(
        result$results$global_group_effect, require_finite = TRUE
    )
})

test_that("DE_timecourse models repeated Subject measurements", {
    se <- example_obj[, colData(example_obj)$Group == "IFNbeta"]
    SummarizedExperiment::colData(se)$Subject <- rep(c("S1", "S2"), 6)

    result <- DE_timecourse(
        se, assay = 2, group = "IFNbeta",
        group_mode = "per_group", time_mode = "factor"
    )

    expect_true(is.finite(result$fits[[1]]$correlation))

    SummarizedExperiment::colData(se)$Subject <- "S1"
    expect_error(
        DE_timecourse(
            se, assay = 2, group = "IFNbeta",
            group_mode = "per_group", time_mode = "factor"
        ),
        "at least two unique subjects"
    )
})
