test_that("data loading works", {
    data(tutorial_data)
    data(tutorial_sample_info)
    data_obj <- create_input(tutorial_data, tutorial_sample_info)

    expect_s4_class(data_obj, "SummarizedExperiment")
})

test_that("color palette works", {
    custom_palette <- c(
        "untreated" = "#1b9e77", "IFNbeta" = "#d95f02",
        "IFNgamma" = "#7570b3", "LPS" = "#e7298a"
    )
    set_custom_palette(custom_palette)
    pal <- get_custom_palette(c("untreated", "IFNbeta", "IFNgamma", "LPS"))

    expect_equal(c("untreated", "IFNbeta", "IFNgamma", "LPS"), names(pal))

    # incomplete palette should be filled with default colors
    custom_palette <- c("untreated" = "#1b9e77")
    set_custom_palette(custom_palette)
    pal <- get_custom_palette(c("untreated", "IFNbeta", "IFNgamma", "LPS"))

    expect_equal(c("untreated", "IFNbeta", "IFNgamma", "LPS"), names(pal))
    expect_equal(pal["untreated"], custom_palette)

    default_pal <- scales::pal_hue()(3)
    names(default_pal) <- c("IFNbeta", "IFNgamma", "LPS")
    expect_equal(pal[c("IFNbeta", "IFNgamma", "LPS")], default_pal)
})

test_that("normalise to start works", {
    data(tutorial_data)
    data(tutorial_sample_info)
    data_obj <- create_input(tutorial_data, tutorial_sample_info)
    data_obj_norm <- normalise_to_start(data_obj)

    expect_s4_class(data_obj_norm, "SummarizedExperiment")
    # Confirm that means at Time 0 is 0
    t0_means <- assays(data_obj_norm[, data_obj_norm$Time == 0])[[2]] |>
        rowMeans()
    expect_true(all(abs(t0_means) < 1e-10))
})

test_that("split_groups works", {
    data(tutorial_data)
    data(tutorial_sample_info)
    data_obj <- create_input(tutorial_data, tutorial_sample_info)
    data_obj_list <- split_groups(data_obj)

    expect_true(is.list(data_obj_list))
    expect_equal(names(data_obj_list), as.character(unique(data_obj$Group)))

    for (group in names(data_obj_list)) {
        se_obj <- data_obj_list[[group]]
        expect_s4_class(se_obj, "SummarizedExperiment")
        expect_true(unique(se_obj$Group) == group)
    }
})

test_that("merge_replicates works", {
    data(tutorial_data)
    data(tutorial_sample_info)
    data_obj <- create_input(tutorial_data, tutorial_sample_info)
    data_obj_list <- split_groups(data_obj)
    data_obj_merged_list <- merge_replicates(data_obj_list)

    expect_true(is.list(data_obj_merged_list))
    expect_equal(names(data_obj_merged_list), names(data_obj_list))

    for (group in names(data_obj_merged_list)) {
        se_obj <- data_obj_merged_list[[group]]
        expect_s4_class(se_obj, "SummarizedExperiment")
        expect_true(unique(se_obj$Group) == group)
        expect_equal(ncol(se_obj), length(unique(se_obj$Time)))
    }
})

test_that("impute_groups works", {
    data(tutorial_data)
    data(tutorial_sample_info)
    data_obj <- create_input(tutorial_data, tutorial_sample_info)

    # Introduce some missing values
    assay(data_obj)[1:10, 1:5] <- NA

    data_obj_list <- split_groups(data_obj)
    data_obj_merged_list <- merge_replicates(data_obj_list)
    data_obj_imputed_list <- impute_groups(data_obj_merged_list)

    for (group in names(data_obj_imputed_list)) {
        data_obj_imputed <- data_obj_imputed_list[[group]]
        expect_s4_class(data_obj_imputed, "SummarizedExperiment")
        expect_false(any(is.na(assay(data_obj_imputed))))

        # Check that non-missing values are unchanged
        data_obj_original <- data_obj_merged_list[[group]]
        expect_equal(assay(data_obj_imputed)[!is.na(assay(data_obj_original))],
                assay(data_obj_original)[!is.na(assay(data_obj_original))])
    }
})

# ---- calc_mean_sd ----

test_that("calc_mean_sd returns list of data.frames", {
    data("example_obj")
    example_obj <- normalise_to_start(example_obj)
    res <- calc_mean_sd(example_obj)
    expect_type(res, "list")
    expect_true("norm" %in% names(res) || length(res) > 0)
    for (i in seq_along(res)) {
        expect_s3_class(res[[i]], "data.frame")
    }
})

test_that("calc_mean_sd data.frame has expected columns", {
    data("example_obj")
    example_obj <- normalise_to_start(example_obj)
    res <- calc_mean_sd(example_obj)
    tb <- res[[1]]
    expect_true(all(c("Group", "Time", "Feature") %in% colnames(tb)))
})

# ---- WGCNA_module with numeric labels ----

test_that("WGCNA_module works with numeric label net", {
    data("example_net")
    mod <- WGCNA_module(example_net)
    expect_s3_class(mod, "data.frame")
    expect_true(nrow(mod) > 0)
})

test_that("prepare_tide returns expected list structure", {
    data("tutorial_data")
    data("tutorial_sample_info")
    tide <- prepare_tide(tutorial_data, tutorial_sample_info,
        keep = "threshold", residual_threshold = 100)

    expect_type(tide, "list")
    expected <- c("se", "se_filtered", "merged_list", "merged_list_filtered",
                  "merged_se", "merged_se_filtered", "variance",
                  "feature_property", "filter_summary",
                  "group_specific_filter",
                  "DE", "enrichment", "WGCNA")
    expect_true(all(expected %in% names(tide)))
    expect_true(is.list(tide$DE))
    expect_true(is.list(tide$enrichment))
    expect_true(is.list(tide$WGCNA))
})

test_that("prepare_tide validates feature_threshold", {
    data("tutorial_data")
    data("tutorial_sample_info")
    expect_error(
        prepare_tide(tutorial_data, tutorial_sample_info,
            keep = "threshold", residual_threshold = 100,
            feature_threshold = "low"),
        "must be numeric"
    )
})

test_that("prepare_tide keep = 'below_quantile' validates residual_threshold", {
    data("tutorial_data")
    data("tutorial_sample_info")
    expect_error(
        prepare_tide(tutorial_data, tutorial_sample_info,
            keep = "below_quantile", residual_threshold = 2)
    )
})

test_that("prepare_tide keep = 'threshold' validates residual_threshold range", {
    data("tutorial_data")
    data("tutorial_sample_info")
    expect_error(
        prepare_tide(tutorial_data, tutorial_sample_info,
            keep = "threshold", residual_threshold = 150),
        "must be 0-100 when keep = 'threshold'"
    )
})

test_that("prepare_tide auto-assigns residual_threshold when NULL", {
    data("tutorial_data")
    data("tutorial_sample_info")
    res_bq <- prepare_tide(tutorial_data, tutorial_sample_info,
        keep = "below_quantile", residual_threshold = NULL,
        min_groups = 2, n_keep = 100)
    expect_match(res_bq$filter_summary$threshold[3], "below_quantile = ")

    res_th <- prepare_tide(tutorial_data, tutorial_sample_info,
        keep = "threshold", residual_threshold = NULL,
        min_groups = 2, n_keep = 100)
    expect_match(res_th$filter_summary$threshold[3], "skipped")
})

test_that("prepare_tide errors when min_groups > n_groups", {
    data("tutorial_data")
    data("tutorial_sample_info")
    expect_error(
        prepare_tide(tutorial_data, tutorial_sample_info,
            keep = "threshold", residual_threshold = 100,
            min_groups = 100),
        "min_groups.*exceeds"
    )
})

test_that("prepare_tide keep = 'below_quantile' runs end-to-end", {
    data("tutorial_data")
    data("tutorial_sample_info")
    res <- prepare_tide(tutorial_data, tutorial_sample_info,
        keep = "below_quantile", residual_threshold = 0.9,
        min_groups = 2, n_keep = 100)
    expect_type(res, "list")
    expect_true("se" %in% names(res))
})

test_that("prepare_tide keep = 'top_n' runs end-to-end", {
    data("tutorial_data")
    data("tutorial_sample_info")
    res <- prepare_tide(tutorial_data, tutorial_sample_info,
        keep = "top_n", n_keep = 50,
        min_groups = 2)
    expect_type(res, "list")
    expect_true("se" %in% names(res))
})

# ---- create_input edge cases ----

test_that("create_input subject_col validates existence", {
    data("tutorial_data")
    data("tutorial_sample_info")

    expect_error(
        create_input(tutorial_data, tutorial_sample_info,
            subject_col = "nonexistent"),
        "not found in sample_ann"
    )
})

test_that("create_input replicate_col validates existence", {
    data("tutorial_data")
    data("tutorial_sample_info")

    expect_error(
        create_input(tutorial_data, tutorial_sample_info,
            replicate_col = "nonexistent"),
        "not found in sample_ann"
    )
})

test_that("create_input batch_col validates existence", {
    data("tutorial_data")
    data("tutorial_sample_info")

    expect_error(
        create_input(tutorial_data, tutorial_sample_info,
            batch_col = "nonexistent"),
        "not found in sample_ann"
    )
})

test_that("create_input errors on missing Feature column", {
    expect_error(
        create_input(data.frame(x = 1:3), data.frame(Sample = "S1",
            Group = "A", Time = 0)),
        "must be named .Feature."
    )
})

test_that("create_input errors on missing required sample_ann columns", {
    data("tutorial_data")
    expect_error(
        create_input(tutorial_data, data.frame(Sample = "S1")),
        "missing one or more required columns"
    )
})

test_that("split_groups errors when Group column missing", {
    data("example_obj")
    se <- example_obj
    SummarizedExperiment::colData(se)$Group <- NULL
    expect_error(
        split_groups(se),
        "must have a 'Group' column"
    )
})

# ---- impute_groups with subject ----


# ---- create_input edge cases ----

test_that("create_input errors when Replicate and user replicate_col conflict", {
    data("tutorial_data")
    ann <- data.frame(Sample = c("S1", "S2"),
        Group = c("A", "A"), Time = c(0, 2),
        Replicate = c("R1", "R2"))
    ann$Rep <- c("A", "B")
    expect_error(
        create_input(data.frame(Feature = c("A", "B"), S1 = 1:2, S2 = 3:4,
            check.names = FALSE), ann, replicate_col = "Rep"),
        "must not be duplicated")
})

test_that("create_input errors when subjects cross groups", {
    ann <- data.frame(Sample = c("S1", "S2", "S3", "S4"),
        Group = c("A", "A", "B", "B"), Time = c(0, 2, 0, 2),
        Replicate = c("R1", "R2", "R1", "R2"))
    ann$Subject <- c("Sub1", "Sub2", "Sub1", "Sub3")
    d <- data.frame(Feature = c("A", "B"), S1 = 1:2, S2 = 3:4,
        S3 = 5:6, S4 = 7:8, check.names = FALSE)
    expect_error(
        create_input(d, ann, subject_col = "Subject"),
        "appear in multiple groups")
})

test_that("create_input errors when Time cannot be converted to numeric", {
    ann <- data.frame(Sample = c("S1", "S2"),
        Group = c("A", "A"), Time = c("abc", "def"))
    d <- data.frame(Feature = "A", S1 = 1, S2 = 2, check.names = FALSE)
    expect_error(
        create_input(d, ann),
        "numeric")
})

test_that("create_input errors on duplicate Replicate labels", {
    d <- data.frame(Feature = c("A", "B"), S1 = 1:2, S2 = 3:4,
        check.names = FALSE)
    ann <- data.frame(Sample = c("S1", "S2"),
        Group = c("A", "A"), Time = c(0, 0),
        Replicate = c("R1", "R1"))
    expect_error(create_input(d, ann), "Duplicate")
})

test_that("create_input errors when samples mismatch annotation", {
    ann <- data.frame(Sample = c("S1", "S2", "Missing"),
        Group = c("A", "A", "A"), Time = c(0, 0, 2))
    d <- data.frame(Feature = c("A", "B"), S1 = 1:2, S2 = 3:4,
        check.names = FALSE)
    expect_error(create_input(d, ann),
        "Some samples in sample annotation are not present")
})

test_that("create_input errors when data has extra samples", {
    ann <- data.frame(Sample = "S1", Group = "A", Time = 0)
    d <- data.frame(Feature = "A", S1 = 1, Extra = 2,
        check.names = FALSE)
    expect_error(create_input(d, ann),
        "Some samples in the data are not present")
})
