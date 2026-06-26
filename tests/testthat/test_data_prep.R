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
    expect_all_true(t0_means == 0)
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
    data("example")
    res <- calc_mean_sd(example_obj)
    expect_type(res, "list")
    expect_true("norm" %in% names(res) || length(res) > 0)
    for (i in seq_along(res)) {
        expect_s3_class(res[[i]], "data.frame")
    }
})

test_that("calc_mean_sd data.frame has expected columns", {
    data("example")
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
