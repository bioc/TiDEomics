# Tests for plot_segments

test_that("plot_segments validates feature argument", {
    data("example_res_list")
    trendy_summary <- summarise_Trendy(example_res_list)

    data("example_obj")
    example_obj <- normalise_to_start(example_obj)
    se_list <- split_groups(example_obj)

    expect_error(
        plot_segments(se_list, example_res_list, feature = 123),
        "must be a character"
    )
    expect_error(
        plot_segments(se_list, example_res_list, feature = "nonexistent"),
        "merged"
    )
})

test_that("plot_segments validates se_obj_imp_list", {
    data("example_res_list")
    expect_error(
        plot_segments(list(), example_res_list, feature = "Actb"),
        "non-empty list"
    )
    expect_error(
        plot_segments(NULL, example_res_list, feature = "Actb"),
        "must be a non-empty list"
    )
})

test_that("plot_segments validates res_list", {
    data("example_obj")
    example_obj <- normalise_to_start(example_obj)
    se_list <- split_groups(example_obj)

    expect_error(
        plot_segments(se_list, list(), feature = "Actb"),
        "merged"
    )
})

test_that("plot_segments returns ggplot for valid feature", {
    data("example_res_list")
    trendy_summary <- summarise_Trendy(example_res_list)

    data("example_obj")
    example_obj <- normalise_to_start(example_obj)
    se_list <- split_groups(example_obj)
    se_list <- merge_replicates(se_list)

    # Pick a feature present in the data
    feat <- rownames(se_list[[1]])[1]

    expect_no_error(
        plot_segments(
            se_obj_imp_list = se_list,
            res_list = example_res_list,
            feature = feat
        )
    )
})

test_that("plot_segments runs with specify_range", {
    data("example_res_list")
    data("example_obj")
    example_obj <- normalise_to_start(example_obj)
    se_list <- split_groups(example_obj)
    se_list <- merge_replicates(se_list)

    feat <- rownames(se_list[[1]])[1]

    # Trendy::plotFeature() dispatches on feature and renders via base graphics
    expect_no_error(
        plot_segments(
            se_obj_imp_list = se_list,
            res_list = example_res_list,
            feature = feat
        )
    )
})

# ---- plot_segments edge cases ----

test_that("plot_segments handles multiple features", {
    data("example_obj")
    example_obj <- normalise_to_start(example_obj)
    example_obj_list <- split_groups(example_obj)
    example_obj_merged_list <- merge_replicates(example_obj_list)
    example_obj_merged_imp_list <- impute_groups(example_obj_merged_list)
    feat <- utils::head(rownames(example_obj_merged_imp_list[[1]]), 3)
    res_list <- run_Trendy(example_obj_merged_imp_list,
        feature = feat, minExp = 0.5, maxK = 1, minNumInSeg = 2,
        meanCut = 0, NCores = 1)
    # plot_segments uses base R graphics and does not return a value
    expect_no_error(
        plot_segments(example_obj_merged_imp_list, res_list,
            feature = feat, nrow = 1)
    )
})

test_that("plot_segments errors on invalid group", {
    data("example_obj")
    example_obj <- normalise_to_start(example_obj)
    se_list <- split_groups(example_obj)
    se_list <- merge_replicates(se_list)
    data("example_res_list")
    feat <- rownames(se_list[[1]])[1]
    expect_error(
        plot_segments(se_list, example_res_list,
            feature = feat, group = "nonexistent"),
        "not found")
})
