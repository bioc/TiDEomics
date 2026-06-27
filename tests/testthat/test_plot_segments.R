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
        "not found"
    )
})

test_that("plot_segments validates se_obj_imp_list", {
    data("example_res_list")
    expect_error(
        plot_segments(list(), example_res_list, feature = "Actb"),
        "SummarizedExperiment"
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
        "run_Trendy"
    )
})

test_that("plot_segments returns ggplot for valid feature", {
    data("example_res_list")
    trendy_summary <- summarise_Trendy(example_res_list)

    data("example_obj")
    example_obj <- normalise_to_start(example_obj)
    se_list <- split_groups(example_obj)

    # Pick a feature present in the data
    feat <- rownames(se_list[[1]])[1]

    p <- plot_segments(
        se_obj_imp_list = se_list,
        res_list = example_res_list,
        feature = feat
    )

    expect_s3_class(p, "ggplot")
})

test_that("plot_segments runs with specify_range", {
    data("example_res_list")
    data("example_obj")
    example_obj <- normalise_to_start(example_obj)
    se_list <- split_groups(example_obj)

    feat <- rownames(se_list[[1]])[1]

    # with specify_range = FALSE
    p <- plot_segments(
        se_obj_imp_list = se_list,
        res_list = example_res_list,
        feature = feat,
        specify_range = FALSE
    )
    expect_s3_class(p, "ggplot")
})
