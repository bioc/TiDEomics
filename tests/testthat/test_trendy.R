# Tests for Trendy pipeline (segmented regression)

test_that("run_Trendy works on imputed merged data", {
    data("example_obj")
    example_obj <- normalise_to_start(example_obj)
    example_obj_list <- split_groups(example_obj)
    example_obj_merged_list <- merge_replicates(example_obj_list)
    example_obj_merged_imp_list <- impute_groups(example_obj_merged_list)

    # Small subset for speed
    test_features <- utils::head(rownames(example_obj_merged_imp_list[[1]]), 10)

    res <- run_Trendy(example_obj_merged_imp_list,
        feature = test_features,
        minExp = 0.5, maxK = 1, minNumInSeg = 2,
        meanCut = 0, NCores = 1)

    expect_true(is.list(res))
    # Results may be a subset (groups with too few points are skipped)
    expect_true(all(names(res) %in% names(example_obj_merged_imp_list)))
})

test_that("summarise_Trendy produces expected columns", {
    data("example_obj")
    example_obj <- normalise_to_start(example_obj)
    example_obj_list <- split_groups(example_obj)
    example_obj_merged_list <- merge_replicates(example_obj_list)
    example_obj_merged_imp_list <- impute_groups(example_obj_merged_list)

    test_features <- utils::head(rownames(example_obj_merged_imp_list[[1]]), 10)
    res <- run_Trendy(example_obj_merged_imp_list,
        feature = test_features,
        minExp = 0.5, maxK = 1, minNumInSeg = 2,
        meanCut = 0, NCores = 1)

    ts <- summarise_Trendy(res)

    expect_s3_class(ts, "data.frame")
    expect_true("Feature" %in% colnames(ts))
    expect_true("Group" %in% colnames(ts))
})

test_that("extract_segment_trends works", {
    data("example_obj")
    example_obj <- normalise_to_start(example_obj)
    example_obj_list <- split_groups(example_obj)
    example_obj_merged_list <- merge_replicates(example_obj_list)
    example_obj_merged_imp_list <- impute_groups(example_obj_merged_list)

    test_features <- utils::head(rownames(example_obj_merged_imp_list[[1]]), 10)
    res <- run_Trendy(example_obj_merged_imp_list,
        feature = test_features,
        minExp = 0.5, maxK = 1, minNumInSeg = 2,
        meanCut = 0, NCores = 1)
    ts <- summarise_Trendy(res)

    tl <- extract_segment_trends(ts)
    expect_true(is.list(tl))
    # Should have one entry per group
    groups <- unique(ts$Group)
    expect_true(all(groups %in% names(tl)))
})

test_that("plot_breakpoints returns ggplot", {
    data("example_obj")
    example_obj <- normalise_to_start(example_obj)
    example_obj_list <- split_groups(example_obj)
    example_obj_merged_list <- merge_replicates(example_obj_list)
    example_obj_merged_imp_list <- impute_groups(example_obj_merged_list)

    test_features <- utils::head(rownames(example_obj_merged_imp_list[[1]]), 10)
    res <- run_Trendy(example_obj_merged_imp_list,
        feature = test_features,
        minExp = 0.5, maxK = 1, minNumInSeg = 2,
        meanCut = 0, NCores = 1)

    p <- plot_breakpoints(res)
    expect_s3_class(p, "ggplot")
})

test_that("run_Trendy validates feature argument", {
    data("example_obj")
    example_obj <- normalise_to_start(example_obj)
    example_obj_list <- split_groups(example_obj)
    example_obj_merged_list <- merge_replicates(example_obj_list)
    example_obj_merged_imp_list <- impute_groups(example_obj_merged_list)

    # feature must be present in data
    expect_error(
        run_Trendy(example_obj_merged_imp_list, feature = "nonexistent",
            maxK = 1, minNumInSeg = 2, NCores = 1),
        "None of the specified features"
    )
})

test_that("run_Trendy validates minNumInSeg", {
    data("example_obj")
    example_obj <- normalise_to_start(example_obj)
    example_obj_list <- split_groups(example_obj)
    example_obj_merged_list <- merge_replicates(example_obj_list)
    example_obj_merged_imp_list <- impute_groups(example_obj_merged_list)
    feat <- rownames(example_obj_merged_imp_list[[1]])[1]

    expect_error(
        run_Trendy(example_obj_merged_imp_list, feature = feat,
            maxK = 1, minNumInSeg = 0, NCores = 1),
        "positive"
    )
})

test_that("summarise_Trendy returns NULL when all groups skipped", {
    res <- list(GroupA = NULL, GroupB = NULL)
    result <- suppressMessages(summarise_Trendy(res))
    expect_null(result)
})

test_that("run_Trendy handles feature = NULL", {
    data("example_obj")
    example_obj <- normalise_to_start(example_obj)
    example_obj_list <- split_groups(example_obj)
    example_obj_merged_list <- merge_replicates(example_obj_list)
    example_obj_merged_list <- calc_feature_property(example_obj_merged_list,
        threshold = 0)
    example_obj_merged_imp_list <- impute_groups(example_obj_merged_list)
    expect_no_error(
        run_Trendy(example_obj_merged_imp_list,
            feature = NULL, maxK = 1, minNumInSeg = 2, NCores = 1)
    )
})

test_that("plot_breakpoints errors on empty res_list", {
    expect_error(plot_breakpoints(list()), "contains no groups")
})

test_that("plot_breakpoints works with a specific valid group name", {
    data("example_res_list")
    p <- plot_breakpoints(example_res_list, group = "IFNbeta")
    expect_s3_class(p, "ggplot")
})

test_that("plot_breakpoints errors on invalid group name", {
    data("example_res_list")
    expect_error(
        plot_breakpoints(example_res_list, group = "nonexistent"),
        "not found"
    )
})

test_that("plot_breakpoints skips NULL groups with message", {
    data("example_res_list")
    res_list <- list(
        IFNbeta = NULL,
        IFNgamma = example_res_list[["IFNgamma"]],
        LPS = example_res_list[["LPS"]]
    )
    expect_message(
        plot_breakpoints(res_list),
        "Skipping group.*IFNbeta"
    )
})

test_that("run_Trendy errors on invalid group name", {
    data("example_obj")
    example_obj <- normalise_to_start(example_obj)
    example_obj_list <- split_groups(example_obj)
    example_obj_merged_list <- merge_replicates(example_obj_list)
    example_obj_merged_imp_list <- impute_groups(example_obj_merged_list)
    feat <- rownames(example_obj_merged_imp_list[[1]])[1]

    expect_error(
        run_Trendy(example_obj_merged_imp_list, feature = feat,
            group = "nonexistent"),
        "not found"
    )
})

