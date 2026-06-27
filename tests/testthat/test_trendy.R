# Tests for Trendy pipeline (segmented regression)

test_that("run_Trendy works on imputed merged data", {
    data("example_obj")
    example_obj <- normalise_to_start(example_obj)
    example_obj_list <- split_groups(example_obj)
    example_obj_merged_list <- merge_replicates(example_obj_list)
    example_obj_merged_imp_list <- impute_groups(example_obj_merged_list)

    # Small subset for speed
    test_features <- head(rownames(example_obj_merged_imp_list[[1]]), 10)

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

    test_features <- head(rownames(example_obj_merged_imp_list[[1]]), 10)
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

    test_features <- head(rownames(example_obj_merged_imp_list[[1]]), 10)
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

    test_features <- head(rownames(example_obj_merged_imp_list[[1]]), 10)
    res <- run_Trendy(example_obj_merged_imp_list,
        feature = test_features,
        minExp = 0.5, maxK = 1, minNumInSeg = 2,
        meanCut = 0, NCores = 1)

    p <- plot_breakpoints(res)
    expect_s3_class(p, "ggplot")
})

