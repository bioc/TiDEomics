# ---- Shared Trendy result ----
data("example_obj")
example_obj <- normalise_to_start(example_obj)
example_obj_list <- split_groups(example_obj)
example_obj_merged_list <- merge_replicates(example_obj_list)
example_obj_merged_imp_list <- impute_groups(example_obj_merged_list)
test_features <- utils::head(rownames(example_obj_merged_imp_list[[1]]), 10)

trendy_res <- run_Trendy(example_obj_merged_imp_list,
    feature = test_features,
    minExp = 0.5, maxK = 1, minNumInSeg = 2,
    meanCut = 0, NCores = 1)

trendy_summary <- summarise_Trendy(trendy_res)

# ---- Tests ----

test_that("run_Trendy works", {
    expect_true(is.list(trendy_res))
    expect_true(all(names(trendy_res) %in% names(example_obj_merged_imp_list)))
})

test_that("summarise_Trendy produces expected columns", {
    expect_s3_class(trendy_summary, "data.frame")
    expect_true("Feature" %in% colnames(trendy_summary))
    expect_true("Group" %in% colnames(trendy_summary))
})

test_that("extract_segment_trends works", {
    tl <- extract_segment_trends(trendy_summary)
    expect_true(is.list(tl))
    groups <- unique(trendy_summary$Group)
    expect_true(all(groups %in% names(tl)))
})

test_that("plot_breakpoints returns ggplot", {
    p <- plot_breakpoints(trendy_res)
    expect_s3_class(p, "ggplot")
})

test_that("run_Trendy validates feature argument", {
    expect_error(
        run_Trendy(example_obj_merged_imp_list, feature = "nonexistent",
            maxK = 1, minNumInSeg = 2, NCores = 1),
        "None of the specified features"
    )
})

test_that("run_Trendy validates minNumInSeg", {
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
    merged_list_prop <- calc_feature_property(example_obj_merged_list,
        threshold = 0)
    imp_list <- impute_groups(merged_list_prop)
    sub_feat <- rownames(imp_list[[1]])[1:10]
    imp_list_sub <- lapply(imp_list, function(x) {
        sub <- x[sub_feat, , drop = FALSE]
        SummarizedExperiment::rowData(sub)[["Feature"]] <- rownames(sub)
        sub
    })
    expect_no_error(
        run_Trendy(imp_list_sub, feature = NULL, maxK = 1, minNumInSeg = 2, NCores = 1)
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
    feat <- rownames(example_obj_merged_imp_list[[1]])[1]
    expect_error(
        run_Trendy(example_obj_merged_imp_list, feature = feat,
            group = "nonexistent"),
        "not found"
    )
})
