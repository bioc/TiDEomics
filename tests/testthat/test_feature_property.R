# Tests for feature property pipeline

test_that("calc_feature_property works on merged data", {
    data("example")
    example_obj <- normalise_to_start(example_obj)
    example_obj_list <- split_groups(example_obj)
    example_obj_merged_list <- merge_replicates(example_obj_list)

    res <- calc_feature_property(example_obj_merged_list, threshold = 0)

    expect_true(is.list(res))
    for (nm in names(res)) {
        rd <- rowData(res[[nm]])
        expect_true("P_trend" %in% colnames(rd))
        expect_true("Max_FC" %in% colnames(rd))
        expect_true("Exp_ratio" %in% colnames(rd))
        expect_s4_class(res[[nm]], "SummarizedExperiment")
    }
})

test_that("calc_feature_property without threshold still works", {
    data("example")
    example_obj <- normalise_to_start(example_obj)
    example_obj_list <- split_groups(example_obj)
    example_obj_merged_list <- merge_replicates(example_obj_list)

    res <- calc_feature_property(example_obj_merged_list, threshold = NULL)
    expect_true(is.list(res))
})

test_that("summarise_feature_property works", {
    data("example")
    example_obj <- normalise_to_start(example_obj)
    example_obj_list <- split_groups(example_obj)
    example_obj_merged_list <- merge_replicates(example_obj_list)
    res <- calc_feature_property(example_obj_merged_list)

    prop_tb <- summarise_feature_property(res)
    expect_s3_class(prop_tb, "data.frame")
    expect_true("Feature" %in% colnames(prop_tb))
    expect_true("Exp_ratio" %in% colnames(prop_tb))
    expect_true("P_trend" %in% colnames(prop_tb))
    expect_true("Max_FC" %in% colnames(prop_tb))
    expect_true("Group" %in% colnames(prop_tb))
})

test_that("group_specific_features filters correctly", {
    data("example")
    example_obj <- normalise_to_start(example_obj)
    example_obj_list <- split_groups(example_obj)
    example_obj_merged_list <- merge_replicates(example_obj_list)
    res <- calc_feature_property(example_obj_merged_list, threshold = 0)
    prop_tb <- summarise_feature_property(res)

    gsf <- group_specific_features(prop_tb, groups = "untreated",
        genename = FALSE, GO = FALSE)

    expect_true(is.character(gsf))
    expect_true(length(gsf) > 0)
    # All returned features should be in the input
    expect_true(all(gsf %in% prop_tb$Feature))
})

test_that("group_specific_features respects filter_ratio", {
    data("example")
    example_obj <- normalise_to_start(example_obj)
    example_obj_list <- split_groups(example_obj)
    example_obj_merged_list <- merge_replicates(example_obj_list)
    res <- calc_feature_property(example_obj_merged_list, threshold = 0)
    prop_tb <- summarise_feature_property(res)

    # Strict filter: features must be present in >= 90% of time points
    gsf_strict <- group_specific_features(prop_tb, groups = "untreated",
        filter_ratio = 0.9, genename = FALSE, GO = FALSE)
    # Lenient filter: 10%
    gsf_lenient <- group_specific_features(prop_tb, groups = "untreated",
        filter_ratio = 0.1, genename = FALSE, GO = FALSE)

    # Strict should return <= lenient (use length, returns character vector)
    expect_true(length(gsf_strict) <= length(gsf_lenient))
})

test_that("group_specific_features with group_pct works", {
    data("example")
    example_obj <- normalise_to_start(example_obj)
    example_obj_list <- split_groups(example_obj)
    example_obj_merged_list <- merge_replicates(example_obj_list)
    res <- calc_feature_property(example_obj_merged_list, threshold = 0)
    prop_tb <- summarise_feature_property(res)

    # Require presence in only 1 out of 2 specified groups
    gsf <- group_specific_features(prop_tb,
        groups = c("untreated", "IFNbeta"),
        filter_ratio = 0.5, group_pct = 0.5,
        genename = FALSE, GO = FALSE)
    expect_true(is.character(gsf) || is.null(gsf))
})

test_that("group_specific_features with NULL groups uses all", {
    data("example")
    example_obj <- normalise_to_start(example_obj)
    example_obj_list <- split_groups(example_obj)
    example_obj_merged_list <- merge_replicates(example_obj_list)
    res <- calc_feature_property(example_obj_merged_list, threshold = 0)
    prop_tb <- summarise_feature_property(res)

    gsf <- group_specific_features(prop_tb, groups = NULL,
        genename = FALSE, GO = FALSE)
    expect_true(is.character(gsf) || is.null(gsf))
})

test_that("summarise_feature_property has expected columns", {
    data("example")
    example_obj <- normalise_to_start(example_obj)
    example_obj_list <- split_groups(example_obj)
    example_obj_merged_list <- merge_replicates(example_obj_list)
    res <- calc_feature_property(example_obj_merged_list, threshold = 0)
    prop_tb <- summarise_feature_property(res)

    expect_true("T_total" %in% colnames(prop_tb))
    expect_true("T_exp" %in% colnames(prop_tb))
    expect_true("Max_FC_time" %in% colnames(prop_tb))
})
