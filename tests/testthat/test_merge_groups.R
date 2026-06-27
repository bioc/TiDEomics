# Tests for merge_groups

test_that("merge_groups preserves features and combines groups", {
    data("example_obj")
    example_obj <- normalise_to_start(example_obj)
    example_obj_list <- split_groups(example_obj)
    example_obj_merged_list <- merge_replicates(example_obj_list)
    example_obj_merged <- merge_groups(example_obj_merged_list)

    expect_s4_class(example_obj_merged, "SummarizedExperiment")
    # Should have all features
    expect_equal(nrow(example_obj_merged), nrow(example_obj))
    # colData should have Group column with all groups
    expect_true("Group" %in% colnames(colData(example_obj_merged)))
    expect_equal(
        sort(levels(colData(example_obj_merged)$Group)),
        sort(names(example_obj_list))
    )
    # Sample names should be prefixed with group names
    expect_true(all(grepl("_", colData(example_obj_merged)$Sample)))
})

test_that("merge_groups handles single group", {
    data("example_obj")
    example_obj <- normalise_to_start(example_obj)
    example_obj_list <- split_groups(example_obj)
    # Take just one group
    single_list <- example_obj_list[1]
    merged <- merge_groups(single_list)

    expect_s4_class(merged, "SummarizedExperiment")
    expect_equal(length(unique(colData(merged)$Group)), 1)
})
