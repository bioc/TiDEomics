# Tests for summarise_module_pattern

test_that("summarise_module_pattern returns list of data.frames", {
    data("example_res_list")
    trendy_summary <- summarise_Trendy(example_res_list)

    data("example_net")
    example_module <- WGCNA_module(example_net)

    res <- summarise_module_pattern(example_module, trendy_summary)

    expect_true(is.list(res))
    expect_true(length(res) > 0)
    for (df in res) {
        expect_s3_class(df, "data.frame")
        expect_true(nrow(df) > 0)
        expect_true("Pattern" %in% colnames(df))
        expect_true("Count" %in% colnames(df))
    }
})

test_that("summarise_module_pattern validates input", {
    expect_error(summarise_module_pattern(NULL, data.frame()),
        "data.frame")
    data("example_res_list")
    trendy_summary <- summarise_Trendy(example_res_list)
    expect_error(summarise_module_pattern(data.frame(x = 1), trendy_summary),
        "Module")
})

test_that("summarise_module_pattern excludes grey module", {
    data("example_res_list")
    trendy_summary <- summarise_Trendy(example_res_list)
    data("example_net")
    example_module <- WGCNA_module(example_net)

    res <- summarise_module_pattern(example_module, trendy_summary)
    # Grey module should not be in result names
    expect_false(any(c("0", "grey", "gray", "M0") %in% names(res)))
})

test_that("summarise_module_pattern result columns are 'Pattern' and 'Count'", {
    data("example_res_list")
    trendy_summary <- summarise_Trendy(example_res_list)
    data("example_net")
    example_module <- WGCNA_module(example_net)

    res <- summarise_module_pattern(example_module, trendy_summary)
    for (df in res) {
        expect_equal(colnames(df), c("Pattern", "Count"))
    }
})
