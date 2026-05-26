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
    }
})
