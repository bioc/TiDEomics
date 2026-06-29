# Tests for summarise_module_metrics

test_that("summarise_module_metrics returns expected columns", {
    data("example_net")
    res <- summarise_module_metrics(example_net)

    expect_s3_class(res, "data.frame")
    expected_cols <- c("Module", "Size", "Proportion",
        "MeanKME", "MeanKME2", "MedianKME", "SDKME", "MinKME", "MaxKME")
    expect_true(all(expected_cols %in% colnames(res)))
    expect_true(nrow(res) > 0)
})

test_that("summarise_module_metrics excludes grey module", {
    data("example_net")
    res <- summarise_module_metrics(example_net)

    # Grey module (M0) should not appear
    expect_false(any(res$Module %in% c("M0", "grey", "gray")))
})

test_that("summarise_module_metrics modules ordered by size", {
    data("example_net")
    res <- summarise_module_metrics(example_net)

    # Module sizes should be decreasing (or at least first is largest)
    expect_true(res$Size[1] >= res$Size[nrow(res)])
})

test_that("summarise_module_metrics validates input", {
    expect_error(summarise_module_metrics(NULL), "must be a list")
    expect_error(summarise_module_metrics(list()), "colors")
})

test_that("summarise_module_metrics kME values are in [-1, 1]", {
    data("example_net")
    res <- summarise_module_metrics(example_net)

    # Signed kME can be negative; should be in [-1, 1]
    expect_true(all(res$MeanKME >= -1 & res$MeanKME <= 1, na.rm = TRUE))
    expect_true(all(res$MinKME >= -1, na.rm = TRUE))
    expect_true(all(res$MaxKME <= 1, na.rm = TRUE))
})

test_that("summarise_module_metrics Proportion sums to <= 1", {
    data("example_net")
    res <- summarise_module_metrics(example_net)

    expect_true(sum(res$Proportion) <= 1)
})

test_that("summarise_module_metrics errors without input_data", {
    data("example_net")
    net_no_input <- example_net
    net_no_input$input_data <- NULL
    expect_error(
        summarise_module_metrics(net_no_input),
        "must contain 'input_data'"
    )
})

test_that("summarise_module_metrics errors when no non-grey modules", {
    data("example_net")
    net_all_grey <- example_net
    net_all_grey$colors <- rep("grey", length(example_net$colors))
    expect_error(
        summarise_module_metrics(net_all_grey),
        "No non-grey modules"
    )
})

test_that("summarise_module_metrics works with unsigned network", {
    data("example_net")
    net_unsigned <- example_net
    net_unsigned$input_data <- example_net$input_data
    if (is.null(net_unsigned$parameters))
        net_unsigned$parameters <- list()
    net_unsigned$parameters$networkType <- "unsigned"

    res <- summarise_module_metrics(net_unsigned)
    expect_s3_class(res, "data.frame")
    # Unsigned kME values should be >= 0
    expect_true(all(res$MinKME >= 0, na.rm = TRUE))
})
