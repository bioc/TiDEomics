# Tests for .module_labels

test_that(".module_labels handles numeric input with prefix", {
    res <- .module_labels(c(0, 1, 2, 3), prefix = TRUE)
    expect_s3_class(res, "factor")
    expect_equal(levels(res), c("M0", "M1", "M2", "M3"))
})

test_that(".module_labels handles numeric input without prefix", {
    res <- .module_labels(c(0, 1, 2, 3), prefix = FALSE)
    expect_s3_class(res, "factor")
    expect_equal(levels(res), c("0", "1", "2", "3"))
})

test_that(".module_labels handles digit string input", {
    res <- .module_labels(c("1", "2", "10"), prefix = TRUE)
    expect_s3_class(res, "factor")
    expect_equal(levels(res), c("M1", "M2", "M10"))
})

test_that(".module_labels handles M-prefixed input", {
    res <- .module_labels(c("M2", "M1", "M3"), prefix = TRUE)
    expect_equal(levels(res), c("M1", "M2", "M3"))

    res2 <- .module_labels(c("M2", "M1", "M3"), prefix = FALSE)
    expect_equal(levels(res2), c("1", "2", "3"))
})

test_that(".module_labels preserves WGCNA colour name order", {
    colours <- c("turquoise", "blue", "brown")
    res <- .module_labels(colours, prefix = TRUE)
    expect_equal(levels(res), c("turquoise", "blue", "brown"))

    # Different order should be preserved
    colours2 <- c("blue", "turquoise", "brown")
    res2 <- .module_labels(colours2, prefix = TRUE)
    expect_equal(levels(res2), c("blue", "turquoise", "brown"))
})

test_that(".module_labels returns correct numeric order", {
    # Even with unordered input
    res <- .module_labels(c(3, 1, 5, 2), prefix = TRUE)
    expect_equal(levels(res), c("M1", "M2", "M3", "M5"))
})

test_that("WGCNA_module preserves feature names as row.names", {
    data("example_net")
    mod <- WGCNA_module(example_net, exclude_grey = FALSE)

    expect_true("Feature" %in% colnames(mod))
    # Feature should contain gene names, not row numbers
    expect_false(any(grepl("^[0-9]+$", mod$Feature[1:5])))
})

test_that(".module_labels handles single-element input", {
    res <- .module_labels(c(5L), prefix = TRUE)
    expect_equal(levels(res), "M5")
    expect_equal(length(res), 1)

    res2 <- .module_labels(c("turquoise"), prefix = TRUE)
    expect_equal(levels(res2), "turquoise")
})

test_that(".module_labels preserves value order for colour names", {
    # Verify factor values match input order, not alphabetical
    colours <- c("yellow", "blue", "turquoise")
    res <- .module_labels(colours, prefix = TRUE)
    expect_equal(as.character(res), colours)
})

test_that("WGCNA_module orders modules by size", {
    data("example_net")
    mod <- WGCNA_module(example_net, exclude_grey = TRUE)

    # Module factor levels should be in decreasing size order
    sizes <- table(mod$Module)
    level_sizes <- sizes[levels(mod$Module)]
    expect_true(all(diff(level_sizes) <= 0))
})
