# Tests for extract_hubs

test_that("extract_hubs returns top features per module", {
    data("example_net")
    hubs <- extract_hubs(example_net, top_n = 3)

    expect_type(hubs, "character")
    expect_true(length(hubs) > 0)
    # Number of hubs = top_n * n_modules
    modules <- setdiff(unique(example_net$colors), c("grey", "0", 0))
    expect_equal(length(hubs), length(modules) * 3)
})

test_that("extract_hubs excludes grey module by default", {
    data("example_net")
    hubs <- extract_hubs(example_net, top_n = 2, exclude_grey = TRUE)

    # Grey/0 features should not appear
    grey_features <- names(example_net$colors)[
        example_net$colors %in% c("grey", "0", 0)]
    expect_false(any(grey_features %in% hubs))
})

test_that("extract_hubs includes grey when exclude_grey = FALSE", {
    data("example_net")
    hubs <- extract_hubs(example_net, top_n = 2, exclude_grey = FALSE)

    grey_features <- names(example_net$colors)[
        example_net$colors %in% c("grey", "0", 0)]
    if (length(grey_features) > 0) {
        expect_true(any(grey_features %in% hubs))
    }
})

test_that("extract_hubs validates arguments", {
    expect_error(extract_hubs(NULL), "must be a list")
    expect_error(extract_hubs(1:3), "must be a list")
    data("example_net")
    expect_error(extract_hubs(example_net, top_n = 0), "positive")
    expect_error(extract_hubs(example_net, top_n = -1), "positive")
    expect_error(extract_hubs(example_net, exclude_grey = "yes"), "must be TRUE or FALSE")
})

test_that("extract_hubs errors on missing required fields", {
    expect_error(
        extract_hubs(list(parameters = list(networkType = "signed"))),
        "input_data"
    )
    expect_error(
        extract_hubs(list(input_data = matrix(1:4, 2, 2))),
        "parameters"
    )
})

test_that("extract_hubs handles top_n larger than module size", {
    data("example_net")
    hubs <- extract_hubs(example_net, top_n = 1000)
    expect_type(hubs, "character")
    # Should return all available features (no error)
    expect_true(length(hubs) <= length(example_net$colors))
})

test_that("extract_hubs handles unsigned networks", {
    # Build a minimal unsigned network mock to test the abs(kME) path.
    # WGCNA convention: samples in rows, features in columns.
    set.seed(42)
    n_genes <- 20
    n_samples <- 5
    expr <- as.data.frame(
        matrix(rnorm(n_samples * n_genes), nrow = n_samples,
               dimnames = list(NULL, paste0("Gene", seq_len(n_genes)))))
    colors <- stats::setNames(
        c(rep("1", 10), rep("2", 10)),
        colnames(expr))
    MEs <- WGCNA::moduleEigengenes(expr, colors,
        excludeGrey = TRUE)$eigengenes
    unsigned_net <- list(
        input_data = expr,
        colors = colors,
        MEs = MEs,
        parameters = list(networkType = "unsigned")
    )
    hubs <- extract_hubs(unsigned_net, top_n = 2)
    expect_type(hubs, "character")
    expect_equal(length(hubs), 4)  # 2 modules x top_n=2
})
