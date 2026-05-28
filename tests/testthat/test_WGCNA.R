test_that("WGCNA works", {
    data("example")
    example_obj <- normalise_to_start(example_obj)

    wgcna_input <- prepare_WGCNA(example_obj, assay = 2, powers = seq(1, 30),
        networkType = "signed", RsquaredCut = 0.8)

    expect_contains(names(wgcna_input),
        c("data", "sample_info", "networkType", "fitIndices", "powerEstimate"))

    # wgcna_input$fitIndices
    picked_power <- wgcna_input$powerEstimate

    example_net <- run_WGCNA(wgcna_input,
        power = picked_power,
        minModuleSize = 10, # only 100 genes in the example data
        numericLabels = TRUE)

    expect_contains(names(example_net),
        c("dendrograms", "colors", "blockGenes", "input_data",
        "sample_info", "parameters"))
    expect_contains(names(example_net$parameters),
        c("power", "networkType", "minModuleSize", "numericLabels"))

    plot_WGCNA(example_net, fontsize = 8)
}
)

# ---- WGCNA_module ----

test_that("WGCNA_module returns correct structure", {
    data("example_net")
    mod <- WGCNA_module(example_net)
    expect_s3_class(mod, "data.frame")
    expect_named(mod, c("Feature", "Module"))
    expect_type(mod$Feature, "character")
    expect_s3_class(mod$Module, "factor")
})

test_that("WGCNA_module excludes grey when requested", {
    data("example_net")
    mod_all <- WGCNA_module(example_net, exclude_grey = FALSE)
    mod_nogrey <- WGCNA_module(example_net, exclude_grey = TRUE)
    expect_true(any(c("0", "grey", "gray") %in% levels(mod_all$Module)))
    expect_true(nrow(mod_nogrey) <= nrow(mod_all))
})

test_that("WGCNA_module errors on invalid input", {
    expect_error(WGCNA_module(list()),
        "must be the output of run_WGCNA")
    expect_error(WGCNA_module(data.frame()),
        "must be the output of run_WGCNA")
})

test_that("WGCNA_module Module is sorted", {
    data("example_net")
    mod <- WGCNA_module(example_net)
    expect_equal(mod, dplyr::arrange(mod, Module))
})

# ---- prepare_WGCNA ----

test_that("prepare_WGCNA validates assay index", {
    data("example")
    example_obj <- normalise_to_start(example_obj)
    expect_error(
        prepare_WGCNA(example_obj, assay = 99,
            powers = seq(1, 10), RsquaredCut = 0.8),
        "assay"
    )
})

# ---- plot_WGCNA ----

test_that("plot_WGCNA runs without error", {
    data("example_net")
    expect_no_error(plot_WGCNA(example_net, fontsize = 8))
})
