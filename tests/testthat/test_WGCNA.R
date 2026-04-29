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
