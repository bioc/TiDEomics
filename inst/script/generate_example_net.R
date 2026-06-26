library(TiDEomics)

# Load example data
data("example")
example_obj <- normalise_to_start(example_obj)

# Prepare WGCNA input
wgcna_input <- prepare_WGCNA(example_obj,
    assay = 2,
    powers = seq(1, 30),
    networkType = "signed",
    RsquaredCut = 0.8
)

# Pick power
picked_power <- wgcna_input$powerEstimate

# Run WGCNA
example_net <- run_WGCNA(wgcna_input,
    power = picked_power,
    minModuleSize = 10, # only 100 genes in the example data
    numericLabels = TRUE)

# Save the data
usethis::use_data(example_net, overwrite = TRUE)
