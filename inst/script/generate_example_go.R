# Generate data/example_go.rda for use in tests and examples
library(TiDEomics)
library(org.Mm.eg.db)
data(example_net)

# select two modules for demonstration
example_module <- WGCNA_module(example_net) |>
    dplyr::filter(Module %in% c("1", "2"))
example_go <- enrichGO_list(example_module,
    OrgDb = org.Mm.eg.db,
    universe = WGCNA_module(example_net, exclude_grey = FALSE)$Feature,
    pvalueCutoff = 0.5, qvalueCutoff = 0.5,
    category = c("BP", "CC"), simplify = TRUE)

usethis::use_data(example_go, overwrite = TRUE)
