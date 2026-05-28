test_that("plot_modules_h works", {
    skip_if_not_installed("org.Mm.eg.db")
    library(dplyr)
    library(org.Mm.eg.db)

    data(example)
    example_obj <- normalise_to_start(example_obj)
    example_obj_list <- split_groups(example_obj)
    example_obj_merged_list <- merge_replicates(example_obj_list)
    example_obj_merged <- merge_groups(example_obj_merged_list)

    data(example_net)
    example_module <- WGCNA_module(example_net, exclude_grey = TRUE) 

    # set cutoff to 1 to show all results for demonstration
    example_go_list <- enrichGO_list(example_module, OrgDb = org.Mm.eg.db,
        universe = example_module$Feature,
        pvalueCutoff = 1, qvalueCutoff = 1,
        category = "BP", simplify = FALSE)
    plot_GO(example_go_list$all, plot_dotplot = TRUE,
        plot_emapplot = FALSE, plot_cnetplot = FALSE)

    expect_null(plot_modules_h(example_module %>% dplyr::filter(Module != "0"),
        example_obj_merged, scale = TRUE,
        ylabel = "Z-score of log2 expression",
        enrich_list = example_go_list$all, enrich_category = "BP",
        heatmap_width = 6, heatmap_height = 4))

    expect_warning(plot_modules_h(example_module %>% dplyr::filter(Module != "0"),
        example_obj_merged, scale = TRUE,
        ylabel = "Z-score of log2 expression",
        enrich_list = example_go_list$all, enrich_category = "CC",
        heatmap_width = 6, heatmap_height = 4))
}
)
