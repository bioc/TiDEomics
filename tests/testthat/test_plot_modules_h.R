test_that("plot_modules_h works", {
    skip_if_not_installed("org.Mm.eg.db")
    library(dplyr)
    library(org.Mm.eg.db)

    data(example_obj)
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

    expect_s4_class(plot_modules_h(example_module |> dplyr::filter(Module != "0"),
        example_obj_merged, scale = TRUE,
        ylabel = "Z-score of log2 expression",
        enrich_list = example_go_list$all, enrich_category = "BP",
        heatmap_width = 6, heatmap_height = 4), "HeatmapList")

    expect_warning(plot_modules_h(example_module |> dplyr::filter(Module != "0"),
        example_obj_merged, scale = TRUE,
        ylabel = "Z-score of log2 expression",
        enrich_list = example_go_list$all, enrich_category = "CC",
        heatmap_width = 6, heatmap_height = 4))
})

test_that("plot_modules_h supports multi-category enrichment", {
    skip_if_not_installed("org.Mm.eg.db")
    library(dplyr)
    library(org.Mm.eg.db)

    data(example_obj)
    example_obj <- normalise_to_start(example_obj)
    example_obj_list <- split_groups(example_obj)
    example_obj_merged_list <- merge_replicates(example_obj_list)
    example_obj_merged <- merge_groups(example_obj_merged_list)

    data(example_net)
    example_module <- WGCNA_module(example_net, exclude_grey = TRUE)

    example_go_list <- enrichGO_list(example_module, OrgDb = org.Mm.eg.db,
        universe = example_module$Feature,
        pvalueCutoff = 1, qvalueCutoff = 1,
        category = c("BP", "MF", "CC"), simplify = FALSE)

    # Multi-category enrichment
    expect_s4_class(plot_modules_h(example_module,
        example_obj_merged, scale = TRUE,
        ylabel = "Z-score of log2 expression",
        enrich_list = example_go_list$all,
        enrich_category = c("BP", "MF"),
        enrich_p_threshold = c(0.05, 0.05),
        heatmap_width = 6, heatmap_height = 4), "HeatmapList")
})

test_that("plot_modules_h handles mark_features with Hub features", {
    skip_if_not_installed("org.Mm.eg.db")
    library(dplyr)
    library(org.Mm.eg.db)

    data(example_obj)
    example_obj <- normalise_to_start(example_obj)
    example_obj_list <- split_groups(example_obj)
    example_obj_merged_list <- merge_replicates(example_obj_list)
    example_obj_merged <- merge_groups(example_obj_merged_list)

    data(example_net)
    example_module <- WGCNA_module(example_net, exclude_grey = TRUE)

    example_go_list <- enrichGO_list(example_module, OrgDb = org.Mm.eg.db,
        universe = example_module$Feature,
        pvalueCutoff = 1, qvalueCutoff = 1,
        category = "BP", simplify = FALSE)

    # Character vector mark_features with Hub features
    hubs <- head(example_module$Feature, 3)
    expect_s4_class(plot_modules_h(example_module,
        example_obj_merged, scale = TRUE,
        ylabel = "Z-score of log2 expression",
        enrich_list = example_go_list$all,
        enrich_category = c("BP", "Hub features"),
        mark_features = hubs,
        heatmap_width = 6, heatmap_height = 4), "HeatmapList")

    # Named list mark_features with Hub features
    hub_list <- list(
        "Hub1" = head(example_module$Feature, 2),
        "Hub2" = tail(example_module$Feature, 2)
    )
    expect_s4_class(plot_modules_h(example_module,
        example_obj_merged, scale = TRUE,
        ylabel = "Z-score of log2 expression",
        enrich_list = example_go_list$all,
        enrich_category = c("BP", "Hub features"),
        mark_features = hub_list,
        heatmap_width = 6, heatmap_height = 4), "HeatmapList")
})

test_that("plot_modules_h handles per-category enrich_top_n and enrich_rank_by", {
    skip_if_not_installed("org.Mm.eg.db")
    library(dplyr)
    library(org.Mm.eg.db)

    data(example_obj)
    example_obj <- normalise_to_start(example_obj)
    example_obj_list <- split_groups(example_obj)
    example_obj_merged_list <- merge_replicates(example_obj_list)
    example_obj_merged <- merge_groups(example_obj_merged_list)

    data(example_net)
    example_module <- WGCNA_module(example_net, exclude_grey = TRUE)

    example_go_list <- enrichGO_list(example_module, OrgDb = org.Mm.eg.db,
        universe = example_module$Feature,
        pvalueCutoff = 1, qvalueCutoff = 1,
        category = c("BP", "MF"), simplify = FALSE)

    # Per-category top_n and rank_by
    expect_s4_class(plot_modules_h(example_module,
        example_obj_merged, scale = TRUE,
        ylabel = "Z-score of log2 expression",
        enrich_list = example_go_list$all,
        enrich_category = c("BP", "MF"),
        enrich_top_n = c(3, 5),
        enrich_rank_by = c("p.adjust", "pvalue"),
        enrich_p_threshold = c(0.05, 0.01),
        heatmap_width = 6, heatmap_height = 4), "HeatmapList")
})

test_that("plot_modules_h validates module argument", {
    expect_error(plot_modules_h(NULL, data.frame()), "data.frame")
    expect_error(plot_modules_h(data.frame(), data.frame()),
        "SummarizedExperiment")
    data(example_obj)
    example_obj <- normalise_to_start(example_obj)
    example_obj_list <- split_groups(example_obj)
    example_obj_merged_list <- merge_replicates(example_obj_list)
    example_obj_merged <- merge_groups(example_obj_merged_list)
    data(example_net)
    mod <- WGCNA_module(example_net, exclude_grey = TRUE)
    expect_error(plot_modules_h(mod, example_obj_merged,
        enrich_category = 123), "character")
    expect_error(plot_modules_h(mod, example_obj_merged,
        enrich_top_n = -1), "positive")
    expect_error(plot_modules_h(mod, example_obj_merged,
        fontsize = -1), "positive")
})

test_that("plot_modules_h runs without enrich_list", {
    library(dplyr)
    data(example_obj)
    example_obj <- normalise_to_start(example_obj)
    example_obj_list <- split_groups(example_obj)
    example_obj_merged_list <- merge_replicates(example_obj_list)
    example_obj_merged <- merge_groups(example_obj_merged_list)

    data(example_net)
    example_module <- WGCNA_module(example_net, exclude_grey = TRUE)

    expect_s4_class(plot_modules_h(example_module,
        example_obj_merged, scale = TRUE,
        ylabel = "Z-score of log2 expression",
        heatmap_width = 6, heatmap_height = 4), "HeatmapList")
})

test_that("plot_modules_h handles NULL mark_features safely", {
    skip_if_not_installed("org.Mm.eg.db")
    library(dplyr)
    library(org.Mm.eg.db)
    data(example_obj)
    example_obj <- normalise_to_start(example_obj)
    example_obj_list <- split_groups(example_obj)
    example_obj_merged_list <- merge_replicates(example_obj_list)
    example_obj_merged <- merge_groups(example_obj_merged_list)
    data(example_net)
    example_module <- WGCNA_module(example_net, exclude_grey = TRUE)
    example_go_list <- enrichGO_list(example_module, OrgDb = org.Mm.eg.db,
        universe = example_module$Feature,
        pvalueCutoff = 1, qvalueCutoff = 1,
        category = "BP", simplify = FALSE)

    # mark_features = NULL should not error
    expect_s4_class(plot_modules_h(example_module,
        example_obj_merged, scale = TRUE,
        ylabel = "Z-score of log2 expression",
        enrich_list = example_go_list$all,
        enrich_category = "BP",
        mark_features = NULL,
        heatmap_width = 6, heatmap_height = 4), "HeatmapList")
})

test_that("plot_modules_h handles enrich_p_threshold = NULL and NA", {
    skip_if_not_installed("org.Mm.eg.db")
    library(dplyr)
    library(org.Mm.eg.db)

    data(example_obj)
    example_obj <- normalise_to_start(example_obj)
    example_obj_list <- split_groups(example_obj)
    example_obj_merged_list <- merge_replicates(example_obj_list)
    example_obj_merged <- merge_groups(example_obj_merged_list)

    data(example_net)
    example_module <- WGCNA_module(example_net, exclude_grey = TRUE)

    example_go_list <- enrichGO_list(example_module, OrgDb = org.Mm.eg.db,
        universe = example_module$Feature,
        pvalueCutoff = 1, qvalueCutoff = 1,
        category = "BP", simplify = FALSE)

    # enrich_p_threshold = NULL (skip colour-coding)
    expect_s4_class(plot_modules_h(example_module,
        example_obj_merged, scale = TRUE,
        ylabel = "Z-score of log2 expression",
        enrich_list = example_go_list$all,
        enrich_category = "BP",
        enrich_p_threshold = NULL,
        heatmap_width = 6, heatmap_height = 4), "HeatmapList")
})
