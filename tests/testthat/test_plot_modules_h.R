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
        pvalueCutoff = 0.9, qvalueCutoff = 0.9,
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
        pvalueCutoff = 0.9, qvalueCutoff = 0.9,
        category = c("BP", "MF", "CC"), simplify = FALSE)

    # Multi-category enrichment
    expect_s4_class(plot_modules_h(example_module,
        example_obj_merged, scale = TRUE,
        ylabel = "Z-score of log2 expression",
        enrich_list = example_go_list$all,
        enrich_category = c("BP", "MF"),
        enrich_p_threshold = NULL,
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
        pvalueCutoff = 0.9, qvalueCutoff = 0.9,
        category = "BP", simplify = FALSE)

    # Character vector mark_features with Hub features
    hubs <- utils::head(example_module$Feature, 3)
    expect_s4_class(plot_modules_h(example_module,
        example_obj_merged, scale = TRUE,
        ylabel = "Z-score of log2 expression",
        enrich_list = example_go_list$all,
        enrich_category = c("BP", "Hub features"),
        mark_features = hubs,
        heatmap_width = 6, heatmap_height = 4), "HeatmapList")

    # Named list mark_features with Hub features
    hub_list <- list(
        "Hub1" = utils::head(example_module$Feature, 2),
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
        pvalueCutoff = 0.9, qvalueCutoff = 0.9,
        category = c("BP", "MF"), simplify = FALSE)

    # Per-category top_n and rank_by
    expect_s4_class(plot_modules_h(example_module,
        example_obj_merged, scale = TRUE,
        ylabel = "Z-score of log2 expression",
        enrich_list = example_go_list$all,
        enrich_category = c("BP", "MF"),
        enrich_top_n = 3,
        enrich_rank_by = c("p.adjust", "pvalue"),
        enrich_p_threshold = NULL,
        heatmap_width = 6, heatmap_height = 4), "HeatmapList")
})

test_that("plot_modules_h validates module argument", {
    expect_error(plot_modules_h(NULL, data.frame()), "data.frame")
    expect_error(plot_modules_h(data.frame(), data.frame()),
        "at least 1 row")
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
        pvalueCutoff = 0.9, qvalueCutoff = 0.9,
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
        pvalueCutoff = 0.9, qvalueCutoff = 0.9,
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

test_that(".reformat_cat_terms pivots category to module", {
    cat_terms <- list(
        BP = list(M1 = data.frame(text = "term1", col = "black",
                   fontsize = 8, stringsAsFactors = FALSE),
                  M2 = data.frame(text = "term2", col = "red",
                   fontsize = 8, stringsAsFactors = FALSE)),
        MF = list(M1 = data.frame(text = "term3", col = "blue",
                   fontsize = 8, stringsAsFactors = FALSE))
    )
    res <- TiDEomics:::.reformat_cat_terms(cat_terms)
    expect_named(res, c("M1", "M2"))
    expect_named(res$M1, c("BP", "MF"))
    expect_equal(res$M1$BP$text, "term1")
    # M2 missing MF gets empty placeholder
    expect_equal(nrow(res$M2$MF), 0)
})

# ---- .textbox_grob ----

test_that(".textbox_grob handles word_wrap = TRUE, first_text_from = 'top'", {
    pdf(file = NULL)
    on.exit(dev.off())

    res <- TiDEomics:::.textbox_grob(
        text = c("hello world", "foo bar baz"),
        word_wrap = TRUE,
        first_text_from = "top"
    )
    expect_s3_class(res, "gTree")
    expect_true("textbox" %in% class(res))
})

test_that(".textbox_grob handles word_wrap = TRUE, first_text_from = 'bottom'", {
    pdf(file = NULL)
    on.exit(dev.off())

    res <- TiDEomics:::.textbox_grob(
        text = c("hello world", "foo bar baz"),
        word_wrap = TRUE,
        first_text_from = "bottom"
    )
    expect_s3_class(res, "gTree")
    expect_true("textbox" %in% class(res))
})

test_that(".textbox_grob handles first_text_from = 'bottom' without word_wrap", {
    pdf(file = NULL)
    on.exit(dev.off())

    res <- TiDEomics:::.textbox_grob(
        text = c("hello", "world"),
        word_wrap = FALSE,
        first_text_from = "bottom"
    )
    expect_s3_class(res, "gTree")
    expect_true("textbox" %in% class(res))
})

test_that(".textbox_grob errors on invalid first_text_from", {
    pdf(file = NULL)
    on.exit(dev.off())

    expect_error(
        TiDEomics:::.textbox_grob(
            text = "hello",
            first_text_from = "invalid"
        ),
        "`first_text_from` can be 'top' or 'bottom'"
    )
})

test_that(".textbox_grob handles padding of length 2", {
    pdf(file = NULL)
    on.exit(dev.off())

    res <- TiDEomics:::.textbox_grob(
        text = "hello",
        padding = grid::unit(c(2, 4), "mm")
    )
    expect_s3_class(res, "gTree")
    expect_true("textbox" %in% class(res))
})

# ---- .multicol_textbox_grob ----

test_that(".multicol_textbox_grob returns nullGrob for empty cats", {
    pdf(file = NULL)
    on.exit(dev.off())

    res <- TiDEomics:::.multicol_textbox_grob(x = list())
    expect_s3_class(res, "null")
})

# ---- .measure_multicol_widths ----

test_that(".measure_multicol_widths handles empty cats", {
    pdf(file = NULL)
    on.exit(dev.off())

    res <- TiDEomics:::.measure_multicol_widths(text = list())
    expect_named(res, character(0))
})

# ---- .anno_multicol_textbox ----

test_that(".anno_multicol_textbox errors on which = 'column'", {
    expect_error(
        TiDEomics:::.anno_multicol_textbox(
            align_to = 1:3,
            text = data.frame(text = "a", col = "black", fontsize = 8),
            which = "column"
        ),
        "only be used as row annotation"
    )
})

test_that(".anno_multicol_textbox handles numeric align_to with data.frame text", {
    res <- TiDEomics:::.anno_multicol_textbox(
        align_to = 1:3,
        text = data.frame(
            text = "term", col = "black", fontsize = 8,
            stringsAsFactors = FALSE
        )
    )
    expect_s4_class(res, "AnnotationFunction")
})

test_that(".anno_multicol_textbox errors when names have no overlap", {
    expect_error(
        TiDEomics:::.anno_multicol_textbox(
            align_to = factor(c("A", "B", "C")),
            text = list(D = list(
                BP = data.frame(
                    text = "t", col = "black", fontsize = 8,
                    stringsAsFactors = FALSE
                )
            ))
        ),
        "names of `text` should have overlap to levels in `align_to`"
    )
})

test_that(".anno_multicol_textbox handles list align_to", {
    text_list <- list(
        M1 = list(BP = data.frame(text = "t1", col = "black",
                                   fontsize = 8, stringsAsFactors = FALSE)),
        M2 = list(BP = data.frame(text = "t2", col = "red",
                                   fontsize = 8, stringsAsFactors = FALSE))
    )
    align_list <- list(M1 = 1:3, M2 = 4:6)

    res <- TiDEomics:::.anno_multicol_textbox(
        align_to = align_list,
        text = text_list
    )
    expect_s4_class(res, "AnnotationFunction")
})

test_that(".anno_multicol_textbox sets default background_gp fill/col/lty/lwd", {
    text_list <- list(
        M1 = list(BP = data.frame(text = "t1", col = "black",
                                   fontsize = 8, stringsAsFactors = FALSE)),
        M2 = list(BP = data.frame(text = "t2", col = "red",
                                   fontsize = 8, stringsAsFactors = FALSE))
    )
    align_list <- list(M1 = 1:3, M2 = 4:6)

    # gpar() without fill/col/lty/lwd -- function fills defaults
    res <- TiDEomics:::.anno_multicol_textbox(
        align_to = align_list,
        text = text_list,
        background_gp = grid::gpar()
    )
    expect_s4_class(res, "AnnotationFunction")
})

test_that(".anno_multicol_textbox sets default header_gp and fontsize", {
    text_list <- list(
        M1 = list(BP = data.frame(text = "t1", col = "black",
                                   fontsize = 10, stringsAsFactors = FALSE))
    )
    align_list <- list(M1 = 1:3)

    res <- TiDEomics:::.anno_multicol_textbox(
        align_to = align_list,
        text = text_list,
        show_headers = TRUE
    )
    expect_s4_class(res, "AnnotationFunction")
})

test_that(".anno_multicol_textbox handles by = 'anno_block'", {
    text_list <- list(
        M1 = list(BP = data.frame(text = "t1", col = "black",
                                   fontsize = 8, stringsAsFactors = FALSE)),
        M2 = list(BP = data.frame(text = "t2", col = "red",
                                   fontsize = 8, stringsAsFactors = FALSE))
    )
    align_list <- list(M1 = 1:3, M2 = 4:6)

    res <- TiDEomics:::.anno_multicol_textbox(
        align_to = align_list,
        text = text_list,
        by = "anno_block"
    )
    expect_s4_class(res, "AnnotationFunction")
})

# ---- plot_modules_h enrichment / save branches ----

test_that("plot_modules_h with pre-existing col in enrichment", {
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

    enrich_mock <- list(
        BP = data.frame(
            Cluster = example_module$Module[1:5],
            Description = paste0("term_", 1:5),
            p.adjust = c(0.001, 0.002, 0.01, 0.02, 0.04),
            col = c("red", "blue", "green", "orange", "purple"),
            stringsAsFactors = FALSE
        )
    )

    expect_s4_class(
        plot_modules_h(example_module, example_obj_merged,
            scale = TRUE, enrich_list = enrich_mock,
            enrich_category = "BP", enrich_top_n = 3,
            heatmap_width = 6, heatmap_height = 4),
        "HeatmapList"
    )
})

test_that("plot_modules_h .check_len length mismatch message", {
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
        pvalueCutoff = 0.9, qvalueCutoff = 0.9,
        category = c("BP", "MF"), simplify = FALSE)

    expect_message(
        plot_modules_h(example_module, example_obj_merged,
            scale = TRUE, enrich_list = example_go_list$all,
            enrich_category = c("BP", "MF"),
            enrich_rank_by = c("p.adjust", "pvalue", "qvalue"),
            enrich_top_n = 3,
            heatmap_width = 6, heatmap_height = 4),
        "First 2 values used"
    )
})

test_that("plot_modules_h mark_features without Hub creates anno_mark", {
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
        pvalueCutoff = 0.9, qvalueCutoff = 0.9,
        category = "BP", simplify = FALSE)

    hubs <- utils::head(example_module$Feature, 3)
    expect_s4_class(
        plot_modules_h(example_module, example_obj_merged,
            scale = TRUE,
            enrich_list = example_go_list$all,
            enrich_category = "BP",
            mark_features = hubs,
            heatmap_width = 6, heatmap_height = 4),
        "HeatmapList"
    )
})

test_that("plot_modules_h with save writes files", {
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
        pvalueCutoff = 0.9, qvalueCutoff = 0.9,
        category = "BP", simplify = FALSE)

    tmpdir <- tempfile()
    dir.create(tmpdir)
    on.exit(unlink(tmpdir, recursive = TRUE))

    expect_s4_class(
        plot_modules_h(example_module, example_obj_merged,
            scale = TRUE,
            enrich_list = example_go_list$all,
            enrich_category = "BP",
            save = tmpdir, device = "png",
            heatmap_width = 6, heatmap_height = 4),
        "HeatmapList"
    )
    png_files <- list.files(tmpdir, pattern = "\\.png$")
    expect_true(length(png_files) > 0)
})

test_that("plot_modules_h skips category with no matching terms", {
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

    # BP has valid terms; MF has module IDs that don't match
    enrich_mock <- list(
        BP = data.frame(
            Cluster = example_module$Module[1:3],
            Description = paste0("term_bp_", 1:3),
            p.adjust = c(0.001, 0.002, 0.01),
            stringsAsFactors = FALSE
        ),
        MF = data.frame(
            Cluster = rep("nonexistent", 2),
            Description = c("term_mf_1", "term_mf_2"),
            p.adjust = c(0.01, 0.02),
            stringsAsFactors = FALSE
        )
    )

    expect_s4_class(
        plot_modules_h(example_module, example_obj_merged,
            scale = TRUE, enrich_list = enrich_mock,
            enrich_category = c("BP", "MF"), enrich_top_n = 3,
            heatmap_width = 6, heatmap_height = 4),
        "HeatmapList"
    )
})
