# ---- Shared data and enrichment precomputation ----
library(dplyr)
data(example_obj)
data(example_net)
data(example_go)
se_list <- split_groups(example_obj)
merged_list <- merge_replicates(se_list)
merged <- merge_groups(merged_list)
example_module <- WGCNA_module(example_net, exclude_grey = TRUE)

# ---- Tests ----

test_that("plot_modules_h works", {
    pq <- plot_modules_h(example_module |> dplyr::filter(Module != '0'),
        merged, assay = 2, scale = TRUE,
        ylabel = "Z-score of log2 (expression)",
        enrich_list = example_go$all, enrich_category = "BP",
        fontsize = 6)

    expect_s4_class(pq, "HeatmapList")
    pq <- plot_modules_h(example_module |> dplyr::filter(Module != '0'),
        merged, assay = 2, scale = FALSE,
        enrich_list = example_go$all, enrich_category = "BP",
        fontsize = 6)
    expect_s4_class(pq, "HeatmapList")
})

test_that("plot_modules_h validates module argument", {
    expect_error(plot_modules_h(NULL, data.frame(), assay = 2), "data.frame")
    expect_error(plot_modules_h(data.frame(), data.frame(), assay = 2), "at least 1 row")
})

test_that("plot_modules_h handles single module without enrichment", {
    mod_one <- example_module |>
        dplyr::filter(Module == levels(Module)[1])

    pq <- plot_modules_h(mod_one, merged,
        assay = 2, scale = TRUE, enrich_category = NULL, fontsize = 6)
    expect_s4_class(pq, "HeatmapList")
})

test_that("plot_modules_h handles Features without enrichment", {
    hubs <- utils::head(example_module$Feature, 3)

    expect_s4_class(plot_modules_h(example_module, merged,
        assay = 2, enrich_list = NULL,
        enrich_category = "Features", mark_features = hubs,
        heatmap_width = 6, heatmap_height = 4), "HeatmapList")
})

test_that("plot_modules_h excludes Features from enrichment parameters", {
    hubs <- utils::head(example_module$Feature, 3)

    expect_s4_class(plot_modules_h(example_module, merged,
        assay = 2, enrich_list = example_go$all,
        enrich_category = c("BP", "Features", "CC"),
        enrich_rank_by = c("p.adjust", "pvalue"),
        enrich_top_n = c(3, 2),
        enrich_threshold = c(0.05, 0.1), mark_features = hubs,
        heatmap_width = 6, heatmap_height = 4), "HeatmapList")
})

test_that("plot_modules_h supports multi-category enrichment and parameters", {
    expect_s4_class(plot_modules_h(example_module,
        merged, assay = 2, scale = TRUE,
        ylabel = "Z-score of log2 expression",
        enrich_list = example_go$all,
        enrich_category = c("BP", "CC"),
        enrich_top_n = c(3, 2),
        enrich_rank_by = c("p.adjust", "pvalue"),
        enrich_threshold = NULL,
        heatmap_width = 6, heatmap_height = 4), "HeatmapList")
})

test_that("plot_modules_h handles mark_features with Features", {
    hubs <- utils::head(example_module$Feature, 3)

    # Character vector mark_features
    expect_s4_class(plot_modules_h(example_module,
        merged, assay = 2, scale = TRUE,
        ylabel = "Z-score of log2 expression",
        enrich_list = example_go$all,
        enrich_category = c("BP", "Features"),
        mark_features = hubs,
        heatmap_width = 6, heatmap_height = 4), "HeatmapList")

    # Named list mark_features
    hub_list <- list(
        "Hub1" = utils::head(example_module$Feature, 2),
        "Hub2" = tail(example_module$Feature, 2)
    )
    expect_s4_class(plot_modules_h(example_module,
        merged, assay = 2, scale = TRUE,
        ylabel = "Z-score of log2 expression",
        enrich_list = example_go$all,
        enrich_category = c("BP", "Features"),
        mark_features = hub_list,
        heatmap_width = 6, heatmap_height = 4), "HeatmapList")
})

test_that("plot_modules_h accepts feature annotation tables", {
    feature_table <- data.frame(
        Feature = utils::head(example_module$Feature, 3),
        Set = c("Consensus", "Exp1", "Consensus"),
        stringsAsFactors = FALSE
    )

    expect_s4_class(plot_modules_h(example_module,
        merged, assay = 2, scale = TRUE,
        enrich_list = example_go$all,
        enrich_category = c("BP", "Features"),
        mark_features = feature_table,
        heatmap_width = 6, heatmap_height = 4), "HeatmapList")
})

test_that("plot_modules_h rejects categories absent from enrich_list", {
    expect_error(plot_modules_h(example_module, merged, assay = 2,
        enrich_category = "Hub features", mark_features =
            utils::head(example_module$Feature, 3)),
        "not present in 'enrich_list'")
    expect_error(plot_modules_h(example_module, merged, assay = 2,
        enrich_list = example_go$all, enrich_category = "Not a category"),
        "not present in 'enrich_list'")
})

test_that("plot_modules_h validates enrichment parameters", {
    expect_error(plot_modules_h(example_module, merged, assay = 2,
        enrich_list = example_go$all, enrich_category = 123), "character")
    expect_error(plot_modules_h(example_module, merged, assay = 2,
        enrich_list = example_go$all, enrich_category = NA_character_),
        "missing or empty")
    expect_error(plot_modules_h(example_module, merged, assay = 2,
        enrich_list = example_go$all,
        enrich_category = c("BP", "BP")), "must be unique")
    expect_error(plot_modules_h(example_module, merged, assay = 2,
        enrich_list = example_go$all, enrich_top_n = -1), "positive")
    expect_error(plot_modules_h(example_module, merged, assay = 2,
        enrich_list = example_go$all, enrich_rank_by = NA_character_),
        "missing or empty")
    expect_error(plot_modules_h(example_module, merged, assay = 2,
        enrich_list = example_go$all, enrich_rank_by = ""),
        "missing or empty")
    expect_error(plot_modules_h(example_module, merged, assay = 2,
        enrich_list = example_go$all, enrich_rank_by = "not_a_column"),
        "not found")
    enrich_non_numeric <- example_go$all
    enrich_non_numeric$BP <- as.data.frame(enrich_non_numeric$BP)
    enrich_non_numeric$BP$p.adjust <-
        as.character(enrich_non_numeric$BP$p.adjust)
    expect_error(plot_modules_h(example_module, merged, assay = 2,
        enrich_list = enrich_non_numeric, enrich_rank_by = "p.adjust"),
        "must be numeric")
    expect_error(plot_modules_h(example_module, merged, assay = 2,
        fontsize = -1), "positive")
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
    expect_equal(nrow(res$M2$MF), 0)
})

test_that(".textbox_grob handles supported text layout options", {
    pdf(file = NULL)
    on.exit(dev.off())

    res <- TiDEomics:::.textbox_grob(
        text = c("hello world", "foo bar baz"),
        word_wrap = TRUE,
        first_text_from = "top"
    )
    expect_s3_class(res, "gTree")
    expect_true("textbox" %in% class(res))
    res_bottom <- TiDEomics:::.textbox_grob(
        text = c("hello world", "foo bar baz"),
        word_wrap = TRUE,
        first_text_from = "bottom"
    )
    expect_s3_class(res_bottom, "gTree")
    expect_true("textbox" %in% class(res_bottom))

    res_no_wrap <- TiDEomics:::.textbox_grob(
        text = c("hello", "world"),
        word_wrap = FALSE,
        first_text_from = "bottom"
    )
    expect_s3_class(res_no_wrap, "gTree")
    expect_true("textbox" %in% class(res_no_wrap))

    res_padding <- TiDEomics:::.textbox_grob(
        text = "hello",
        padding = grid::unit(c(2, 4), "mm")
    )
    expect_s3_class(res_padding, "gTree")
    expect_true("textbox" %in% class(res_padding))
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

test_that(".multicol_textbox_grob returns nullGrob for empty cats", {
    pdf(file = NULL)
    on.exit(dev.off())

    res <- TiDEomics:::.multicol_textbox_grob(x = list())
    expect_s3_class(res, "null")
})

test_that(".measure_multicol_widths handles empty cats", {
    pdf(file = NULL)
    on.exit(dev.off())

    res <- TiDEomics:::.measure_multicol_widths(text = list())
    expect_named(res, character(0))
})

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

test_that("plot_modules_h with pre-existing col in enrichment", {
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
        plot_modules_h(example_module, merged,
            assay = 2, scale = TRUE, enrich_list = enrich_mock,
            enrich_category = "BP", enrich_top_n = 3,
            heatmap_width = 6, heatmap_height = 4),
        "HeatmapList"
    )
})

test_that("plot_modules_h rejects enrichment parameter length mismatch", {
    expect_error(
        plot_modules_h(example_module, merged,
            assay = 2, scale = TRUE, enrich_list = example_go$all,
            enrich_category = c("BP", "CC"),
            enrich_rank_by = c("p.adjust", "pvalue", "qvalue"),
            enrich_top_n = 3,
            heatmap_width = 6, heatmap_height = 4),
        "length 1 or length 2"
    )
})

test_that("plot_modules_h mark_features without Hub creates anno_mark", {
    hubs <- utils::head(example_module$Feature, 3)
    expect_s4_class(
        plot_modules_h(example_module, merged,
            assay = 2, scale = TRUE,
            enrich_list = example_go$all,
            enrich_category = "BP",
            mark_features = hubs,
            heatmap_width = 6, heatmap_height = 4),
        "HeatmapList"
    )
})

test_that("plot_modules_h with save writes files", {
    tmpdir <- tempfile()
    dir.create(tmpdir)
    on.exit(unlink(tmpdir, recursive = TRUE))

    expect_s4_class(
        plot_modules_h(example_module, merged,
            assay = 2, scale = TRUE,
            enrich_list = example_go$all,
            enrich_category = "BP",
            save = tmpdir, device = "png",
            heatmap_width = 6, heatmap_height = 4),
        "HeatmapList"
    )
    png_files <- list.files(tmpdir, pattern = "\\.png$")
    expect_true(length(png_files) > 0)
})

test_that("plot_modules_h skips category with no matching terms", {
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
        plot_modules_h(example_module, merged,
            assay = 2, scale = TRUE, enrich_list = enrich_mock,
            enrich_category = c("BP", "MF"), enrich_top_n = 3,
            heatmap_width = 6, heatmap_height = 4),
        "HeatmapList"
    )
})

test_that("plot_modules_h warns on features missing from assay", {
    mod_with_fake <- rbind(
        example_module,
        data.frame(Feature = "FAKE_NOT_IN_ASSAY", Module = "1",
                   stringsAsFactors = FALSE)
    )

    expect_warning(
        plot_modules_h(mod_with_fake, merged,
            assay = 2, scale = FALSE, heatmap_width = 6, heatmap_height = 4),
        "feature\\(s\\) in the module are missing from the assay"
    )
})

test_that("plot_modules_h stops when no module features are in the SE", {
    mod_wrong <- data.frame(Feature = "NOT_IN_SE", Module = "1",
        stringsAsFactors = FALSE)

    expect_error(
        plot_modules_h(mod_wrong, merged,
            assay = 2, scale = FALSE, heatmap_width = 6, heatmap_height = 4),
        "derived from the same data"
    )
})
