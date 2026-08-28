# ---- Shared data ----
library(dplyr)

data(example_obj)
example_obj_list <- split_groups(example_obj)
example_obj_merged_list <- merge_replicates(example_obj_list)
example_obj_merged <- merge_groups(example_obj_merged_list)

data(example_net)
example_module <- WGCNA_module(example_net, exclude_grey = TRUE)

se_list <- list(Exp1 = example_obj_merged, Exp2 = example_obj_merged)

# ---- Tests ----

test_that("plot_modules_h multiple datasets: works with per-dataset panels", {
    pq <- plot_modules_h(example_module, se_list, assay = 2,
        heatmap_width = 6, heatmap_height = 4)
    expect_s4_class(pq, "HeatmapList")
    expect_equal(base::names(pq@ht_list), c("Exp1", "Exp2"))
    expect_equal(pq@ht_list[[1]]@column_title, "Exp1")
    expect_equal(pq@ht_list[[2]]@column_title, "Exp2")
    expect_s4_class(pq@ht_list[[1]]@top_annotation, "HeatmapAnnotation")
    expect_true("Group" %in%
        base::names(pq@ht_list[[1]]@top_annotation@anno_list))
})

test_that("plot_modules_h supports group annotation modes", {
    pq_name <- plot_modules_h(example_module, se_list, assay = 2,
        group_annotation = "name", heatmap_width = 6,
        heatmap_height = 4)
    expect_true("Group" %in%
        base::names(pq_name@ht_list[[1]]@top_annotation@anno_list))

    pq_colour <- plot_modules_h(example_module, se_list, assay = 2,
        group_annotation = "colour", heatmap_width = 6,
        heatmap_height = 4)
    expect_true("Group" %in%
        base::names(pq_colour@ht_list[[1]]@top_annotation@anno_list))

    pq_both <- plot_modules_h(example_module, se_list, assay = 2,
        group_annotation = "both", heatmap_width = 6,
        heatmap_height = 4)
    expect_equal(
        base::names(pq_both@ht_list[[1]]@top_annotation@anno_list),
        c("Group_name", "Group_colour")
    )

    pq_none <- plot_modules_h(example_module, se_list, assay = 2,
        group_annotation = "none", heatmap_width = 6,
        heatmap_height = 4)
    expect_true(is.null(pq_none@ht_list[[1]]@top_annotation) ||
        length(pq_none@ht_list[[1]]@top_annotation@anno_list) == 0L)
})

test_that("plot_modules_h supports group annotations for one dataset", {
    pq_both <- plot_modules_h(example_module, example_obj_merged, assay = 2,
        group_annotation = "both", heatmap_width = 6,
        heatmap_height = 4)
    expect_equal(
        base::names(pq_both@top_annotation@anno_list),
        c("Group_name", "Group_colour")
    )

    pq_none <- plot_modules_h(example_module, example_obj_merged, assay = 2,
        group_annotation = "none", heatmap_width = 6,
        heatmap_height = 4)
    expect_true(is.null(pq_none@top_annotation) ||
        length(pq_none@top_annotation@anno_list) == 0L)
})

test_that("plot_modules_h can return native legends", {
    result <- plot_modules_h(example_module, se_list, assay = 2,
        return_legends = TRUE, heatmap_width = 6, heatmap_height = 4)
    expect_named(result, c("heatmap", "legends"))
    expect_s4_class(result$heatmap, "HeatmapList")
    expect_length(result$legends, 1L)
})

test_that("plot_modules_h validates group annotation mode", {
    expect_error(
        plot_modules_h(example_module, se_list, assay = 2,
            group_annotation = 1),
        "must be a character vector"
    )
    expect_error(
        plot_modules_h(example_module, se_list, assay = 2,
            group_annotation = NA_character_),
        "should be one of"
    )
    expect_error(
        plot_modules_h(example_module, se_list, assay = 2,
            group_annotation = "invalid"),
        "should be one of"
    )
    expect_error(
        plot_modules_h(example_module, se_list, assay = 2,
            return_legends = 1),
        "must be TRUE or FALSE"
    )
})

test_that("plot_modules_h multiple datasets: names override titles", {
    pq <- plot_modules_h(example_module, se_list, assay = 2,
        names = c("A", "B"),
        heatmap_width = 6, heatmap_height = 4)
    expect_equal(base::names(pq@ht_list), c("A", "B"))
})

test_that("plot_modules_h multiple datasets: saves with multi prefix", {
    tmp_dir <- tempdir()
    plot_modules_h(example_module, se_list, assay = 2, suffix = "test_sfx",
        save = tmp_dir, device = "png",
        heatmap_width = 6, heatmap_height = 4)
    expected_file <- file.path(tmp_dir,
        paste0("WGCNA_h_multi_test_sfx_", Sys.Date(), ".png"))
    expect_true(file.exists(expected_file))
    unlink(expected_file)
})

test_that("plot_modules_h multiple datasets: errors", {
    expect_error(plot_modules_h(example_module,
        list(example_obj_merged, example_obj_merged), assay = 2),
        "must have names")
    expect_error(plot_modules_h(example_module,
        list(Exp1 = example_obj_merged, Bad = list(1)), assay = 2),
        "not a SummarizedExperiment")
    expect_error(plot_modules_h(example_module, se_list, assay = 2,
        names = "OnlyOne"), "one entry per dataset")
    expect_error(plot_modules_h(example_module, se_list, assay = 2,
        ylabel = c("Y one", "Y two", "Y three")),
        "single character string")
    expect_error(plot_modules_h(example_module,
        list(Exp1 = example_obj_merged, Exp2 = example_obj), assay = 2),
        "'Exp2' must contain merged data")
    se_diff <- example_obj_merged
    SummarizedExperiment::assayNames(se_diff) <- c("orig", "other")
    expect_error(plot_modules_h(example_module,
        list(Exp1 = example_obj_merged, Exp2 = se_diff), assay = 2),
        "must be the same across all datasets")
})

test_that("plot_modules_h multi-input warns with per-dataset counts (grey excluded)", {
    mod_with_fake <- rbind(
        example_module,
        data.frame(Feature = "FAKE_NOT_IN_ASSAY", Module = "1",
                   stringsAsFactors = FALSE)
    )

    expect_warning(
        plot_modules_h(mod_with_fake,
            list(Exp1 = example_obj_merged, Exp2 = example_obj_merged),
            assay = 2, scale = FALSE,
            heatmap_width = 6, heatmap_height = 4),
        "1 feature \\(Exp1\\) and 1 feature \\(Exp2\\) are not shared "
    )

    mod_all_with_fake <- rbind(
        WGCNA_module(example_net),
        data.frame(Feature = "FAKE_NOT_IN_ASSAY", Module = "1",
                   stringsAsFactors = FALSE)
    )
    expect_warning(
        plot_modules_h(mod_all_with_fake,
            list(Exp1 = example_obj_merged, Exp2 = example_obj_merged),
            assay = 2, scale = FALSE,
            heatmap_width = 6, heatmap_height = 4),
        "1 feature \\(Exp1\\) and 1 feature \\(Exp2\\) are not shared "
    )
})

test_that("plot_modules_h multi-input stops when no common features", {
    merged_wrong <- example_obj_merged
    rownames(merged_wrong) <- paste0("X_", rownames(merged_wrong))

    expect_error(
        plot_modules_h(example_module,
            list(Exp1 = merged_wrong, Exp2 = merged_wrong),
            assay = 2, scale = FALSE,
            heatmap_width = 6, heatmap_height = 4),
        "No common features across the provided SE objects"
    )
})
