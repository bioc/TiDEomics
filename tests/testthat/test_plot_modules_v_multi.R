# ---- Shared data ----
library(dplyr)

data(example_obj)
example_obj_list <- split_groups(example_obj)
example_obj_merged_list <- merge_replicates(example_obj_list)
example_obj_merged <- merge_groups(example_obj_merged_list)

data(example_net)
module_df <- WGCNA_module(example_net) |>
    dplyr::filter(Module != "0")

se_list <- list(Exp1 = example_obj_merged, Exp2 = example_obj_merged)

# ---- Tests ----

test_that("plot_modules_v multiple datasets: titles and subtitles", {
    cases <- list(
        list(names = NULL, title = "Module profiles: Exp1, Exp2",
            subs = c("Exp1", "Exp2")),
        list(names = c("A", "B"), title = "Module profiles: A, B",
            subs = c("A", "B"))
    )
    for (case in cases) {
        pq <- plot_modules_v(module_df, se_list, assay = 2,
            names = case$names)
        expect_true("patchwork" %in% class(pq))
        expect_match(pq$meta$patches$annotation$title, case$title)
        panels <- as.list(pq)
        expect_equal(length(panels), 4)
        expect_equal(vapply(panels[c(1, 3)], function(p) p$labels$subtitle,
            character(1)), case$subs)
        expect_equal(vapply(panels[c(1, 3)], function(p) p$labels$y,
            character(1)), c("Z-score of log2 abundance",
            "Z-score of log2 abundance"))
        expect_null(panels[[2]]$labels$subtitle)
        expect_null(panels[[4]]$labels$subtitle)
    }
})

test_that("plot_modules_v single SE with names: no multi titles", {
    pq <- plot_modules_v(module_df, example_obj_merged, assay = 2,
        names = "OnlyOne")
    expect_true("patchwork" %in% class(pq))
    expect_null(pq$meta$patches$annotation$title)
    panels <- as.list(pq)
    expect_equal(length(panels), 2)
    expect_null(panels[[1]]$labels$subtitle)
})

test_that("plot_modules_v multiple datasets: per-dataset ylabel", {
    pq <- plot_modules_v(module_df, se_list, assay = 2,
        ylabel = c("Y one", "Y two"))
    panels <- as.list(pq)
    expect_equal(vapply(panels[c(1, 3)], function(p) p$labels$y,
        character(1)), c("Y one", "Y two"))
    expect_null(panels[[2]]$labels$y)
    expect_null(panels[[4]]$labels$y)
})

test_that("plot_modules_v multiple datasets: saves with multi prefix", {
    tmp_dir <- tempdir()
    plot_modules_v(module_df, se_list, assay = 2, suffix = "test_sfx",
        save = tmp_dir, device = "png")
    expected_file <- file.path(tmp_dir,
        paste0("WGCNA_v_multi_test_sfx_", Sys.Date(), ".png"))
    expect_true(file.exists(expected_file))
    unlink(expected_file)
})

test_that("plot_modules_v multiple datasets: errors", {
    expect_error(plot_modules_v(module_df, list(example_obj_merged,
        example_obj_merged), assay = 2), "must have names")
    expect_error(plot_modules_v(module_df,
        list(Exp1 = example_obj_merged, Bad = list(1)), assay = 2),
        "not a SummarizedExperiment")
    expect_error(plot_modules_v(module_df, se_list, assay = 2,
        names = "OnlyOne"), "one entry per dataset")
    expect_error(plot_modules_v(module_df, se_list, assay = 2,
        ylabel = c("Y one", "Y two", "Y three")),
        "length equal to the number of datasets \\(2\\)")
    expect_error(plot_modules_v(module_df,
        list(Exp1 = example_obj_merged, Exp2 = example_obj), assay = 2),
        "'Exp2' must contain merged data")
    se_diff <- example_obj_merged
    SummarizedExperiment::assayNames(se_diff) <- c("orig", "other")
    expect_error(plot_modules_v(module_df,
        list(Exp1 = example_obj_merged, Exp2 = se_diff), assay = 2),
        "must be the same across all datasets")
})

test_that("plot_modules_v multiple datasets: stops when no shared features", {
    mod_wrong <- data.frame(Feature = "NOT_IN_SE", Module = "1",
        stringsAsFactors = FALSE)

    expect_error(
        plot_modules_v(mod_wrong, se_list, assay = 2),
        "No common features across the provided SE objects and module table"
    )
})
