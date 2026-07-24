# ---- Shared data ----
library(dplyr)

data(example_obj)
example_obj <- normalise_to_start(example_obj)
example_obj_list <- split_groups(example_obj)
example_obj_merged_list <- merge_replicates(example_obj_list)
example_obj_merged <- merge_groups(example_obj_merged_list)

data(example_net)
example_module <- WGCNA_module(example_net)

# ---- Tests ----

test_that("plot_modules_v works, handles variants, and saves", {
    pq <- plot_modules_v(example_module |> dplyr::filter(Module != '0'),
        example_obj_merged, scale = TRUE,
        ylabel = "Z-score of log2 (expression)",
        height_ratio = 2,
        fontsize = 6)

    expect_true("patchwork" %in% class(pq))
    expect_true("gg" %in% class(pq))

    pq <- plot_modules_v(example_module |> dplyr::filter(Module != '0'),
        example_obj_merged, scale = FALSE, height_ratio = 2, fontsize = 6)
    expect_true("patchwork" %in% class(pq))

    mod_one <- example_module |>
        dplyr::filter(Module == levels(Module)[1])
    pq <- plot_modules_v(mod_one, example_obj_merged,
        scale = TRUE, height_ratio = 2, fontsize = 6)
    expect_true("patchwork" %in% class(pq))

    tmp_dir <- tempdir()
    pq <- plot_modules_v(example_module |> dplyr::filter(Module != '0'),
        example_obj_merged, scale = TRUE, ylabel = "Abundance",
        suffix = "test_sfx", save = tmp_dir,
        device = "png", height_ratio = 2, fontsize = 6)
    expect_true("patchwork" %in% class(pq))
    expected_file <- file.path(tmp_dir,
        paste0("WGCNA_v_test_sfx_", Sys.Date(), ".png"))
    expect_true(file.exists(expected_file))
    unlink(expected_file)
})

test_that("plot_modules_v validates module argument", {
    expect_error(plot_modules_v(NULL, data.frame()), "data.frame")
    expect_error(plot_modules_v(data.frame(), data.frame()), "at least 1 row")
})
