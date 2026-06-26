test_that("plot_modules_v works", {
    library(dplyr)

    data(example)
    example_obj <- normalise_to_start(example_obj)
    example_obj_list <- split_groups(example_obj)
    example_obj_merged_list <- merge_replicates(example_obj_list)
    example_obj_merged <- merge_groups(example_obj_merged_list)

    data(example_net)
    example_module <- WGCNA_module(example_net)

    pq <- plot_modules_v(example_module |> dplyr::filter(Module != '0'),
        example_obj_merged, scale = TRUE,
        ylabel = "Z-score of log2 (expression)",
        height_ratio = 2,
        fontsize = 6)

    expect_true("patchwork" %in% class(pq))
    expect_true("gg" %in% class(pq))
}
)
