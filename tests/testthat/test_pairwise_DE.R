test_that("DE_between_group works", {
    data("example_obj")
    example_obj <- normalise_to_start(example_obj)

    DE_between_group_out <- DE_between_group(example_obj, assay = 2)

    expect_true(is.list(DE_between_group_out))
    expect_true(all(c("all_list", "de_list", "fit_list",
        "ref_groups", "all_groups") %in% names(DE_between_group_out)))
    expect_all_true(names(DE_between_group_out$all_list) ==
        names(DE_between_group_out$de_list))
    expect_true(paste0("IFNbeta-untreated") %in%
        names(DE_between_group_out$all_list))
})

test_that("DE_between_group with trend=TRUE and filter works", {
    data("example_obj")
    example_obj <- normalise_to_start(example_obj)

    DE_res <- DE_between_group(example_obj, assay = 2,
        filter = 1, trend = TRUE)

    expect_true(is.list(DE_res))
    expect_true("all_list" %in% names(DE_res))
})

test_that("DE_between_group validates group argument", {
    data("example_obj")
    example_obj <- normalise_to_start(example_obj)

    expect_error(
        DE_between_group(example_obj, assay = 2, group = "nonexistent"),
        "not found"
    )
})

test_that("DE_between_group de_list contains significant features", {
    data("example_obj")
    example_obj <- normalise_to_start(example_obj)

    DE_res <- DE_between_group(example_obj, assay = 2,
        adjP_thres = 0.05, logFC_thres = 1)

    # de_list should have features passing thresholds
    for (contrast in names(DE_res$de_list)) {
        df <- DE_res$de_list[[contrast]]
        if (is.data.frame(df) && nrow(df) > 0) {
            expect_true(all(df$adj.P.Val < 0.05))
            expect_true(all(abs(df$logFC) > 1))
        }
    }
})

test_that("DE_between_time validates parameters", {
    data("example_obj")
    example_obj <- normalise_to_start(example_obj)

    expect_error(
        DE_between_time(example_obj, assay = 1, group = "nonexistent"),
        "not found"
    )
})

test_that("DE_between_time works", {
    data("example_obj")
    example_obj <- normalise_to_start(example_obj)

    DE_between_time_out <- DE_between_time(example_obj, assay = 1)

    plot_DE_between_time(example_obj,
        de_list = DE_between_time_out$de_list,
        fontsize = 8, value = TRUE, nrow = 1, heatmap_width = 3)

    expect_true(is.list(DE_between_time_out))
    expect_true(all(c("all_list", "de_list") %in% names(DE_between_time_out)))
    expect_all_true(names(DE_between_time_out$all_list) ==
                        names(DE_between_time_out$de_list))
    expect_true(paste0("IFNbeta") %in%
                    names(DE_between_time_out$all_list))
    expect_true("t2-t0" %in% names(DE_between_time_out$all_list$IFNbeta))
})

test_that("plot_volcano works", {
    data("example_obj")
    example_obj <- normalise_to_start(example_obj)

    DE_between_group_out <- DE_between_group(example_obj, assay = 2)

    plot_volcano(DE_between_group_out, group1 = "untreated",
        group2 = "IFNbeta", time = 24,
        logFC_thres = 0.5, adjP_thres = 0.05, label = TRUE)

    plot_volcano(DE_between_group_out, group1 = "untreated",
        group2 = "IFNbeta", time = 24, label = TRUE)

    expect_error(plot_volcano(DE_between_group_out, group1 = "untreated",
        group2 = "IFNbeta", time = 30),
        "No DE results found")

    expect_error(plot_volcano(DE_between_group_out, group1 = "untreated",
        group2 = "Nonexisiting", time = 24),
        "No DE results found")

    expect_error(plot_volcano(DE_between_group_out, group = "untreated",
        group2 = "IFNbeta", time = 24),
        "provide either")

    expect_error(plot_volcano(DE_between_group_out, group1 = "untreated",
        group2 = "IFNbeta", time1 = 24),
        "provide either")

    DE_between_time_out <- DE_between_time(example_obj, assay = 1)

    plot_volcano(DE_between_time_out, group = "IFNbeta", time1 = 0,
        time2 = 24, logFC_thres = 0.5, adjP_thres = 0.05, label = TRUE)

    plot_volcano(DE_between_time_out, group = "IFNbeta", time1 = 0,
        time2 = 24)

    expect_error(plot_volcano(DE_between_time_out, group = "IFNbeta", time1 = 0,
        time2 = 30),
        "No DE results found")

    expect_error(plot_volcano(DE_between_time_out, group = "Nonexisting",
        time1 = 0, time2 = 24),
        "No DE results found")

    expect_error(plot_volcano(DE_between_time_out, group1 = "IFNbeta",
        time1 = 0, time2 = 24),
        "provide either")

    expect_error(plot_volcano(DE_between_time_out, group = "IFNbeta",
        time = 0, time2 = 24),
        "provide either")
})
