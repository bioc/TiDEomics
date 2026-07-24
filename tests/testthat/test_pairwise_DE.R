# ---- Shared data and DE precomputation ----
data("example_obj")
example_obj <- normalise_to_start(example_obj)

DE_between_group_out <- DE_between_group(example_obj, assay = 2,
    adjP_thres = 0.05, logFC_thres = 1)

DE_between_time_out <- DE_between_time(example_obj, assay = 1)

# ---- Tests ----

test_that("DE_between_group works and de_list contains significant features", {
    expect_true(is.list(DE_between_group_out))
    expect_true(all(c("all_list", "de_list", "fit_list",
        "ref_groups", "all_groups") %in% names(DE_between_group_out)))
    expect_all_true(names(DE_between_group_out$all_list) ==
        names(DE_between_group_out$de_list))
    expect_true(paste0("IFNbeta-untreated") %in%
        names(DE_between_group_out$all_list))

    # de_list should have features passing thresholds
    for (contrast in names(DE_between_group_out$de_list)) {
        df <- DE_between_group_out$de_list[[contrast]]
        if (is.data.frame(df) && nrow(df) > 0) {
            expect_true(all(df$adj.P.Val < 0.05))
            expect_true(all(abs(df$logFC) > 1))
        }
    }
})

test_that("DE_between_group with trend=TRUE and filter works", {
    DE_res <- DE_between_group(example_obj, assay = 2,
        filter = 1, trend = TRUE)

    expect_true(is.list(DE_res))
    expect_true("all_list" %in% names(DE_res))
})

test_that("DE_between_group validates group argument", {
    expect_error(
        DE_between_group(example_obj, assay = 2, group = "nonexistent"),
        "not found"
    )
})

test_that("DE_between_time works", {
    expect_true(is.list(DE_between_time_out))
    expect_true(all(c("all_list", "de_list") %in% names(DE_between_time_out)))
    expect_all_true(names(DE_between_time_out$all_list) ==
                        names(DE_between_time_out$de_list))
    expect_true(paste0("IFNbeta") %in%
                    names(DE_between_time_out$all_list))
    expect_true("t2-t0" %in% names(DE_between_time_out$all_list$IFNbeta))
})

test_that("DE_between_time validates parameters", {
    expect_error(
        DE_between_time(example_obj, assay = 1, group = "nonexistent"),
        "not found"
    )
})

test_that("DE_between_time errors when filter exceeds replicates", {
    expect_error(
        DE_between_time(example_obj, assay = 1, filter = 100),
        "Filter value is larger than the number of replicates"
    )
})

test_that("plot_volcano works", {
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

test_that("plot_DE_between_time handles value=TRUE and value=FALSE", {
    DE_between_time_out_a2 <- DE_between_time(example_obj, assay = 2)

    # value=TRUE with thresholds
    expect_error(
        plot_DE_between_time(DE_between_time_out_a2, fontsize = 8,
            adjP_thres = 0.05, logFC_thres = 1,
            value = TRUE, nrow = 1, heatmap_width = 3),
        NA
    )
    # value=FALSE returns ggplot
    p <- plot_DE_between_time(DE_between_time_out_a2, fontsize = 8,
        value = FALSE, nrow = 1, heatmap_width = 3)
    expect_true(inherits(p, "gg"))
})

test_that("plot_DE_between_group filters by adjP and logFC thresholds", {
    DE_res <- DE_between_group(example_obj, assay = 2, filter = 1)

    expect_error(
        plot_DE_between_group(DE_res, group = "untreated",
            adjP_thres = 0.05, logFC_thres = 1),
        NA
    )
})

test_that("DE_between_time uses duplicateCorrelation when Subject present", {
    set.seed(42)
    n_genes <- 20
    times <- c(0, 2, 4)
    groups <- c("A")
    subjects <- paste0("S", 1:4)
    n_total <- length(times) * length(subjects)
    mat <- matrix(rnorm(n_genes * n_total), nrow = n_genes)
    rownames(mat) <- paste0("Gene", 1:n_genes)
    colnames(mat) <- paste0("T", times, "_", rep(subjects, each = length(times)))

    ann <- data.frame(
        Group = rep(groups, n_total),
        Time = rep(times, length(subjects)),
        Subject = rep(subjects, each = length(times)),
        Replicate = rep(1L, n_total),
        stringsAsFactors = FALSE
    )
    ann$Sample <- colnames(mat)

    se <- SummarizedExperiment::SummarizedExperiment(
        assays = list(orig = mat),
        colData = ann
    )
    se <- normalise_to_start(se)

    res <- DE_between_time(se, assay = "norm")
    expect_true(is.list(res))
    expect_true("all_list" %in% names(res))
})
