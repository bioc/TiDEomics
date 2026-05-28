# Tests for plotting functions - validate returns and no-error behavior

test_that("plot_variance returns ggplot", {
    data("example")
    se <- normalise_to_start(example_obj)
    vd <- decomp_variance(se, features = rownames(se)[1:20],
        fixed_effect_var = NULL, core = 1)

    p <- plot_variance(vd, rank = "Time", top_n = 5)
    expect_s3_class(p, "ggplot")

    p2 <- plot_variance(vd, rank = "Group", top_n = 10)
    expect_s3_class(p2, "ggplot")
})

test_that("plot_trend with errorbar prints without error", {
    data("example")
    se <- normalise_to_start(example_obj)
    tbl <- calc_mean_sd(se)

    expect_no_error(
        plot_trend(tbl$orig,
            features = head(unique(tbl$orig$Feature), 3))
    )
})

test_that("plot_trend without errorbar prints without error", {
    data("example")
    se <- normalise_to_start(example_obj)
    tbl <- calc_mean_sd(se)

    expect_no_error(
        plot_trend(tbl$orig,
            features = head(unique(tbl$orig$Feature), 3),
            errorbar = FALSE)
    )
})

test_that("plot_trend with specified groups and title", {
    data("example")
    se <- normalise_to_start(example_obj)
    tbl <- calc_mean_sd(se)
    grp <- as.character(unique(tbl$orig$Group)[1])

    expect_no_error(
        plot_trend(tbl$orig,
            features = head(unique(tbl$orig$Feature), 2),
            groups = grp, title = "Test", ylab = "Value")
    )
})

test_that("plot_cor_matrix runs without error", {
    data("example")
    se <- normalise_to_start(example_obj)
    expect_no_error(
        plot_cor_matrix(se, method = "spearman")
    )
})

test_that("plot_distribution with facet_by='Group' returns ggplot", {
    data("example")
    se <- normalise_to_start(example_obj)
    p <- plot_distribution(se, facet_by = "Group")
    expect_s3_class(p, "ggplot")
})

test_that("plot_distribution with facet_by='Time' returns ggplot", {
    data("example")
    se <- normalise_to_start(example_obj)
    p <- plot_distribution(se, facet_by = "Time")
    expect_s3_class(p, "ggplot")
})

test_that("plot_distribution without facet_by returns ggplot", {
    data("example")
    se <- normalise_to_start(example_obj)
    p <- plot_distribution(se)
    expect_s3_class(p, "ggplot")
})

test_that("plot_cv returns ggplot", {
    data("example")
    se <- normalise_to_start(example_obj)
    p <- plot_cv(se)
    expect_s3_class(p, "ggplot")
})

test_that("plot_ID returns ggplot", {
    data("example")
    se <- normalise_to_start(example_obj)
    p <- plot_ID(se)
    expect_s3_class(p, "ggplot")
})

test_that("plot_missing returns ggplot", {
    data("example")
    se <- normalise_to_start(example_obj)
    p <- plot_missing(se)
    expect_s3_class(p, "ggplot")
})

test_that("plot_pca with plot=FALSE returns data", {
    data("example")
    se <- normalise_to_start(example_obj)
    res <- plot_pca(se, plot = FALSE, plot_screeplot = FALSE,
        plot_loadings = FALSE, plot_morepc = FALSE)
    expect_true(is.list(res))
    expect_true("rotated" %in% names(res))
    expect_true("variance" %in% names(res))
})

test_that("plot_pca with plot=TRUE runs", {
    data("example")
    se <- normalise_to_start(example_obj)
    expect_no_error(
        plot_pca(se, plot = TRUE, plot_screeplot = FALSE,
            plot_loadings = FALSE, plot_morepc = FALSE,
            circle = FALSE)
    )
})

test_that("plot_pca with circle=TRUE runs", {
    data("example")
    se <- normalise_to_start(example_obj)
    expect_no_error(
        plot_pca(se, plot = TRUE, plot_screeplot = FALSE,
            plot_loadings = FALSE, plot_morepc = FALSE,
            circle = TRUE)
    )
})

test_that("plot_pca with screeplot and loadings", {
    data("example")
    se <- normalise_to_start(example_obj)
    expect_no_error(
        plot_pca(se, plot = TRUE, plot_screeplot = TRUE,
            plot_loadings = TRUE, plot_morepc = FALSE,
            circle = FALSE)
    )
})


test_that("plot_umap with plot=FALSE returns data", {
    data("example")
    se <- normalise_to_start(example_obj)
    res <- plot_umap(se, plot = FALSE, seed = 42)
    expect_s3_class(res, "data.frame")
    expect_true("V1" %in% colnames(res))
    expect_true("V2" %in% colnames(res))
})

test_that("plot_umap with circle=TRUE runs without error", {
    data("example")
    se <- normalise_to_start(example_obj)
    expect_no_error(
        plot_umap(se, circle = TRUE, seed = 42)
    )
})

test_that("plot_umap with plot_ID=TRUE runs without error", {
    data("example")
    se <- normalise_to_start(example_obj)
    expect_no_error(
        plot_umap(se, plot_ID = TRUE, seed = 42)
    )
})

test_that("plot_pca_3D runs without error", {
    data("example")
    se <- normalise_to_start(example_obj)
    pca_res <- plot_pca(se, plot = FALSE, plot_screeplot = FALSE,
        plot_loadings = FALSE, plot_morepc = FALSE)
    expect_no_error(
        plot_pca_3D(pca_res, pcs = 1:3)
    )
})

test_that("plot_volcano runs with DE_between_group output", {
    data("example")
    se <- normalise_to_start(example_obj)
    de_out <- suppressMessages(DE_between_group(se, assay = 2, filter = 1))
    groups <- levels(se$Group)
    p <- plot_volcano(de_out,
        group1 = groups[1], group2 = groups[2],
        time = unique(se$Time)[2],
        logFC_thres = 0.5, adjP_thres = 0.99, label = FALSE)
    expect_s3_class(p, "ggplot")
})

test_that("plot_DE_between_time runs", {
    data("example")
    se <- normalise_to_start(example_obj)
    de_out <- suppressMessages(DE_between_time(se, assay = 1, filter = 1))
    expect_no_error(
        plot_DE_between_time(se, de_list = de_out$de_list, fontsize = 8)
    )
})

test_that("plot_pca_by_group returns ggplot", {
    data("example")
    se <- normalise_to_start(example_obj)
    p <- plot_pca_by_group(se, circle = FALSE, arrow = FALSE,
        legend_pos = "none")
    expect_s3_class(p, "ggplot")
})

test_that("plot_umap_by_group returns ggplot", {
    data("example")
    se <- normalise_to_start(example_obj)
    p <- plot_umap_by_group(se, seed = 42, legend_pos = "none")
    expect_s3_class(p, "ggplot")
})

test_that("plot_pca_arrows returns ggplot", {
    data("example")
    se <- normalise_to_start(example_obj)
    se_one <- se[, se$Group == unique(se$Group)[1]]
    p <- plot_pca_arrows(se_one, circle = FALSE, arrow = FALSE)
    expect_s3_class(p, "ggplot")
})
