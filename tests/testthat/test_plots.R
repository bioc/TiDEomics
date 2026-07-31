# ---- Shared data and precomputation ----
data(example_obj)
se <- normalise_to_start(example_obj)
vd <- decomp_variance(se, features = rownames(se)[1:20],
    fixed_effect_var = NULL, assay = "orig", core = 1)

# ---- Tests ----

test_that("plot_variance works, validates, and handles edge cases", {
    expect_s3_class(plot_variance(vd, rank = "Time", top_n = 5), "ggplot")
    expect_s3_class(plot_variance(vd, rank = "Group", top_n = 10), "ggplot")
    expect_error(plot_variance(NULL), "data.frame")
    expect_message(plot_variance(vd, rank = "Group", top_n = 100),
        "Plotting all")
    expect_error(plot_variance(vd, features = c("NotAFeature", "AlsoNot")),
        "not found")
    expect_s3_class(plot_variance(vd, rank = "Group", top_n = 5,
        show_ylab = FALSE), "ggplot")
})

test_that("plot_trend works with SE, table, errorbar, groups, and title", {
    tbl <- calc_mean_sd(se)
    feats <- utils::head(unique(tbl$orig$Feature), 3)
    grp <- as.character(unique(tbl$orig$Group)[1])

    # table input with errorbar
    expect_no_error(plot_trend(tbl$orig, features = feats))
    # table input without errorbar
    expect_no_error(plot_trend(tbl$orig, features = feats, errorbar = FALSE))
    # with groups and title
    expect_no_error(plot_trend(tbl$orig, features = utils::head(feats, 2),
        groups = grp, title = "Test", ylabel = "Value"))
    # SE directly
    expect_no_error(plot_trend(se, assay = "orig",
        features = utils::head(rownames(se), 2)))
})

test_that("plot_cor_matrix runs and validates", {
    expect_no_error(plot_cor_matrix(se, method = "spearman"))
    expect_error(plot_cor_matrix(NULL), "SummarizedExperiment")
})

test_that("plot_distribution works with and without facet_by and validates", {
    expect_s3_class(plot_distribution(se, facet_by = "Group"), "ggplot")
    expect_s3_class(plot_distribution(se, facet_by = "Time"), "ggplot")
    expect_s3_class(plot_distribution(se), "ggplot")
    expect_error(plot_distribution(se, facet_by = "invalid"),
        "should be one of")
})

test_that("plot_cv returns ggplot and validates", {
    p <- plot_cv(se)
    expect_s3_class(p, "ggplot")
    expect_error(plot_cv(NULL), "SummarizedExperiment")
})

test_that("plot_ID and plot_missing return list, validate, and handle signif", {
    p <- plot_ID(se)
    expect_type(p, "list")
    expect_s3_class(p$comparison, "ggplot")
    expect_error(plot_ID(NULL), "SummarizedExperiment")
    expect_no_error(plot_ID(se, signif = TRUE))

    p <- plot_missing(se)
    expect_type(p, "list")
    expect_s3_class(p$comparison, "ggplot")
    expect_error(plot_missing(NULL), "SummarizedExperiment")
    res <- plot_missing(example_obj, signif = TRUE)
    expect_s3_class(res$comparison, "ggplot")
})

test_that("plot_pca handles circle, screeplot, and loadings variants", {
    res <- plot_pca(se, plot_screeplot = FALSE,
        plot_loadings = FALSE, plot_morepc = FALSE,
        circle = FALSE)
    expect_type(res, "list")
    expect_true("p_list" %in% names(res))
    expect_true("pca" %in% names(res))
    expect_true("rotated" %in% names(res$pca))
    expect_true("variance" %in% names(res$pca))

    expect_no_error(
        plot_pca(se, plot_screeplot = FALSE,
            plot_loadings = FALSE, plot_morepc = FALSE,
            circle = TRUE)
    )
    expect_no_error(
        plot_pca(se, plot_screeplot = TRUE,
            plot_loadings = TRUE, plot_morepc = FALSE,
            circle = FALSE)
    )
})

test_that("plot_umap returns layout, plots, and handles circle/plot_ID", {
    res <- plot_umap(se, seed = 42)
    expect_type(res, "list")
    expect_true(all(c("p_list", "umap_layout") %in% names(res)))
    expect_s3_class(res$umap_layout, "data.frame")
    expect_true(all(c("V1", "V2") %in% colnames(res$umap_layout)))
    expect_no_error(plot_umap(se, circle = TRUE, seed = 42))
    expect_no_error(plot_umap(se, plot_ID = TRUE, seed = 42))
})

pca_res <- plot_pca(se, plot_screeplot = FALSE,
    plot_loadings = FALSE, plot_morepc = FALSE)$pca

test_that("plot_pca_3D runs and validates", {
    expect_no_error(
        plot_pca_3D(pca_res, pcs = 1:3)
    )
    expect_error(
        plot_pca_3D(NULL),
        "must be a PCA result"
    )
    expect_error(
        plot_pca_3D(list()),
        "must be a PCA result"
    )
    expect_error(
        plot_pca_3D(list(rotated = NULL, variance = NULL)),
        "must be a PCA result"
    )

    # Wrong length
    expect_error(plot_pca_3D(pca_res, pcs = 1:2),
        "Please specify three principal components")
    expect_error(plot_pca_3D(pca_res, pcs = 1:4),
        "Please specify three principal components")
    expect_error(plot_pca_3D(pca_res, pcs = c(1, 2)),
        "Please specify three principal components")
    # Out of range
    expect_error(plot_pca_3D(pca_res, pcs = c(1, 2, 999)),
        "Invalid PCs specified")
    expect_error(plot_pca_3D(pca_res, pcs = c(0, 1, 2)),
        "Invalid PCs specified")
})

test_that("plot_DE_between_group runs and validates input", {
    expect_error(
        plot_DE_between_group(list()),
        "must be the output of DE_between_group"
    )

    de_out <- suppressMessages(DE_between_group(se, assay = 1, filter = 1))

    expect_no_error(
        plot_DE_between_group(de_out, fontsize = 8)
    )
    expect_error(
        plot_DE_between_group(de_out, group = "nonexistent"),
        "not found"
    )
})

test_that("plot_pca_by_group works with SE, list, and circle/arrow variants", {
    se_list <- split_groups(se)

    expect_s3_class(plot_pca_by_group(se, circle = FALSE, arrow = FALSE,
        legend_pos = "none"), "ggplot")
    expect_s3_class(plot_pca_by_group(se, circle = TRUE, arrow = FALSE,
        legend_pos = "none"), "ggplot")
    expect_s3_class(plot_pca_by_group(se, circle = FALSE, arrow = TRUE,
        legend_pos = "none"), "ggplot")
    expect_s3_class(plot_pca_by_group(se, circle = TRUE, arrow = TRUE,
        legend_pos = "none"), "ggplot")
    # also accepts a list from split_groups()
    expect_s3_class(plot_pca_by_group(se_list, circle = FALSE, arrow = FALSE),
        "ggplot")
})

test_that("plot_umap_by_group works with SE and list", {
    se_list <- split_groups(se)

    expect_s3_class(plot_umap_by_group(se, seed = 42, legend_pos = "none"),
        "ggplot")
    expect_s3_class(plot_umap_by_group(se_list, seed = 42, legend_pos = "none"),
        "ggplot")
})

test_that("plot_pca_arrows returns ggplot and warns on multiple groups", {
    expect_warning(
        plot_pca_arrows(se, circle = FALSE, arrow = FALSE),
        "multiple groups"
    )

    se_one <- se[, se$Group == unique(se$Group)[1]]
    expect_s3_class(plot_pca_arrows(se_one, circle = FALSE, arrow = FALSE),
        "ggplot")
})

# ---- Validation tests for plot functions ----

test_that("plot_WGCNA validates net argument", {
    expect_error(plot_WGCNA(NULL), "must be a list")
    expect_error(plot_WGCNA(list()), "'net' must contain a 'colors' element")
})

test_that("plot_breakpoints validates arguments", {
    expect_error(plot_breakpoints(NULL), "must be a list")
    expect_error(plot_breakpoints(list()), "contains no groups")
})

test_that("plot_pca warns on out-of-range morepc", {
    expect_warning(
        plot_pca(se, plot_screeplot = FALSE,
            plot_loadings = FALSE,
            plot_morepc = TRUE, morepc = c(1, 2, 100),
            circle = FALSE),
        "out of range"
    )
})

# ---- plot_trend stop conditions ----

test_that("plot_trend errors on invalid inputs", {
    expect_error(plot_trend(list(), features = "Gene1", title = "test",
        ylabel = "Abundance"), "must be a SummarizedExperiment or a data.frame")
    expect_error(plot_trend(example_obj, features = character(0),
        title = "test"), "must be non-empty")
    expect_error(plot_trend(structure(list(), class = "x"), features = "Gene1",
        title = "test"), "must be a SummarizedExperiment")
    expect_error(plot_trend(example_obj,
        features = rownames(example_obj)[1:4],
        groups = "nonexistent_group", title = "test"), "not found")
})

# ---- plot_variance edge branches ----

test_that("plot_variance edge cases", {
    vd_small <- decomp_variance(se, features = rownames(se)[1:10],
        fixed_effect_var = NULL, assay = "orig", core = 1)

    expect_message(
        plot_variance(vd_small, rank = "Group", top_n = 100),
        "Plotting all"
    )
    expect_error(
        plot_variance(vd_small, features = c("NotAFeature", "AlsoNot")),
        "not found"
    )
    p <- plot_variance(vd_small, rank = "Group", top_n = 5, show_ylab = FALSE)
    expect_s3_class(p, "ggplot")
})

# ---- plot_pca edge branches ----

test_that("plot_pca validates pc1 and pc2 via missing()", {
    # Explicit pc1/pc2 triggers !missing() validation
    expect_error(
        plot_pca(se, pc1 = 0, plot_screeplot = FALSE,
            plot_loadings = FALSE, plot_morepc = FALSE, circle = FALSE),
        "positive integer"
    )
    expect_error(
        plot_pca(se, pc2 = -1, plot_screeplot = FALSE,
            plot_loadings = FALSE, plot_morepc = FALSE, circle = FALSE),
        "positive integer"
    )
})

test_that("plot_pca validates xlim and ylim", {
    expect_error(
        plot_pca(se, xlim_min = "not_numeric", plot_screeplot = FALSE,
            plot_loadings = FALSE, plot_morepc = FALSE, circle = FALSE),
        "numeric"
    )
    expect_error(
        plot_pca(se, ylim_max = "bad", plot_screeplot = FALSE,
            plot_loadings = FALSE, plot_morepc = FALSE, circle = FALSE),
        "numeric"
    )
})

test_that("plot_pca messages on missing values", {
    # Inject NA into first feature
    assay(se, 1)[1, 1] <- NA

    expect_message(
        plot_pca(se, plot_screeplot = FALSE, plot_loadings = FALSE,
            plot_morepc = FALSE, circle = FALSE),
        "missing values"
    )
})

test_that("plot_pca warns when fewer than 2 valid PCs for pairs plot", {
    # Only 1 valid PC: triggers length(morepc) < 2
    expect_warning(
        plot_pca(se, plot_screeplot = FALSE, plot_loadings = FALSE,
            plot_morepc = TRUE, morepc = 1, circle = FALSE),
        "At least two valid PCs"
    )
})
