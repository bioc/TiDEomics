# ---- Shared data and precomputation ----
data(example_obj)
se <- example_obj
feats <- utils::head(rownames(se), 3)
se_list <- list(Exp1 = se, Exp2 = se)
se_disjoint <- se
rownames(se_disjoint) <- paste0("other_", seq_len(nrow(se_disjoint)))

# ---- Tests for multiple datasets (named list of SummarizedExperiment) ----

test_that("plot_trend multiple datasets: layouts, titles, and y labels", {
    p_default <- plot_trend(se_list, features = feats)
    expect_s3_class(p_default, "ggplot")
    expect_equal(p_default$labels$title, paste0("Trends across datasets in ",
        paste(unique(se$Group), collapse = ", ")))
    expect_equal(p_default$labels$y, "Log2 abundance")

    expect_s3_class(plot_trend(se_list, features = feats,
        multi_layout = "horizontal"), "ggplot")
    expect_s3_class(plot_trend(se_list, features = feats,
        multi_layout = "wrap"), "ggplot")
    expect_s3_class(plot_trend(se_list, features = feats,
        multi_layout = "overlay"), "ggplot")
})

test_that("plot_trend multiple datasets: options", {
    # custom y-axis label
    p <- plot_trend(se_list, features = feats, ylabel = "Intensity")
    expect_equal(p$labels$y, "Intensity")

    # groups subset and unknown-group error
    grp <- as.character(unique(se$Group)[1])
    expect_no_error(plot_trend(se_list, features = feats, groups = grp))
    expect_error(plot_trend(se_list, features = feats,
        groups = "nonexistent_group"), "not found")

    # errorbar off, ncol, shared y-axes, explicit title
    expect_no_error(plot_trend(se_list, features = feats, errorbar = FALSE))
    expect_no_error(plot_trend(se_list, features = feats,
        multi_layout = "wrap", ncol = 2))
    expect_no_error(plot_trend(se_list, features = feats,
        independent_y = FALSE))
    p2 <- plot_trend(se_list, features = feats, title = "My title")
    expect_equal(p2$labels$title, paste0("My title in ",
        paste(unique(se$Group), collapse = ", ")))
})

test_that("plot_trend multiple datasets: errors", {
    expect_error(plot_trend(list(se, se), features = feats),
        "must have names")
    expect_error(plot_trend(list(Exp1 = se, Bad = list(1)),
        features = feats), "not a SummarizedExperiment")
    expect_error(plot_trend(list(Exp1 = se, Exp2 = se),
        features = "NotAFeature"), "None of the requested features found")
    expect_error(plot_trend(list(Exp1 = se, Exp2 = se),
        features = character(0)), "must be non-empty")
    expect_error(plot_trend(list(Exp1 = se, Exp2 = se),
        features = feats, assay = 99), "out of range")
    expect_error(plot_trend(list(Exp1 = se, Exp2 = se),
        features = feats, multi_layout = "bogus"), "should be one of")
    expect_error(plot_trend(list(Exp1 = se, Exp2 = se),
        features = feats, ylabel = c("A", "B")),
        "must be a single character")
})

# ---- Tests for multiple datasets: PCA ----

se_with_extra <- function() {
    extra <- matrix(stats::rnorm(3 * ncol(se)), nrow = 3,
        dimnames = list(paste0("extra_", seq_len(3)), colnames(se)))
    assays <- lapply(SummarizedExperiment::assays(se), rbind, extra)
    SummarizedExperiment::SummarizedExperiment(
        assays = assays,
        colData = SummarizedExperiment::colData(se))
}

test_that("plot_pca multiple datasets: structure, titles, annotation", {
    res <- plot_pca(list(Exp1 = se, Exp2 = se),
        plot_screeplot = FALSE, plot_loadings = FALSE, plot_morepc = FALSE)
    expect_equal(names(res), c("p_list", "pca", "plot"))
    expect_equal(names(res$p_list), c("Exp1", "Exp2"))
    expect_equal(names(res$pca), c("Exp1", "Exp2"))
    expect_s3_class(res$plot, "patchwork")
    expect_match(res$plot$meta$patches$annotation$title,
        "PCA across datasets \\(all features without missing values\\)")
    expect_match(res$p_list$Exp1$p1$labels$title, "^Exp1 \\([0-9]+ features\\)$")
    expect_match(res$p_list$Exp2$p1$labels$title,
        "^Exp2 \\([0-9]+ features\\)$")
})

test_that("plot_pca multiple datasets: shared_only subsets features", {
    se2 <- se_with_extra()
    expect_message(res <- plot_pca(list(Exp1 = se, Exp2 = se2),
        shared_only = TRUE,
        plot_screeplot = FALSE, plot_loadings = FALSE, plot_morepc = FALSE),
        "Subsetting features: Exp1=")
    expect_match(res$plot$meta$patches$annotation$title,
        "PCA across datasets \\(shared features without missing values\\)")
    expect_equal(nrow(res$pca$Exp1$loadings), nrow(res$pca$Exp2$loadings))
})

test_that("plot_pca multiple datasets: names override", {
    res <- plot_pca(list(Exp1 = se, Exp2 = se), names = c("A", "B"),
        plot_screeplot = FALSE, plot_loadings = FALSE, plot_morepc = FALSE)
    expect_equal(names(res$p_list), c("Exp1", "Exp2"))
    expect_match(res$p_list$Exp1$p1$labels$title, "^A \\([0-9]+ features\\)$")
    expect_match(res$p_list$Exp2$p1$labels$title,
        "^B \\([0-9]+ features\\)$")
})

test_that("plot_pca single SE with names: single return, multi-style title", {
    res <- plot_pca(se, names = "MyData",
        plot_screeplot = FALSE, plot_loadings = FALSE, plot_morepc = FALSE)
    expect_equal(names(res), c("p_list", "pca"))
    expect_match(res$p_list$p1$labels$title, "^MyData \\([0-9]+ features\\)$")
})

test_that("plot_pca multiple datasets: errors", {
    expect_error(plot_pca(list(se, se)), "must have names")
    expect_error(plot_pca(list(Exp1 = se, Bad = list(1))),
        "not a SummarizedExperiment")
    expect_error(plot_pca(list(Exp1 = se, Exp2 = se), names = "OnlyOne"),
        "one entry per dataset")
    expect_error(plot_pca(list(Exp1 = se, Exp2 = se), shared_only = "yes"),
        "TRUE or FALSE")
    expect_error(plot_pca(list(Exp1 = se, Other = se_disjoint),
        shared_only = TRUE), "No features are shared by all datasets")
    expect_error(plot_pca(se, pc1 = 99,
        plot_screeplot = FALSE, plot_loadings = FALSE, plot_morepc = FALSE),
        "'pc1' and 'pc2' must be within the available PCs")
    expect_error(plot_pca(se, pc2 = 99,
        plot_screeplot = FALSE, plot_loadings = FALSE, plot_morepc = FALSE),
        "must be within the available PCs")
})

test_that("plot_pca single SE shared_only: no-op without message", {
    msgs <- testthat::capture_messages(plot_pca(se, shared_only = TRUE,
        plot_screeplot = FALSE, plot_loadings = FALSE, plot_morepc = FALSE))
    expect_false(any(grepl("Subsetting features", msgs)))
})

# ---- Tests for multiple datasets: UMAP ----

test_that("plot_umap multiple datasets: structure, titles, annotation", {
    res <- plot_umap(list(Exp1 = se, Exp2 = se), seed = 42)
    expect_equal(names(res), c("p_list", "umap_layout", "plot"))
    expect_equal(names(res$umap_layout), c("Exp1", "Exp2"))
    expect_s3_class(res$plot, "patchwork")
    expect_match(res$plot$meta$patches$annotation$title,
        "UMAP across datasets \\(all features without missing values\\)")
    expect_true(all(c("V1", "V2") %in% colnames(res$umap_layout$Exp1)))
    expect_match(res$p_list$Exp1$p1$labels$title, "^Exp1 \\([0-9]+ features\\)$")
    expect_match(res$p_list$Exp2$p1$labels$title,
        "^Exp2 \\([0-9]+ features\\)$")
})

test_that("plot_umap multiple datasets: shared_only and names override", {
    se2 <- se_with_extra()
    expect_message(res <- plot_umap(list(Exp1 = se, Exp2 = se2),
        seed = 42, shared_only = TRUE), "Subsetting features: Exp1=")
    expect_match(res$plot$meta$patches$annotation$title,
        "UMAP across datasets \\(shared features without missing values\\)")
    res2 <- plot_umap(list(Exp1 = se, Exp2 = se), seed = 42,
        names = c("A", "B"))
    expect_match(res2$p_list$Exp1$p1$labels$title, "^A \\([0-9]+ features\\)$")
    expect_equal(names(res2$umap_layout), c("Exp1", "Exp2"))
})

test_that("plot_umap multiple datasets: errors", {
    expect_error(plot_umap(list(se, se)), "must have names")
    expect_error(plot_umap(list(Exp1 = se, Exp2 = se), names = "OnlyOne"),
        "one entry per dataset")
    expect_error(plot_umap(list(Exp1 = se, Other = se_disjoint),
        shared_only = TRUE), "No features are shared by all datasets")
})

test_that("plot_pca/plot_umap multi: messages identify the dataset", {
    se_na <- se
    SummarizedExperiment::assay(se_na)[1, 1] <- NA
    expect_message(plot_pca(list(Exp1 = se_na, Exp2 = se),
        plot_screeplot = FALSE, plot_loadings = FALSE, plot_morepc = FALSE),
        "missing values in 'Exp1'")
    expect_message(plot_umap(list(Exp1 = se_na, Exp2 = se), seed = 42),
        "missing values in 'Exp1'")
    expect_message(plot_umap(list(Exp1 = se, Exp2 = se), seed = 42),
        "Exp1: using n_neighbors")
})

test_that("plot_pca/plot_umap multi: circle panels combined", {
    res_pca <- plot_pca(list(Exp1 = se, Exp2 = se), circle = TRUE,
        plot_screeplot = FALSE, plot_loadings = FALSE, plot_morepc = FALSE)
    panels_pca <- as.list(res_pca$plot)
    expect_length(panels_pca, 2)
    expect_equal(length(panels_pca[[1]]$layers),
        length(res_pca$p_list$Exp1$p1$layers) + 1)
    expect_equal(length(panels_pca[[2]]$layers),
        length(res_pca$p_list$Exp2$p1$layers) + 1)

    res_umap <- plot_umap(list(Exp1 = se, Exp2 = se), seed = 42,
        circle = TRUE)
    panels_umap <- as.list(res_umap$plot)
    expect_length(panels_umap, 2)
    expect_equal(length(panels_umap[[1]]$layers),
        length(res_umap$p_list$Exp1$p1$layers) + 1)
})

# ---- Tests for multiple datasets: correlation matrix ----

test_that("plot_cor_matrix multiple datasets: structure and titles", {
    res <- suppressMessages(plot_cor_matrix(list(Exp1 = se, Exp2 = se)))
    expect_equal(names(res), c("Exp1", "Exp2"))
    expect_s4_class(res[[1]], "Heatmap")
    expect_match(res[[1]]@column_title, "^Exp1 \\([0-9]+\\)$")
    expect_match(res[[2]]@column_title, "^Exp2 \\([0-9]+\\)$")
})

test_that("plot_cor_matrix multiple datasets: shared_only subsets features", {
    se2 <- se_with_extra()
    expect_message(res <- plot_cor_matrix(list(Exp1 = se, Exp2 = se2),
        shared_only = TRUE), "Subsetting features: Exp1=")
    expect_equal(res[[2]]@column_title,
        paste0("Exp2 (", nrow(se), ")"))
})

test_that("plot_cor_matrix multiple datasets: cophenetic similarity message", {
    expect_message(plot_cor_matrix(list(Exp1 = se, Exp2 = se)),
        "Cophenetic correlation between datasets \\(Pearson\\): 1.00")
    expect_message(
        plot_cor_matrix(list(Exp1 = se, Exp2 = se, Exp3 = se)),
        paste0("Cophenetic correlation between datasets \\(Pearson\\): ",
            "Exp1-Exp2 = 1.00, Exp1-Exp3 = 1.00, ",
            "Exp2-Exp3 = 1.00"))
    expect_message(
        plot_cor_matrix(list(A = se, B = se, C = se, D = se)),
        paste0("Cophenetic correlation between datasets \\(Pearson\\): ",
            "A-B = 1.00, A-C = 1.00, A-D = 1.00",
            "\\nB-C = 1.00, B-D = 1.00, C-D = 1.00"))
})

test_that("plot_cor_matrix multiple datasets: same samples any column
    order", {
    se_rev <- se[, rev(seq_len(ncol(se)))]
    expect_message(plot_cor_matrix(list(Exp1 = se, Exp2 = se_rev)),
        "Cophenetic correlation between datasets \\(Pearson\\): 1.00")
})

test_that("plot_cor_matrix multiple datasets: unmatched samples skip
    similarity", {
    se_sub <- se[, -ncol(se)]
    msgs <- testthat::capture_messages(
        plot_cor_matrix(list(Exp1 = se, Exp2 = se_sub)))
    expect_true(any(grepl("Skipping similarity statistics", msgs)))
    expect_false(any(grepl("Cophenetic", msgs)))
})

test_that("plot_cor_matrix multiple datasets: names override and errors", {
    res <- suppressMessages(plot_cor_matrix(list(Exp1 = se, Exp2 = se),
        names = c("A", "B")))
    expect_equal(names(res), c("Exp1", "Exp2"))
    expect_match(res[[1]]@column_title, "^A \\([0-9]+\\)$")
    expect_error(plot_cor_matrix(list(se, se)), "must have names")
    expect_error(plot_cor_matrix(list(Exp1 = se, Bad = list(1))),
        "not a SummarizedExperiment")
    expect_error(plot_cor_matrix(list(Exp1 = se, Exp2 = se),
        names = "OnlyOne"), "one entry per dataset")
    expect_error(plot_cor_matrix(list(Exp1 = se, Other = se_disjoint),
        shared_only = TRUE), "No features are shared by all datasets")
})
