# ---- Shared data and precomputation ----
data(example_obj)
se <- normalise_to_start(example_obj)

# ---- Tests for multiple datasets (named list of SummarizedExperiment) ----

test_that("plot_trend multiple datasets: layouts, titles, and y labels", {
    feats <- utils::head(rownames(se), 3)
    se_list <- list(RNA = se, Protein = se)

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
    feats <- utils::head(rownames(se), 3)
    se_list <- list(RNA = se, Protein = se)

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
    feats <- utils::head(rownames(se), 3)
    expect_error(plot_trend(list(se, se), features = feats),
        "must have names")
    expect_error(plot_trend(list(RNA = se, Bad = list(1)),
        features = feats), "not a SummarizedExperiment")
    expect_error(plot_trend(list(RNA = se, Protein = se),
        features = "NotAFeature"), "None of the requested features found")
    expect_error(plot_trend(list(RNA = se, Protein = se),
        features = character(0)), "must be non-empty")
    expect_error(plot_trend(list(RNA = se, Protein = se),
        features = feats, assay = 99), "out of range")
    expect_error(plot_trend(list(RNA = se, Protein = se),
        features = feats, multi_layout = "bogus"), "should be one of")
    expect_error(plot_trend(list(RNA = se, Protein = se),
        features = feats, ylabel = c("A", "B")),
        "must be a single character")
})
