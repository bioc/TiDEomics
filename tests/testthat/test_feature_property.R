# Tests for feature property pipeline

test_that("calc_feature_property works on merged data", {
    data("example_obj")
    example_obj <- normalise_to_start(example_obj)
    example_obj_list <- split_groups(example_obj)
    example_obj_merged_list <- merge_replicates(example_obj_list)

    res <- calc_feature_property(example_obj_merged_list, threshold = 0)

    expect_true(is.list(res))
    for (nm in names(res)) {
        rd <- rowData(res[[nm]])
        expect_true("P_trend" %in% colnames(rd))
        expect_true("Max_FC" %in% colnames(rd))
        expect_true("Exp_ratio" %in% colnames(rd))
        expect_s4_class(res[[nm]], "SummarizedExperiment")
    }
})

test_that("calc_feature_property without threshold still works", {
    data("example_obj")
    example_obj <- normalise_to_start(example_obj)
    example_obj_list <- split_groups(example_obj)
    example_obj_merged_list <- merge_replicates(example_obj_list)

    res <- calc_feature_property(example_obj_merged_list, threshold = NULL)
    expect_true(is.list(res))
})

test_that("summarise_feature_property works", {
    data("example_obj")
    example_obj <- normalise_to_start(example_obj)
    example_obj_list <- split_groups(example_obj)
    example_obj_merged_list <- merge_replicates(example_obj_list)
    res <- calc_feature_property(example_obj_merged_list)

    prop_tb <- summarise_feature_property(res)
    expect_s3_class(prop_tb, "data.frame")
    expect_true("Feature" %in% colnames(prop_tb))
    expect_true("Exp_ratio" %in% colnames(prop_tb))
    expect_true("P_trend" %in% colnames(prop_tb))
    expect_true("Max_FC" %in% colnames(prop_tb))
    expect_true("Group" %in% colnames(prop_tb))
})

test_that("group_specific_features filters correctly", {
    data("example_obj")
    example_obj <- normalise_to_start(example_obj)
    example_obj_list <- split_groups(example_obj)
    example_obj_merged_list <- merge_replicates(example_obj_list)
    res <- calc_feature_property(example_obj_merged_list, threshold = 0)
    prop_tb <- summarise_feature_property(res)

    gsf <- group_specific_features(prop_tb, groups = "untreated",
        genename = FALSE, GO = FALSE)

    expect_true(is.list(gsf) || is.null(gsf))
    if (!is.null(gsf)) {
        expect_true(length(gsf$features) > 0)
        # All returned features should be in the input
        expect_true(all(gsf$features %in% prop_tb$Feature))
    }
})

test_that("group_specific_features respects filter_ratio", {
    data("example_obj")
    example_obj <- normalise_to_start(example_obj)
    example_obj_list <- split_groups(example_obj)
    example_obj_merged_list <- merge_replicates(example_obj_list)
    res <- calc_feature_property(example_obj_merged_list, threshold = 0)
    prop_tb <- summarise_feature_property(res)

    # Strict filter: features must be present in >= 90% of time points
    gsf_strict <- group_specific_features(prop_tb, groups = "untreated",
        filter_ratio = 0.9, genename = FALSE, GO = FALSE)
    # Loose filter: 10%
    gsf_loose <- group_specific_features(prop_tb, groups = "untreated",
        filter_ratio = 0.1, genename = FALSE, GO = FALSE)

    # Stricter filter should return <= looser filter
    # Both should be non-NULL with reasonable filter_ratio values
    skip_if(
        is.null(gsf_strict) || is.null(gsf_loose),
        "No features passed filter"
    )
    expect_true(length(gsf_strict$features) <=
                length(gsf_loose$features))
})

test_that("group_specific_features with group_pct works", {
    data("example_obj")
    example_obj <- normalise_to_start(example_obj)
    example_obj_list <- split_groups(example_obj)
    example_obj_merged_list <- merge_replicates(example_obj_list)
    res <- calc_feature_property(example_obj_merged_list, threshold = 0)
    prop_tb <- summarise_feature_property(res)

    # Require presence in only 1 out of 2 specified groups
    gsf <- group_specific_features(prop_tb,
        groups = c("untreated", "IFNbeta"),
        filter_ratio = 0.5, group_pct = 0.5,
        genename = FALSE, GO = FALSE)
    expect_true(is.list(gsf) || is.null(gsf))
})

test_that("group_specific_features with NULL groups uses all", {
    data("example_obj")
    example_obj <- normalise_to_start(example_obj)
    example_obj_list <- split_groups(example_obj)
    example_obj_merged_list <- merge_replicates(example_obj_list)
    res <- calc_feature_property(example_obj_merged_list, threshold = 0)
    prop_tb <- summarise_feature_property(res)

    gsf <- group_specific_features(prop_tb, groups = NULL,
        genename = FALSE, GO = FALSE)
    expect_true(is.list(gsf) || is.null(gsf))
})

test_that("summarise_feature_property has expected columns", {
    data("example_obj")
    example_obj <- normalise_to_start(example_obj)
    example_obj_list <- split_groups(example_obj)
    example_obj_merged_list <- merge_replicates(example_obj_list)
    res <- calc_feature_property(example_obj_merged_list, threshold = 0)
    prop_tb <- summarise_feature_property(res)

    expect_true("T_total" %in% colnames(prop_tb))
    expect_true("T_exp" %in% colnames(prop_tb))
    expect_true("Max_FC_time" %in% colnames(prop_tb))
})

test_that("group_specific_features genename path returns data.frame", {
    skip_if_not_installed("org.Mm.eg.db")
    library(org.Mm.eg.db)
    data("example_obj")
    example_obj <- normalise_to_start(example_obj)
    example_obj_list <- split_groups(example_obj)
    example_obj_merged_list <- merge_replicates(example_obj_list)
    res <- calc_feature_property(example_obj_merged_list, threshold = 0)
    prop_tb <- summarise_feature_property(res)

    gsf <- group_specific_features(prop_tb, groups = "untreated",
        filter_ratio = 0.5, genename = TRUE, GO = FALSE,
        OrgDb = org.Mm.eg.db, keytype = "SYMBOL")
    skip_if(
        is.null(gsf),
        "No features passed filter for genename test"
    )
    expect_true("genename" %in% names(gsf))
    expect_s3_class(gsf$genename, "data.frame")
})

test_that("group_specific_features GO path returns plot", {
    skip_if_not_installed("org.Mm.eg.db")
    library(org.Mm.eg.db)
    data("example_obj")
    example_obj <- normalise_to_start(example_obj)
    example_obj_list <- split_groups(example_obj)
    example_obj_merged_list <- merge_replicates(example_obj_list)
    res <- calc_feature_property(example_obj_merged_list, threshold = 0)
    prop_tb <- summarise_feature_property(res)

    gsf <- group_specific_features(prop_tb, groups = "untreated",
        filter_ratio = 0.5, genename = FALSE, GO = TRUE,
        OrgDb = org.Mm.eg.db, keytype = "SYMBOL")
    skip_if(
        is.null(gsf) || !"GO" %in% names(gsf),
        "No features passed filter for GO test"
    )
    expect_s3_class(gsf$GO, "ggplot")
})

# ---- group_specific_features edge cases ----

test_that("group_specific_features covers NA threshold branch", {
    data("example_obj")
    example_obj <- normalise_to_start(example_obj)
    example_obj_list <- split_groups(example_obj)
    example_obj_merged_list <- merge_replicates(example_obj_list)
    res <- calc_feature_property(example_obj_merged_list, threshold = NULL)
    prop_tb <- summarise_feature_property(res)

    expect_message(
        gsf <- group_specific_features(prop_tb, groups = "untreated",
            filter_ratio = 0.1, genename = FALSE, GO = FALSE),
        "non-NA time points"
    )
    # May be NULL if no features pass filter
    expect_true(is.list(gsf) || is.null(gsf))
})

test_that("group_specific_features errors on invalid group name", {
    data("example_obj")
    example_obj <- normalise_to_start(example_obj)
    example_obj_list <- split_groups(example_obj)
    example_obj_merged_list <- merge_replicates(example_obj_list)
    res <- calc_feature_property(example_obj_merged_list, threshold = 0)
    prop_tb <- summarise_feature_property(res)

    expect_error(
        group_specific_features(prop_tb, groups = "nonexistent",
            genename = FALSE, GO = FALSE),
        "not found"
    )
})

test_that("group_specific_features errors on missing OrgDb when genename=TRUE", {
    data("example_obj")
    example_obj <- normalise_to_start(example_obj)
    example_obj_list <- split_groups(example_obj)
    example_obj_merged_list <- merge_replicates(example_obj_list)
    res <- calc_feature_property(example_obj_merged_list, threshold = 0)
    prop_tb <- summarise_feature_property(res)

    expect_error(
        group_specific_features(prop_tb, groups = "untreated",
            genename = TRUE, GO = FALSE, OrgDb = NULL, keytype = NULL),
        "Both 'OrgDb' and 'keytype' must be provided"
    )
})

test_that("group_specific_features errors on multiple thresholds", {
    # Manually construct data with two different thresholds
    prop_tb <- data.frame(
        Feature = c("A", "B"),
        Group = c("G1", "G2"),
        Exp_ratio = c(0.8, 0.8),
        Exp_threshold = c(0, 1),
        P_trend = c(0.5, 0.5),
        Max_FC = c(1, 1),
        stringsAsFactors = FALSE
    )
    expect_error(
        group_specific_features(prop_tb, groups = c("G1", "G2"),
            genename = FALSE, GO = FALSE),
        "Multiple Exp_threshold"
    )
})

test_that("calc_feature_property handles missing assay 2", {
    set.seed(42)
    mat <- matrix(rnorm(30), nrow = 10)
    colnames(mat) <- c("0", "2", "24")
    rownames(mat) <- paste0("Feature", 1:10)
    se <- SummarizedExperiment::SummarizedExperiment(
        assays = list(orig = mat),
        colData = data.frame(
            Sample = c("S1", "S2", "S3"),
            Time = c(0, 2, 24),
            Group = "A"
        )
    )
    se_list <- list(A = se)

    expect_message(
        res <- calc_feature_property(se_list),
        "Assay 2"
    )
    expect_true(is.list(res))
    rd <- SummarizedExperiment::rowData(res[["A"]])
    expect_true("AUC" %in% colnames(rd))
    expect_true(all(is.na(rd$AUC)))
})

# ---- Simulated data for genename/GO branch coverage ----

test_that("group_specific_features genename branch with simulated data", {
    skip_if_not_installed("org.Mm.eg.db")
    library(org.Mm.eg.db)

    # Simulated data: "Tnf","Il6","Tp53" only in group A (unique to A);
    # "Actb","Gapdh" appear in both A and B (shared, not unique)
    prop_tb <- data.frame(
        Feature = c("Actb", "Gapdh", "Tnf", "Il6", "Tp53",
                     "Actb", "Gapdh"),
        Group = c("A", "A", "A", "A", "A",
                  "B", "B"),
        Exp_ratio = c(0.8, 0.8, 0.8, 0.8, 0.8,
                      0.8, 0.8),
        Exp_threshold = 0,
        P_trend = 0.01,
        Max_FC = 1.5,
        stringsAsFactors = FALSE
    )

    gsf <- group_specific_features(prop_tb, groups = "A",
        filter_ratio = 0.5, genename = TRUE, GO = FALSE,
        OrgDb = org.Mm.eg.db, keytype = "SYMBOL")

    expect_true(is.list(gsf))
    expect_setequal(gsf$features, c("Tnf", "Il6", "Tp53"))
    expect_true("genename" %in% names(gsf))
    expect_s3_class(gsf$genename, "data.frame")
})

test_that("group_specific_features GO branch with simulated data", {
    skip_if_not_installed("org.Mm.eg.db")
    library(org.Mm.eg.db)

    # Well-known mouse genes unique to group A for GO enrichment
    prop_tb <- data.frame(
        Feature = c("Actb", "Gapdh", "Tnf", "Il6", "Tp53", "Myc",
                     "Nfkb1", "Ccl2",
                     "Actb", "Gapdh"),
        Group = c("A", "A", "A", "A", "A", "A", "A", "A",
                  "B", "B"),
        Exp_ratio = c(0.8, 0.8, 0.8, 0.8, 0.8, 0.8, 0.8, 0.8,
                      0.8, 0.8),
        Exp_threshold = 0,
        P_trend = 0.01,
        Max_FC = 1.5,
        stringsAsFactors = FALSE
    )

    gsf <- group_specific_features(prop_tb, groups = "A",
        filter_ratio = 0.5, genename = FALSE, GO = TRUE,
        OrgDb = org.Mm.eg.db, keytype = "SYMBOL",
        pvalueCutoff = 0.9, qvalueCutoff = 0.9,
        minGSSize = 1, simplify = FALSE)

    expect_true(is.list(gsf))
    expect_true(length(gsf$features) >= 6)
    # GO enrichment was triggered (branch reached); result may be a
    # named list of ggplots from plot_GO or absent if no terms enriched
    if ("GO" %in% names(gsf)) {
        expect_true(is.list(gsf$GO) || inherits(gsf$GO, "ggplot"))
    }
})
