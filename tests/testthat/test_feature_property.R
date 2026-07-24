# ---- Shared pipeline output ----
data("example_obj")
example_obj <- normalise_to_start(example_obj)
example_obj_list <- split_groups(example_obj)
example_obj_merged_list <- merge_replicates(example_obj_list)

se_list_prop <- calc_feature_property(example_obj_merged_list, threshold = 0)
prop_tb <- summarise_feature_property(se_list_prop)

# ---- Tests ----

test_that("calc_feature_property works on merged data", {
    expect_true(is.list(se_list_prop))
    for (nm in names(se_list_prop)) {
        rd <- rowData(se_list_prop[[nm]])
        expect_true("P_trend" %in% colnames(rd))
        expect_true("Max_FC" %in% colnames(rd))
        expect_true("Exp_ratio" %in% colnames(rd))
        expect_s4_class(se_list_prop[[nm]], "SummarizedExperiment")
    }
})

test_that("calc_feature_property without threshold still works", {
    res <- calc_feature_property(example_obj_merged_list, threshold = NULL)
    expect_true(is.list(res))
})

test_that("summarise_feature_property works", {
    expect_s3_class(prop_tb, "data.frame")
    expect_true("Feature" %in% colnames(prop_tb))
    expect_true("Exp_ratio" %in% colnames(prop_tb))
    expect_true("P_trend" %in% colnames(prop_tb))
    expect_true("Max_FC" %in% colnames(prop_tb))
    expect_true("Group" %in% colnames(prop_tb))
})

test_that("summarise_feature_property has expected columns", {
    expect_true("T_total" %in% colnames(prop_tb))
    expect_true("T_exp" %in% colnames(prop_tb))
    expect_true("Max_FC_time" %in% colnames(prop_tb))
})

test_that("group_specific_features filters correctly", {
    gsf <- group_specific_features(prop_tb, groups = "untreated",
        genename = FALSE, GO = FALSE)
    expect_true(is.list(gsf) || is.null(gsf))
    if (!is.null(gsf)) {
        expect_true(length(gsf$features) > 0)
        expect_true(all(gsf$features %in% prop_tb$Feature))
    }
})

test_that("group_specific_features respects filter_ratio", {
    gsf_strict <- group_specific_features(prop_tb, groups = "untreated",
        filter_ratio = 0.9, genename = FALSE, GO = FALSE)
    gsf_loose <- group_specific_features(prop_tb, groups = "untreated",
        filter_ratio = 0.1, genename = FALSE, GO = FALSE)

    skip_if(is.null(gsf_strict) || is.null(gsf_loose),
        "No features passed filter")
    expect_true(length(gsf_strict$features) <= length(gsf_loose$features))
})

test_that("group_specific_features with group_pct works", {
    gsf <- group_specific_features(prop_tb,
        groups = c("untreated", "IFNbeta"),
        filter_ratio = 0.5, group_pct = 0.5,
        genename = FALSE, GO = FALSE)
    expect_true(is.list(gsf) || is.null(gsf))
})

test_that("group_specific_features with NULL groups uses all", {
    gsf <- group_specific_features(prop_tb, groups = NULL,
        genename = FALSE, GO = FALSE)
    expect_true(is.list(gsf) || is.null(gsf))
})

test_that("group_specific_features genename path returns data.frame", {
    skip_if_not_installed("org.Mm.eg.db")
    library(org.Mm.eg.db)

    gsf <- group_specific_features(prop_tb, groups = "untreated",
        filter_ratio = 0.5, genename = TRUE, GO = FALSE,
        OrgDb = org.Mm.eg.db, keytype = "SYMBOL")
    skip_if(is.null(gsf), "No features passed filter for genename test")
    expect_true("genename" %in% names(gsf))
    expect_s3_class(gsf$genename, "data.frame")
})

test_that("group_specific_features GO path returns plot", {
    skip_if_not_installed("org.Mm.eg.db")
    library(org.Mm.eg.db)

    gsf <- group_specific_features(prop_tb, groups = "untreated",
        filter_ratio = 0.5, genename = FALSE, GO = TRUE,
        OrgDb = org.Mm.eg.db, keytype = "SYMBOL")
    skip_if(is.null(gsf) || !"GO" %in% names(gsf),
        "No features passed filter for GO test")
    expect_s3_class(gsf$GO, "ggplot")
})

test_that("group_specific_features covers NA threshold branch", {
    se_list <- calc_feature_property(example_obj_merged_list, threshold = NULL)
    prop_tb_na <- summarise_feature_property(se_list)

    expect_message(
        gsf <- group_specific_features(prop_tb_na, groups = "untreated",
            filter_ratio = 0.1, genename = FALSE, GO = FALSE),
        "non-NA time points"
    )
    expect_true(is.list(gsf) || is.null(gsf))
})

test_that("group_specific_features errors on invalid group name", {
    expect_error(
        group_specific_features(prop_tb, groups = "nonexistent",
            genename = FALSE, GO = FALSE),
        "not found"
    )
})

test_that("group_specific_features errors on missing OrgDb when genename=TRUE", {
    expect_error(
        group_specific_features(prop_tb, groups = "untreated",
            genename = TRUE, GO = FALSE, OrgDb = NULL, keytype = NULL),
        "Both 'OrgDb' and 'keytype' must be provided"
    )
})

test_that("group_specific_features errors on multiple thresholds", {
    prop_tb_multi <- data.frame(
        Feature = c("A", "B"),
        Group = c("G1", "G2"),
        Exp_ratio = c(0.8, 0.8),
        Exp_threshold = c(0, 1),
        P_trend = c(0.5, 0.5),
        Max_FC = c(1, 1),
        stringsAsFactors = FALSE
    )
    expect_error(
        group_specific_features(prop_tb_multi, groups = c("G1", "G2"),
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

test_that("group_specific_features genename branch with simulated data", {
    skip_if_not_installed("org.Mm.eg.db")
    library(org.Mm.eg.db)

    prop_tb_sim <- data.frame(
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

    gsf <- group_specific_features(prop_tb_sim, groups = "A",
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

    prop_tb_sim <- data.frame(
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

    gsf <- group_specific_features(prop_tb_sim, groups = "A",
        filter_ratio = 0.5, genename = FALSE, GO = TRUE,
        OrgDb = org.Mm.eg.db, keytype = "SYMBOL",
        pvalueCutoff = 0.9, qvalueCutoff = 0.9,
        minGSSize = 1, simplify = FALSE)

    expect_true(is.list(gsf))
    expect_true(length(gsf$features) >= 6)
    if ("GO" %in% names(gsf)) {
        expect_true(is.list(gsf$GO) || inherits(gsf$GO, "ggplot"))
    }
})
