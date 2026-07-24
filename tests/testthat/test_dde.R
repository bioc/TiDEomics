# Tests for flatten_DE, flatten_enrich, and enrichment format compatibility

# ---- Shared setup ----
data("example_net")
data("example_go")
data("tutorial_data")
data("tutorial_sample_info")
example_module_all <- WGCNA_module(example_net, exclude_grey = TRUE)
example_module_filt <- example_module_all |>
    dplyr::filter(Module %in% c("1", "2"))

tide <- prepare_tide(tutorial_data, tutorial_sample_info,
    keep = "threshold", residual_threshold = 100)

if (requireNamespace("msigdbr", quietly = TRUE)) {
    msigdb_res_dde <- enrich_msigdb(example_module_all,
        species = "Mus musculus", db_species = "MM",
        category = "MH", universe = example_module_all$Feature,
        pvalueCutoff = 0.9, minGSSize = 1, maxGSSize = 1000)
}

# ---- flatten_DE ----

test_that("flatten_DE validates input", {
    expect_error(flatten_DE(NULL), "must be a list")
    expect_equal(flatten_DE(list()), list())
})

test_that("flatten_DE handles empty DE within prepare_tide output", {
    result <- flatten_DE(tide$DE)
    expect_type(result, "list")
})

test_that("flatten_DE flattens and renames DE_between_group output", {
    tide_de <- tide
    tide_de$DE <- DE_between_group(tide$se, assay = "norm",
        filter = 1, trend = TRUE)

    result <- flatten_DE(tide_de$DE)
    expect_type(result, "list")
    expect_true(length(result) > 0)
    for (nm in names(result)) {
        df <- result[[nm]]
        expect_s3_class(df, "data.frame")
        expect_true("log2FoldChange" %in% colnames(df))
        expect_true("pvalue" %in% colnames(df))
        expect_true("padj" %in% colnames(df))
        expect_false("logFC" %in% colnames(df))
        expect_false("P.Value" %in% colnames(df))
        expect_false("adj.P.Val" %in% colnames(df))
    }
})

test_that("flatten_DE handles DE_between_time naming without double T prefix", {
    tide_de <- tide
    tide_de$DE <- DE_between_time(tide$se, assay = "norm", filter = 1)

    result <- flatten_DE(tide_de$DE)
    expect_type(result, "list")
    expect_true(length(result) > 0)
    time_labels <- grep("Tt", names(result), value = TRUE)
    expect_length(time_labels, 0)
})

test_that("flatten_DE handles nested DE structure", {
    data(example_obj)
    example_obj <- normalise_to_start(example_obj)

    make_de_df <- function(feats) {
        data.frame(
            Feature = feats,
            logFC = rnorm(length(feats), 0, 1),
            P.Value = runif(length(feats), 0, 0.05),
            adj.P.Val = runif(length(feats), 0, 0.1),
            stringsAsFactors = FALSE
        )
    }

    de_list <- list(
        untreated = list(all_list = list(
            T2 = make_de_df(rownames(example_obj)[1:20]),
            T6 = make_de_df(rownames(example_obj)[21:40])
        )),
        IFNbeta = list(all_list = list(
            T2 = make_de_df(rownames(example_obj)[41:60])
        ))
    )

    result <- flatten_DE(de_list)
    expect_type(result, "list")
    expect_true(length(result) >= 3)
    for (nm in names(result)) {
        expect_true("log2FoldChange" %in% colnames(result[[nm]]))
    }
})

# ---- flatten_enrich ----

test_that("flatten_enrich validates input", {
    expect_error(flatten_enrich(NULL), "must be a list")
    expect_equal(flatten_enrich(list()), list())
})

test_that("flatten_enrich handles enrichGO_list with DDE-compatible output", {
    enrich_flat <- flatten_enrich(example_go)
    expect_type(enrich_flat, "list")
    expect_true(length(enrich_flat) > 0)
    expect_true(all(grepl("^clusterProfiler_", names(enrich_flat))))
    for (nm in names(enrich_flat)) {
        expect_true(is.data.frame(enrich_flat[[nm]]) ||
            inherits(enrich_flat[[nm]], "enrichResult"))
    }
})

test_that("flatten_enrich handles multiple enrichment sources and flat msigdb", {
    skip_if_not_installed("msigdbr")

    # Already-flat enrich_msigdb
    enrich_flat_flat <- flatten_enrich(msigdb_res_dde)
    expect_type(enrich_flat_flat, "list")
    expect_true(length(enrich_flat_flat) > 0)
    expect_true(all(grepl("^clusterProfiler_", names(enrich_flat_flat))))

    # Multiple sources (GO + MSigDB)
    enrich_list <- list(
        GO = example_go,
        MSigDB = msigdb_res_dde
    )
    enrich_flat <- flatten_enrich(enrich_list)
    expect_type(enrich_flat, "list")
    expect_true(length(enrich_flat) >= 2)
    has_prefix <- grepl("^clusterProfiler_", names(enrich_flat))
    expect_true(all(has_prefix))
})

test_that("flatten_enrich handles already-flat enrichR_list", {
    skip_if_not_installed("enrichR")

    enrichR_res <- enrichR_list(example_module_filt,
        databases = "DSigDB",
        pvalueCutoff = 0.9)

    enrich_flat <- flatten_enrich(enrichR_res)
    expect_type(enrich_flat, "list")
    expect_true("enrichr_DSigDB" %in% names(enrich_flat))
    df <- enrich_flat[["enrichr_DSigDB"]]
    expect_true("Term" %in% colnames(df))
    expect_true("Adjusted.P.value" %in% colnames(df))
})

# ---- Enrichment output format compatibility ----

test_that("enrichGO_list output format is compatible with DeeDeeExperiment", {
    enrich <- example_go$all
    expect_type(enrich, "list")
    expect_true(length(enrich) > 0)
    for (nm in names(enrich)) {
        df <- enrich[[nm]]
        if (!is.null(df) && is.data.frame(df) && nrow(df) > 0) {
            expect_true("Description" %in% colnames(df))
            expect_true("p.adjust" %in% colnames(df))
            expect_true("Cluster" %in% colnames(df))
        }
    }
})

test_that("enrichR_list output format is compatible with DeeDeeExperiment", {
    skip_if_not_installed("enrichR")

    enrichR_res <- enrichR_list(example_module_filt,
        databases = "DSigDB",
        pvalueCutoff = 0.9)

    expect_type(enrichR_res, "list")
    expect_true("DSigDB" %in% names(enrichR_res))
    expect_true(is.data.frame(enrichR_res[["DSigDB"]]))

    df <- enrichR_res[["DSigDB"]]
    expect_true("Description" %in% colnames(df))
    expect_true("Cluster" %in% colnames(df))
})
