# Tests for as_DeeDeeExperiment and enrichment format compatibility

test_that("as_DeeDeeExperiment validates input", {
    expect_error(as_DeeDeeExperiment(NULL), "must be a list")
    expect_error(as_DeeDeeExperiment(list()), "Missing elements")
})

test_that("as_DeeDeeExperiment handles empty results and WGCNA metadata", {
    skip_if_not_installed("DeeDeeExperiment")
    data("tutorial_data")
    data("tutorial_sample_info")
    tide <- prepare_tide(tutorial_data, tutorial_sample_info,
        keep = "threshold", residual_threshold = 100)

    # Empty DE and enrichment
    dde <- as_DeeDeeExperiment(tide)
    expect_s4_class(dde, "DeeDeeExperiment")

    # WGCNA stored in metadata
    tide$WGCNA <- list(modules = "test")
    dde2 <- as_DeeDeeExperiment(tide)
    expect_equal(dde2@metadata$WGCNA, list(modules = "test"))
})

test_that("as_DeeDeeExperiment correctly renames DE columns", {
    skip_if_not_installed("DeeDeeExperiment")
    data("tutorial_data")
    data("tutorial_sample_info")
    tide <- prepare_tide(tutorial_data, tutorial_sample_info,
        keep = "threshold", residual_threshold = 100)

    # Run a simple DE
    tide$DE <- DE_between_group(tide$se, assay = "norm",
        filter = 1, trend = TRUE)

    dde <- as_DeeDeeExperiment(tide)
    expect_s4_class(dde, "DeeDeeExperiment")
})

test_that("enrichGO_list output format is compatible with DeeDeeExperiment", {
    skip_if_not_installed("org.Mm.eg.db")
    library(org.Mm.eg.db)

    data("example_net")
    example_module <- WGCNA_module(example_net, exclude_grey = TRUE) |>
        dplyr::filter(Module %in% c("1", "2"))

    go_list <- enrichGO_list(example_module, OrgDb = org.Mm.eg.db,
        universe = example_module$Feature,
        pvalueCutoff = 0.9, qvalueCutoff = 0.9,
        category = c("BP", "MF"), simplify = FALSE)

    # $all is a named list of data.frames
    enrich <- go_list$all
    expect_type(enrich, "list")
    expect_true(length(enrich) > 0)

    # Check non-null data.frames have clusterProfiler columns
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
    data("example_net")
    example_module <- WGCNA_module(example_net, exclude_grey = TRUE)

    enrichR_res <- enrichR_list(example_module,
        databases = "DSigDB",
        pvalueCutoff = 0.9)

    # Should be a named list of data.frames
    expect_type(enrichR_res, "list")
    expect_true("DSigDB" %in% names(enrichR_res))
    expect_true(is.data.frame(enrichR_res[["DSigDB"]]))

    # Has required columns
    df <- enrichR_res[["DSigDB"]]
    expect_true("Description" %in% colnames(df))
    expect_true("Cluster" %in% colnames(df))
})

test_that("as_DeeDeeExperiment handles enrichGO_list unmerged enrichResults", {
    skip_if_not_installed("DeeDeeExperiment")
    skip_if_not_installed("org.Mm.eg.db")
    library(org.Mm.eg.db)
    data("tutorial_data")
    data("tutorial_sample_info")
    tide <- prepare_tide(tutorial_data, tutorial_sample_info,
        keep = "threshold", residual_threshold = 100)

    data("example_net")
    example_module <- WGCNA_module(example_net, exclude_grey = TRUE)
    tide$enrichment$GO <- enrichGO_list(example_module, OrgDb = org.Mm.eg.db,
        universe = example_module$Feature,
        pvalueCutoff = 0.9, qvalueCutoff = 0.9,
        category = "BP", simplify = FALSE)

    dde <- as_DeeDeeExperiment(tide)
    expect_s4_class(dde, "DeeDeeExperiment")
    # Should have fea entries from unmerged enrichResult objects
    fea_names <- DeeDeeExperiment::getFEANames(dde)
    expect_true(length(fea_names) > 0)
})

test_that("as_DeeDeeExperiment handles multiple enrichment sources", {
    skip_if_not_installed("DeeDeeExperiment")
    skip_if_not_installed("org.Mm.eg.db")
    skip_if_not_installed("msigdbr")
    library(org.Mm.eg.db)
    data("tutorial_data")
    data("tutorial_sample_info")
    tide <- prepare_tide(tutorial_data, tutorial_sample_info,
        keep = "threshold", residual_threshold = 100)

    data("example_net")
    example_module <- WGCNA_module(example_net, exclude_grey = TRUE)

    tide$enrichment$GO <- enrichGO_list(example_module, OrgDb = org.Mm.eg.db,
        universe = example_module$Feature,
        pvalueCutoff = 0.9, qvalueCutoff = 0.9,
        category = "BP", simplify = FALSE)
    tide$enrichment$MSigDB <- enrich_msigdb(example_module,
        species = "Mus musculus", db_species = "MM",
        category = "MH", universe = example_module$Feature,
        pvalueCutoff = 0.9, minGSSize = 1, maxGSSize = 1000)

    dde <- as_DeeDeeExperiment(tide)
    expect_s4_class(dde, "DeeDeeExperiment")
    fea_names <- DeeDeeExperiment::getFEANames(dde)
    # Should have entries from both sources
    expect_true(length(fea_names) >= 2)
})

test_that("as_DeeDeeExperiment handles DE_between_time naming", {
    skip_if_not_installed("DeeDeeExperiment")
    data("tutorial_data")
    data("tutorial_sample_info")
    tide <- prepare_tide(tutorial_data, tutorial_sample_info,
        keep = "threshold", residual_threshold = 100)

    tide$DE <- DE_between_time(tide$se, assay = "norm", filter = 1)

    dde <- as_DeeDeeExperiment(tide)
    expect_s4_class(dde, "DeeDeeExperiment")
    dea_names <- DeeDeeExperiment::getDEANames(dde)
    # Time labels (t2-t0) should NOT get double T-prefix
    time_labels <- grep("Tt", dea_names, value = TRUE)
    expect_length(time_labels, 0)
    # Should have entries
    expect_true(length(dea_names) > 0)
})

test_that("as_DeeDeeExperiment handles nested DE (Case 2: DE_between_time)", {
    skip_if_not_installed("DeeDeeExperiment")
    library(DeeDeeExperiment)
    library(SingleCellExperiment)

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

    tide <- list(
        se = example_obj,
        DE = list(
            untreated = list(all_list = list(
                T2 = make_de_df(rownames(example_obj)[1:20]),
                T6 = make_de_df(rownames(example_obj)[21:40])
            )),
            IFNbeta = list(all_list = list(
                T2 = make_de_df(rownames(example_obj)[41:60])
            ))
        ),
        enrichment = list(go = list(), enrichr = list(), msigdb = list()),
        variance = list(),
        feature_property = list(),
        filter_summary = list(),
        WGCNA = list()
    )

    dde <- as_DeeDeeExperiment(tide)
    expect_s4_class(dde, "DeeDeeExperiment")
})
