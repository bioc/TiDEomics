# Tests for enrich_msigdb

# ---- Shared setup ----
data("example_net")
example_module <- WGCNA_module(example_net, exclude_grey = TRUE)
example_module_filt <- WGCNA_module(example_net) |>
    dplyr::filter(Module %in% c("1", "2"))

if (requireNamespace("msigdbr", quietly = TRUE)) {
    msigdb_res <- enrich_msigdb(example_module, species = "Mus musculus",
        db_species = "MM", category = "MH",
        universe = example_module$Feature,
        minGSSize = 1, maxGSSize = 1000, pvalueCutoff = 0.99)
}

test_that("enrich_msigdb validates gene_list format", {
    skip_if_not_installed("msigdbr")

    # Data.frame input (from WGCNA_module)
    expect_error(
        enrich_msigdb(example_module, species = "Mus musculus", db_species = "MM",
            category = "MH", universe = example_module$Feature),
        NA
    )

    # Named list input
    gene_list <- list(A = example_module$Feature[1:10],
                      B = example_module$Feature[11:20])
    expect_error(
        enrich_msigdb(gene_list, category = "MH", species = "Mus musculus", db_species = "MM",
            universe = example_module$Feature),
        NA
    )
})

test_that("enrich_msigdb validates parameters (no msigdbr required)", {
    # pvalueCutoff validation
    expect_error(
        enrich_msigdb(example_module, universe = example_module$Feature,
            pvalueCutoff = -0.1),
        "between 0 and 1"
    )
    expect_error(
        enrich_msigdb(example_module, universe = example_module$Feature,
            pvalueCutoff = 1.5),
        "between 0 and 1"
    )
    # minGSSize / maxGSSize validation
    expect_error(
        enrich_msigdb(example_module, universe = example_module$Feature,
            minGSSize = -1),
        "positive"
    )
    expect_error(
        enrich_msigdb(example_module, universe = example_module$Feature,
            minGSSize = 0),
        "positive"
    )
    expect_error(
        enrich_msigdb(example_module, universe = example_module$Feature,
            maxGSSize = 0),
        "positive"
    )
    # pAdjustMethod must be character
    expect_error(
        enrich_msigdb(example_module, universe = example_module$Feature,
            pAdjustMethod = 1),
        "character"
    )
    # species must be character
    expect_error(
        enrich_msigdb(example_module, universe = example_module$Feature,
            species = 123),
        "character"
    )
    # db_species must be character
    expect_error(
        enrich_msigdb(example_module, universe = example_module$Feature,
            db_species = 123),
        "character"
    )
    # db_species rejects invalid values
    expect_error(
        enrich_msigdb(example_module, universe = example_module$Feature,
            db_species = "XX"),
        "should be one of"
    )
    # name must be character
    expect_error(
        enrich_msigdb(example_module_filt, universe = example_module_filt$Feature,
            category = "MH", species = "Mus musculus", db_species = "MM",
            minGSSize = 1, pvalueCutoff = 0.99, name = 42),
        "character"
    )
})

test_that("enrich_msigdb requires category or gene_sets", {
    expect_error(
        enrich_msigdb(example_module,
            universe = example_module$Feature),
        "Either .category. or .gene_sets. must be specified"
    )
})

test_that("enrich_msigdb returns valid output with expected columns", {
    skip_if_not_installed("msigdbr")

    # Returns a named list
    expect_type(msigdb_res, "list")
    expect_true(length(msigdb_res) > 0)
    expect_true(length(names(msigdb_res)) > 0)

    df <- msigdb_res[[1]]
    if (is.data.frame(df) && nrow(df) > 0) {
        # Expected columns
        expect_true("Description" %in% colnames(df))
        expect_true("Cluster" %in% colnames(df))

        expected <- c("Cluster", "ID", "Description", "GeneRatio",
                      "BgRatio", "pvalue", "p.adjust", "Count", "geneID")
        expect_true(all(expected %in% colnames(df)))
        expect_true(all(nchar(df$geneID) > 0))
        expect_equal(df$Count, lengths(strsplit(df$geneID, ";")))
    }
})

# ---- .prepare_gene_list and .hypergeometric_test ----

test_that(".prepare_gene_list converts data.frame to named list", {
    df <- data.frame(Feature = c("A", "B", "C", "D"),
                     Module  = c("M1", "M1", "M2", "M2"))
    res <- TiDEomics:::.prepare_gene_list(df)
    expect_type(res, "list")
    expect_equal(names(res), c("M1", "M2"))
    expect_equal(res$M1, c("A", "B"))
    expect_equal(res$M2, c("C", "D"))
})

test_that(".prepare_gene_list passes through named list", {
    x <- list(G1 = c("a", "b"), G2 = c("c"))
    res <- TiDEomics:::.prepare_gene_list(x)
    expect_equal(res, x)
})

test_that(".prepare_gene_list errors on unnamed list", {
    expect_error(
        TiDEomics:::.prepare_gene_list(list(c("a", "b"))),
        "must be a named list"
    )
})

test_that(".prepare_gene_list errors on bad data.frame", {
    expect_error(
        TiDEomics:::.prepare_gene_list(data.frame(x = 1:3)),
        "must have .Feature. and .Module. columns"
    )
})

test_that(".hypergeometric_test returns correct structure", {
    gene_sets <- list(GS1 = c("A", "B", "C", "D", "E"),
                      GS2 = c("F", "G", "H"))
    universe <- LETTERS[1:10]
    query <- c("A", "B", "F")

    res <- TiDEomics:::.hypergeometric_test(query, gene_sets, universe)
    expect_s3_class(res, "data.frame")
    expect_true(all(c("Term", "N", "m", "k", "q", "pvalue") %in% colnames(res)))
    expect_true(nrow(res) >= 1)
})

test_that(".hypergeometric_test respects min_overlap", {
    gene_sets <- list(GS1 = c("A", "B"), GS2 = c("C"))
    universe <- LETTERS[1:10]
    query <- c("A", "C")

    res <- TiDEomics:::.hypergeometric_test(query, gene_sets, universe,
        min_overlap = 2L)
    expect_null(res)
})

test_that(".hypergeometric_test returns NULL for empty query-universe intersect", {
    gene_sets <- list(GS1 = c("A", "B"))
    universe <- LETTERS[1:5]
    query <- c("X", "Y", "Z")
    res <- TiDEomics:::.hypergeometric_test(query, gene_sets, universe)
    expect_null(res)
})

test_that(".hypergeometric_test returns NULL for no passing gene sets", {
    gene_sets <- list()
    universe <- LETTERS[1:10]
    query <- c("A", "B")
    res <- TiDEomics:::.hypergeometric_test(query, gene_sets, universe)
    expect_null(res)
})

test_that("enrich_msigdb handles gene_sets parameter", {
    skip_if_not_installed("msigdbr")

    res <- enrich_msigdb(example_module, species = "Mus musculus", db_species = "MM",
        gene_sets = "HALLMARK_APOPTOSIS",
        universe = example_module$Feature)
    expect_type(res, "list")
})

test_that("enrich_msigdb skips empty gene lists", {
    skip_if_not_installed("msigdbr")

    gene_list <- list(M1 = character(0),
                      M2 = example_module$Feature[example_module$Module == "2"])
    res <- enrich_msigdb(gene_list, species = "Mus musculus",
        db_species = "MM", category = "MH",
        universe = example_module$Feature)
    expect_type(res, "list")
})

test_that("enrich_msigdb warns when no gene sets match size range", {
    skip_if_not_installed("msigdbr")

    expect_warning(
        res <- enrich_msigdb(example_module, species = "Mus musculus",
            db_species = "MM", category = "MH",
            universe = example_module$Feature,
            minGSSize = 1e6, maxGSSize = 1e7),
        "No gene sets with size"
    )
    expect_null(res[[1]])
})

test_that("enrich_msigdb warns on missing gene_sets", {
    skip_if_not_installed("msigdbr")

    expect_warning(
        enrich_msigdb(example_module,
            species = "Mus musculus", db_species = "MM",
            gene_sets = c("HALLMARK_APOPTOSIS", "NONEXISTENT_SET_XYZ"),
            universe = example_module$Feature),
        "Gene set")
})

test_that("enrich_msigdb errors when gene_sets finds nothing", {
    skip_if_not_installed("msigdbr")

    expect_error(
        enrich_msigdb(example_module,
            species = "Mus musculus", db_species = "MM",
            gene_sets = "NONEXISTENT_SET_XYZ123",
            universe = example_module$Feature),
        "None of the requested gene sets")
})

test_that("enrich_msigdb errors on unknown collection", {
    skip_if_not_installed("msigdbr")

    expect_error(
        enrich_msigdb(example_module,
            species = "Mus musculus", db_species = "MM",
            category = "C5", subcategory = "NonExistent_XYZ",
            universe = example_module$Feature,
            pvalueCutoff = 0.99, minGSSize = 1),
        "Unknown collection"
    )
})

test_that("enrich_msigdb warns when no terms pass pvalue cutoff", {
    skip_if_not_installed("msigdbr")

    expect_warning(
        enrich_msigdb(example_module_filt, universe = example_module_filt$Feature,
            category = "MH", species = "Mus musculus", db_species = "MM",
            minGSSize = 1, pvalueCutoff = 1e-10),
        "No enriched terms pass cutoffs"
    )
})
