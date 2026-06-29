# Tests for enrichment functions (GO rank, plot_GO)

test_that("enrichGO_rank works with variance decomposition output", {
    data("example_obj")
    se <- normalise_to_start(example_obj)
    vd <- decomp_variance(se, features = rownames(se)[1:20],
        fixed_effect_var = NULL, assay = "orig", core = 1)

    skip_if_not_installed("org.Mm.eg.db")
    library(org.Mm.eg.db)

    res <- enrichGO_rank(vd, gene_rank_by = "Time",
        OrgDb = org.Mm.eg.db, keyType = "SYMBOL",
        category = "BP", pvalueCutoff = 0.99)

    # May return NULL if no significant terms
    if (!is.null(res)) {
        expect_type(res, "list")
        expect_s4_class(res$BP, "gseaResult")
    }
})

test_that("enrichGO_rank validates gene_rank_by", {
    skip_if_not_installed("org.Mm.eg.db")
    library(org.Mm.eg.db)
    data("example_obj")
    se <- normalise_to_start(example_obj)
    vd <- decomp_variance(se, features = rownames(se)[1:10],
        fixed_effect_var = NULL, assay = "orig", core = 1)

    expect_error(
        enrichGO_rank(vd, gene_rank_by = "NotAColumn",
            OrgDb = org.Mm.eg.db, keyType = "SYMBOL"),
        "Gene ranking variable not found"
    )
})

test_that("enrichGO_rank requires exactly one ranking variable", {
    data("example_obj")
    se <- normalise_to_start(example_obj)
    vd <- decomp_variance(se, features = rownames(se)[1:10],
        fixed_effect_var = NULL, assay = "orig", core = 1)
    rank_table <- as.data.frame(vd)

    # Empty gene_rank_by should error (caught by .check_character)
    expect_error(
        enrichGO_rank(rank_table, gene_rank_by = character(0),
            OrgDb = NULL, keyType = "SYMBOL"),
        "non-empty"
    )
})

test_that("enrichGO_rank errors with multiple gene_rank_by values (line 64)", {
    skip_if_not_installed("org.Mm.eg.db")
    library(org.Mm.eg.db)
    data("example_obj")
    se <- normalise_to_start(example_obj)
    vd <- decomp_variance(se, features = rownames(se)[1:10],
        fixed_effect_var = NULL, assay = "orig", core = 1)
    rank_table <- as.data.frame(vd)

    # c("a","b") passes .check_character, then hits length>1 at lines 63-66
    expect_error(
        enrichGO_rank(rank_table, gene_rank_by = c("a", "b"),
            OrgDb = org.Mm.eg.db, keyType = "SYMBOL"),
        "Please specify one variable for ranking the genes."
    )
})

test_that("enrichGO_rank errors with NULL gene_rank_by (caught by .check_character)", {
    data("example_obj")
    se <- normalise_to_start(example_obj)
    vd <- decomp_variance(se, features = rownames(se)[1:10],
        fixed_effect_var = NULL, assay = "orig", core = 1)
    rank_table <- as.data.frame(vd)

    # NULL is caught by .check_character() before reaching lines 63-66
    expect_error(
        enrichGO_rank(rank_table, gene_rank_by = NULL,
            OrgDb = NULL, keyType = "SYMBOL"),
        "gene_rank_by.*must be a character"
    )
})

test_that("enrichGO_rank messages when go_rank_by column missing (lines 112-115)", {
    skip_if_not_installed("org.Mm.eg.db")
    library(org.Mm.eg.db)

    data("example_obj")
    se <- normalise_to_start(example_obj)
    vd <- decomp_variance(se, features = rownames(se)[1:100],
        fixed_effect_var = NULL, assay = "orig", core = 1)

    # Use a non-existent column name as go_rank_by; lenient cutoff ensures terms exist
    expect_message(
        res <- enrichGO_rank(vd, gene_rank_by = "Time",
            OrgDb = org.Mm.eg.db, keyType = "SYMBOL",
            category = "BP", pvalueCutoff = 0.99,
            go_rank_by = "nonexistent_column"),
        "GO term ranking variable not found in the result"
    )
    # Result should still be a valid list since terms were found
    if (!is.null(res)) {
        expect_type(res, "list")
    }
})

test_that("enrichGO_rank returns NULL with strict pvalueCutoff (lines 121-123)", {
    skip_if_not_installed("org.Mm.eg.db")
    library(org.Mm.eg.db)

    data("example_obj")
    se <- normalise_to_start(example_obj)
    vd <- decomp_variance(se, features = rownames(se)[1:10],
        fixed_effect_var = NULL, assay = "orig", core = 1)

    # Mock gseGO to return an empty gseaResult, triggering the NULL
    # return path at lines 120-123 and the associated message
    local_mocked_bindings(
        gseGO = function(...) methods::new("gseaResult"),
        .package = "clusterProfiler"
    )

    expect_message(
        res <- enrichGO_rank(vd, gene_rank_by = "Time",
            OrgDb = org.Mm.eg.db, keyType = "SYMBOL",
            category = "BP", pvalueCutoff = 0.05),
        "No significant GO terms found in any category."
    )
    expect_null(res)
})

test_that("plot_GO runs with GO enrichment output", {
    skip_if_not_installed("org.Mm.eg.db")
    library(org.Mm.eg.db)

    data("example_net")
    example_module <- WGCNA_module(example_net) |>
        dplyr::filter(Module %in% c("1", "2"))

    go_res <- enrichGO_list(example_module, OrgDb = org.Mm.eg.db,
        universe = example_module$Feature,
        pvalueCutoff = 0.5, qvalueCutoff = 0.5,
        category = "BP", simplify = FALSE)

    # Dotplot
    expect_error(
        plot_GO(go_res$all, plot_dotplot = TRUE,
            plot_emapplot = FALSE, plot_cnetplot = FALSE,
            showCategory_dotplot = 3),
        NA
    )
    # Cnetplot
    expect_error(
        plot_GO(go_res$all, plot_dotplot = FALSE,
            plot_cnetplot = TRUE, plot_emapplot = FALSE,
            showCategory_cnetplot = 3),
        NA
    )
    # Emaplot
    p <- plot_GO(go_res$all, plot_dotplot = FALSE,
        plot_emapplot = TRUE, plot_cnetplot = FALSE,
        showCategory_emapplot = 3)
    expect_type(p, "list")
    emap_names <- grep("^emapplot_", names(p), value = TRUE)
    if (length(emap_names) > 0)
        expect_s3_class(p[[emap_names[1]]], "ggplot")
})

test_that("plot_GO messages when no plot types specified", {
    expect_message(
        plot_GO(list(BP = data.frame()),
            plot_dotplot = FALSE, plot_cnetplot = FALSE,
            plot_emapplot = FALSE),
        "specify at least one"
    )
})

test_that("plot_GO errors when all/simplified wrapper passed", {
    expect_error(
        plot_GO(list(all = list(), simplified = list())),
        "all.*simplified"
    )
})

test_that("plot_GO errors when no valid categories found", {
    expect_error(
        plot_GO(list(x = 1), plot_dotplot = TRUE),
        "No valid"
    )
})

test_that("enrichR_list errors on mutual universe args", {
    skip_if_not_installed("enrichR")
    gene_list <- list(A = c("TP53", "BRCA1"), B = c("EGFR", "MYC"))
    expect_error(
        enrichR_list(gene_list,
            databases = "DSigDB",
            universe = c("TP53", "BRCA1", "EGFR", "MYC"),
            universe_list = list(A = c("TP53"), B = c("EGFR"))),
        "only one of universe"
    )
})

test_that("enrichR_list errors when no databases available", {
    skip_if_not_installed("enrichR")
    gene_list <- list(A = c("TP53", "BRCA1"))
    expect_error(
        enrichR_list(gene_list,
            databases = "NonExistentDB_XYZ123",
            universe = c("TP53", "BRCA1")),
        "None of the requested databases"
    )
})

# ---- .prepare_gene_list ----

test_that(".prepare_gene_list converts data.frame to named list", {
    df <- data.frame(Feature = c("A", "B", "C", "D"),
                     Module = factor(c("1", "1", "2", "2")),
                     stringsAsFactors = FALSE)
    res <- TiDEomics:::.prepare_gene_list(df)
    expect_type(res, "list")
    expect_named(res, c("1", "2"))
    expect_equal(res[["1"]], c("A", "B"))
    expect_equal(res[["2"]], c("C", "D"))
})

test_that(".prepare_gene_list passes through named list unchanged", {
    x <- list(A = c("TP53", "BRCA1"), B = c("EGFR"))
    res <- TiDEomics:::.prepare_gene_list(x)
    expect_equal(res, x)
})

test_that(".prepare_gene_list errors on unnamed list", {
    expect_error(
        TiDEomics:::.prepare_gene_list(list(c("A", "B"), c("C"))),
        "must be a named list"
    )
})

test_that(".prepare_gene_list errors on bad data.frame", {
    df <- data.frame(Gene = c("A", "B"), Cluster = c("1", "1"))
    expect_error(
        TiDEomics:::.prepare_gene_list(df),
        "must have 'Feature' and 'Module' columns"
    )
})

# ---- enrichGO_list data.frame input ----

test_that("enrichGO_list accepts data.frame input", {
    skip_if_not_installed("org.Mm.eg.db")
    library(org.Mm.eg.db)
    data("example_net")
    example_module <- WGCNA_module(example_net) |>
        dplyr::filter(Module %in% c("1", "2"))

    expect_error(
        enrichGO_list(example_module, OrgDb = org.Mm.eg.db,
            category = "BP", pvalueCutoff = 0.5, qvalueCutoff = 0.5),
        NA
    )
})

test_that("enrichGO_list simplify=TRUE returns unmerged_simplified", {
    skip_if_not_installed("org.Mm.eg.db")
    library(org.Mm.eg.db)
    data("example_net")
    example_module <- WGCNA_module(example_net) |>
        dplyr::filter(Module %in% c("1", "2"))

    go_simplified <- enrichGO_list(example_module, OrgDb = org.Mm.eg.db,
        category = "BP", pvalueCutoff = 0.5, qvalueCutoff = 0.5,
        simplify = TRUE)

    expect_true("all" %in% names(go_simplified))
    expect_true("simplified" %in% names(go_simplified))
    expect_true("unmerged_all" %in% names(go_simplified))
    expect_true("unmerged_simplified" %in% names(go_simplified))
    # Check unmerged_simplified structure
    expect_type(go_simplified$unmerged_simplified, "list")
    expect_true("BP" %in% names(go_simplified$unmerged_simplified))
})

# ---- enrichR_list site parameter ----

test_that("enrichR_list validates site parameter", {
    gene_list <- list(A = c("TP53", "BRCA1"))
    expect_error(
        enrichR_list(gene_list, databases = "DSigDB",
            site = "NotASite")
    )
})

test_that("enrichR_list validates universe_list names match", {
    skip_if_not_installed("enrichR")
    gene_list <- list(A = c("TP53", "BRCA1"), B = c("EGFR"))
    expect_error(
        enrichR_list(gene_list, databases = "DSigDB",
            universe_list = list(X = c("TP53"))),
        "must match"
    )
})


test_that("enrichR_list rejects non-list, non-data.frame input", {
    expect_error(
        enrichR_list(gene_list = c("TP53", "BRCA1")),
        "must be a data.frame"
    )
})

# ---- enrichR_list coverage tests ----

test_that("enrichR_list does not print background message when universe is provided", {
    skip_if_not_installed("enrichR")

    # Mock enrichR functions to avoid internet dependency
    local_mocked_bindings(
        setEnrichrSite = function(site) NULL,
        listEnrichrDbs = function() data.frame(libraryName = "DSigDB",
            stringsAsFactors = FALSE),
        enrichr = function(genes, databases, ...) {
            stats::setNames(lapply(databases, function(db) {
                data.frame(
                    Term = "Test pathway",
                    Overlap = "1/100",
                    P.value = 0.001,
                    Adjusted.P.value = 0.01,
                    Old.P.value = 0.001,
                    Old.Adjusted.P.value = 0.01,
                    Z.score = 2.0,
                    Combined.Score = 10.0,
                    Genes = paste(genes, collapse = ","),
                    stringsAsFactors = FALSE
                )
            }), databases)
        },
        .package = "enrichR"
    )

    gene_list <- list(A = c("TP53", "BRCA1"))

    msgs <- capture_messages(
        enrichR_list(gene_list, databases = "DSigDB",
            universe = c("TP53", "BRCA1", "EGFR", "MYC"))
    )

    # "Background genes not specified" should NOT appear because universe
    # is provided (non-NULL), so the else-branch at line 115 is taken
    # instead of the message-branch at line 111.
    expect_false(any(grepl("Background genes not specified", msgs)))
})

test_that("enrichR_list errors when universe_list is unnamed", {
    skip_if_not_installed("enrichR")

    # Need to mock enrichR functions to reach the universe_list check
    # which is after the database availability check
    local_mocked_bindings(
        setEnrichrSite = function(site) NULL,
        listEnrichrDbs = function() data.frame(libraryName = "DSigDB",
            stringsAsFactors = FALSE),
        .package = "enrichR"
    )

    gene_list <- list(A = c("TP53", "BRCA1"))
    expect_error(
        enrichR_list(gene_list, databases = "DSigDB",
            universe_list = list(c("TP53", "BRCA1"))),
        "must be a named list"
    )
})

test_that("enrichR_list skips empty gene lists during enrichment", {
    skip_if_not_installed("enrichR")

    local_mocked_bindings(
        setEnrichrSite = function(site) NULL,
        listEnrichrDbs = function() data.frame(libraryName = "DSigDB",
            stringsAsFactors = FALSE),
        enrichr = function(genes, databases, ...) {
            stats::setNames(lapply(databases, function(db) {
                data.frame(
                    Term = "Test pathway",
                    Overlap = "1/100",
                    P.value = 0.001,
                    Adjusted.P.value = 0.01,
                    Old.P.value = 0.001,
                    Old.Adjusted.P.value = 0.01,
                    Z.score = 2.0,
                    Combined.Score = 10.0,
                    Genes = paste(genes, collapse = ","),
                    stringsAsFactors = FALSE
                )
            }), databases)
        },
        .package = "enrichR"
    )

    # Group B is empty; the function should skip it via the
    # `if (length(gene_list[[i]]) == 0) next` guard at line 131
    # and only process group A.
    gene_list <- list(A = c("TP53", "BRCA1"), B = character(0))
    result <- suppressMessages(
        enrichR_list(gene_list, databases = "DSigDB")
    )

    expect_type(result, "list")
    expect_true("DSigDB" %in% names(result))
    expect_gt(nrow(result[["DSigDB"]]), 0)
})

# Test 4: length(db_tables) == 0 -> "No enriched terms found."
# This is hard to trigger without a real web API call, as it requires
# enrichR::enrichr() to return empty results for all databases.

test_that("enrichR_list messages when no enriched terms found (db_tables empty)", {
    skip_if_not_installed("enrichR")

    # Mock enrichR to return 0-row data.frames so nrow(enr[[db]]) == 0
    # is TRUE and the skip at line 138 is taken, leaving db_tables empty.
    local_mocked_bindings(
        setEnrichrSite = function(site) NULL,
        listEnrichrDbs = function() data.frame(libraryName = "DSigDB",
            stringsAsFactors = FALSE),
        enrichr = function(genes, databases, ...) {
            stats::setNames(lapply(databases, function(db) {
                data.frame(
                    Term = character(0),
                    Overlap = character(0),
                    P.value = numeric(0),
                    Adjusted.P.value = numeric(0),
                    Old.P.value = numeric(0),
                    Old.Adjusted.P.value = numeric(0),
                    Z.score = numeric(0),
                    Combined.Score = numeric(0),
                    Genes = character(0),
                    stringsAsFactors = FALSE
                )
            }), databases)
        },
        .package = "enrichR"
    )

    gene_list <- list(A = c("TP53", "BRCA1"))

    expect_message(
        result <- enrichR_list(gene_list, databases = "DSigDB"),
        "No enriched terms found"
    )
    expect_type(result, "list")
    expect_length(result, 0)
})

# ---- enrichGO_list edge branches ----

test_that("enrichGO_list defaults category to all three with message", {
    skip_if_not_installed("org.Mm.eg.db")
    library(org.Mm.eg.db)
    data("example_net")
    example_module <- WGCNA_module(example_net) |>
        dplyr::filter(Module %in% c("1", "2"))

    go_res <- enrichGO_list(example_module, OrgDb = org.Mm.eg.db,
        pvalueCutoff = 0.9, qvalueCutoff = 0.9, category = NULL)
    expect_true(all(c("BP", "MF", "CC") %in% names(go_res$all)))
})

test_that("enrichGO_list errors when both universe and universe_list given", {
    skip_if_not_installed("org.Mm.eg.db")
    library(org.Mm.eg.db)
    data("example_net")
    example_module <- WGCNA_module(example_net) |>
        dplyr::filter(Module %in% c("1", "2"))

    gene_list <- list(a = example_module$Feature[1:10])
    expect_error(
        enrichGO_list(gene_list, OrgDb = org.Mm.eg.db,
            category = "BP",
            universe = example_module$Feature,
            universe_list = list(a = example_module$Feature),
            pvalueCutoff = 0.5, qvalueCutoff = 0.5),
        "only one of universe"
    )
})

test_that("enrichGO_list errors on unnamed universe_list", {
    skip_if_not_installed("org.Mm.eg.db")
    library(org.Mm.eg.db)
    data("example_net")
    example_module <- WGCNA_module(example_net) |>
        dplyr::filter(Module %in% c("1", "2"))

    gene_list <- list(a = example_module$Feature[1:10])
    expect_error(
        enrichGO_list(gene_list, OrgDb = org.Mm.eg.db,
            category = "BP",
            universe_list = list(example_module$Feature),
            pvalueCutoff = 0.5, qvalueCutoff = 0.5),
        "universe_list must be a named list"
    )
})

test_that("enrichGO_list errors when gene_list and universe_list names mismatch", {
    skip_if_not_installed("org.Mm.eg.db")
    library(org.Mm.eg.db)
    data("example_net")
    example_module <- WGCNA_module(example_net) |>
        dplyr::filter(Module %in% c("1", "2"))

    gene_list <- list(a = example_module$Feature[1:10],
                      b = example_module$Feature[11:20])
    expect_error(
        enrichGO_list(gene_list, OrgDb = org.Mm.eg.db,
            category = "BP",
            universe_list = list(a = example_module$Feature,
                                 c = example_module$Feature),
            pvalueCutoff = 0.5, qvalueCutoff = 0.5),
        "Names of gene_list and universe_list must match"
    )
})

test_that("enrichGO_list skips empty gene list element", {
    skip_if_not_installed("org.Mm.eg.db")
    library(org.Mm.eg.db)
    data("example_net")
    example_module <- WGCNA_module(example_net) |>
        dplyr::filter(Module %in% c("1", "2"))

    gene_list <- list(empty = character(0),
                      a = example_module$Feature[1:10])
    expect_message(
        go_res <- enrichGO_list(gene_list, OrgDb = org.Mm.eg.db,
            category = "BP",
            pvalueCutoff = 0.5, qvalueCutoff = 0.5),
        "empty"
    )
    expect_false("empty" %in% names(go_res$unmerged_all$BP))
})
