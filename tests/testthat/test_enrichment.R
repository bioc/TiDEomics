# Tests for enrichment functions (GO rank, plot_GO)

data("example_net")
data("example_go")

test_that("enrichGO_rank works and validates parameters", {
    data("example_obj")
    se <- normalise_to_start(example_obj)

    vd_20 <- decomp_variance(se, features = rownames(se)[1:20],
        fixed_effect_var = NULL, assay = "orig", core = 1)
    vd_10 <- decomp_variance(se, features = rownames(se)[1:10],
        fixed_effect_var = NULL, assay = "orig", core = 1)
    rank_table_10 <- as.data.frame(vd_10)

    skip_if_not_installed("org.Mm.eg.db")
    library(org.Mm.eg.db)

    # Basic output
    res <- enrichGO_rank(vd_20, gene_rank_by = "Time",
        OrgDb = org.Mm.eg.db, keyType = "SYMBOL",
        category = "BP", pvalueCutoff = 0.99)
    if (!is.null(res)) {
        expect_type(res, "list")
        expect_s4_class(res$BP, "gseaResult")
    }

    # Invalid gene_rank_by column
    expect_error(
        enrichGO_rank(vd_10, gene_rank_by = "NotAColumn",
            OrgDb = org.Mm.eg.db, keyType = "SYMBOL"),
        "Gene ranking variable not found"
    )

    # Empty gene_rank_by
    expect_error(
        enrichGO_rank(rank_table_10, gene_rank_by = character(0),
            OrgDb = NULL, keyType = "SYMBOL"),
        "non-empty"
    )

    # Multiple gene_rank_by values
    expect_error(
        enrichGO_rank(vd_10, gene_rank_by = c("a", "b"),
            OrgDb = org.Mm.eg.db, keyType = "SYMBOL"),
        "Please specify one variable for ranking the genes."
    )

    # NULL gene_rank_by
    expect_error(
        enrichGO_rank(rank_table_10, gene_rank_by = NULL,
            OrgDb = NULL, keyType = "SYMBOL"),
        "gene_rank_by.*must be a character"
    )
})

test_that("enrichGO_rank messages on invalid go_rank_by and returns NULL on strict cutoff", {
    skip_if_not_installed("org.Mm.eg.db")
    library(org.Mm.eg.db)

    data("example_obj")
    se <- normalise_to_start(example_obj)

    vd <- decomp_variance(se, features = rownames(se)[1:100],
        fixed_effect_var = NULL, assay = "orig", core = 1)
    vd_small <- decomp_variance(se, features = rownames(se)[1:10],
        fixed_effect_var = NULL, assay = "orig", core = 1)

    # go_rank_by not found in result
    expect_message(
        res <- enrichGO_rank(vd, gene_rank_by = "Time",
            OrgDb = org.Mm.eg.db, keyType = "SYMBOL",
            category = "BP", pvalueCutoff = 0.99,
            go_rank_by = "nonexistent_column"),
        "GO term ranking variable not found in the result"
    )
    if (!is.null(res)) {
        expect_type(res, "list")
    }

    # Mock gseGO to test NULL return with strict pvalueCutoff
    local_mocked_bindings(
        gseGO = function(...) methods::new("gseaResult"),
        .package = "clusterProfiler"
    )

    expect_message(
        res <- enrichGO_rank(vd_small, gene_rank_by = "Time",
            OrgDb = org.Mm.eg.db, keyType = "SYMBOL",
            category = "BP", pvalueCutoff = 0.05),
        "No significant GO terms found in any category."
    )
    expect_null(res)
})

test_that("enrichGO_list and plot_GO work", {
    skip_if_not_installed("org.Mm.eg.db")
    library(org.Mm.eg.db)

    example_module <- WGCNA_module(example_net) |>
        dplyr::filter(Module %in% c("1", "2")) |>
        dplyr::slice_head(n = 25, by = Module)

    # Dotplot
    expect_error(
        plot_GO(example_go$all, plot_dotplot = TRUE,
            plot_emapplot = FALSE, plot_cnetplot = FALSE,
            showCategory_dotplot = 3),
        NA
    )
    # Cnetplot
    expect_error(
        plot_GO(example_go$all, plot_dotplot = FALSE,
            plot_cnetplot = TRUE, plot_emapplot = FALSE,
            showCategory_cnetplot = 3),
        NA
    )
    # Emaplot
    p <- plot_GO(example_go$all, plot_dotplot = FALSE,
        plot_emapplot = TRUE, plot_cnetplot = FALSE,
        showCategory_emapplot = 3)
    expect_type(p, "list")
    emap_names <- grep("^emapplot_", names(p), value = TRUE)
    if (length(emap_names) > 0)
        expect_s3_class(p[[emap_names[1]]], "ggplot")

    # enrichGO_list accepts data.frame input (no universe)
    expect_error(
        enrichGO_list(example_module, OrgDb = org.Mm.eg.db,
            category = "BP", pvalueCutoff = 0.5, qvalueCutoff = 0.5),
        NA
    )

    # enrichGO_list supports multiple categories with simplify
    example_module <- WGCNA_module(example_net) |>
        dplyr::filter(Module %in% c("1", "2"))
    go_multi <- enrichGO_list(example_module,
        OrgDb = org.Mm.eg.db,
        universe = WGCNA_module(example_net, exclude_grey = FALSE)$Feature,
        pvalueCutoff = 0.5, qvalueCutoff = 0.5,
        category = c("BP", "CC"), simplify = TRUE)
    expect_true(all(c("BP", "CC") %in% names(go_multi$all)))
    expect_true(all(c("BP", "CC") %in% names(go_multi$simplified)))
    expect_true("unmerged_all" %in% names(go_multi))
    expect_true("unmerged_simplified" %in% names(go_multi))
    expect_type(go_multi$unmerged_simplified, "list")
    expect_true("BP" %in% names(go_multi$unmerged_simplified))
    expect_true("CC" %in% names(go_multi$unmerged_simplified))
})

test_that("enrichGO_list validates universe and universe_list", {
    skip_if_not_installed("org.Mm.eg.db")
    library(org.Mm.eg.db)

    data("example_net")
    example_module <- WGCNA_module(example_net) |>
        dplyr::filter(Module %in% c("1", "2")) |>
        dplyr::slice_head(n = 25, by = Module)
    gene_list <- list(a = example_module$Feature[1:10],
                      b = example_module$Feature[11:20])

    # Both universe and universe_list given
    expect_error(
        enrichGO_list(gene_list, OrgDb = org.Mm.eg.db,
            category = "BP",
            universe = example_module$Feature,
            universe_list = list(a = example_module$Feature),
            pvalueCutoff = 0.5, qvalueCutoff = 0.5),
        "only one of universe"
    )

    # Unnamed universe_list
    expect_error(
        enrichGO_list(gene_list, OrgDb = org.Mm.eg.db,
            category = "BP",
            universe_list = list(example_module$Feature),
            pvalueCutoff = 0.5, qvalueCutoff = 0.5),
        "universe_list must be a named list"
    )

    # Name mismatch between gene_list and universe_list
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
        dplyr::filter(Module %in% c("1", "2")) |>
        dplyr::slice_head(n = 25, by = Module)

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

test_that("plot_GO messages and errors on invalid input", {
    expect_message(
        plot_GO(list(BP = data.frame()),
            plot_dotplot = FALSE, plot_cnetplot = FALSE,
            plot_emapplot = FALSE),
        "specify at least one"
    )

    expect_error(
        plot_GO(list(all = list(), simplified = list())),
        "all.*simplified"
    )

    expect_error(
        plot_GO(list(x = 1), plot_dotplot = TRUE),
        "No valid"
    )
})

test_that("enrichR_list validates parameters", {
    skip_if_not_installed("enrichR")

    gene_list <- list(A = c("TP53", "BRCA1"), B = c("EGFR", "MYC"))

    # Mutual universe args
    expect_error(
        enrichR_list(gene_list,
            databases = "DSigDB",
            universe = c("TP53", "BRCA1", "EGFR", "MYC"),
            universe_list = list(A = c("TP53"), B = c("EGFR"))),
        "only one of universe"
    )

    # No databases available
    expect_error(
        enrichR_list(list(A = c("TP53", "BRCA1")),
            databases = "NonExistentDB_XYZ123",
            universe = c("TP53", "BRCA1")),
        "None of the requested databases"
    )

    # Site parameter validation
    gene_list2 <- list(A = c("TP53", "BRCA1"))
    expect_error(
        enrichR_list(gene_list2, databases = "DSigDB",
            site = "NotASite")
    )

    # universe_list names mismatch
    gene_list3 <- list(A = c("TP53", "BRCA1"), B = c("EGFR"))
    expect_error(
        enrichR_list(gene_list3, databases = "DSigDB",
            universe_list = list(X = c("TP53"))),
        "must match"
    )

    # Non-list, non-data.frame input
    expect_error(
        enrichR_list(gene_list = c("TP53", "BRCA1")),
        "must be a data.frame"
    )
})

test_that(".prepare_gene_list works", {
    # data.frame to named list
    df <- data.frame(Feature = c("A", "B", "C", "D"),
                     Module = factor(c("1", "1", "2", "2")),
                     stringsAsFactors = FALSE)
    res <- TiDEomics:::.prepare_gene_list(df)
    expect_type(res, "list")
    expect_named(res, c("1", "2"))
    expect_equal(res[["1"]], c("A", "B"))
    expect_equal(res[["2"]], c("C", "D"))

    # Passes through named list unchanged
    x <- list(A = c("TP53", "BRCA1"), B = c("EGFR"))
    res2 <- TiDEomics:::.prepare_gene_list(x)
    expect_equal(res2, x)

    # Errors on unnamed list
    expect_error(
        TiDEomics:::.prepare_gene_list(list(c("A", "B"), c("C"))),
        "must be a named list"
    )

    # Errors on bad data.frame
    df2 <- data.frame(Gene = c("A", "B"), Cluster = c("1", "1"))
    expect_error(
        TiDEomics:::.prepare_gene_list(df2),
        "must have 'Feature' and 'Module' columns"
    )
})

# enrichR_list mocked tests (unchanged)
test_that("enrichR_list does not print background message when universe is provided", {
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

    gene_list <- list(A = c("TP53", "BRCA1"))

    msgs <- capture_messages(
        enrichR_list(gene_list, databases = "DSigDB",
            universe = c("TP53", "BRCA1", "EGFR", "MYC"))
    )

    expect_false(any(grepl("Background genes not specified", msgs)))
})

test_that("enrichR_list errors when universe_list is unnamed", {
    skip_if_not_installed("enrichR")

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

    gene_list <- list(A = c("TP53", "BRCA1"), B = character(0))
    result <- suppressMessages(
        enrichR_list(gene_list, databases = "DSigDB")
    )

    expect_type(result, "list")
    expect_true("DSigDB" %in% names(result))
    expect_gt(nrow(result[["DSigDB"]]), 0)
})

test_that("enrichR_list messages when no enriched terms found", {
    skip_if_not_installed("enrichR")

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
