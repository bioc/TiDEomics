# Tests for enrichment functions (GO rank, plot_GO)

test_that("enrichGO_rank works with variance decomposition output", {
    data("example")
    se <- normalise_to_start(example_obj)
    vd <- decomp_variance(se, features = rownames(se)[1:20],
        fixed_effect_var = NULL, core = 1)

    skip_if_not_installed("org.Mm.eg.db")
    library(org.Mm.eg.db)

    res <- enrichGO_rank(vd, gene_rank_by = "Time",
        OrgDb = org.Mm.eg.db, keyType = "SYMBOL",
        category = "BP", pvalueCutoff = 0.99)

    # May return NULL if no significant terms
    if (!is.null(res)) {
        expect_s4_class(res, "gseaResult")
    }
})

test_that("enrichGO_rank validates gene_rank_by", {
    data("example")
    se <- normalise_to_start(example_obj)
    vd <- decomp_variance(se, features = rownames(se)[1:10],
        fixed_effect_var = NULL, core = 1)

    expect_error(
        enrichGO_rank(vd, gene_rank_by = "NotAColumn",
            OrgDb = NULL, keyType = "SYMBOL"),
        "Gene ranking variable not found"
    )
})

test_that("plot_GO dotplot runs with GO enrichment output", {
    skip_if_not_installed("org.Mm.eg.db")
    library(org.Mm.eg.db)

    data("example_net")
    example_module <- WGCNA_module(example_net) 

    go_res <- enrichGO_list(example_module, OrgDb = org.Mm.eg.db,
        universe = example_module$Feature,
        pvalueCutoff = 1, qvalueCutoff = 1,
        category = "BP", simplify = FALSE)

    expect_error(
        plot_GO(go_res$all, plot_dotplot = TRUE,
            plot_emapplot = FALSE, plot_cnetplot = FALSE,
            showCategory_dotplot = 3),
        NA
    )
})

test_that("plot_GO cnetplot runs with GO enrichment output", {
    skip_if_not_installed("org.Mm.eg.db")
    library(org.Mm.eg.db)

    data("example_net")
    example_module <- WGCNA_module(example_net) 

    go_res <- enrichGO_list(example_module, OrgDb = org.Mm.eg.db,
        universe = example_module$Feature,
        pvalueCutoff = 1, qvalueCutoff = 1,
        category = "BP", simplify = FALSE)

    expect_error(
        plot_GO(go_res$all, plot_dotplot = FALSE,
            plot_cnetplot = TRUE, plot_emapplot = FALSE,
            showCategory_cnetplot = 3),
        NA
    )
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
    example_module <- WGCNA_module(example_net) %>%
        dplyr::filter(Module %in% c("1", "2"))

    expect_error(
        enrichGO_list(example_module, OrgDb = org.Mm.eg.db,
            category = "BP", pvalueCutoff = 1, qvalueCutoff = 1),
        NA
    )
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
    gene_list <- list(A = c("TP53", "BRCA1"), B = c("EGFR"))
    expect_error(
        enrichR_list(gene_list, databases = "DSigDB",
            universe_list = list(X = c("TP53"))),
        "must match"
    )
})

