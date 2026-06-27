# Tests for argument validation

test_that(".check_numeric accepts numeric", {
    expect_silent(.check_numeric(1, "x"))
    expect_silent(.check_numeric(1.5, "x"))
    expect_silent(.check_numeric(c(1, 2, 3), "x"))
})

test_that(".check_numeric handles edge cases", {
    expect_error(.check_numeric(NA, "x"), "must be numeric")    # NA (logical)
    expect_silent(.check_numeric(NA_real_, "x"))                 # NA_real_ is numeric
    expect_error(.check_numeric(NULL, "x"), "must be numeric")
    expect_silent(.check_numeric(Inf, "x"))                      # Inf is numeric
})

# DE_between_group validation
test_that("DE_between_group validates se_obj", {
    expect_error(DE_between_group(NULL, assay = 2), "SummarizedExperiment")
    expect_error(DE_between_group(data.frame(), assay = 2), "SummarizedExperiment")
})

test_that("DE_between_group validates pval thresholds", {
    data("example_obj")
    example_obj <- normalise_to_start(example_obj)
    expect_error(DE_between_group(example_obj, assay = 2, adjP_thres = -0.1), "between 0 and 1")
    expect_error(DE_between_group(example_obj, assay = 2, adjP_thres = 1.5), "between 0 and 1")
    expect_error(DE_between_group(example_obj, assay = 2, logFC_thres = -1), "non-negative")
})

test_that("DE_between_group returns ref_groups and all_groups", {
    data("example_obj")
    example_obj <- normalise_to_start(example_obj)
    res <- DE_between_group(example_obj, assay = 2)

    expect_true("ref_groups" %in% names(res))
    expect_true("all_groups" %in% names(res))
    expect_type(res$ref_groups, "character")
    expect_type(res$all_groups, "character")
})

# plot_volcano validation
test_that("plot_volcano validates group/time params conditionally", {
    data("example_obj")
    example_obj <- normalise_to_start(example_obj)
    DE_res <- DE_between_group(example_obj, assay = 2)

    # Missing required params
    expect_error(
        plot_volcano(DE_res, group1 = "untreated"),
        "provide either"
    )

    # Valid between-group call
    expect_error(
        plot_volcano(DE_res, group1 = "nonexistent", group2 = "IFNbeta",
            time = 24),
        "No DE results found"
    )
})

# extract_hubs validation
test_that("extract_hubs validates arguments", {
    expect_error(extract_hubs(NULL), "must be a list")
    data("example_net")
    expect_error(extract_hubs(example_net, top_n = 0), "positive")
    expect_error(extract_hubs(example_net, top_n = -1), "positive")
    expect_error(extract_hubs(example_net, exclude_grey = "no"), "logical")
})

# WGCNA_module validation
test_that("WGCNA_module validates net argument", {
    expect_error(WGCNA_module(NULL), "must be a list")
    expect_error(WGCNA_module(list(colors = 1:3)), "colors")
    data("example_net")
    expect_error(WGCNA_module(example_net, exclude_grey = "yes"), "logical")
})

# create_input validation
test_that("create_input validates data and sample_ann", {
    expect_error(create_input(NULL, data.frame()), "data.frame")
    expect_error(create_input(data.frame(a = 1), NULL), "data.frame")
})

# impute_groups validation
test_that("impute_groups validates fun argument", {
    data("example_obj")
    example_obj <- normalise_to_start(example_obj)
    se_list <- split_groups(example_obj)

    expect_error(
        impute_groups(se_list, fun = "min"),
        "must be a function"
    )
})

# enrichGO_rank pAdjustMethod validation
test_that("enrichGO_rank validates pAdjustMethod", {
    skip_if_not_installed("org.Mm.eg.db")
    library(org.Mm.eg.db)

    data("example_obj")
    example_obj <- normalise_to_start(example_obj)
    # Create a dummy rank table
    rank_table <- data.frame(
        Feature = rownames(example_obj)[1:50],
        pvalue = runif(50),
        stringsAsFactors = FALSE
    )

    expect_error(
        enrichGO_rank(rank_table, gene_rank_by = "pvalue",
            OrgDb = org.Mm.eg.db,
            pAdjustMethod = "invalid_method"),
        "should be one of"
    )
})
