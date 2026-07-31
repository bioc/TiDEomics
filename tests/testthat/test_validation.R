# Tests for argument validation

data("example_obj")
se_norm <- normalise_to_start(example_obj)
DE_res <- DE_between_group(se_norm, assay = 2)

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
    expect_error(DE_between_group(se_norm, assay = 2, adjP_thres = -0.1), "between 0 and 1")
    expect_error(DE_between_group(se_norm, assay = 2, adjP_thres = 1.5), "between 0 and 1")
    expect_error(DE_between_group(se_norm, assay = 2, logFC_thres = -1), "non-negative")
})

test_that("DE_between_group returns ref_groups and all_groups", {
    expect_true("ref_groups" %in% names(DE_res))
    expect_true("all_groups" %in% names(DE_res))
    expect_s3_class(DE_res$ref_groups, "factor")
    expect_type(DE_res$all_groups, "character")
})

# plot_volcano validation
test_that("plot_volcano validates group/time params conditionally", {
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
    expect_error(extract_hubs(example_net, exclude_grey = "no"), "TRUE or FALSE")
})

# WGCNA_module validation
test_that("WGCNA_module validates net argument", {
    expect_error(WGCNA_module(NULL), "must be a list")
    expect_error(WGCNA_module(list()), "colors")
    data("example_net")
    expect_error(WGCNA_module(example_net, exclude_grey = "yes"), "TRUE or FALSE")
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
        "merged"
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

# ---- check_inputs edge cases ----

test_that(".check_numeric rejects character and factors", {
    expect_error(TiDEomics:::.check_numeric("a", "x"), "must be numeric")
    expect_error(TiDEomics:::.check_numeric(factor(1), "x"), "must be numeric")
})

test_that(".check_positive_int rejects non-integer", {
    expect_error(TiDEomics:::.check_positive_int(1.5, "x"), "positive integer")
    expect_error(TiDEomics:::.check_positive_int(-1, "x"), "positive integer")
})

test_that(".check_df validates data.frame with no rows", {
    expect_error(TiDEomics:::.check_df(data.frame(), "x"), "at least 1 row")
})

test_that(".check_character validates length-0 input", {
    expect_error(TiDEomics:::.check_character(character(0), "x"), "non-empty")
})

test_that(".check_pval rejects NA", {
    expect_error(TiDEomics:::.check_pval(NA_real_, "x"), "between 0 and 1")
})


test_that(".check_se_list rejects empty list", {
    expect_error(TiDEomics:::.check_se_list(list(), "x"), "non-empty")
})

test_that(".check_se_list rejects non-SE elements", {
    expect_error(TiDEomics:::.check_se_list(list(1), "x"), "SummarizedExperiment")
})

test_that(".check_se_list rejects unnamed lists", {
    expect_error(TiDEomics:::.check_se_list(list(se_norm), "x"), "names")
    expect_error(
        TiDEomics:::.check_se_list(
            stats::setNames(list(se_norm, se_norm), c("A", "")), "x"),
        "names")
    expect_error(
        TiDEomics:::.check_se_list(
            stats::setNames(list(se_norm, se_norm), c("A", NA)), "x"),
        "names")
})

test_that(".check_se_list accepts named lists", {
    expect_silent(TiDEomics:::.check_se_list(list(A = se_norm), "x"))
})

# ---- .check_se_list_no_na ----

test_that(".check_se_list_no_na detects missing values", {
    data("example_obj")
    example_obj <- normalise_to_start(example_obj)
    example_obj_list <- split_groups(example_obj)
    example_obj_merged_list <- merge_replicates(example_obj_list)
    assay(example_obj_merged_list[[1]])[1, 1] <- NA
    expect_error(
        TiDEomics:::.check_se_list_no_na(example_obj_merged_list, "x"),
        "missing values")
})

# ---- .check_se_has_properties ----

test_that(".check_se_has_properties errors on missing rowData columns", {
    se <- SummarizedExperiment::SummarizedExperiment(
        assays = list(counts = matrix(1:2, 2, 1)),
        colData = data.frame(Group = "A", Time = 0))
    expect_error(
        TiDEomics:::.check_se_has_properties(se, "x"),
        "missing rowData columns")
})

# ---- .check_df_trendy_summary ----

test_that(".check_df_trendy_summary errors on missing columns", {
    expect_error(
        TiDEomics:::.check_df_trendy_summary(data.frame(x = 1), "x"),
        "summarise_Trendy")
})

# ---- .check_df_feature_property ----

test_that(".check_df_feature_property errors on missing columns", {
    expect_error(
        TiDEomics:::.check_df_feature_property(data.frame(x = 1), "x"),
        "summarise_feature_property")
})

# ---- .check_limma_input ----

test_that(".check_limma_input warns on raw-count-like data", {
    mat <- matrix(round(2^rnorm(100, 10, 2)), nrow = 10)
    expect_warning(
        TiDEomics:::.check_limma_input(mat, "test"),
        "raw counts")
})
