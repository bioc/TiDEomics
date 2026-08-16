# ---- Shared data ----
set.seed(1)
df_exp1 <- data.frame(
    Feature = paste0("Gene", 1:100),
    Signed_Max_FC = stats::rnorm(100, 0, 2),
    stringsAsFactors = FALSE
)
df_exp2 <- data.frame(
    Feature = paste0("Gene", 1:100),
    Signed_Max_FC = stats::rnorm(100, 0, 2),
    stringsAsFactors = FALSE
)

# ---- Tests ----

test_that("plot_scatter works", {
    p <- plot_scatter(df_exp1, df_exp2,
        col1 = "Signed_Max_FC", col2 = "Signed_Max_FC",
        name1 = "Exp1", name2 = "Exp2")
    expect_s3_class(p, "ggplot")
    expect_equal(p$labels$title, "Exp1 vs Exp2")
    expect_match(p$labels$subtitle, "100 shared features")
    expect_match(p$labels$subtitle, "\\+\\+:\\d+")
})

test_that("plot_scatter derives default axis labels", {
    p <- plot_scatter(df_exp1, df_exp2,
        col1 = "Signed_Max_FC", col2 = "Signed_Max_FC")
    expect_equal(p$labels$x, "Dataset 1: Signed_Max_FC")
    expect_equal(p$labels$y, "Dataset 2: Signed_Max_FC")
    expect_equal(p$labels$title,
        "Dataset 1: Signed_Max_FC vs Dataset 2: Signed_Max_FC")
})

test_that("plot_scatter compares two properties in one data frame", {
    one_table <- df_exp1
    one_table$Variance <- abs(one_table$Signed_Max_FC)
    p <- plot_scatter(one_table,
        col1 = "Signed_Max_FC", col2 = "Variance")
    expect_s3_class(p, "ggplot")
    expect_equal(p$labels$x, "Dataset 1: Signed_Max_FC")
    expect_equal(p$labels$y, "Dataset 1: Variance")
})

test_that("plot_scatter threshold and label_top_n", {
    p <- plot_scatter(df_exp1, df_exp2,
        col1 = "Signed_Max_FC", col2 = "Signed_Max_FC",
        colour_threshold = 1, label_top_n = 5)
    expect_s3_class(p, "ggplot")
    expect_match(p$labels$subtitle, "below colour threshold")
    expect_no_error(plot_scatter(df_exp1, df_exp2,
        col1 = "Signed_Max_FC", col2 = "Signed_Max_FC",
        colour_threshold = 1,
        label_top_n = c("++" = 5, "+-" = 2, "-+" = 2, "--" = 5)))
    p_asym <- plot_scatter(df_exp1, df_exp2,
        col1 = "Signed_Max_FC", col2 = "Signed_Max_FC",
        colour_threshold = c(x = 1, y = 0.5),
        label_threshold = c(x = 1, y = 0.5), label_top_n = 3)
    expect_equal(
        p_asym$data$colored,
        abs(p_asym$data$x) >= 1 & abs(p_asym$data$y) >= 0.5
    )
})

test_that("plot_scatter highlight_features, correlation options", {
    p <- plot_scatter(df_exp1, df_exp2,
        col1 = "Signed_Max_FC", col2 = "Signed_Max_FC",
        highlight_features = c("Gene1", "Gene2"), label_top_n = 0)
    expect_s3_class(p, "ggplot")
    label_layer <- p$layers[[length(p$layers)]]
    expect_s3_class(label_layer$geom, "GeomTextRepel")
    expect_setequal(
        label_layer$data$Feature[!is.na(label_layer$data$Feature)],
        c("Gene1", "Gene2")
    )
    expect_true(any(label_layer$data$label_type == "correlation"))
    expect_match(
        label_layer$data$label[label_layer$data$label_type == "correlation"],
        "Spearman r ="
    )
    expect_match(
        label_layer$data$label[label_layer$data$label_type == "correlation"],
        "p = [0-9]\\.[0-9]{2}e[+-][0-9]{2}"
    )
    expect_equal(
        unique(label_layer$data$segment_colour[
            label_layer$data$label_type == "feature"]),
        "black"
    )
    expect_true(all(is.na(label_layer$data$segment_colour[
        label_layer$data$label_type == "correlation"])))
    expect_no_error(plot_scatter(df_exp1, df_exp2,
        col1 = "Signed_Max_FC", col2 = "Signed_Max_FC",
        show_cor = FALSE))
    expect_no_error(plot_scatter(df_exp1, df_exp2,
        col1 = "Signed_Max_FC", col2 = "Signed_Max_FC",
        cor_method = "pearson"))
    expect_no_error(plot_scatter(df_exp1, df_exp2,
        col1 = "Signed_Max_FC", col2 = "Signed_Max_FC",
        cor_method = "kendall"))
    expect_error(plot_scatter(df_exp1, df_exp2,
        col1 = "Signed_Max_FC", col2 = "Signed_Max_FC",
        cor_method = "bogus"), "should be one of")
})

test_that("highlighted features are labelled and marked", {
    p <- plot_scatter(df_exp1, df_exp2,
        col1 = "Signed_Max_FC", col2 = "Signed_Max_FC",
        colour_threshold = 100, highlight_features = "Gene1")
    gene1 <- p$data[p$data$Feature == "Gene1", , drop = FALSE]
    expect_true(gene1$highlighted)
    expect_false(gene1$colored)
    expect_match(p$labels$subtitle,
        "Highlighted features are labelled and marked in black")
})

test_that("plot_scatter reports unshared highlighted features", {
    expect_message(
        plot_scatter(df_exp1, df_exp2,
            col1 = "Signed_Max_FC", col2 = "Signed_Max_FC",
            highlight_features = c("Gene1", "NotShared")),
        "Ignoring 1 highlighted feature")
})

test_that("plot_scatter errors", {
    df_bad <- df_exp1
    df_bad$Feature <- NULL
    expect_error(plot_scatter(df_bad, df_exp2,
        col1 = "Signed_Max_FC", col2 = "Signed_Max_FC"),
        "df1 must have a 'Feature' column")
    expect_error(plot_scatter(df_exp1, df_bad,
        col1 = "Signed_Max_FC", col2 = "Signed_Max_FC"),
        "df2 must have a 'Feature' column")
    expect_error(plot_scatter(df_exp1, df_exp2,
        col1 = "Missing", col2 = "Signed_Max_FC"),
        "Column 'Missing' not found in df1")
    expect_error(plot_scatter(df_exp1, df_exp2,
        col1 = "Signed_Max_FC", col2 = "Missing"),
        "Column 'Missing' not found in df2")
    df_disjoint <- data.frame(
        Feature = paste0("Other", 1:100),
        Signed_Max_FC = stats::rnorm(100),
        stringsAsFactors = FALSE
    )
    expect_error(plot_scatter(df_exp1, df_disjoint,
        col1 = "Signed_Max_FC", col2 = "Signed_Max_FC"),
        "No shared features between df1 and df2")
    df_na <- df_exp2
    df_na$Signed_Max_FC[seq_len(99)] <- NA_real_
    expect_error(plot_scatter(df_exp1, df_na,
        col1 = "Signed_Max_FC", col2 = "Signed_Max_FC"),
        "Fewer than 2 features with non-missing")
    expect_error(plot_scatter(df_exp1, df_exp2,
        col1 = "Signed_Max_FC", col2 = "Signed_Max_FC",
        label_top_n = c(1, 2, 3)),
        "non-negative integer")
    expect_error(plot_scatter("not_a_df", df_exp2,
        col1 = "Signed_Max_FC", col2 = "Signed_Max_FC"),
        "'df1' must be a data.frame")
    expect_error(plot_scatter(as.matrix(df_exp1), df_exp2,
        col1 = "Signed_Max_FC", col2 = "Signed_Max_FC"),
        "'df1' must be a data.frame")
    expect_error(plot_scatter(df_exp1, df_exp2,
        col1 = "Signed_Max_FC", col2 = "Signed_Max_FC",
        fontsize = -1), "positive")
    expect_error(plot_scatter(df_exp1, df_exp2,
        col1 = "Signed_Max_FC", col2 = "Signed_Max_FC",
        colour_threshold = c(x = 1, z = 0.5)), "name")
    expect_error(plot_scatter(df_exp1, df_exp2,
        col1 = "Signed_Max_FC", col2 = "Signed_Max_FC",
        colour_threshold = -1), "non-negative")
    expect_error(plot_scatter(df_exp1, df_exp2,
        col1 = c("Signed_Max_FC", "Other"), col2 = "Signed_Max_FC"),
        "must each contain one value")
    expect_error(plot_scatter(df_exp1, df_exp2,
        col1 = "Signed_Max_FC", col2 = "Signed_Max_FC",
        label_top_n = c("++" = 1, "bad" = 2, "-+" = 3, "--" = 4)),
        "quadrant names")
    expect_error(plot_scatter(df_exp1, df_exp2,
        col1 = "Signed_Max_FC", col2 = "Signed_Max_FC",
        show_cor = NA), "TRUE or FALSE")
    expect_error(plot_scatter(df_exp1, df_exp2,
        col1 = "Signed_Max_FC", col2 = "Signed_Max_FC",
        show_cor = 1), "TRUE or FALSE")
    expect_error(plot_scatter(df_exp1, df_exp2,
        col1 = "Signed_Max_FC", col2 = "Signed_Max_FC",
        show_cor = "TRUE"), "TRUE or FALSE")
    df_duplicate <- df_exp1
    df_duplicate$Feature[2] <- df_duplicate$Feature[1]
    expect_error(plot_scatter(df_duplicate, df_exp2,
        col1 = "Signed_Max_FC", col2 = "Signed_Max_FC"),
        "Feature.*unique")
    df_character <- df_exp1
    df_character$Signed_Max_FC <- as.character(df_character$Signed_Max_FC)
    expect_error(plot_scatter(df_character, df_exp2,
        col1 = "Signed_Max_FC", col2 = "Signed_Max_FC"), "must be numeric")
})

test_that("plot_scatter handles non-finite and constant properties", {
    df_nonfinite <- df_exp2
    df_nonfinite$Signed_Max_FC[1:2] <- c(Inf, -Inf)
    expect_message(
        p <- plot_scatter(df_exp1, df_nonfinite,
            col1 = "Signed_Max_FC", col2 = "Signed_Max_FC"),
        "missing or infinite"
    )
    expect_equal(nrow(p$data), 98)

    df_two <- df_exp2
    df_two$Signed_Max_FC[3:100] <- NA_real_
    expect_message(
        p_two <- plot_scatter(df_exp1, df_two,
            col1 = "Signed_Max_FC", col2 = "Signed_Max_FC"),
        "fewer than 3 feature pairs"
    )
    expect_equal(nrow(p_two$data), 2)

    df_constant <- df_exp2
    df_constant$Signed_Max_FC <- 1
    expect_message(
        expect_no_error(plot_scatter(df_exp1, df_constant,
            col1 = "Signed_Max_FC", col2 = "Signed_Max_FC")),
        "one property is constant"
    )
})
