# Tests for vertically aligned UpSet plots.

upset_v_sets <- list(
    RNA = paste0("G", 1:6),
    Protein = paste0("G", 4:9)
)
upset_v_comb <- ComplexHeatmap::comb_name(
    ComplexHeatmap::make_comb_mat(upset_v_sets, set_on_rows = TRUE)
)
upset_v_labels <- stats::setNames(
    paste0("Intersection_", seq_along(upset_v_comb)),
    upset_v_comb
)

test_that("plot_upset_v returns an UpSet panel without enrichment", {
    ht <- suppressMessages(plot_upset_v(upset_v_sets))

    expect_s4_class(ht, "Heatmap")
    expect_warning(
        suppressMessages(plot_upset_v(upset_v_sets)),
        NA
    )
    expect_invisible(suppressMessages(plot_upset_v(upset_v_sets)))
})

test_that("plot_upset_v skips enrichment checks without an enrichment list", {
    expect_no_error(
        suppressMessages(plot_upset_v(
            upset_v_sets,
            enrich_list = NULL,
            enrich_category = NA_character_,
            enrich_top_n = 0L,
            enrich_rank_by = NA_character_,
            include_all = "invalid",
            dot_min_size = 0,
            dot_max_size = -1,
            dot_outline = "not-a-colour",
            dot_outline_width = 0,
            row_gap_mm = -1,
            term_wrap_width = 0L,
            show_dot_column_names = "invalid",
            dot_size_col = NA_character_,
            dot_colour_col = ""
        ))
    )
})

test_that("plot_upset_v supports custom numeric enrichment columns", {
    enrich <- data.frame(
        Cluster = c(upset_v_comb[1], upset_v_labels[2], "unmatched"),
        Description = c("term A", "term B", "term C"),
        size = c(3, 2, 1),
        score = c(0.01, 0.02, 0.03),
        stringsAsFactors = FALSE
    )

    messages <- capture_messages(
        result <- plot_upset_v(
            upset_v_sets,
            enrich_list = list(BP = enrich),
            enrich_category = "BP",
            enrich_rank_by = "score",
            dot_size_col = "size",
            dot_colour_col = "score",
            intersection_labels = upset_v_labels,
            intersection_order = rev(upset_v_comb)
        )
    )

    expect_s4_class(result, "HeatmapList")
    expect_true(any(grepl("not matched", messages)))
})

test_that("plot_upset_v accepts bare NA for dot outlines", {
    expect_no_error(
        suppressMessages(plot_upset_v(upset_v_sets, dot_outline = NA))
    )
})

test_that("plot_upset_v validates numeric enrichment columns", {
    enrich <- data.frame(
        Cluster = upset_v_comb[1],
        Description = "term A",
        Count = 3,
        p.adjust = 0.01,
        stringsAsFactors = FALSE
    )

    bad_rank <- enrich
    bad_rank$p.adjust <- "0.01"
    expect_error(
        suppressMessages(plot_upset_v(
            upset_v_sets, enrich_list = list(BP = bad_rank),
            enrich_category = "BP"
        )),
        "must be numeric"
    )

    bad_size <- enrich
    bad_size$Count <- "3"
    expect_error(
        suppressMessages(plot_upset_v(
            upset_v_sets, enrich_list = list(BP = bad_size),
            enrich_category = "BP"
        )),
        "must be numeric"
    )

    bad_colour <- enrich
    bad_colour$p.adjust <- "0.01"
    expect_error(
        suppressMessages(plot_upset_v(
            upset_v_sets, enrich_list = list(BP = bad_colour),
            enrich_category = "BP", enrich_rank_by = "Count"
        )),
        "must be numeric"
    )
})

test_that("plot_upset_v messages when no enrichment cluster matches", {
    enrich <- data.frame(
        Cluster = "unmatched",
        Description = "term A",
        Count = 3,
        p.adjust = 0.01,
        stringsAsFactors = FALSE
    )

    messages <- capture_messages(
        result <- plot_upset_v(
            upset_v_sets,
            enrich_list = list(BP = enrich),
            enrich_category = "BP"
        )
    )

    expect_s4_class(result, "Heatmap")
    expect_true(any(grepl("No Cluster values", messages)))
    expect_true(any(grepl("No matching enrichment terms", messages)))
})

test_that("plot_upset_v validates enrichment list names", {
    enrich <- data.frame(
        Cluster = upset_v_comb[1],
        Description = "term A",
        Count = 3,
        p.adjust = 0.01,
        stringsAsFactors = FALSE
    )

    expect_error(
        suppressMessages(plot_upset_v(
            upset_v_sets,
            enrich_list = list(enrich),
            enrich_category = "BP"
        )),
        "must be a named list"
    )
})
