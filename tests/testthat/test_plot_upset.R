# ---- Shared data ----
genes <- paste0("G", 1:200)
set_list <- list(
    Exp1_Up = genes[1:60],
    Exp1_Down = genes[61:100],
    Exp2_Up = genes[c(1:40, 101:115)],
    Exp2_Down = genes[c(61:80, 116:130)]
)
enrich_mock <- list(BP = data.frame(
    Cluster = c("Exp1_Up_Exp2_Up", "Exp1_Down", "Exp2_Down"),
    Description = c("term A", "term B", "term C"),
    p.adjust = c(0.001, 0.4, 0.9),
    stringsAsFactors = FALSE
))
enrich_mock$CC <- data.frame(
    Cluster = c("Exp1_Up_Exp2_Up", "Exp1_Down"),
    Description = c("term CC A", "term CC B"),
    p.adjust = c(0.002, 0.03),
    stringsAsFactors = FALSE
)

# ---- Tests ----

test_that("plot_upset handles basic and empty sets", {
    pq <- plot_upset(set_list)
    expect_s4_class(pq, "Heatmap")
    expect_invisible(plot_upset(set_list))

    set_empty <- set_list
    set_empty$Empty <- character(0)
    expect_warning(
        pq <- plot_upset(set_empty),
        "Removing empty set\\(s\\): Empty")
    expect_s4_class(pq, "Heatmap")
})

test_that("plot_upset errors", {
    expect_error(plot_upset(genes), "must be a list")
    expect_error(plot_upset(list(genes, genes)),
        "must be a named list")
    expect_error(plot_upset(list(A = 1:3)),
        "must contain character vectors")
    expect_error(plot_upset(list(A = c("G1", NA_character_))),
        "non-missing, non-empty")
    expect_error(plot_upset(list(A = c("G1", ""))),
        "non-missing, non-empty")
    expect_error(plot_upset(stats::setNames(list(genes), NA_character_)),
        "must be a named list")
    expect_error(suppressWarnings(
        plot_upset(list(A = character(0)))),
        "Need at least 1 non-empty set")
})

test_that("plot_upset enrichment textbox and legend", {
    pq <- plot_upset(set_list,
        enrich_list = enrich_mock,
        enrich_category = "BP",
        enrich_top_n = 3,
        enrich_rank_by = "p.adjust",
        enrich_threshold = 0.05)
    expect_s4_class(pq, "Heatmap")
})

test_that("plot_upset accepts per-category enrichment parameters", {
    expect_s4_class(
        plot_upset(set_list,
            enrich_list = enrich_mock,
            enrich_category = c("BP", "CC"),
            enrich_top_n = c(3, 2),
            enrich_rank_by = c("p.adjust", "p.adjust"),
            enrich_threshold = c(0.05, 0.1)),
        "Heatmap"
    )
    expect_error(
        plot_upset(set_list,
            enrich_list = enrich_mock,
            enrich_category = c("BP", "CC"),
            enrich_top_n = c(1, 2, 3)),
        "length 1 or length 2"
    )
})

test_that("plot_upset warns on unmatched clusters", {
    enrich_bad <- enrich_mock
    enrich_bad$BP <- rbind(enrich_bad$BP, data.frame(
        Cluster = "NoSuchIntersection",
        Description = "term D", p.adjust = 0.01))
    expect_warning(
        pq <- plot_upset(set_list,
            enrich_list = enrich_bad,
            enrich_category = "BP"),
        "not found in UpSet intersections")
    expect_s4_class(pq, "Heatmap")
})

test_that("plot_upset warns on missing rank column", {
    enrich_norank <- list(BP = data.frame(
        Cluster = "Exp1_Up_Exp2_Up",
        Description = "term A",
        stringsAsFactors = FALSE))
    expect_warning(
        pq <- plot_upset(set_list,
            enrich_list = enrich_norank,
            enrich_category = "BP",
            enrich_rank_by = "p.adjust"),
        "enrich_rank_by column 'p.adjust' not found")
    expect_s4_class(pq, "Heatmap")
})

test_that("plot_upset intersection labels, colours, and order", {
    expect_no_error(plot_upset(set_list,
        intersection_labels = c("1010" = "Concordant_Up",
            "0101" = "Concordant_Down"),
        intersection_colours = c("1010" = "#0072B5",
            "0101" = "#20854E")))
    expect_no_error(plot_upset(set_list,
        set_order = names(set_list),
        intersection_order = order(
            ComplexHeatmap::comb_size(
                ComplexHeatmap::make_comb_mat(
                    set_list, set_on_rows = FALSE)), decreasing = TRUE)))
    expect_s4_class(plot_upset(set_list,
        intersection_order = c(
            "1010", "1000", "0101",
            "0100", "0010", "0001")), "Heatmap")

    m <- ComplexHeatmap::make_comb_mat(set_list, set_on_rows = FALSE)
    comb_names <- ComplexHeatmap::comb_name(m)
    intersection_colours <- stats::setNames("#0072B5", comb_names[1])
    expect_s4_class(
        plot_upset(
            set_list, intersection_colours = intersection_colours),
        "Heatmap"
    )
})

test_that("plot_upset rejects duplicate intersection labels", {
    m <- ComplexHeatmap::make_comb_mat(set_list, set_on_rows = FALSE)
    comb_names <- ComplexHeatmap::comb_name(m)
    expect_gte(length(comb_names), 2)
    duplicate_labels <- stats::setNames(
        c("Duplicate", "Duplicate"), comb_names[seq_len(2)]
    )
    expect_error(
        plot_upset(set_list,
            intersection_labels = duplicate_labels),
        "must be unique"
    )
})

test_that("plot_upset rejects duplicate mapping names", {
    m <- ComplexHeatmap::make_comb_mat(set_list, set_on_rows = FALSE)
    comb_names <- ComplexHeatmap::comb_name(m)
    expect_gte(length(comb_names), 2)

    duplicate_labels <- stats::setNames(
        c("Label 1", "Label 2"), comb_names[seq_len(2)]
    )
    names(duplicate_labels)[2] <- names(duplicate_labels)[1]
    expect_error(
        plot_upset(set_list, intersection_labels = duplicate_labels),
        "names must be unique"
    )

    duplicate_colours <- stats::setNames(
        c("#0072B5", "#D55E00"), comb_names[seq_len(2)]
    )
    names(duplicate_colours)[2] <- names(duplicate_colours)[1]
    expect_error(
        plot_upset(set_list, intersection_colours = duplicate_colours),
        "names must be unique"
    )
})
