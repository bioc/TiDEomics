# Tests for the added plot_modules_h selected module layout, term wrapping,
# and legend options.

data(example_obj)
data(example_net)

added_se <- example_obj |>
    split_groups() |>
    merge_replicates() |>
    merge_groups()
added_se_list <- list(Exp1 = added_se, Exp2 = added_se)
added_module <- WGCNA_module(example_net, exclude_grey = TRUE)
added_module_levels <- levels(added_module$Module)
added_selected <- added_module_levels[seq_len(2L)]

added_terms <- function(prefix) {
    data.frame(
        Cluster = rep(added_selected, each = 3L),
        Description = rep(
            paste0(
                prefix,
                c(
                    " long biological process term for wrapping",
                    " another extended annotation term for wrapping",
                    " third extended annotation term for wrapping"
                )
            ),
            length(added_selected)
        ),
        p.adjust = rep(c(0.001, 0.02, 0.2), length(added_selected)),
        stringsAsFactors = FALSE
    )
}

added_enrich <- list(
    BP = added_terms("BP"),
    CC = added_terms("CC")
)

.draw_added_plot <- function(...) {
    grDevices::pdf(file = NULL)
    on.exit(grDevices::dev.off(), add = TRUE)
    plot_modules_h(...)
}

test_that("selected layouts draw with automatic and fixed heights", {
    auto_height <- .draw_added_plot(
        added_module, added_se, assay = 2, scale = TRUE,
        selected_modules = added_selected,
        enrich_list = added_enrich,
        enrich_category = c("BP", "CC"),
        enrich_wrap_width = 1.2,
        profile_width = 2.5, heatmap_width = 6,
        heatmap_height = NULL, fontsize = 6
    )
    fixed_height <- .draw_added_plot(
        added_module, added_se, assay = 2, scale = TRUE,
        selected_modules = added_selected,
        enrich_list = added_enrich,
        enrich_category = c("BP", "CC"),
        enrich_wrap_width = 1.2,
        profile_width = 2.5, heatmap_width = 6,
        heatmap_height = 12, fontsize = 6
    )

    expect_s4_class(auto_height, "Heatmap")
    expect_s4_class(fixed_height, "Heatmap")
})

test_that("selected layouts work for multiple datasets without wrapping", {
    result <- .draw_added_plot(
        added_module, added_se_list, assay = 2, scale = TRUE,
        selected_modules = added_selected,
        enrich_list = added_enrich,
        enrich_category = c("BP", "CC"),
        profile_width = 2.5, heatmap_width = 6,
        heatmap_height = 12, fontsize = 6
    )

    expect_s4_class(result, "HeatmapList")
    expect_length(result@ht_list, 2L)
})

test_that("all-module textboxes retain wrapped and unwrapped paths", {
    unwrapped <- .draw_added_plot(
        added_module, added_se, assay = 2, scale = TRUE,
        enrich_list = added_enrich,
        enrich_category = c("BP", "CC"),
        profile_width = 2.5, heatmap_width = 6,
        heatmap_height = 8, fontsize = 6
    )
    wrapped <- .draw_added_plot(
        added_module, added_se, assay = 2, scale = TRUE,
        enrich_list = added_enrich,
        enrich_category = c("BP", "CC"),
        enrich_wrap_width = 1.2,
        legend_position = "bottom", return_legends = TRUE,
        profile_width = 2.5, heatmap_width = 6,
        heatmap_height = 8, fontsize = 6
    )

    expect_s4_class(unwrapped, "Heatmap")
    expect_named(wrapped, c("heatmap", "legends"))
    expect_s4_class(wrapped$heatmap, "Heatmap")
    expect_true(length(wrapped$legends) > 0L)
})

test_that("textbox helpers support custom term and line spacing", {
    grDevices::pdf(file = NULL)
    on.exit(grDevices::dev.off(), add = TRUE)
    df <- data.frame(
        text = c("first long term", "second long term"),
        col = c("black", "grey50"), fontsize = c(8, 8),
        stringsAsFactors = FALSE
    )
    spacing_args <- list(
        line_space = grid::unit(0.15 * 8, "pt"),
        term_space = grid::unit(0.75 * 8, "pt")
    )
    no_wrap <- do.call(
        TiDEomics:::.textbox_grob,
        c(list(text = df$text, gp = grid::gpar(col = df$col,
            fontsize = df$fontsize), add_new_line = TRUE,
            word_wrap = FALSE), spacing_args)
    )
    wrapped <- do.call(
        TiDEomics:::.textbox_grob,
        c(list(text = df$text, gp = grid::gpar(col = df$col,
            fontsize = df$fontsize), add_new_line = TRUE,
            word_wrap = TRUE, max_width = grid::unit(12, "mm")),
            spacing_args)
    )

    text <- list(
        M1 = list(BP = df, CC = df[1, , drop = FALSE]),
        M2 = list(BP = df[2, , drop = FALSE], CC = df[1, , drop = FALSE])
    )
    layout <- TiDEomics:::.multicol_textbox_layout(
        text,
        textbox_args = c(
            list(
                word_wrap = TRUE,
                max_width = grid::unit(12, "mm"),
                text_padding = grid::unit(0.65 * 8, "pt")
            ),
            spacing_args
        )
    )
    panel_fun <- TiDEomics:::.multicol_textbox_panel_fun(
        layout,
        background_gp = grid::gpar(
            fill = "grey98", col = "grey50", lty = 1, lwd = 0.5
        ),
        show_headers = TRUE,
        header_gp = grid::gpar(fontsize = 8, fontface = "bold")
    )
    grid::grid.newpage()
    panel_fun(1L, "M1")

    expect_s3_class(no_wrap, "textbox")
    expect_s3_class(wrapped, "textbox")
    expect_true(all(is.finite(layout$heights_mm)))
    expect_true(is.function(panel_fun))
})

test_that("new arguments validate only where applicable", {
    expect_error(
        plot_modules_h(added_module, added_se, assay = 2,
            selected_modules = TRUE),
        "character or numeric"
    )
    expect_error(
        plot_modules_h(added_module, added_se, assay = 2,
            selected_modules = character(0)),
        "at least one module"
    )
    expect_error(
        plot_modules_h(added_module, added_se, assay = 2,
            selected_modules = "not_a_module"),
        "unknown module"
    )
    expect_error(
        plot_modules_h(added_module, added_se, assay = 2,
            enrich_list = added_enrich,
            enrich_category = "BP", enrich_wrap_width = 0),
        "positive number"
    )
    expect_error(
        plot_modules_h(added_module, added_se, assay = 2,
            enrich_list = added_enrich,
            enrich_category = "BP", enrich_wrap_width = c(1, 2)),
        "positive number"
    )
    expect_error(
        plot_modules_h(added_module, added_se, assay = 2,
            legend_position = "left"),
        "one of"
    )

    no_enrichment <- .draw_added_plot(
        added_module, added_se, assay = 2,
        enrich_wrap_width = "ignored_without_enrichment",
        heatmap_width = 6, heatmap_height = 8, fontsize = 6
    )
    expect_s4_class(no_enrichment, "Heatmap")
})
