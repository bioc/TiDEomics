comb_helper_sets <- list(
    A = c("a", "b", "c"),
    B = c("a", "b", "d")
)

test_that("local UpSet matrix helpers build combinations", {
    m <- TiDEomics:::.make_comb_mat(
        comb_helper_sets, set_on_rows = TRUE
    )

    expect_equal(dim(m), c(2L, 3L))
    expect_equal(
        TiDEomics:::.comb_name(m),
        c("11", "10", "01")
    )
    expect_equal(
        TiDEomics:::.comb_size(m),
        c("11" = 2L, "10" = 1L, "01" = 1L)
    )
    expect_equal(attr(m, "set_size"), c(3L, 3L))
    expect_equal(rownames(m), names(comb_helper_sets))
    expect_equal(
        unname(as.vector(m)),
        c(1L, 1L, 1L, 0L, 0L, 1L)
    )
})

test_that("local UpSet matrix orientations", {
    m_rows <- TiDEomics:::.make_comb_mat(
        comb_helper_sets, set_on_rows = TRUE
    )
    m_columns <- TiDEomics:::.make_comb_mat(
        comb_helper_sets, set_on_rows = FALSE
    )

    expect_equal(
        unname(as.vector(m_rows)),
        unname(as.vector(t(m_columns)))
    )
    expect_equal(
        TiDEomics:::.comb_name(m_rows),
        TiDEomics:::.comb_name(m_columns)
    )
    expect_equal(
        TiDEomics:::.comb_size(m_rows),
        TiDEomics:::.comb_size(m_columns)
    )
})

test_that("local combination helpers handle a singleton set", {
    m <- TiDEomics:::.make_comb_mat(list(A = "a"))

    expect_equal(dim(m), c(1L, 1L))
    expect_equal(TiDEomics:::.comb_name(m), "1")
    expect_equal(TiDEomics:::.comb_size(m), c("1" = 1L))
    expect_equal(unname(as.vector(m)), 1L)
})

test_that("UpSet annotations and plots accept both orientations", {
    for (set_on_rows in c(TRUE, FALSE)) {
        m <- TiDEomics:::.make_comb_mat(
            comb_helper_sets, set_on_rows = set_on_rows
        )
        annotation_gp <- grid::gpar(fontsize = 8, fontface = "bold")
        axis_param <- list(gp = grid::gpar(fontsize = 8))
        top_annotation <- TiDEomics:::.comb_annotation(
            m,
            position = "top",
            gp = grid::gpar(fill = "#8C8C8C"),
            annotation_name_gp = annotation_gp,
            axis_param = axis_param
        )
        right_annotation <- TiDEomics:::.comb_annotation(
            m,
            position = "right",
            gp = grid::gpar(fill = "#8C8C8C"),
            annotation_name_gp = annotation_gp,
            axis_param = axis_param
        )
        n_sets <- length(comb_helper_sets)
        n_combinations <- length(TiDEomics:::.comb_name(m))
        set_order <- seq_len(n_sets)
        comb_order <- seq_len(n_combinations)

        expect_s4_class(top_annotation, "HeatmapAnnotation")
        expect_s4_class(right_annotation, "HeatmapAnnotation")
        expect_no_error(ComplexHeatmap::UpSet(
            m,
            comb_col = rep("#666666", n_combinations),
            set_order = set_order,
            comb_order = comb_order,
            top_annotation = top_annotation,
            right_annotation = right_annotation
        ))
    }
})
