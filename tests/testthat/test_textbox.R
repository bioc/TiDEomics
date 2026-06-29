# Tests for internal textbox / annotation functions in anno_multicol_textbox.R

# ---- .textbox_grob ----

test_that(".textbox_grob with word_wrap = TRUE", {
    gb <- TiDEomics:::.textbox_grob(
        text = c("term one", "term two with more text"),
        word_wrap = TRUE, add_new_line = FALSE
    )
    expect_s3_class(gb, "textbox")
})

test_that(".textbox_grob with first_text_from = 'bottom'", {
    gb <- TiDEomics:::.textbox_grob(
        text = c("top term", "middle term", "bottom term"),
        first_text_from = "bottom"
    )
    expect_s3_class(gb, "textbox")
})

test_that(".textbox_grob word_wrap with first_text_from = 'bottom'", {
    gb <- TiDEomics:::.textbox_grob(
        text = c("longer first term", "second"),
        word_wrap = TRUE, add_new_line = FALSE,
        first_text_from = "bottom"
    )
    expect_s3_class(gb, "textbox")
})

test_that(".textbox_grob with length(padding) == 2", {
    gb <- TiDEomics:::.textbox_grob(
        text = c("a", "b"),
        padding = unit(c(1, 2), "mm")
    )
    expect_s3_class(gb, "textbox")
})

test_that(".textbox_grob with round_corners = TRUE", {
    gb <- TiDEomics:::.textbox_grob(
        text = c("a", "b"),
        round_corners = TRUE, r = unit(0.1, "snpc")
    )
    expect_s3_class(gb, "textbox")
})

# ---- .measure_multicol_widths ----

test_that(".measure_multicol_widths with empty cats", {
    res <- TiDEomics:::.measure_multicol_widths(list())
    expect_equal(res, stats::setNames(numeric(0), character(0)))
})

# ---- .multicol_textbox_grob ----

test_that(".multicol_textbox_grob with length(cats) == 0", {
    gb <- TiDEomics:::.multicol_textbox_grob(stats::setNames(list(), character(0)))
    expect_s3_class(gb, "null")
})

# ---- .anno_multicol_textbox input validation ----

test_that(".anno_multicol_textbox errors when align_to is list but text is not", {
    df <- data.frame(text = "x", col = "black", fontsize = 8,
                     stringsAsFactors = FALSE)
    expect_error(
        TiDEomics:::.anno_multicol_textbox(
            align_to = list(a = 1:2), text = "not_a_list"
        ),
        "text.*should be a list"
    )
})

test_that(".anno_multicol_textbox errors when align_to length != text length", {
    df <- data.frame(text = "x", col = "black", fontsize = 8,
                     stringsAsFactors = FALSE)
    expect_error(
        TiDEomics:::.anno_multicol_textbox(
            align_to = list(a = 1:2),
            text = list(a = list(BP = df), b = list(BP = df),
                         c = list(BP = df))
        ),
        "should have the same length"
    )
})

test_that(".anno_multicol_textbox errors when names do not match", {
    df <- data.frame(text = "x", col = "black", fontsize = 8,
                     stringsAsFactors = FALSE)
    expect_error(
        TiDEomics:::.anno_multicol_textbox(
            align_to = stats::setNames(list(1:2, 3:4), c("mod1", "mod2")),
            text = stats::setNames(list(list(BP = df), list(BP = df)),
                            c("mod1", "mod3"))
        ),
        "should have the same names"
    )
})

test_that(".anno_multicol_textbox errors on wrong format", {
    expect_error(
        TiDEomics:::.anno_multicol_textbox(
            align_to = factor("x"), text = 42
        ),
        "Format of.*align_to.*or.*text.*is wrong"
    )
})

# ---- .anno_multicol_textbox side / by branches ----

.make_test_text <- function() {
    df <- data.frame(text = c("term1", "term2"),
                     col = c("black", "red"),
                     fontsize = c(8, 8),
                     stringsAsFactors = FALSE)
    list(
        mod1 = list(BP = df, MF = df[1, , drop = FALSE]),
        mod2 = list(BP = df)
    )
}

test_that(".anno_multicol_textbox with side = 'left'", {
    text <- .make_test_text()
    anno <- TiDEomics:::.anno_multicol_textbox(
        align_to = list(mod1 = 1:5, mod2 = 6:10),
        text = text, side = "left", by = "anno_link"
    )
    expect_s4_class(anno, "AnnotationFunction")
})

test_that(".anno_multicol_textbox with by = 'anno_block'", {
    text <- .make_test_text()
    anno <- TiDEomics:::.anno_multicol_textbox(
        align_to = list(mod1 = 1:5, mod2 = 6:10),
        text = text, by = "anno_block"
    )
    expect_s4_class(anno, "AnnotationFunction")
})
