# Tests for utility functions (palette, theme)

# ---- set_custom_palette / get_custom_palette ----

test_that("set_custom_palette stores palette in options", {
    pal <- c("A" = "red", "B" = "blue")
    set_custom_palette(pal)
    expect_equal(getOption("custom.palette"), pal)
    # clean up
    options(custom.palette = NULL)
})

test_that("set_custom_palette errors on unnamed vector", {
    expect_error(set_custom_palette(c("red", "blue")),
        "named character")
})

test_that("set_custom_palette errors on non-character input", {
    expect_error(set_custom_palette(1:3),
        "named character")
})

test_that("get_custom_palette returns requested groups in order", {
    pal <- c("A" = "red", "B" = "blue", "C" = "green")
    set_custom_palette(pal)
    res <- get_custom_palette(c("C", "A"))
    expect_named(res, c("C", "A"))
    expect_equal(res[["A"]], "red")
    # clean up
    options(custom.palette = NULL)
})

test_that("get_custom_palette generates default palette when unset", {
    options(custom.palette = NULL)
    res <- get_custom_palette(c("X", "Y", "Z"))
    expect_named(res, c("X", "Y", "Z"))
    expect_equal(length(res), 3)
})

test_that("get_custom_palette messages on missing groups", {
    pal <- c("A" = "red")
    set_custom_palette(pal)
    expect_message(get_custom_palette(c("A", "B")), "missing")
    options(custom.palette = NULL)
})

# ---- theme_custom ----

test_that("theme_custom returns a ggplot2 theme", {
    th <- theme_custom()
    expect_s3_class(th, "theme")
})

test_that("theme_custom with panel border uses different settings", {
    th_border <- theme_custom(panel_border = TRUE)
    th_noborder <- theme_custom(panel_border = FALSE)
    # theme elements differ
    expect_false(identical(th_border, th_noborder))
})

test_that("theme_custom respects base_size", {
    th <- theme_custom(base_size = 12)
    expect_equal(th$text$size, 12)
})
