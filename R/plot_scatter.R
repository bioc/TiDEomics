#' Feature property scatter plot between two datasets
#'
#' Compares two per-feature numeric properties. Features present in both
#' input tables are matched by `Feature` and plotted using the selected numeric
#' column from each table on the x- and y-axes. When `df2` is NULL, two columns
#' from `df1` are compared.
#'
#' @param df1 Data frame with columns `Feature` and a numeric value
#'   column (e.g., logFC, variance proportion).
#' @param df2 Optional data frame for the second property. When NULL, `df1` is
#'   used for both properties.
#' @param col1 Name of the numeric column in `df1` to plot on the
#'   x-axis (e.g., `"logFC"`).
#' @param col2 Name of the numeric column in `df2` to plot on the
#'   y-axis.
#' @param name1 Label for dataset 1 (x-axis). Default: `"Dataset 1: <col1>"`.
#' @param name2 Label for dataset 2 (y-axis). Default: `"Dataset 2: <col2>"`.
#' @param highlight_features Optional: character vector of feature names to
#'   label and mark in black. Features not shared by both datasets are reported
#'   and ignored.
#' @param label_top_n Number of top features to label per quadrant. Can be:
#'   - A single integer: same number of labels for all four quadrants.
#'   - A vector of 4 integers: labels for quadrants (++, +-, -+, --)
#'     respectively. Use `c(5, 3, 3, 5)` or named
#'     `c("++"=5, "+-"=3, "-+"=3, "--"=5)`.
#'   Within each quadrant, features are ranked by distance from origin
#'   (sqrt(x^2 + y^2)). These features are labelled along with
#'   threshold-qualified and highlighted features.
#'   Default: 0 (no additional labels).
#' @param label_threshold Minimum absolute value for labelling.
#'   A single value applies to both axes, or use a named numeric
#'   vector, e.g., `c(x = 1, y = 0.5)`. Both the absolute x and absolute y
#'   values must meet their respective thresholds. Qualifying features are
#'   labelled along with the top-N and highlighted features.
#'   Default: NULL.
#' @param show_cor If TRUE (default), annotate the plot with the correlation
#'   coefficient and p-value between the two properties.
#' @param cor_method Correlation method: `"spearman"` (default),
#'   `"pearson"`, or `"kendall"`.
#' @param colour_threshold Minimum absolute value for point colouring. A single
#'   value applies to both axes, or use a named numeric vector, e.g.,
#'   `c(x = 1, y = 0.5)`. Points are coloured only when both absolute x and y
#'   values meet their respective thresholds; points below either threshold
#'   are shown as grey. Default: NULL (colour all points).
#' @param fontsize Base font size (default: 8).
#' @param show_diagonal If NULL (default), show the dashed `y = x`
#'   line when comparing the same property across two data frames.
#'   Set TRUE or FALSE to override this behaviour.
#' @param ... Additional arguments passed to `ggrepel::geom_text_repel()`
#'   for controlling label appearance.
#'
#' @returns A ggplot object of scatter plot.
#' @export
#'
#' @examples
#' # Simulated per-feature numeric properties
#' df_1 <- data.frame(
#'     Feature = paste0("Gene", 1:100),
#'     Property_1 = rnorm(100, 0, 2),
#'     Property_2 = rnorm(100, 0, 1),
#'     stringsAsFactors = FALSE
#' )
#' df_2 <- data.frame(
#'     Feature = paste0("Gene", 1:100),
#'     Property_1 = rnorm(100, 0, 2),
#'     stringsAsFactors = FALSE
#' )
#' plot_scatter(df_1, df_2,
#'     col1 = "Property_1", col2 = "Property_1",
#'     colour_threshold = 1, label_top_n = 5,
#'     highlight_features = c("Gene1", "Gene2", "Gene3")
#' )
#' # Per-quadrant control: 5 labels in concordant, 2 in discordant
#' plot_scatter(df_1, df_2,
#'     col1 = "Property_1", col2 = "Property_1",
#'     colour_threshold = 1,
#'     label_top_n = c("++" = 5, "+-" = 2, "-+" = 2, "--" = 5)
#' )
#' # Compare two properties from one data frame
#' plot_scatter(df_1, col1 = "Property_1", col2 = "Property_2",
#'     colour_threshold = 1, label_top_n = 5)
plot_scatter <- function(
    df1, df2 = NULL,
    col1, col2,
    name1 = NULL, name2 = NULL,
    highlight_features = NULL,
    label_top_n = 0,
    label_threshold = NULL,
    show_cor = TRUE,
    cor_method = c("spearman", "pearson", "kendall"),
    colour_threshold = NULL,
    fontsize = 8,
    show_diagonal = NULL,
    ...
) {
    .check_df(df1, "df1")
    if (is.null(df2)) {
        df2 <- df1
        one_df <- TRUE
    } else {
        .check_df(df2, "df2")
        one_df <- FALSE
    }
    .check_character(col1, "col1")
    .check_character(col2, "col2")
    if (length(col1) != 1L || length(col2) != 1L) {
        stop("'col1' and 'col2' must each contain one value.")
    }
    if (is.null(name1)) name1 <- paste0("Dataset 1: ", col1)
    if (is.null(name2)) name2 <-
        if (one_df) {
            paste0("Dataset 1: ", col2)
        } else {
            paste0("Dataset 2: ", col2)
        }
    .check_character(name1, "name1")
    .check_character(name2, "name2")
    if (length(name1) != 1L || length(name2) != 1L) {
        stop("'name1' and 'name2' must each contain one value.")
    }
    .check_positive(fontsize, "fontsize")
    .check_logical(show_cor, "show_cor")
    if (!is.null(show_diagonal)) {
        .check_logical(show_diagonal, "show_diagonal")
    }
    if (!is.null(highlight_features)) {
        .check_character(highlight_features, "highlight_features")
    }
    if (!is.numeric(label_top_n) ||
        !length(label_top_n) %in% c(1L, 4L) || anyNA(label_top_n) ||
        any(!is.finite(label_top_n)) || any(label_top_n < 0) ||
        any(label_top_n != floor(label_top_n))) {
        stop("'label_top_n' must contain non-negative integer counts, with ",
            "length 1 or 4.")
    }
    if (length(label_top_n) == 4L && !is.null(names(label_top_n)) &&
        !setequal(names(label_top_n), c("++", "+-", "-+", "--"))) {
        stop("Named 'label_top_n' must use the quadrant names ++, +-, ",
            "-+, and --.")
    }
    cor_method <- match.arg(cor_method)
    resolve_threshold <- function(value, arg) {
        if (is.null(value)) return(NULL)
        if (!is.numeric(value) || length(value) < 1L || length(value) > 2L ||
            anyNA(value) || any(!is.finite(value)) || any(value < 0)) {
            stop("'", arg, "' must be NULL, a non-negative number, ",
                "or a numeric vector of length 2 named 'x' and 'y'.")
        }
        if (length(value) == 1) {
            return(stats::setNames(rep(value, 2), c("x", "y")))
        }
        if (is.null(names(value)) ||
            !setequal(names(value), c("x", "y"))) {
            stop("'", arg, "' must name its two values 'x' and 'y'.")
        }
        value[c("x", "y")]
    }
    colour_threshold <- resolve_threshold(
        colour_threshold, "colour_threshold")
    label_threshold <- resolve_threshold(label_threshold, "label_threshold")
    if (!"Feature" %in% colnames(df1)) {
        stop("df1 must have a 'Feature' column.")
    }
    if (!"Feature" %in% colnames(df2)) {
        stop("df2 must have a 'Feature' column.")
    }
    if (!col1 %in% colnames(df1)) {
        stop("Column '", col1, "' not found in df1.")
    }
    if (!col2 %in% colnames(df2)) {
        stop("Column '", col2, "' not found in df2.")
    }
    if (!is.numeric(df1[[col1]])) {
        stop("Column '", col1, "' in df1 must be numeric.")
    }
    if (!is.numeric(df2[[col2]])) {
        stop("Column '", col2, "' in df2 must be numeric.")
    }
    if (is.null(show_diagonal)) {
        show_diagonal <- !one_df && identical(col1, col2)
    }
    if (anyDuplicated(df1$Feature) || anyDuplicated(df2$Feature)) {
        stop("'Feature' values must be unique within each data frame.")
    }

    # ---- Merge on Feature ----
    shared <- intersect(df1$Feature, df2$Feature)
    if (length(shared) == 0) {
        stop("No shared features between df1 and df2.")
    }
    if (!is.null(highlight_features)) {
        not_shared <- setdiff(highlight_features, shared)
        if (length(not_shared) > 0) {
            if (one_df) {
                message(
                    "Ignoring ", length(not_shared),
                    " highlighted feature(s) not present in df1."
                )
            } else {
                message(
                    "Ignoring ", length(not_shared),
                    " highlighted feature(s) not present in both df1 and df2."
                )
            }
        }
    }

    x_vals <- as.numeric(df1[[col1]][match(shared, df1$Feature)])
    y_vals <- as.numeric(df2[[col2]][match(shared, df2$Feature)])

    # Remove incomplete and non-finite pairs.
    keep <- is.finite(x_vals) & is.finite(y_vals)
    n_removed <- sum(!keep)
    if (n_removed > 0L) {
        message("Ignoring ", n_removed,
            " feature pair(s) with missing or infinite values.")
    }
    x_vals <- x_vals[keep]
    y_vals <- y_vals[keep]
    n <- length(x_vals)
    if (n < 2) {
        stop("Fewer than 2 features with non-missing values in both datasets.")
    }

    # ---- Build plot data ----
    df <- data.frame(
        Feature = shared[keep], x = x_vals, y = y_vals,
        stringsAsFactors = FALSE
    )

    # Quadrant relative to origin
    df$quadrant <- ifelse(df$x >= 0 & df$y >= 0, "++",
        ifelse(df$x >= 0 & df$y < 0, "+-",
            ifelse(df$x < 0 & df$y >= 0, "-+", "--")
        )
    )

    # Highlighting
    if (!is.null(highlight_features)) {
        df$highlighted <- df$Feature %in% highlight_features
    } else {
        df$highlighted <- FALSE
    }

    # Colour only points above both axis-specific thresholds.
    if (!is.null(colour_threshold)) {
        df$colored <- abs(df$x) >= colour_threshold["x"] &
            abs(df$y) >= colour_threshold["y"]
        n_grey <- sum(!df$colored)
        n_col <- sum(df$colored)
        sub_msg <- paste0(
            n, " shared features (", n_col, " coloured, ",
            n_grey, " below colour threshold)"
        )
    } else {
        df$colored <- TRUE
        sub_msg <- paste0(n, " shared features")
    }

    # ---- Quadrant counts ----
    # All-points counts (for subtitle)
    qtab_all <- table(factor(df$quadrant,
        levels = c("++", "+-", "-+", "--")))
    # Coloured-only counts (for legend labels)
    qtab_col <- table(factor(df$quadrant[df$colored],
        levels = c("++", "+-", "-+", "--")))
    quadrant_labels <- c(
        sprintf("++ (n=%d)", qtab_col["++"]),
        sprintf("+- (n=%d)", qtab_col["+-"]),
        sprintf("-+ (n=%d)", qtab_col["-+"]),
        sprintf("-- (n=%d)", qtab_col["--"])
    )
    names(quadrant_labels) <- c("++", "+-", "-+", "--")

    # Append all-feature quadrant counts to subtitle
    qstr <- paste0(
        "++:", qtab_all["++"], "  ",
        "+-:", qtab_all["+-"], "  ",
        "-+:", qtab_all["-+"], "  ",
        "--:", qtab_all["--"]
    )
    sub_msg <- paste0(sub_msg, "\n", qstr)
    if (any(df$highlighted)) {
        sub_msg <- paste0(
            sub_msg,
            "\nHighlighted features are labelled and marked in black"
        )
    }

    # ---- Plot ----
    p <- ggplot(df, aes(x = x, y = y))
    if (show_diagonal) {
        p <- p + geom_abline(
            slope = 1, intercept = 0, linewidth = 0.3,
            linetype = "dashed", color = "grey50"
)
    }
    p <- p +
        geom_hline(yintercept = 0, linewidth = 0.2, color = "grey80") +
        geom_vline(xintercept = 0, linewidth = 0.2, color = "grey80")

    # Below-threshold: grey points
    if (any(!df$colored)) {
        p <- p + geom_point(
            data = df[!df$colored & !df$highlighted, , drop = FALSE],
            color = "grey85", size = 1.2, alpha = 0.4
        )
    }

    # Above-threshold: coloured by quadrant
    p <- p +
        geom_point(
            data = df[df$colored & !df$highlighted, , drop = FALSE],
            aes(color = quadrant), size = 1.2, alpha = 0.8
        ) +
        geom_point(
            data = df[df$highlighted, , drop = FALSE],
            color = "black", size = 1.2, alpha = 1
        ) +
        scale_color_manual(
            values = c(
                "++" = "#0072B5",
                "+-" = "#E18727",
                "-+" = "#BC3C29",
                "--" = "#20854E"
            ),
            labels = quadrant_labels,
            name = "Quadrant"
        ) +
        labs(
            x = name1, y = name2, title = paste0(name1, " vs ", name2),
            subtitle = sub_msg
        ) +
        theme_custom(base_size = fontsize)

    # ---- Correlation annotation ----
    cor_data <- NULL
    if (show_cor) {
        if (n < 3L) {
            message("Correlation cannot be calculated because fewer than ",
                "3 feature pairs remain.")
        } else if (length(unique(df$x)) < 2L ||
            length(unique(df$y)) < 2L) {
            message("Correlation cannot be calculated ",
            "because at least one property is constant.")
        } else {
            rho <- stats::cor(df$x, df$y, method = cor_method,
                use = "complete.obs")
            p_val <- stats::cor.test(df$x, df$y,
                method = cor_method, exact = FALSE
            )$p.value
            cor_method_label <- paste0(
                toupper(substr(cor_method, 1, 1)),
                substr(cor_method, 2, nchar(cor_method)))
            p_text <- if (is.na(p_val)) {
                "p = NA"
            } else if (p_val == 0) {
                sprintf("p < %.2e", .Machine$double.eps)
            } else {
                sprintf("p = %.2e", p_val)
            }
            cor_label <- sprintf("%s r = %.3f\n%s",
                cor_method_label, rho, p_text)
            # Place annotation in the quadrant with fewest points.
            x_rng <- range(df$x, na.rm = TRUE)
            y_rng <- range(df$y, na.rm = TRUE)
            margin_x <- 0.05 * diff(x_rng)
            margin_y <- 0.05 * diff(y_rng)
            # Count all points per quadrant.
            q_counts <- c(
                sum(df$quadrant == "++"),
                sum(df$quadrant == "+-"),
                sum(df$quadrant == "-+"),
                sum(df$quadrant == "--")
            )
            names(q_counts) <- c("++", "+-", "-+", "--")
            # Place annotation in the emptiest quadrant.
            best_q <- names(which.min(q_counts))
            pos <- list(
                "++" = list(x = x_rng[2] - margin_x,
                    y = y_rng[2] - margin_y),
                "+-" = list(x = x_rng[2] - margin_x,
                    y = y_rng[1] + margin_y),
                "-+" = list(x = x_rng[1] + margin_x,
                    y = y_rng[2] - margin_y),
                "--" = list(x = x_rng[1] + margin_x,
                    y = y_rng[1] + margin_y)
            )[[best_q]]
            cor_data <- df[1, , drop = FALSE]
            cor_data$Feature <- NA_character_
            cor_data$x <- pos$x
            cor_data$y <- pos$y
            cor_data$quadrant <- NA_character_
            cor_data$highlighted <- FALSE
            cor_data$colored <- TRUE
            cor_data$label <- cor_label
            cor_data$label_type <- "correlation"
        }
    }

    # ---- Label top N per quadrant ----
    if (any(label_top_n > 0) || any(df$highlighted) ||
        !is.null(label_threshold) || !is.null(cor_data)) {
        quadrants <- c("++", "+-", "-+", "--")
        if (length(label_top_n) == 1) {
            label_counts <- stats::setNames(rep(label_top_n, 4), quadrants)
        } else {
            label_counts <- stats::setNames(as.integer(label_top_n), quadrants)
            if (!is.null(names(label_top_n))) {
                label_counts <- stats::setNames(
                    as.integer(label_top_n[quadrants]), quadrants)
            }
        }

        # Threshold-qualified and top-N labels
        label_data <- df[0, , drop = FALSE]
        if (any(df$highlighted)) {
            label_data <- rbind(label_data,
                df[df$highlighted, , drop = FALSE])
        }
        if (!is.null(label_threshold)) {
            df_threshold <- df[abs(df$x) >= label_threshold["x"] &
                abs(df$y) >= label_threshold["y"], , drop = FALSE]
            label_data <- rbind(label_data, df_threshold)
        }
        for (q in quadrants) {
            n_q <- label_counts[q]
            if (n_q > 0) {
                df_q <- df[df$quadrant == q, ,
                    drop = FALSE]
                if (nrow(df_q) > 0) {
                    dist_origin <- sqrt(df_q$x^2 + df_q$y^2)
                    top_q <- utils::head(df_q[order(-dist_origin), ,
                        drop = FALSE], n_q)
                    label_data <- rbind(label_data, top_q)
                }
            }
        }

        label_data <- label_data[!duplicated(label_data$Feature), ]
        label_data$label <- label_data$Feature
        label_data$label_type <- rep("feature", nrow(label_data))
        if (!is.null(cor_data)) {
            label_data <- rbind(label_data, cor_data)
        }
        label_data$segment_colour <- ifelse(
            label_data$label_type == "correlation", NA_character_, "black")

        if (nrow(label_data) > 0) {
            p <- p + ggrepel::geom_text_repel(
                data = label_data,
                aes(label = label, segment.colour = segment_colour),
                size = fontsize * 0.35, max.overlaps = Inf, ...
            )
        }
    }

    p
}
