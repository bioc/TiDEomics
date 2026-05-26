#' Plot feature abundance over time
#'
#' @description Plot trend of feature abundances / expression over time by
#' mean and standard deviation (SD) for each group.
#'
#' @param table_mean_sd Output by `calc_mean_sd()`, a dataframe with columns:
#' Feature, Time, Group, Mean, SD
#' @param groups Groups to be plotted, if NULL, all groups will be used
#' (default is NULL)
#' @param features Features to be plotted, if NULL, an error will be raised
#' @param title Title of the plot (default is "Feature")
#' @param ylab Y axis label of the plot (default is "Abundance")
#' @param errorbar Whether to plot error bars (default is TRUE)
#' @param fontsize Font size for the plot (default is 8)
#'
#' @import ggplot2
#' @returns Plot of feature abundances over time by mean and SD
#' @export
#' @examples
#' data("example")
#' table_mean_sd_list <- calc_mean_sd(example_obj)
#' table_mean_sd <- table_mean_sd_list$norm0
#' plot_trend(table_mean_sd,
#'     features = sample(unique(table_mean_sd$Feature), 4))
plot_trend <- function(table_mean_sd, groups = NULL, features,
    title = "Feature", ylab = "Abundance", errorbar = TRUE, fontsize = 8) {
    # table_mean_sd = table_mean_sd_orig or table_mean_sd_norm

    if (is.null(features)) {
        stop("Please specify features to be plotted.")
    }

    if (is.null(groups)) {
        groups <- unique(table_mean_sd$Group)
        message(sprintf("Group not specified. Plotting all groups: %s",
            paste(groups, collapse = ", ")))
    } else if (!all(groups %in% unique(table_mean_sd$Group))) {
        stop("At least one of the specified groups is not found in the input.")
    }

    if (errorbar) {
        p <- table_mean_sd %>%
            dplyr::filter(Group %in% groups & Feature %in% features) %>%
            # dplyr::filter(!is.na(Mean)) %>% 
            # filtering out will make non-NA points directly connected
            # by the input order
            dplyr::mutate(Feature = factor(Feature, levels = features)) %>%
            ggplot(aes(x = Time, y = Mean, group = Group, color = Group)) +
            geom_line(linewidth = 0.7) +
            geom_point(size = 0.3) +
            geom_errorbar(aes(ymin = Mean - SD, ymax = Mean + SD),
                width = 0.5,
                linewidth = 0.7,
                position = position_dodge(0.05),
            ) +
            scale_color_manual(values = get_custom_palette(groups)) +
            theme_custom(
                panel_border = TRUE, legend_position = "bottom",
                base_size = fontsize
            ) +
            ggtitle(paste0(title, " in ", paste(groups, collapse = ", "))) +
            xlab("Time") +
            ylab(ylab) +
            facet_wrap(~Feature, scales = "free_y")
    } else {
        p <- table_mean_sd %>%
            dplyr::filter(Group %in% groups & Feature %in% features) %>%
            # dplyr::filter(!is.na(Mean)) %>% 
            # filtering out will make non-NA points directly connected
            # by the input order
            dplyr::mutate(Feature = factor(Feature, levels = features)) %>%
            ggplot(aes(x = Time, y = Mean, group = Group, color = Group)) +
            geom_line(linewidth = 0.7) +
            geom_point(size = 0.3) +
            scale_color_manual(values = get_custom_palette(groups)) +
            theme_custom(
                panel_border = TRUE, legend_position = "bottom",
                base_size = fontsize
            ) +
            ggtitle(paste0(title, " in ", paste(groups, collapse = ", "))) +
            xlab("Time") +
            ylab(ylab) +
            facet_wrap(~Feature, scales = "free_y")
    }
    print(p)
}
