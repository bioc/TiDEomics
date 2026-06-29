#' Set custom color palette
#' @description Set a custom color palette for groups to use in all subsequent
#' plotting functions that support custom palettes.
#'
#' @param palette A named vector where names correspond to group names
#' and values are the corresponding colors (e.g., `c("Group1" = "#FF0000",
#' "Group2" = "#00FF00")`).
#' @returns The function does not return a value but sets the custom palette as
#' an option that can be accessed by other plotting functions.
#' The palette will be stored in the R session options under the
#' name "custom.palette".
#' @export
#' @examples
#' custom_palette <- c("untreated" = "#1b9e77", "IFNbeta" = "#d95f02",
#'     "IFNgamma" = "#7570b3", "LPS" = "#e7298a")
#' set_custom_palette(custom_palette)
set_custom_palette <- function(palette) {
    if (!is.character(palette) || is.null(names(palette))) {
        stop("Palette must be a named character vector with group names ",
        "as names and colors as values.")
    }
    options(custom.palette = palette)
    message(sprintf("Custom palette has been set successfully for groups %s.",
        paste(names(palette), collapse = ", ")))
}

#' Get custom color palette
#' @description Get the custom color palette set by `set_custom_palette`.
#' If no custom palette has been set, the function will return a default
#' color palette based on the number of groups specified in the
#' `groups` argument.
#'
#' @param groups A character vector of group names for which to retrieve
#' the colors from the custom palette.
#'
#' @returns A named vector of colors corresponding to the specified group names.
#' @export
#' @examples
#' get_custom_palette(c("untreated", "IFNbeta"))
get_custom_palette <- function(groups) {
    groups <- as.character(groups)
    if (length(groups) == 0) stop("'groups' must not be empty.")
    pal <- getOption("custom.palette")

    if (is.null(pal)) { # no custom palette set, return default palette
        pal <- scales::pal_hue()(length(groups))
        names(pal) <- groups
    } else {
        missing <- setdiff(groups, names(pal))
        if (length(missing) > 0) {
            missing_groups <- paste(missing, collapse = ", ")
            message("The following groups are missing in the custom palette: ",
                missing_groups)
            # generate default colors for missing groups
            pal_2 <- scales::pal_hue()(length(missing))
            names(pal_2) <- missing
            pal <- c(pal, pal_2) # combine custom and default palettes
        }
    }

    # return colors in the order of the specified groups
    pal <- pal[groups |> as.character()]

    return(pal)
}
