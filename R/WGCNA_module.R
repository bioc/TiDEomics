#' Convert WGCNA output to feature-module data frame
#'
#' @description Converts the module assignments from `run_WGCNA()` output
#'   into two-column data.frame with `Feature` and `Module` columns,
#'   expected by `plot_modules_v()`,
#'   `plot_modules_h()`, `summarise_module_pattern()`, and all enrichment
#'   functions.
#'
#' @param net The output of `run_WGCNA()` (a list containing `$colors`)
#' @param exclude_grey Logical. If `TRUE`, features assigned to
#'   module `0` (grey / unassigned) are removed. Default: FALSE
#'
#' @returns A data.frame with columns `Feature` (character) and `Module`
#'   (factor ordered by decreasing module size). When `exclude_grey = TRUE`,
#'   grey/unassigned features are excluded.
#' @export
#'
#' @examples
#' data(example_net)
#' module <- WGCNA_module(example_net)
WGCNA_module <- function(net, exclude_grey = FALSE) {
    .check_list(net, "net", "run_WGCNA")
    .check_logical(exclude_grey, "exclude_grey")
    if (!"colors" %in% names(net)) {
        stop("'net' must contain a 'colors' element. Use run_WGCNA().")
    }
    colors <- net$colors

    gene_module <- data.frame(
        Module = .module_labels(colors, prefix = FALSE),
        row.names = names(colors)
    )

    if (exclude_grey) {
        gene_module <- gene_module[
            !gene_module$Module %in% c("0", "grey", "gray"), , drop = FALSE]
        gene_module$Module <- droplevels(gene_module$Module)
    }

    # Order modules by size (largest first)
    mod_sizes <- table(gene_module$Module)
    gene_module$Module <- factor(gene_module$Module,
        levels = names(sort(mod_sizes, decreasing = TRUE)))

    gene_module <- gene_module |>
        tibble::rownames_to_column("Feature") |>
        dplyr::arrange(Module)

    return(gene_module)
}
