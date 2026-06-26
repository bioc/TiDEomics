#' Summarise module patterns
#'
#' @description Summarise module patterns by integrating WGCNA output with
#' Trendy results
#'
#' @param module A data frame with columns "Feature" and "Module"
#' @param trendy_summary A data frame of trendy results, output of
#' `summarise_Trendy()`
#'
#' @returns A list of data frames, each data frame shows the pattern counts
#' in a module
#' @export
#' @examples
#' data("example_res_list")
#' trendy_summary <- summarise_Trendy(example_res_list)
#'
#' data(example_net)
#' example_module <- WGCNA_module(example_net) 
#' summarise_module_pattern(example_module, trendy_summary)
summarise_module_pattern <- function(module, trendy_summary) {
    .check_df(module, "module")
    .check_df_trendy_summary(trendy_summary, "trendy_summary")
    trendy_tb_wider <- trendy_summary |>
        dplyr::select(-Breakpoint) |>
        tidyr::pivot_wider(names_from = Group, values_from = Pattern,
            id_cols = Feature)

    module <- module |>
        dplyr::filter(!Module %in% c("0", "grey", "gray")) |>
        dplyr::mutate(Module = droplevels(Module))
    module_levels <- levels(module$Module)

    patterns <- paste0(colnames(trendy_tb_wider)[-1], collapse = ", ")
    message(sprintf("Most common pattern in each module (%s):", patterns))

    trendy_genes_pattern_list <- list()
    for (lvl in module_levels) {
        module_genes <- module |>
            dplyr::filter(Module == lvl) |>
            dplyr::pull(Feature)

        trendy_genes <- trendy_tb_wider |>
            dplyr::filter(Feature %in% module_genes)

        trendy_genes_pattern <- trendy_genes |>
            tidyr::unite("Pattern", -Feature, sep = ", ") |>
            dplyr::count(Pattern, name = "Count") |>
            dplyr::arrange(dplyr::desc(Count))

        trendy_genes_pattern_list[[lvl]] <- trendy_genes_pattern
    }

    for (lvl in module_levels) {
        tp <- trendy_genes_pattern_list[[lvl]]
        if (is.null(tp) || nrow(tp) == 0) {
            message("Module ", lvl, ": no patterns found.")
        } else {
            message(sprintf("Module %s: %s", lvl, tp[1, 1]))
        }
    }

    return(trendy_genes_pattern_list)
}
