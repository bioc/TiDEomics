#' Summarise module patterns
#'
#' @description Summarise module patterns by integrating WGCNA output with
#' Trendy results
#'
#' @param module A data frame with columns "Feature" and "Module"
#' @param trendy_summary A data frame of trendy results, output of
#' `summarise_Trendy()`
#' @param print_top_n Whether to output the top n patterns per module as text
#' (default is TRUE)
#' @param top_n Number of top patterns to show per module when `print_top_n` is
#' TRUE (default is 5)
#'
#' @import magrittr
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
summarise_module_pattern <- function(module, trendy_summary,
    print_top_n = TRUE, top_n = 5) {
    trendy_tb_wider <- trendy_summary %>%
        dplyr::select(-Breakpoint) %>%
        tidyr::pivot_wider(names_from = Group, values_from = Pattern,
            id_cols = Feature)

    module <- module %>%
        dplyr::filter(!Module %in% c("0", "grey", "gray")) %>%
        dplyr::mutate(Module = droplevels(Module))
    module_levels <- levels(module$Module)

    patterns <- paste0(colnames(trendy_tb_wider)[-1], collapse = ", ")
    message(sprintf("Most common pattern in each module (%s):", patterns))

    trendy_genes_pattern_list <- list()
    for (lvl in module_levels) {
        module_genes <- module %>%
            dplyr::filter(Module == lvl) %>%
            dplyr::pull(Feature)

        trendy_genes <- trendy_tb_wider %>%
            dplyr::filter(Feature %in% module_genes)

        trendy_genes_pattern <- trendy_genes %>%
            dplyr::select(-Feature) %>%
            apply(1, function(x) paste0(x, collapse = ", ")) %>%
            table() %>%
            data.frame() %>%
            dplyr::arrange(dplyr::desc(Freq))

        colnames(trendy_genes_pattern) <- c(
            paste0(colnames(trendy_tb_wider)[-1], collapse = ", "),
            "Count"
        )

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

    if (print_top_n) {
        for (lvl in module_levels) {
            tp <- trendy_genes_pattern_list[[lvl]]
            if (is.null(tp) || nrow(tp) == 0) next
            message(sprintf("Module %s top %d patterns:", lvl, top_n))
            print(utils::head(tp, n = top_n))
        }
    }

    return(trendy_genes_pattern_list)
}
