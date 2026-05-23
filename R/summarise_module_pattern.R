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
#' library(magrittr)
#' data("example_res_list")
#' trendy_summary <- summarise_Trendy(example_res_list)
#'
#' data(example_net)
#' example_module <- data.frame(Module = as.factor(example_net$colors)) %>%
#'     tibble::rownames_to_column("Feature") %>% dplyr::arrange(Module)
#' summarise_module_pattern(example_module, trendy_summary)
summarise_module_pattern <- function(module, trendy_summary,
    print_top_n = TRUE, top_n = 5) {
    trendy_tb_wider <- trendy_summary %>%
        dplyr::select(-Breakpoint) %>%
        tidyr::pivot_wider(names_from = Group, values_from = Pattern,
            id_cols = Feature)

    trendy_genes_pattern_list <- list()
    n_module <- length(unique(module$Module)) - 1

    patterns <- paste0(colnames(trendy_tb_wider)[-1], collapse = ", ")
    message(sprintf("Most common pattern in each module (%s):", patterns))

    for (i in seq(1, n_module)) {
        module_genes <- module %>%
            dplyr::filter(Module == i) %>%
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

        trendy_genes_pattern_list[[i]] <- trendy_genes_pattern
    }

    for (i in seq(1, n_module)) {
        message(sprintf("Module %d: %s", i,
            trendy_genes_pattern_list[[i]][1, 1]))
    }

    if (print_top_n) {
        for (i in seq(1, n_module)) {
            message(sprintf("Module %d top %d patterns:", i, top_n))
            print(utils::head(trendy_genes_pattern_list[[i]], n = top_n))
        }
    }

    return(trendy_genes_pattern_list)
}
