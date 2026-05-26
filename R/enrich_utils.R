# Utility for enrichment functions

#' Prepare input gene list for enrichment functions
#' @description Convert a WGCNA_module() output data.frame to a named gene 
#' list, or pass through an existing named list unchanged.
#' @param x A data.frame with 'Feature' and 'Module' columns 
#' (as returned by WGCNA_module()), or a named list of gene vectors.
#' @return A named list of gene vectors, where names correspond to module 
#' names.
#' @keywords internal
.prepare_gene_list <- function(x) {
    if (is.data.frame(x)) {
        if (!all(c("Feature", "Module") %in% colnames(x))) {
            stop("Input data.frame must have 'Feature' and 'Module' columns ",
                "(as returned by WGCNA_module()).")
        }
        x <- split(x$Feature, x$Module)
    }
    if (is.null(names(x))) {
        stop("Input gene_list must be a named list or a data.frame ",
            "from WGCNA_module().")
    }
    return(x)
}
