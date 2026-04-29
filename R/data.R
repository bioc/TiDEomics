#' Dataset for TiDEomics tutorial, sample information
#'
#' A subset of GSE263759 data set published in 
#' [Integrated time-series analysis and high-content CRISPR screening delineate the dynamics of macrophage immune regulation](https://doi.org/10.1016/j.cels.2025.101346)
#'
#' Code for preparing the data is available in `data-raw/tutorial_input.R`
#' - Ensembl IDs were mapped to symbols, genes with all zero counts 
#' were excluded.
#' - Use time 0 untreated samples for other groups' time 0.
#' - Include "untreated", "IFNbeta", "IFNgamma", "LPS" groups.
#' - Sample 500 random genes.
#'
#' @format A data.frame with 40 rows and 5 variables:
#' \describe{
#'   \item{Sample}{Sample ID}
#'   \item{Group}{Experimental group (untreated, different treatments)}
#'   \item{Time}{Time point}
#'   \item{Replicate}{Replicate ID for each group and time point}
#'   \item{Batch}{Batch information}
#' }
#' @source GSE263759
#' @usage data(tutorial_sample_info)
"tutorial_sample_info"

#' Dataset for TiDEomics tutorial, expression matrix
#'
#' A subset of GSE263759 data set published in 
#' [Integrated time-series analysis and high-content CRISPR screening delineate the dynamics of macrophage immune regulation](https://doi.org/10.1016/j.cels.2025.101346)
#'
#' Code for preparing the data is available in `data-raw/tutorial_input.R`
#' - Ensembl IDs were mapped to symbols, genes with all zero counts 
#' were excluded.
#' - Use time 0 untreated samples for other groups' time 0.
#' - Include "untreated", "IFNbeta", "IFNgamma", "LPS" groups.
#' - Sample 500 random genes.
#'
#' @format A data.frame with 500 rows and 41 variables:
#' \describe{
#'   \item{Feature}{Feature ID, e.g. gene symbols}
#'   \item{Sample1, Sample2, ...}{Expression values for each sample}
#' }
#' @source GSE263759
#' @usage data(tutorial_data)
"tutorial_data"

#' SummarizedExperiment object for runnable examples
#'
#' A subset of `data_obj <- create_input(data = tutorial_data, 
#' sample_ann = tutorial_sample_info)` for use in runnable examples 
#' in function documentation.
#'
#' Code for preparing the data is available in `data-raw/tutorial_input.R`
#'
#' @format A SummarizedExperiment object with assays of 100 rows and 
#' 40 columns, colData of 40 rows and 5 columns:
#' \describe{
#'   \item{colData}{tutorial_sample_info}
#'   \item{assays}{tutorial_data first 100 rows, "Feature" column as rownames}
#' }
#' @source GSE263759
#' @usage data(example)
"example_obj"

#' `run_WGCNA()` output object for runnable examples
#'
#' Code for producing the data is available in `run_WGCNA()` examples
#'
#' @format A list including WGCNA module assignments, module eigengenes,
#' dendrogram, input data, sample information, and parameters used
#' @source GSE263759
#' @usage data(example_net)
"example_net"

#' `run_Trendy` output object for runnable examples
#'
#' Code for preparing the data is available in `run_Trendy()` examples
#'
#' @format A nested list, each element is a list with Trendy results 
#' for one group
#' @source GSE263759
#' @usage data(example_res_list)
"example_res_list"
