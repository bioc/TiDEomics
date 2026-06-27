#' Group specific features
#'
#' @description Identify features that are unique to selected groups (present
#' in those groups but not in any others) and annotate them with gene names
#' using the `clusterProfiler::bitr()` function. The function also provides an
#' option to perform Gene Ontology (GO) enrichment analysis on the identified
#' unique features using the `enrichGO_list()` function and visualize the
#' results with a dot plot.
#'
#' @param property_random_fc A data frame containing the results of the
#' `calc_feature_property()` function, which includes the feature names, group
#' names, and the proportion of expressed values for each feature
#' in each group. The
#' data frame should have at least the following columns: "Feature", "Group",
#' "Exp_ratio" and "Exp_threshold".
#' @param groups A character vector of group names to be compared. If NULL, all
#' groups in the input data will be used (default is NULL).
#' @param filter_ratio A numeric value between 0 and 1 specifying the minimum
#' proportion of expressed time points (minimum Exp_ratio) for a feature to be
#' considered present (default is 0.5, meaning that a feature must have >=50%
#' values > threshold in a group (>= 0.5 Exp_ratio) to be considered present
#' in that group).
#' @param group_pct A numeric value between 0 and 1 specifying the percentage
#' of groups in which a feature must be present. (default is 1, meaning that a
#' feature must be present in all specified groups).
#' @param OrgDb An OrgDb object from the `AnnotationDbi` package corresponding
#' to the organism of interest (e.g., `org.Hs.eg.db` for human, `org.Mm.eg.db`
#' for mouse). This will be used for gene annotation with the
#' `clusterProfiler::bitr()` function.
#' @param keytype A character string specifying the type of gene identifiers
#' used in the row names of the assay data (e.g., "SYMBOL", "ENTREZID",
#' "ENSEMBL"). This will be used for gene annotation with the
#' `clusterProfiler::bitr()` function. Available key types depend on the
#' `OrgDb` database and can be checked with the `AnnotationDbi::keytypes`
#' function.
#' @param genename A logical value indicating whether to output a table of
#' gene names. If TRUE, the function will use the `clusterProfiler::bitr()`
#' function to annotate the features with gene names based on the specified
#' `OrgDb` and `keytype`. (default is TRUE).
#' @param GO A logical value indicating whether to perform Gene Ontology
#' (GO) enrichment analysis on the identified unique features. If TRUE, the
#' function will use the `enrichGO_list()` function to perform GO enrichment
#' analysis and visualize the results with a dot plot. (default is TRUE).
#' @param ... Additional arguments to be passed to the `enrichGO_list()`
#' function for GO enrichment analysis (e.g., `pvalueCutoff`, `qvalueCutoff`,
#' etc.).
#'
#' @returns A character vector of features that are identified as unique to
#' the specified groups based on the filtering criteria. If `genename` is
#' TRUE, a table of gene names corresponding to the unique features will be
#' included. If `GO` is TRUE, a dot plot of GO enrichment results for the unique
#' features will be included.
#' @export
#' @examples
#' data(example_obj)
#' example_obj <- normalise_to_start(example_obj)
#' example_obj_list <- split_groups(example_obj)
#' example_obj_merged_list <- merge_replicates(example_obj_list)
#'
#' example_obj_merged_list <-
#'     calc_feature_property(example_obj_merged_list, threshold = 0)
#' property_random_fc <- summarise_feature_property(example_obj_merged_list)
#'
#' group_specific_features(property_random_fc, groups = c("untreated"),
#'     genename = FALSE, GO = FALSE)
group_specific_features <- function(
    property_random_fc, groups = NULL,
    filter_ratio = 0.5,
    group_pct = 1,
    genename = TRUE,
    GO = TRUE,
    OrgDb = NULL, keytype = NULL,
    ...
) {
    .check_pval(filter_ratio, "filter_ratio")
    .check_pval(group_pct, "group_pct")
    .check_df_feature_property(property_random_fc, "property_random_fc")
    .check_logical(genename, "genename")
    .check_logical(GO, "GO")
    if (is.null(groups)) {
        groups <- unique(property_random_fc$Group)
    } else {
        if (!all(groups %in% unique(property_random_fc$Group))) {
            stop("At least one of the specified groups is not found ",
            "in the input data.")
        }
    }

    if (genename || GO) {
        if (is.null(OrgDb) || is.null(keytype)) {
            stop("Both 'OrgDb' and 'keytype' must be provided when ",
            "'genename' or 'GO' is TRUE.")
        }
    }

    threshold <- property_random_fc$Exp_threshold |> unique()
    if (length(threshold) > 1) {
        stop("Multiple Exp_threshold values found across groups. ",
            "All groups should use the same threshold in ",
            "calc_feature_property().")
    }


    # count of included expressed groups for each feature
    filter_count_groups <- property_random_fc |>
        dplyr::filter(Exp_ratio >= filter_ratio) |>
        dplyr::filter(Group %in% groups) |>
        dplyr::distinct(Feature, Group) |>
        dplyr::group_by(Feature) |>
        dplyr::summarise(Count = dplyr::n()) |>
        dplyr::ungroup()

    # count of excluded expressed groups for each feature
    filter_count_groups_op <- property_random_fc |>
        dplyr::filter(Exp_ratio >= filter_ratio) |>
        dplyr::filter(!(Group %in% groups)) |>
        dplyr::distinct(Feature, Group) |>
        dplyr::group_by(Feature) |>
        dplyr::summarise(Count = dplyr::n()) |>
        dplyr::ungroup()

    # features that are present in at least 1 excluded group
    other_group_genes <- filter_count_groups_op |>
        dplyr::filter(Count > 0) |>
        dplyr::pull(Feature)

    # features that are present in at least group_pct of the included groups
    # and not present in any excluded group
    group_num <- ceiling(length(groups) * group_pct)
    unique_genes <- filter_count_groups |>
        dplyr::filter(Count >= group_num) |>
        dplyr::pull(Feature) |>
        setdiff(other_group_genes)

    if (is.na(threshold)) {
        message("Filtering criteria: >=", 100 * filter_ratio,
            "% non-NA time points in >=", group_num, " of groups: ",
            paste(groups, collapse = ", "))
    } else {
        message("Filtering criteria: >=", 100 * filter_ratio,
            "% values >", threshold, " in >=", group_num, " of groups: ",
            paste(groups, collapse = ", "))
    }

    if (length(unique_genes) == 0) {
        message("No unique features with the specified ",
        "filter and groups.")
        return(NULL)
    }

    return_list <- list(features = unique_genes)

    if (genename && length(unique_genes) > 0) {
        return_list$genename <- unique_genes |>
            clusterProfiler::bitr(
                fromType = keytype, toType = c(keytype, "GENENAME"),
                OrgDb = OrgDb
            ) |>
            dplyr::arrange(.data[[keytype]]) |>
            DT::datatable(
                options = list(pageLength = 10),
                caption = paste0(
                    "Features with >=", 100 * filter_ratio,
                    "% ", if (is.na(threshold)) "non-NA time points"
                    else paste0("values >", threshold), " in >=",
                    group_num, " of groups: ",
                    paste(groups, collapse = ", ")
                )
            )
    }

    if (GO == TRUE && length(unique_genes) > 0) {
        unique_genes_go <- enrichGO_list(list("Unique" = unique_genes),
            OrgDb = OrgDb,
            universe = property_random_fc$Feature |> unique(),
            keyType = keytype, ...
        )

        if (length(intersect(names(unique_genes_go$all),
                c("BP", "MF", "CC"))) == 0) {
            message("No GO terms enriched for these features. Skipping plot.")
        } else {
            return_list$GO <- plot_GO(unique_genes_go$all,
                plot_dotplot = TRUE,
                showCategory_dotplot = 10,
                label = paste0(
                    "features with >=", 100 * filter_ratio,
                    "% ", if (is.na(threshold)) "non-NA time points"
                    else paste0("values >", threshold), " in >=",
                    group_num, " of groups: ",
                    paste(groups, collapse = ", ")
                )
            )
        }
    }

    return(return_list)
}
