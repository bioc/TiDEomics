# Prepare tutorial data set from GSE263759
if (!requireNamespace("pacman", quietly = TRUE)) install.packages("pacman")
library(pacman)
pacman::p_load(GEOquery, tibble, stringr, dplyr, org.Mm.eg.db, clusterProfiler)

# download sample information from GEO
geo_data <- getGEO("GSE263759", GSEMatrix = TRUE, returnType = 'ExpressionSet')
# getGEO() now returns SummarizedExperiment objects by default.
# Pass returnType = 'ExpressionSet' for the previous behavior.

geo_sample_info <- pData(geo_data[[1]]) |>
    dplyr::select(
        "title", "stimulus:ch1", "time point:ch1",
        "bio-replicate:ch1", "batch:ch1"
    ) |>
    dplyr::rename(
        "Group" = "stimulus:ch1",
        "Time" = "time point:ch1",
        "Replicate" = "bio-replicate:ch1",
        "Batch" = "batch:ch1",
        "Sample" = "title"
    ) |>
    dplyr::mutate(
        Time = Time |> str_replace("h", "") |> as.numeric(),
        Replicate = Replicate |> str_replace("R", "") |> as.numeric(),
        Batch = Batch |> as.numeric()
    )
row.names(geo_sample_info) <- NULL

geo_sample_info |>
    dplyr::group_by(Group) |>
    dplyr::summarise(n = dplyr::n())

# download expression matrix from GEO
tmp_dir <- tempdir()
getGEOSuppFiles("GSE263759", baseDir = tmp_dir, makeDirectory = FALSE)
geo_data_tb_gz <-
    gzfile(file.path(tmp_dir, "GSE263759_RNA_raw_counts.csv.gz"), 'rt')
geo_data_tb <- read.csv(geo_data_tb_gz) |>
    dplyr::rename("Feature" = "gene")

# map ensembl IDs to gene symbols
geo_data_tb_symbol <- geo_data_tb$Feature |>
    bitr(fromType = "ENSEMBL", toType = "SYMBOL", OrgDb = org.Mm.eg.db) |>
    merge(geo_data_tb, by.x = "ENSEMBL", by.y = "Feature") |>
    dplyr::distinct(SYMBOL, .keep_all = TRUE) |>
    dplyr::select(-ENSEMBL) |>
    dplyr::rename("Feature" = "SYMBOL")

# remove genes with all zero counts
gene_keep <- geo_data_tb_symbol[, -1] |>
    apply(1, function(x) sum(x))
gene_keep <- which(gene_keep != 0)
geo_data_tb_filtered <- geo_data_tb_symbol[gene_keep, ]

# apply time 0 of untreated to all groups
# duplicate the time 0 untreated samples to all other samples' time 0
geo_sample_info_new <- geo_sample_info
geo_data_tb_new <- geo_data_tb_filtered

for (i in unique(geo_sample_info$Group)) {
    if (i != "untreated") {
        group_t0 <- geo_data_tb_filtered |>
        dplyr::select(starts_with("RNA_untreated_0h"))
        group_t0_sample <- geo_sample_info |>
        dplyr::filter(Group == "untreated" & Time == 0)

        colnames(group_t0) <- colnames(group_t0) |>
        str_replace("untreated", i)
        geo_data_tb_new <- cbind(geo_data_tb_new, group_t0)

        group_t0_sample <- group_t0_sample |>
        dplyr::mutate(
            Sample = Sample |>
            str_replace("RNA_untreated", paste0("RNA_", i)),
            Group = i
        )
        geo_sample_info_new <- rbind(geo_sample_info_new, group_t0_sample)
    }
}

# Tutorial dataset
tutorial_sample_info <- geo_sample_info_new |>
    dplyr::filter(Group %in% c("untreated", "IFNbeta", "IFNgamma", "LPS"))
gene_random <- geo_data_tb_new$Feature |> sample(500)
tutorial_data <- geo_data_tb_new |>
    dplyr::select(Feature, dplyr::all_of(tutorial_sample_info$Sample)) |>
    dplyr::filter(Feature %in% gene_random)

# Normalise: log2(CPM + 1)
count_mat <- as.matrix(tutorial_data[, -1])
lib_sizes <- colSums(count_mat)
cpm <- t(t(count_mat) / lib_sizes * 1e6)
tutorial_data[, -1] <- log2(cpm + 1)

usethis::use_data(tutorial_sample_info, tutorial_data, overwrite = TRUE)

# Save subset for use in runnable examples
library(TiDEomics)
data_obj <- create_input(data = tutorial_data,
    sample_ann = tutorial_sample_info)
data_obj <- normalise_to_start(data_obj)

example_obj <- data_obj[1:100, ]
usethis::use_data(example_obj, overwrite = TRUE)
