# Prepare tutorial data set from GSE263759
if (!requireNamespace("pacman", quietly = TRUE)) install.packages("pacman")
library(pacman)
pacman::p_load(GEOquery, magrittr, tibble, stringr, dplyr, org.Mm.eg.db)

# download sample information from GEO
geo_data <- getGEO("GSE263759", GSEMatrix = TRUE)

geo_sample_info <- pData(geo_data[[1]]) %>%
  dplyr::select(
    "title", "stimulus:ch1", "time point:ch1",
    "bio-replicate:ch1", "batch:ch1"
  ) %>%
  dplyr::rename(
    "Group" = "stimulus:ch1",
    "Time" = "time point:ch1",
    "Replicate" = "bio-replicate:ch1",
    "Batch" = "batch:ch1",
    "Sample" = "title"
  ) %>%
  mutate(
    Time = Time %>% str_replace("h", "") %>% as.numeric(),
    Replicate = Replicate %>% str_replace("R", "") %>% as.numeric(),
    Batch = Batch %>% as.numeric()
  ) # %>%
# arrange(Group, Time, Replicate)
# row.names(geo_sample_info) = geo_sample_info$Sample
row.names(geo_sample_info) <- NULL

geo_sample_info %>%
  dplyr::group_by(Group) %>%
  dplyr::summarise(n = n())

# download expression matrix from GEO
getGEOSuppFiles("GSE263759", baseDir = "data-raw", makeDirectory = FALSE)#
geo_data_tb_gz <- gzfile("data-raw//GSE263759_RNA_raw_counts.csv.gz", 'rt')
geo_data_tb <- read.csv(geo_data_tb_gz) %>%
  dplyr::rename("Feature" = "gene")

# map ensembl IDs to gene symbols
geo_data_tb_symbol <- geo_data_tb$Feature %>%
  bitr(fromType = "ENSEMBL", toType = "SYMBOL", OrgDb = org.Mm.eg.db) %>%
  merge(geo_data_tb, by.x = "ENSEMBL", by.y = "Feature") %>%
  distinct(SYMBOL, .keep_all = TRUE) %>%
  dplyr::select(-ENSEMBL) %>%
  dplyr::rename("Feature" = "SYMBOL")

# remove genes with all zero counts
gene_keep <- geo_data_tb_symbol[, -1] %>%
  apply(1, function(x) sum(x)) %>%
  `!=`(0) %>%
  which()
# sample_order = geo_sample_info$Sample
# geo_data_tb = geo_data_tb[gene_keep, c("Genes", sample_order)]
geo_data_tb_filtered <- geo_data_tb_symbol[gene_keep, ]

# apply time 0 of untreated to all groups
# duplicate the time 0 untreated samples to all other samples' time 0
geo_sample_info_new <- geo_sample_info
geo_data_tb_new <- geo_data_tb_filtered

for (i in unique(geo_sample_info$Group)) {
  if (i != "untreated") {
    group_t0 <- geo_data_tb_filtered %>%
      dplyr::select(starts_with("RNA_untreated_0h"))
    group_t0_sample <- geo_sample_info %>%
      filter(Group == "untreated" & Time == 0)

    colnames(group_t0) <- colnames(group_t0) %>%
      str_replace("untreated", i)
    geo_data_tb_new <- cbind(geo_data_tb_new, group_t0)

    group_t0_sample <- group_t0_sample %>%
      mutate(
        Sample = Sample %>%
          str_replace("RNA_untreated", paste0("RNA_", i)),
        Group = i
      )
    geo_sample_info_new <- rbind(geo_sample_info_new, group_t0_sample)
  }
}

# Tutorial dataset
tutorial_sample_info <- geo_sample_info_new %>%
  filter(Group %in% c("untreated", "IFNbeta", "IFNgamma", "LPS"))
gene_random <- geo_data_tb_new$Feature %>% sample(500)
tutorial_data <- geo_data_tb_new %>%
  dplyr::select(Feature, example_sample_info$Sample) %>%
  filter(Feature %in% gene_random)

usethis::use_data(tutorial_sample_info, tutorial_data)

# Save subset for use in runnable examples
data_obj <- create_input(data = tutorial_data, sample_ann = tutorial_sample_info)
data_obj <- normalise_to_start(data_obj)
# data_obj_list <- split_groups(data_obj)
# data_obj_merged_list <- merge_replicates(data_obj_list)
# data_obj_merged <- merge_groups(data_obj_merged_list)
# data_obj_merged_list <- calc_feature_property(data_obj_merged_list)
# data_obj_merged_imp_list <- impute_groups(data_obj_merged_list)

example_obj = data_obj[1:100, ]
save(example_obj, file = "../data/example.rda")
