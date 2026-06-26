library(TiDEomics)

# Load example data
data("example")
example_obj <- normalise_to_start(example_obj)

example_obj_list <- split_groups(example_obj)
example_obj_merged_list <- merge_replicates(example_obj_list)

# Calculate feature properties (required for imputation)
example_obj_merged_list <- calc_feature_property(example_obj_merged_list,
    threshold = 0
)

# Impute missing values (no missing values in example, but run for completeness)
example_obj_merged_imp_list <- impute_groups(example_obj_merged_list)

# Run Trendy analysis
# Note: "untreated" group has only 3 time points, so Trendy will not
# run for this group
example_res_list <- run_Trendy(example_obj_merged_imp_list,
    maxK = 1,
    minNumInSeg = 2,
    meanCut = 0
)

# Save the data
usethis::use_data(example_res_list, overwrite = TRUE)
