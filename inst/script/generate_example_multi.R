# Generate example_obj_2, example_net_2 and example_cons for
# TiDEomics multi-dataset examples and tests
#
# Dataset 1 (Exp1) = real expression data (example_obj).
# Dataset 2 (Exp2) is a perturbed copy of Exp1:
# each gene's assay 1 values receive independent
# Gaussian noise scaled to 0.5x its standard deviation.
# The first 5 feature names of Exp2 then receive the suffix
# "_Exp2" to introduce difference.
#
# Saves:
#   example_obj_2: SummarizedExperiment for the simulated Exp2
#   example_net_2: run_WGCNA() output for Exp2
#   example_cons: run_consensus_WGCNA() output (Exp1 + Exp2)

data(example_obj)

# Perturb Exp1 to create Exp2
exp1_mat <- as.matrix(SummarizedExperiment::assays(example_obj)[[1]])
gene_sd <- apply(exp1_mat, 1, sd)

set.seed(5678)
noise <- matrix(rnorm(length(exp1_mat), 0, 1), nrow = nrow(exp1_mat))
exp2_mat <- exp1_mat + noise * 0.5 * gene_sd

exp2_df <- data.frame(Feature = rownames(exp1_mat), exp2_mat,
    stringsAsFactors = FALSE, check.names = FALSE)

# The first 5 feature names of Exp2 get the suffix "_Exp2".
# The remaining features stay identical to Exp1.
exp2_df$Feature[seq_len(5)] <- paste0(
    exp2_df$Feature[seq_len(5)], "_Exp2"
)

example_obj_2 <- normalise_to_start(create_input(exp2_df,
    as.data.frame(SummarizedExperiment::colData(example_obj))))

# WGCNA
inp2 <- prepare_WGCNA(example_obj_2, assay = 2, powers = seq(1, 30),
    networkType = "signed", RsquaredCut = 0.8)

example_net_2 <- run_WGCNA(inp2, power = 8, minModuleSize = 10,
    numericLabels = TRUE)

# Consensus WGCNA is run on the shared features only
shared_features <- intersect(rownames(example_obj), rownames(example_obj_2))
inp1_cons <- prepare_WGCNA(example_obj[shared_features, , drop = FALSE],
    assay = 2, powers = seq(1, 30), networkType = "signed", RsquaredCut = 0.8)
inp2_cons <- prepare_WGCNA(example_obj_2[shared_features, , drop = FALSE],
    assay = 2, powers = seq(1, 30), networkType = "signed", RsquaredCut = 0.8)

example_cons <- run_consensus_WGCNA(
    nets = list(Exp1 = inp1_cons, Exp2 = inp2_cons),
    powers = c(Exp1 = 8, Exp2 = 8),
    minModuleSize = 10
)

# Save
usethis::use_data(example_obj_2, overwrite = TRUE)
usethis::use_data(example_net_2, overwrite = TRUE)
usethis::use_data(example_cons, overwrite = TRUE)
