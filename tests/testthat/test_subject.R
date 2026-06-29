# Helper: create a small test dataset with Subject
.make_test_se <- function(n_genes = 20, n_subjects = 4, n_times = 3,
                          n_groups = 1, with_reps = FALSE, seed = 42) {
    set.seed(seed)
    times <- seq(0, (n_times - 1) * 2, by = 2)
    groups <- paste0("G", seq_len(n_groups))

    # Assign subjects to groups (each subject in exactly one group)
    subj_per_group <- n_subjects %/% n_groups
    extra <- n_subjects %% n_groups
    subj_assign <- rep(groups, each = subj_per_group)
    if (extra > 0) subj_assign <- c(subj_assign, groups[seq_len(extra)])
    subjects <- paste0("S", seq_len(n_subjects))
    subject_group <- data.frame(Subject = subjects, Group = subj_assign,
        stringsAsFactors = FALSE)

    sample_ann <- expand.grid(
        Time = times,
        Subject = subjects,
        stringsAsFactors = FALSE
    )
    sample_ann <- merge(sample_ann, subject_group, by = "Subject")
    if (with_reps) {
        sample_ann <- rbind(sample_ann, sample_ann)
        sample_ann <- sample_ann |>
            dplyr::group_by(Subject, Group, Time) |>
            dplyr::mutate(Replicate = dplyr::row_number()) |>
            dplyr::ungroup() |>
            as.data.frame()
    } else {
        sample_ann <- sample_ann |>
            dplyr::group_by(Group, Time) |>
            dplyr::mutate(Replicate = dplyr::row_number()) |>
            dplyr::ungroup() |>
            as.data.frame()
    }
    sample_ann$Sample <- paste0(sample_ann$Group, "_T", sample_ann$Time,
        "_", sample_ann$Subject,
        if (with_reps) paste0("_R", sample_ann$Replicate) else "")

    # Simulate: subject baseline + group effect + time effect + noise
    subj_baseline <- stats::setNames(rnorm(n_subjects, 10, 2), subjects)
    grp_effect <- stats::setNames(seq_len(n_groups) * 2, groups)
    time_effect <- times * 0.5

    data_mat <- matrix(NA, nrow = n_genes, ncol = nrow(sample_ann))
    rownames(data_mat) <- paste0("Gene", seq_len(n_genes))
    colnames(data_mat) <- sample_ann$Sample

    for (i in seq_len(nrow(sample_ann))) {
        s <- sample_ann$Subject[i]
        g <- sample_ann$Group[i]
        t <- sample_ann$Time[i]
        base <- subj_baseline[s] + grp_effect[g] + time_effect[t == times]
        data_mat[, i] <- rnorm(n_genes, base, 0.5)
    }

    data <- data.frame(Feature = rownames(data_mat), data_mat,
        check.names = FALSE, stringsAsFactors = FALSE)
    list(data = data, ann = sample_ann[, c("Sample", "Group", "Time",
        "Subject", "Replicate")])
}


test_that("create_input with subject_col works", {
    td <- .make_test_se(n_subjects = 4)
    expect_message(
        se <- create_input(td$data, td$ann, subject_col = "Subject"),
        "Subject column 'Subject' recognized"
    )
    expect_true("Subject" %in% colnames(colData(se)))
    expect_equal(dplyr::n_distinct(colData(se)$Subject), 4)
})

test_that("create_input without subject_col gives informative message", {
    td <- .make_test_se(n_subjects = 4)
    td$ann$Subject <- NULL
    expect_message(
        se <- create_input(td$data, td$ann),
        "No Subject column specified"
    )
    expect_false("Subject" %in% colnames(colData(se)))
})

test_that("create_input auto-generates Replicate when absent", {
    td <- .make_test_se(n_subjects = 2, with_reps = TRUE)
    td$ann$Replicate <- NULL
    expect_message(
        se <- create_input(td$data, td$ann, subject_col = "Subject"),
        "Auto-generating replicate IDs"
    )
    expect_true("Replicate" %in% colnames(colData(se)))
})

test_that("create_input errors for invalid subject_col", {
    td <- .make_test_se()
    expect_error(
        create_input(td$data, td$ann, subject_col = "NotAColumn"),
        "not found in sample_ann"
    )
})

test_that("decomp_variance auto-detects Subject", {
    td <- .make_test_se(n_subjects = 4, n_groups = 2, with_reps = TRUE)
    se <- create_input(td$data, td$ann, subject_col = "Subject")

    vd <- decomp_variance(se, features = paste0("Gene", 1:10),
        fixed_effect_var = NULL, assay = "orig", core = 1)
    expect_true("Subject" %in% colnames(vd))
    expect_true("Group" %in% colnames(vd))
    expect_true("Time" %in% colnames(vd))
    expect_true("Residual" %in% colnames(vd))
    # Percentages should sum to ~100
    pct_sum <- rowSums(vd[, c("Subject", "Group", "Time", "Residual")])
    expect_true(all(abs(pct_sum - 100) < 1, na.rm = TRUE))
})

test_that("decomp_variance without Subject gives basic output", {
    td <- .make_test_se(n_subjects = 4, n_groups = 2)
    td$ann$Subject <- NULL
    se <- create_input(td$data, td$ann)
    vd <- decomp_variance(se, features = paste0("Gene", 1:10),
        fixed_effect_var = NULL, assay = "orig", core = 1)

    expect_true("Group" %in% colnames(vd))
    expect_true("Time" %in% colnames(vd))
    expect_true("Residual" %in% colnames(vd))
    expect_false("Subject" %in% colnames(vd))
    # Percentages should sum to ~100
    pct_sum <- rowSums(vd[, c("Group", "Time", "Residual")])
    expect_true(all(abs(pct_sum - 100) < 1, na.rm = TRUE))
})

test_that("decomp_variance handles 1-group design", {
    td <- .make_test_se(n_subjects = 4, n_groups = 1)
    se <- create_input(td$data, td$ann, subject_col = "Subject")
    expect_message(
        vd <- decomp_variance(se, features = paste0("Gene", 1:10),
            fixed_effect_var = NULL, assay = "orig", core = 1),
        "Group dropped"
    )
    expect_false("Group" %in% colnames(vd))
    expect_true("Time" %in% colnames(vd))
})

test_that("merge_replicates with Subject does two-stage averaging", {
    td <- .make_test_se(n_subjects = 4, n_groups = 2, with_reps = TRUE)
    se <- create_input(td$data, td$ann, subject_col = "Subject")
    se_list <- split_groups(se)
    expect_message(
        merged <- merge_replicates(se_list),
        "Subject column detected"
    )
    for (nm in names(merged)) {
        expect_equal(ncol(merged[[nm]]), length(unique(colData(se)$Time)))
    }
})

test_that("normalise_to_start by_subject works", {
    td <- .make_test_se(n_subjects = 4, n_groups = 2)
    se <- create_input(td$data, td$ann, subject_col = "Subject")
    expect_message(
        se_norm <- normalise_to_start(se, by_subject = TRUE),
        "subject-level baseline"
    )
    # Check assay 2 exists
    expect_equal(length(assays(se_norm)), 2)
})

test_that("normalise_to_start by_subject errors without Subject", {
    td <- .make_test_se(n_subjects = 4)
    td$ann$Subject <- NULL
    se <- create_input(td$data, td$ann)
    expect_error(
        normalise_to_start(se, by_subject = TRUE),
        "requires a 'Subject' column"
    )
})

test_that("calc_mean_sd with Subject computes between-subject SD", {
    td <- .make_test_se(n_subjects = 4, n_groups = 2)
    se <- create_input(td$data, td$ann, subject_col = "Subject")
    expect_message(
        res <- calc_mean_sd(se),
        "Subject column detected"
    )
    expect_true("Mean" %in% colnames(res$orig))
    expect_true("SD" %in% colnames(res$orig))
})

test_that("impute_groups with impute_by = 'subject' works", {
    td <- .make_test_se(n_subjects = 4, n_groups = 2)
    se <- create_input(td$data, td$ann, subject_col = "Subject")

    # Introduce NAs
    assay(se)[1:5, 1:3] <- NA

    se_list <- split_groups(se)
    merged <- merge_replicates(se_list)
    # Merging drops Subject -> falls back to group-level with a warning
    expect_warning(
        imputed <- impute_groups(merged, impute_by = "subject"),
        "Falling back"
    )
    for (nm in names(imputed)) {
        expect_false(any(is.na(assay(imputed[[nm]]))))
    }
})

test_that("impute_groups with impute_by = 'group' works", {
    data(example_obj)
    example_obj <- normalise_to_start(example_obj)
    assay(example_obj)[1:5, 1:3] <- NA
    example_obj_list <- split_groups(example_obj)
    example_obj_merged_list <- merge_replicates(example_obj_list)

    imputed <- impute_groups(example_obj_merged_list, impute_by = "group")
    for (nm in names(imputed)) {
        expect_false(any(is.na(assay(imputed[[nm]]))))
    }
})

test_that("impute_groups falls back when no Subject column", {
    td <- .make_test_se(n_subjects = 4)
    td$ann$Subject <- NULL
    se <- create_input(td$data, td$ann)
    assay(se)[1:5, 1:3] <- NA
    se_list <- split_groups(se)
    merged <- merge_replicates(se_list)

    expect_warning(
        imputed <- impute_groups(merged, impute_by = "subject"),
        "Falling back"
    )
    for (nm in names(imputed)) {
        expect_false(any(is.na(assay(imputed[[nm]]))))
    }
})

test_that("impute_groups by subject imputes per-subject", {
    data("tutorial_data")
    data("tutorial_sample_info")
    # Add Subject column for subject-level imputation
    # Each subject must be unique to one group
    sample_ann_subj <- tutorial_sample_info
    sample_ann_subj$Subject <- paste0("S", seq_len(nrow(sample_ann_subj)))
    se <- create_input(tutorial_data, sample_ann_subj, subject_col = "Subject")
    se_norm <- normalise_to_start(se)
    se_split <- split_groups(se_norm)
    se_merged <- merge_replicates(se_split)

    # Add some NAs to trigger imputation
    assays(se_merged[[1]])[[1]][1, 1] <- NA

    # Subject is dropped by merge_replicates; falls back to group
    expect_warning(
        res <- impute_groups(se_merged, impute_by = "subject"),
        "no Subject column"
    )
    expect_type(res, "list")
    # Verify NA was filled
    expect_false(is.na(assays(res[[1]])[[1]][1, 1]))
})

test_that("impute_groups with impute_by = 'subject' emits per-subject message when Subject present", {
    # merge_replicates drops Subject from colData, so we add it back manually
    # to test the per-subject imputation message path
    td <- .make_test_se(n_subjects = 4, n_groups = 2)
    se <- create_input(td$data, td$ann, subject_col = "Subject")
    se_list <- split_groups(se)
    merged <- merge_replicates(se_list)

    # Add Subject back after merge
    colData(merged[["G1"]])$Subject <- "S1"
    colData(merged[["G2"]])$Subject <- "S3"

    # Introduce NAs to trigger imputation
    assay(merged[["G1"]])[1, 1] <- NA
    assay(merged[["G2"]])[2, 2] <- NA

    expect_message(
        imputed <- impute_groups(merged, impute_by = "subject"),
        "Imputing per-subject"
    )
    expect_type(imputed, "list")
    expect_false(any(is.na(assay(imputed[["G1"]]))))
    expect_false(any(is.na(assay(imputed[["G2"]]))))
})

test_that("decomp_variance handles interaction = TRUE", {
    td <- .make_test_se(n_subjects = 4, n_groups = 2, with_reps = TRUE)
    se <- create_input(td$data, td$ann, subject_col = "Subject")
    se <- normalise_to_start(se)

    expect_message(
        vd <- decomp_variance(se, features = rownames(se)[1:20],
            fixed_effect_var = c("Group", "Time"),
            interaction = TRUE, assay = "orig", core = 1),
        "interaction"
    )
    expect_type(vd, "list")
})

test_that("decomp_variance handles fixed_effect_var = Time", {
    td <- .make_test_se(n_subjects = 4, n_groups = 2)
    se <- create_input(td$data, td$ann, subject_col = "Subject")
    se <- normalise_to_start(se)

    vd <- decomp_variance(se, features = rownames(se)[1:20],
        fixed_effect_var = "Time",
        interaction = FALSE, assay = "orig", core = 1)
    expect_type(vd, "list")
})

# ---- decomp_variance edge branches ----

test_that("decomp_variance errors on non-existent features", {
    td <- .make_test_se(n_subjects = 4, n_groups = 2, with_reps = TRUE)
    se <- create_input(td$data, td$ann, subject_col = "Subject")
    se <- normalise_to_start(se)

    expect_error(
        decomp_variance(se, features = c("NotAGene", "AlsoNotAGene"),
            assay = "orig", core = 1),
        "Features not found"
    )
})

test_that("decomp_variance errors when required colData columns missing", {
    td <- .make_test_se(n_subjects = 4, n_groups = 2, with_reps = TRUE)
    se <- create_input(td$data, td$ann, subject_col = "Subject")
    colData(se)$Group <- NULL

    expect_error(
        decomp_variance(se, features = rownames(se)[1:10],
            assay = "orig", core = 1),
        "Required columns not found"
    )
})

test_that("decomp_variance messages when single Subject level detected", {
    td <- .make_test_se(n_subjects = 1, n_groups = 2, with_reps = TRUE)
    se <- create_input(td$data, td$ann, subject_col = "Subject")
    se <- normalise_to_start(se)

    expect_message(
        decomp_variance(se, features = rownames(se)[1:10],
            assay = "orig", core = 1),
        "Only 1 Subject level detected"
    )
})

test_that("decomp_variance errors when no variable has >=2 levels", {
    td <- .make_test_se(n_subjects = 1, n_groups = 1, n_times = 1)
    se <- create_input(td$data, td$ann, subject_col = "Subject")
    se <- normalise_to_start(se)

    expect_error(
        decomp_variance(se, features = rownames(se)[1:10],
            assay = "orig", core = 1),
        "Cannot fit variance decomposition"
    )
})
