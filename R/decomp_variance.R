#' Variance decomposition
#'
#' @description Variance decomposition by linear mixed models (LMM), to estimate
#' the contribution of Group and Time to the variance of each feature.
#' Modified from `PALMO::lmeVariance()` function.
#'
#' Group and Time are always included as random effects. The formula is
#' extended automatically based on colData:
#' - Subject column present: adds `(1|Subject)` to model between-subject
#'   differences
#' - `interaction = TRUE`: adds `(1|Group:Time)` (or `(1|Subject:Time)`
#'   when Subject present) to capture interaction variance. Requires >= 2
#'   replicates per combination. Default: FALSE.
#'
#' Output always includes: Group, Time, Residual. Subject is included when
#' present. All values are percentages of total variance.
#'
#' Allow using original data or data normalised to time point 0.
#'
#' @param se_obj A SummarizedExperiment object created by `create_input()`
#' @param features A vector of features to include. If NULL, all features.
#' @param fixed_effect_var Column names to include as fixed effects
#'   (regressed out, not appearing as variance components). Default is
#'   `NULL` (no fixed effects). Set to `"Batch"` to regress out batch
#'   effects, or provide other columns (e.g., `c("Sex", "Age")`).
#' @param interaction Logical. If TRUE, adds `(1|Group:Time)` (or
#'   `(1|Subject:Time)` when Subject present) to the model to capture
#'   interaction variance. Default: FALSE.
#' @param assay Assay to use: `"orig"` for original data, `"norm"` for
#' normalised-to-start data. Numeric indices (1, 2) are also accepted.
#' No default, must be specified explicitly.
#' @param core Number of cores for parallel processing (default: 1)
#'
#' @import SummarizedExperiment
#' @returns A data frame with variance decomposition results (percentages).
#' @references https://github.com/aifimmunology/PALMO/blob/main/R/lmeVariance.R
#' @export
#' @examples
#' data(example_obj)
#' example_obj <- normalise_to_start(example_obj)
#'
#' var_decomp <- decomp_variance(example_obj, assay = 1)
decomp_variance <- function(
    se_obj, features = NULL,
    fixed_effect_var = NULL,
    interaction = FALSE,
    assay, core = 1
) {
    .check_se(se_obj)
    .check_logical(interaction, "interaction")
    .check_positive_int(core, "core")

    assay <- .match_assay(assay, se_obj)

    if (!is.null(features)) {
        if (!all(features %in% rownames(se_obj))) {
            stop(sprintf("Features not found: %s",
                paste(setdiff(features, rownames(se_obj)), collapse = ", ")))
        }
        se_obj <- se_obj[features, ]
    }

    required <- c("Group", "Time")
    if (!all(required %in% colnames(colData(se_obj)))) {
        stop(sprintf("Required columns not found: %s",
            paste(setdiff(required, colnames(colData(se_obj))),
                collapse = ", ")))
    }

    # ---- Detect data structure ----
    has_subject <- "Subject" %in% colnames(colData(se_obj))

    ann <- colData(se_obj) |>
        as.data.frame() |>
        dplyr::mutate(Sample_new = if (has_subject) {
            paste0(Group, "_", Subject, "_", Time, "_", Replicate)
        } else {
            paste0(Group, "_", Time, "_", Replicate)
        }) |>
        dplyr::mutate(Time = factor(Time))

    mat <- assays(se_obj)[[assay]] |> as.data.frame()
    colnames(mat) <- ann$Sample_new[match(colnames(mat), ann$Sample)]

    ann <- ann |>
        dplyr::select(-Sample) |>
        dplyr::rename(Sample = Sample_new)
    row.names(ann) <- ann$Sample

    # ---- Build model ----
    # Auto-drop terms with only 1 unique level (e.g., 1-group designs)
    n_groups <- dplyr::n_distinct(ann$Group)
    n_times  <- dplyr::n_distinct(ann$Time)
    n_subjects <- if (has_subject) dplyr::n_distinct(ann$Subject) else 0

    featureSet <- character(0)
    rand_terms <- character(0)

    if (has_subject) {
        if (n_subjects >= 2) {
            featureSet <- "Subject"
        } else {
            message("Only 1 Subject level detected - ",
                "Subject dropped from model.")
        }
    }
    if (n_groups >= 2) {
        featureSet <- c(featureSet, "Group")
    } else {
        message("Only 1 Group level detected - Group dropped from model.")
    }
    if (n_times >= 2) {
        featureSet <- c(featureSet, "Time")
    } else {
        message("Only 1 Time level detected - Time dropped from model.")
    }

    if (length(featureSet) == 0) {
        stop("No variable has >=2 levels. Cannot fit variance decomposition.")
    }

    # ---- Optional interaction term ----
    if (interaction) {
        inter_term <- if (has_subject) "Subject:Time" else "Group:Time"
        n_per_comb <- if (has_subject) {
            table(paste(ann$Subject, ann$Group, ann$Time))
        } else {
            table(paste(ann$Group, ann$Time))
        }
        max_reps <- max(n_per_comb)
        if (max_reps < 2) {
            message("No replicates detected. Interaction term skipped.")
        } else {
            featureSet <- c(featureSet, inter_term)
            if (max_reps <= 2) {
                message("Max ", max_reps, " replicate(s): interaction ",
                    "variance estimate may be unstable. ")
            }
        }
    }
    rand_terms <- paste(paste("(1|", featureSet, ")", sep = ""),
        collapse = " + ")

    # ---- Fixed effects ----
    fixed_terms <- character(0)
    if (!is.null(fixed_effect_var)) {
        missing <- setdiff(fixed_effect_var, colnames(ann))
        if (length(missing) > 0) {
            stop("fixed_effect_var not found in colData: ",
                paste(missing, collapse = ", "))
        }
        fixed_terms <- fixed_effect_var
    }

    if (length(fixed_terms) > 0) {
        form <- paste(paste(fixed_terms, collapse = " + "), "+", rand_terms)
        message(sprintf("Fixed effects: %s (regressed out, not in output).",
            paste(fixed_terms, collapse = ", ")))
    } else {
        form <- rand_terms
    }

    # Define featureList
    featureList <- c(featureSet, "Residual")

    message(sprintf(
        "LMM: exp ~ %s  |  Output: %s",
        form, paste(featureList, collapse = ", ")
    ))

    rowN <- row.names(mat)
    op <- pbapply::pboptions(type = "none")
    on.exit(pbapply::pboptions(op))

    lmem_res <- pbapply::pblapply(seq_len(length(rowN)), cl = core,
        function(gn) {
            geneName <- rowN[gn]
            df <- data.frame(
                exp = as.numeric(mat[geneName, ]), ann,
                stringsAsFactors = FALSE
            )

            # Skip features where any grouping factor has only 1 level
            gene_group <- vapply(seq_len(length(featureSet)), function(i) {
                # Interaction terms (e.g. Group:Time): level check is
                # redundant - if component factors pass, the interaction
                # necessarily has >=1 level with >1 observation
                if (grepl(":", featureSet[i])) return(1)
                df1 <- df[!is.na(df$exp) & !is.na(df[, featureSet[i]]), ]
                res <- data.frame(table(df1[, featureSet[i]]))
                sum(res$Freq > 1)
            }, numeric(1))
            if (any(gene_group == 0)) return(NULL)

            # Skip features with near-zero variance (lmer cannot decompose)
            resp_var <- stats::var(df$exp, na.rm = TRUE)
            if (is.na(resp_var) || resp_var < 1e-12) return(NULL)

            form1 <- stats::as.formula(paste("exp ~ ", form, sep = ""))
            lmem <- suppressMessages(
                tryCatch(
                    lme4::lmer(formula = form1, data = df,
                        control = lme4::lmerControl(
                            optimizer = "bobyqa", calc.derivs = FALSE)),
                    error = function(e) NULL
                )
            )
            if (is.null(lmem)) return(NULL)

            lmem_re <- as.data.frame(lme4::VarCorr(lmem))
            row.names(lmem_re) <- lmem_re$grp
            lmem_re <- lmem_re[featureList, ]

            c(
                geneName, mean(df$exp, na.rm = TRUE),
                stats::median(df$exp, na.rm = TRUE),
                stats::sd(df$exp, na.rm = TRUE),
                max(df$exp, na.rm = TRUE),
                (lmem_re$vcov) / sum(lmem_re$vcov)
            )
        }
    )
    pbapply::pboptions(op)

    lmem_res <- do.call(rbind, lmem_res)
    lmem_res <- data.frame(lmem_res, check.names = FALSE,
        stringsAsFactors = FALSE)
    colnames(lmem_res) <- c("Feature", "mean", "median", "sd", "max",
        featureList)
    lmem_res$mean <- as.numeric(lmem_res$mean)
    lmem_res <- lmem_res[!is.na(lmem_res$mean), ]
    row.names(lmem_res) <- lmem_res$Feature

    # Convert to numeric and sort by Group
    temp <- apply(lmem_res[, -1, drop = FALSE], 1, function(x) as.numeric(x))
    if (!is.matrix(temp)) {
        temp <- t(as.matrix(temp))
        rownames(temp) <- rownames(lmem_res)
    }
    row.names(temp) <- colnames(lmem_res)[-1]
    lmem_res <- data.frame(
        Feature = colnames(temp), t(temp), check.names = FALSE,
        stringsAsFactors = FALSE
    )
    if ("Group" %in% colnames(lmem_res)) {
        lmem_res <- lmem_res[order(lmem_res[["Group"]], decreasing = TRUE), ]
    }
    lmem_res[, featureList] <- 100 * lmem_res[, featureList]

    return(lmem_res)
}
