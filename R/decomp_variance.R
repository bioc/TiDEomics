#' Variance decomposition
#'
#' @description Variance decomposition by linear mixed models, to estimate
#' the contribution of different variables to the variance of each feature.
#' Modified from `PALMO::lmeVariance()` function.
#'
#' Allow using original data or data normalised to time point 0
#'
#' @param se_obj A SummarizedExperiment object created by `create_input()`
#' @param features A vector of features (e.g. genes) to include in the
#' analysis. If NULL, all features are used. (default is NULL)
#' @param variables Variables to be included in linear mixed model for variance
#' analysis. (default is c("Group", "Time"))
#' @param fixed_effect_var Fixed effect variables to be included in linear
#' mixed model, variance contribution obtained by adding them as random
#' variables. (default is NULL)
#' @param assay 1 for the original input data, or 2 for the data normalised
#' to time point 0
#' @param core Number of cores to use for parallel processing (default is 2)
#'
#' @import SummarizedExperiment
#' @import magrittr
#' @returns A data frame with variance decomposition results
#' @references https://github.com/aifimmunology/PALMO/blob/main/R/lmeVariance.R
#' @export
#' @examples
#' data("example")
#' example_obj <- normalise_to_start(example_obj)
#'
#' var_decomp <- decomp_variance(example_obj, assay = 1)
#' plot_variance(var_decomp, rank = "Time", top_n = 20)
decomp_variance <- function(
    se_obj, features = NULL, variables = c("Group", "Time"),
    fixed_effect_var = NULL,
    assay = c(1, 2), core = 2
) {
    if (!is.null(features)) {
        if (!all(features %in% rownames(se_obj))) {
            stop(sprintf("Features not found in the input: %s",
                paste(setdiff(features, rownames(se_obj)), collapse = ", ")))
        } else {
            palmo_obj <- se_obj[features, ]
        }
    } else {
        palmo_obj <- se_obj
    }

    if (!all(variables %in% colnames(colData(palmo_obj)))) {
        stop(sprintf("Variables not found in the input: %s",
            paste(setdiff(variables, colnames(colData(palmo_obj))),
                collapse = ", ")))
    }

    ann <- colData(palmo_obj) %>%
        as.data.frame() %>%
        dplyr::mutate(Sample_new = paste0(.data$Group, "_", .data$Time, "_",
            .data$Replicate)) %>%
        dplyr::mutate(Time = factor(.data$Time))

    mat <- assays(palmo_obj)[[assay]] %>% as.data.frame()
    colnames(mat) <- colnames(mat) %>%
        plyr::mapvalues(from = ann$Sample, to = ann$Sample_new)

    ann <- ann %>%
        dplyr::select(-Sample) %>%
        dplyr::rename(Sample = Sample_new)
    row.names(ann) <- ann$Sample

    ## Define formula
    featureSet <- variables
    lmer_control <- TRUE
    form <- paste(paste("(1|", featureSet, ")", sep = ""), collapse = " + ")

    ## check fixed effect variables
    if (!is.null(fixed_effect_var)) {
        check_fixedeffect <- intersect(fixed_effect_var, colnames(ann))
        if (length(check_fixedeffect) > 0) {
            featureSet <- c(featureSet, check_fixedeffect)
            form <- paste(paste("(1|", featureSet, ")", sep = ""),
                collapse = " + ")
            form <- paste(paste(check_fixedeffect, sep = "", collapse = " + "),
                " + ", form,
                sep = "", collapse = " + "
            )
        } else {
            stop("Fixed effect variables not found")
        }
    }
    # data_object@result$var_formula <- form

    ## Define featureList
    featureList <- c(featureSet, "Residual")
    rowN <- row.names(mat)
    op <- pbapply::pboptions(type = "timer") # default
    suppressMessages(
        lmem_res <- pbapply::pblapply(seq(1, length(rowN)), cl = core,
        function(gn) {
            geneName <- rowN[gn]
            df <- data.frame(
                exp = as.numeric(mat[geneName, ]), ann,
                stringsAsFactors = FALSE
            )

            ## remove features with not enough values
            gene_group <- lapply(seq(1, length(featureSet)), function(i) {
                df1 <- df[!is.na(df$exp) & !is.na(df[, featureSet[i]]), ]
                res <- data.frame(table(df1[, featureSet[i]]))
                res <- res[res$Freq > 1, ]
                return(nrow(res))
            })
            gene_group <- as.numeric(gene_group)
            ## Check any of attribute has only one group
            check_gene_group <- sum(gene_group == 1)
            if (check_gene_group == 0) {
                ## Define formula form
                form1 <- stats::as.formula(paste("exp ~ ", form, sep = ""))

                ## linear mixed effect model
                if (lmer_control == TRUE) {
                    lmem <- lme4::lmer(
                        formula = form1, data = df,
                        control = lme4::lmerControl(
                            optimizer = "bobyqa",
                            calc.derivs = FALSE
                        )
                    )
                } else {
                    lmem <- lme4::lmer(formula = form1, data = df)
                }

                lmem_re <- as.data.frame(lme4::VarCorr(lmem))
                row.names(lmem_re) <- lmem_re$grp
                lmem_re <- lmem_re[featureList, ]
                fix_effect <- lme4::fixef(lmem) # get fixed effect
                lmem_re$CV <- lmem_re$sdcor / fix_effect ## Calculate CV
                return(c(
                    geneName, mean(df$exp, na.rm = TRUE),
                    stats::median(df$exp, na.rm = TRUE), 
                    stats::sd(df$exp, na.rm = TRUE),
                    max(df$exp, na.rm = TRUE),
                    (lmem_re$vcov) / sum(lmem_re$vcov)
                ))
            }
        }
    ) , classes = "message")
    pbapply::pboptions(op)
    lmem_res <- do.call(rbind, lmem_res)
    lmem_res <- data.frame(lmem_res, check.names = FALSE,
        stringsAsFactors = FALSE)
    colnames(lmem_res) <- c("Feature", "mean", "median", "sd", "max",
        featureList)
    lmem_res$mean <- as.numeric(lmem_res$mean)
    lmem_res <- lmem_res[!is.na(lmem_res$mean), ]
    row.names(lmem_res) <- lmem_res$Feature
    # Some data converted into character hence convert into numeric
    temp <- apply(lmem_res[, -1], 1, function(x) {
        as.numeric(x)
    })
    row.names(temp) <- colnames(lmem_res)[-1]
    lmem_res <- data.frame(
        Feature = colnames(temp), t(temp), check.names = FALSE,
        stringsAsFactors = FALSE
    )
    lmem_res <- lmem_res[order(lmem_res[, featureList[1]], decreasing = TRUE), ]

    ## Return object
    lmem_res[, featureList] <- 100 * lmem_res[, featureList] # in percentage
    # data_object@result$variance_decomposition <- lmem_res

    var_decomp <- lmem_res

    return(var_decomp)
}
