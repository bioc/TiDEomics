# ---- extract_consensus_hubs ----

data(example_cons)

test_that("extract_consensus_hubs returns consensus and per-dataset hubs", {
    colors <- example_cons$network$colors
    n_mod <- length(unique(colors[colors != 0]))

    hubs <- extract_consensus_hubs(example_cons, top_n = 1)
    expect_s3_class(hubs, "data.frame")
    expect_named(hubs,
        c("Feature", "Module", "Set", "Rank", "Score", "ScoreType"))
    expect_equal(as.integer(table(hubs$Set)),
        c(n_mod, n_mod, n_mod))
    expect_setequal(unique(hubs$Set),
        c("Consensus hub", "Exp1 hub", "Exp2 hub"))
    expect_false(any(colors[hubs$Feature] %in% c("grey", "0", 0)))
    expect_equal(hubs$Rank, rep(1L, nrow(hubs)))

    # top_n = 2 returns more hubs from every set; overlaps are retained
    hubs2 <- extract_consensus_hubs(example_cons, top_n = 2)
    expect_equal(as.integer(table(hubs2$Set)),
        c(2 * n_mod, 2 * n_mod, 2 * n_mod))
    expect_true(all(hubs2$Feature %in% example_cons$common_features))
    expect_true(all(hubs2$Feature[hubs2$Set == "Exp1 hub"] %in%
        colnames(example_cons$multiExpr$Exp1$data)))
    expect_true(all(hubs2$Feature[hubs2$Set == "Exp2 hub"] %in%
        colnames(example_cons$multiExpr$Exp2$data)))
    expect_false(any(colors[hubs2$Feature] %in% c("grey", "0", 0)))
})

test_that("extract_consensus_hubs handles exclude_grey", {
    colors <- example_cons$network$colors
    n_mod <- length(unique(colors[colors != 0]))
    n_grey_mod <- if (any(colors == 0)) 1 else 0

    hubs_grey <- extract_consensus_hubs(example_cons, top_n = 1,
        exclude_grey = FALSE)
    expect_equal(nrow(hubs_grey),
        3 * (n_mod + n_grey_mod))
    expect_equal(sum(hubs_grey$Module == "0"), 3 * n_grey_mod)
})

test_that("extract_consensus_hubs validates input", {
    expect_error(extract_consensus_hubs(list()),
        "'consensus\\$network' must be a list")

    cons <- example_cons
    cons$parameters$networkType <- NULL
    expect_error(extract_consensus_hubs(cons),
        "must contain 'networkType'")
})

test_that("extract_consensus_hubs handles edge-case networks", {
    # multiExpr length must match the number of module eigengene sets
    cons <- example_cons
    cons$multiExpr$Exp3 <- cons$multiExpr$Exp1
    expect_error(extract_consensus_hubs(cons))

    # all-grey network: no modules, empty hubs
    cons <- example_cons
    cons$network$colors[] <- 0
    expect_message(hubs <- extract_consensus_hubs(cons),
        "No modules to extract hubs from")
    expect_s3_class(hubs, "data.frame")
    expect_named(hubs,
        c("Feature", "Module", "Set", "Rank", "Score", "ScoreType"))
    expect_equal(nrow(hubs), 0L)

    # unsigned networks use absolute kME values
    cons <- example_cons
    cons$parameters$networkType <- "unsigned"
    hubs <- extract_consensus_hubs(cons, top_n = 1)
    expect_true(nrow(hubs) >= 1)
    expect_true(all(grepl("abs", hubs$ScoreType)))

    # a module with no signedKME column triggers a missing-column error
    cons <- example_cons
    cons$network$colors[1] <- 99
    expect_error(extract_consensus_hubs(cons))
})

# ---- run_consensus_WGCNA ----

data(example_obj)
data(example_obj_2)
shared_features <- intersect(
    rownames(example_obj), rownames(example_obj_2)
)
nets <- list(
    Exp1 = prepare_WGCNA(example_obj[shared_features, , drop = FALSE],
        powers = seq(1, 20), assay = 2, networkType = "signed"),
    Exp2 = prepare_WGCNA(example_obj_2[shared_features, , drop = FALSE],
        powers = seq(1, 20), assay = 2, networkType = "signed")
)

test_that("run_consensus_WGCNA validates arguments", {
    expect_error(run_consensus_WGCNA(list(), powers = c(Exp1 = 6)),
        "must contain at least 2 networks")
    expect_error(run_consensus_WGCNA(list(nets$Exp1, nets$Exp2),
        powers = c(Exp1 = 6, Exp2 = 6)),
        "non-empty names for all elements")
    expect_error(run_consensus_WGCNA("bad", powers = c(Exp1 = 6)),
        "must be a list")
    expect_error(run_consensus_WGCNA(list(Exp1 = 1, Exp2 = 2),
        powers = c(Exp1 = 6, Exp2 = 6)),
        "must be a list, the output of prepare_WGCNA")
    expect_error(run_consensus_WGCNA(nets, powers = c(6, 6)),
        "one uniquely named value for each dataset")
    expect_error(run_consensus_WGCNA(nets, powers = c(Exp1 = 6)),
        "one uniquely named value for each dataset")
})

test_that("run_consensus_WGCNA validates numeric arguments", {
    nets <- list(Exp1 = list(), Exp2 = list())

    expect_error(run_consensus_WGCNA(nets,
        powers = c(Exp1 = 6, Exp2 = 6), plot = "yes"),
        "TRUE or FALSE")
    expect_error(run_consensus_WGCNA(nets,
        powers = c(Exp1 = 6, Exp2 = 6), numericLabels = "yes"),
        "TRUE or FALSE")
})

test_that("run_consensus_WGCNA validates minModuleSize", {
    nets <- list(
        Exp1 = list(data = data.frame(a = 1:2, b = 3:4),
            networkType = "signed"),
        Exp2 = list(data = data.frame(a = 1:2, b = 3:4),
            networkType = "signed")
    )

    expect_error(run_consensus_WGCNA(nets,
        powers = c(Exp1 = 6, Exp2 = 6), minModuleSize = "x"),
        "positive integer")
    expect_error(run_consensus_WGCNA(nets,
        powers = c(Exp1 = 6, Exp2 = 6), minModuleSize = 1.5),
        "positive integer")
    expect_error(run_consensus_WGCNA(nets,
        powers = c(Exp1 = 6, Exp2 = 6), minModuleSize = 0),
        "positive integer")
})

test_that("run_consensus_WGCNA handles networkType from nets and dots", {
    expect_error(run_consensus_WGCNA(nets,
        powers = c(Exp1 = 8, Exp2 = 8), networkType = "unsigned"),
        "does not match the networkType stored")
    expect_error(run_consensus_WGCNA(nets,
        powers = c(Exp1 = 8, Exp2 = 8), networkType = 123),
        "does not match the networkType stored")

    nets_bad <- nets
    nets_bad$Exp1$networkType <- 123
    nets_bad$Exp2$networkType <- 123
    expect_error(run_consensus_WGCNA(nets_bad,
        powers = c(Exp1 = 8, Exp2 = 8)),
        "same non-empty 'networkType'")

    nets_bad <- nets
    nets_bad$Exp1$data <- data.frame()
    expect_error(run_consensus_WGCNA(nets_bad,
        powers = c(Exp1 = 8, Exp2 = 8), networkType = "signed"),
        "must be a non-empty data.frame or matrix")
})

test_that("run_consensus_WGCNA stops when too few common features", {
    nets_bad <- nets
    colnames(nets_bad$Exp2$data) <- paste0("x_", colnames(nets_bad$Exp2$data))

    expect_error(run_consensus_WGCNA(nets_bad,
        powers = c(Exp1 = 8, Exp2 = 8)),
        "Too few common features")
    expect_error(run_consensus_WGCNA(nets_bad,
        powers = c(Exp1 = 8, Exp2 = 8), minModuleSize = 50),
        "Too few common features")
})

test_that("run_consensus_WGCNA runs and returns consensus network", {
    expect_message(cons <- run_consensus_WGCNA(
        nets = nets,
        powers = c(Exp1 = 8, Exp2 = 8),
        minModuleSize = 10, plot = FALSE,
        useDiskCache = FALSE),
        "Running blockwise consensus WGCNA")

    expect_named(cons, c("network", "multiExpr", "common_features",
        "names", "parameters"))
    expect_equal(cons$names, c("Exp1", "Exp2"))
    expect_equal(names(cons$multiExpr), c("Exp1", "Exp2"))
    expect_true(length(cons$common_features) >= 10)
    expect_equal(ncol(cons$multiExpr$Exp1$data),
        length(cons$common_features))
    expect_equal(length(cons$network$colors),
        length(cons$common_features))
    expect_equal(cons$parameters$networkType, "signed")
    expect_equal(cons$parameters$powers,
        c(Exp1 = 8, Exp2 = 8))
    expect_equal(cons$parameters$minModuleSize, 10)
    expect_equal(cons$parameters$useDiskCache, FALSE)
    expect_equal(cons$parameters$numericLabels, TRUE)
    expect_equal(cons$network$consensusQuantile, 0)

    hubs <- extract_consensus_hubs(cons, top_n = 2)
    expect_s3_class(hubs, "data.frame")
    expect_true(all(hubs$Feature %in% cons$common_features))
})

test_that("run_consensus_WGCNA warns when feature sets differ", {
    nets_bad <- nets
    cn <- colnames(nets_bad$Exp2$data)
    colnames(nets_bad$Exp2$data) <- paste0("x_", cn)

    expect_warning(
        expect_error(run_consensus_WGCNA(nets_bad,
            powers = c(Exp1 = 8, Exp2 = 8), minModuleSize = 100),
            "Too few common features"),
        "Feature sets differ across datasets"
    )
})

test_that("run_consensus_WGCNA plots with numericLabels = FALSE", {
    expect_no_error(run_consensus_WGCNA(
        nets = nets,
        powers = c(Exp1 = 8, Exp2 = 8),
        minModuleSize = 10, plot = TRUE, numericLabels = FALSE,
        useDiskCache = FALSE))
})
