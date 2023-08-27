test_that("`create_noncoding_impact_score_df` -- works as expected", {
    path2annovar_csv <- system.file("extdata/example.hg19_multianno.csv", package = "driveR")
    expect_is(noncoding_df <- create_noncoding_impact_score_df(path2annovar_csv),
        "data.frame")
    expect_equal(ncol(noncoding_df), 2)
})

test_that("`create_gene_level_scna_df` -- works as expected", {
    expect_is(create_gene_level_scna_df(example_scna_table), "data.frame")
    expect_is(create_gene_level_scna_df(example_scna_table, build = "GRCh38"), "data.frame")

})

test_that("`create_gene_level_scna_df` -- argument checks work", {

    expect_error(create_gene_level_scna_df(example_scna_table, build = "INVALID"),
        paste0("`build should be one of ", paste(dQuote(c("GRCh37", "GRCh38")), collapse = ", ")))

    tmp <- example_scna_table[, 1:3]
    nec_cols <- c("chr", "start", "end", "log2ratio")
    expect_error(create_gene_level_scna_df(tmp), paste0("`scna_segs_df` should contain all of: ",
        paste(dQuote(nec_cols), collapse = ", ")))

    expect_error(create_gene_level_scna_df(example_scna_table, gene_overlap_threshold = "INVALID"),
        "`gene_overlap_threshold` should be numberic")

    expect_error(create_gene_level_scna_df(example_scna_table, gene_overlap_threshold = -1),
        "`gene_overlap_threshold` should be between 0-100")
})

gene_SCNA_df_batch <- create_gene_level_scna_df(example_cohort_scna_table)

test_that("`create_SCNA_score_df` -- works as expected", {
    expect_is(SCNA_scores_df <- create_SCNA_score_df(example_gene_scna_table), "data.frame")
    expect_equal(ncol(SCNA_scores_df), 2)

    expect_is(SCNA_scores_df2 <- create_SCNA_score_df(example_gene_scna_table, build = "GRCh38"),
        "data.frame")
    expect_equal(ncol(SCNA_scores_df2), 2)

    # corner case 1
    tmp <- example_gene_scna_table[1:2, ]
    tmp$log2ratio <- 0.1
    expect_warning(SCNA_scores_df <- create_SCNA_score_df(tmp))
    expect_equal(ncol(SCNA_scores_df), 2)
    expect_equal(nrow(SCNA_scores_df), 0)

    # corner case 2
    expect_warning(SCNA_scores_df <- create_SCNA_score_df(example_gene_scna_table[1:2,
        ]))
    expect_equal(ncol(SCNA_scores_df), 2)
    expect_equal(nrow(SCNA_scores_df), 0)
})

test_that("`create_SCNA_score_df` -- argument checks work", {

    expect_error(create_SCNA_score_df(example_gene_scna_table, build = "INVALID"),
        paste0("`build should be one of ", paste(dQuote(c("GRCh37", "GRCh38")), collapse = ", ")))


    expect_error(create_SCNA_score_df(example_gene_scna_table, log2_ratio_threshold = "INVALID"),
        "`log2_ratio_threshold` should be numberic")


    expect_error(create_SCNA_score_df(example_gene_scna_table, MCR_overlap_threshold = "INVALID"),
        "`MCR_overlap_threshold` should be numberic")

    expect_error(create_SCNA_score_df(example_gene_scna_table, MCR_overlap_threshold = -1),
        "`MCR_overlap_threshold` should be between 0-100")
})

test_that("`determine_hotspot_genes` -- works as expected", {
    path2annovar_csv <- system.file("extdata/example.hg19_multianno.csv", package = "driveR")
    expect_is(hotspot_genes <- determine_hotspot_genes(path2annovar_csv), "character")
})

test_that("`determine_hotspot_genes` -- argument check works", {
    expect_error(determine_hotspot_genes(path2annovar_csv, hotspot_threshold = "INVALID"),
        "`hotspot_threshold` should be numeric")
})

test_that("`determine_double_hit_genes` -- works as expected", {
    path2annovar_csv <- system.file("extdata/example.hg19_multianno.csv", package = "driveR")
    expect_is(dhit_genes <- determine_double_hit_genes(path2annovar_csv, example_gene_scna_table),
        "character")

    path2annovar_csv <- system.file("extdata/example_cohort.hg19_multianno.csv",
        package = "driveR")
    expect_is(dhit_genes <- determine_double_hit_genes(path2annovar_csv, gene_SCNA_df_batch,
        batch_analysis = TRUE), "character")
})

test_that("`determine_double_hit_genes` -- argument checks work", {
    path2annovar_csv <- system.file("extdata/example.hg19_multianno.csv", package = "driveR")
    expect_error(determine_double_hit_genes(path2annovar_csv, example_gene_scna_table,
        log2_hom_loss_threshold = "INVALID"), "`log2_hom_loss_threshold` should be numberic")

    expect_error(determine_double_hit_genes(path2annovar_csv, example_gene_scna_table,
        batch_analysis = "INVALID"), "`batch_analysis` should be `TRUE` or `FALSE`")

    expect_error(determine_double_hit_genes(path2annovar_csv, example_gene_scna_table,
        batch_analysis = TRUE), "'tumor id' should be present in both ANNOVAR output and SCNA table if `batch_analysis == TRUE`")

})
