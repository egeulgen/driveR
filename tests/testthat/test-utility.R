##################################################
## Project: driveR
## Script purpose: Testthat testing script for
## utility functions
## Date: Sep 6, 2020
## Author: Ege Ulgen
##################################################

# create_metaprediction_score_df ------------------------------------------
test_that("`create_metaprediction_score_df` works", {
    path2annovar_csv <- system.file("extdata/example.hg19_multianno.csv",
                                    package = "driveR")
    expect_is(metapred_df <- driveR:::create_metaprediction_score_df(path2annovar_csv), "data.frame")
    expect_equal(ncol(metapred_df), 2)
})

# create_noncoding_impact_score_df ----------------------------------------
test_that("`create_noncoding_impact_score_df` works", {
    path2annovar_csv <- system.file("extdata/example.hg19_multianno.csv",
                                    package = "driveR")
    expect_is(noncoding_df <- driveR:::create_noncoding_impact_score_df(path2annovar_csv), "data.frame")
    expect_equal(ncol(noncoding_df), 2)
})

# create_gene_level_scna_df -----------------------------------------------
test_that("`create_gene_level_scna_df` works", {
    expect_is(driveR:::create_gene_level_scna_df(example_scna_table), "data.frame")
})

test_that("`create_gene_level_scna_df` argument checks work", {
    tmp <- example_scna_table[, 1:3]
    nec_cols <- c("chr", "start", "end", "log2ratio")
    expect_error(driveR:::create_gene_level_scna_df(tmp),
                 paste0("`scna_df` should contain all of: ",
                        paste(dQuote(nec_cols), collapse = ", ")))

    expect_error(driveR:::create_gene_level_scna_df(example_scna_table,
                                                    gene_overlap_threshold = "INVALID"),
                 "`gene_overlap_threshold` should be numberic")

    expect_error(driveR:::create_gene_level_scna_df(example_scna_table,
                                                    gene_overlap_threshold = -1),
                 "`gene_overlap_threshold` should be between 0-100")
})

# create_SCNA_score_df ----------------------------------------------------
test_that("`create_SCNA_score_df` works", {
    expect_is(SCNA_scores_df <- driveR:::create_SCNA_score_df(example_scna_table), "data.frame")
    expect_equal(ncol(SCNA_scores_df), 2)

    # corner case
    expect_is(SCNA_scores_df <- driveR:::create_SCNA_score_df(example_scna_table[1:10, ]), "data.frame")
    expect_equal(ncol(SCNA_scores_df), 2)
    expect_equal(nrow(SCNA_scores_df), 0)
})

test_that("`create_SCNA_score_df` argument checks work", {
    expect_error(driveR:::create_SCNA_score_df(example_scna_table,
                                               log2_ratio_threshold = "INVALID"),
                 "`log2_ratio_threshold` should be numberic")


    expect_error(driveR:::create_SCNA_score_df(example_scna_table,
                                               MCR_overlap_threshold = "INVALID"),
                 "`MCR_overlap_threshold` should be numberic")

    expect_error(driveR:::create_SCNA_score_df(example_scna_table,
                                                    MCR_overlap_threshold = -1),
                 "`MCR_overlap_threshold` should be between 0-100")
})

# determine_hotspot_genes -------------------------------------------------
test_that("`determine_hotspot_genes` works", {
    path2annovar_csv <- system.file("extdata/example.hg19_multianno.csv",
                                    package = "driveR")
    expect_is(hotspot_genes <- driveR:::determine_hotspot_genes(path2annovar_csv), "character")
})

test_that("`determine_hotspot_genes` arg check works", {
    expect_error(driveR:::determine_hotspot_genes(path2annovar_csv,
                                                  hotspot_threshold = "INVALID"),
                 "`hotspot_threshold` should be numeric")
})

# determine_double_hit_genes ----------------------------------------------
test_that("`determine_double_hit_genes` works", {
    path2annovar_csv <- system.file("extdata/example.hg19_multianno.csv",
                                    package = "driveR")
    expect_is(dhit_genes <- driveR:::determine_double_hit_genes(path2annovar_csv,
                                                                example_scna_table),
              "character")

    path2annovar_csv <- system.file("extdata/example_cohort.hg19_multianno.csv",
                                    package = "driveR")
    expect_is(dhit_genes <- driveR:::determine_double_hit_genes(path2annovar_csv,
                                                                example_cohort_scna_table),
              "character")
})

test_that("`determine_double_hit_genes` argument checks work", {
    path2annovar_csv <- system.file("extdata/example.hg19_multianno.csv",
                                    package = "driveR")
    expect_error(driveR:::determine_double_hit_genes(path2annovar_csv,
                                                     example_scna_table,
                                                     log2_hom_loss_threshold = "INVALID"),
                 "`log2_hom_loss_threshold` should be numberic")
})
