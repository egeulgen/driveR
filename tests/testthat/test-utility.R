##################################################
## Project: driveR
## Script purpose: Testthat testing script for
## utility functions
## Date: Sep 4, 2020
## Author: Ege Ulgen
##################################################

# create_metaprediction_score_df ------------------------------------------
test_that("`create_metaprediction_score_df` works", {
    path2annovar_csv <- system.file("extdata/imielinski.hg19_multianno.csv",
                                    package = "driveR")
    expect_is(metapred_df <- driveR:::create_metaprediction_score_df(path2annovar_csv), "data.frame")
    expect_equal(ncol(metapred_df), 2)
})


# create_SCNA_score_df ----------------------------------------------------
test_that("`create_SCNA_score_df` works", {
    expect_is(SCNA_scores_df <- driveR:::create_SCNA_score_df(imielinski_scna_table), "data.frame")
    expect_equal(ncol(SCNA_scores_df), 2)
})

test_that("`create_SCNA_score_df` format check works", {
    tmp <- imielinski_scna_table[, 1:3]
    nec_cols <- c("chr", "start", "end", "log2ratio")
    expect_error(driveR:::create_SCNA_score_df(tmp),
                 paste0("`scna_df` should contain all of: ",
                        paste(dQuote(nec_cols), collapse = ", ")))
})

# determine_hotspot_genes -------------------------------------------------
test_that("`determine_hotspot_genes` works", {
    path2annovar_csv <- system.file("extdata/imielinski.hg19_multianno.csv",
                                    package = "driveR")
    expect_is(hotspot_genes <- driveR:::determine_hotspot_genes(path2annovar_csv), "character")
})

test_that("`determine_hotspot_genes` arg check works", {
    expect_error(driveR:::determine_hotspot_genes(path2annovar_csv, threshold = "INVALID"),
                 "`threshold` should be numeric")
})
