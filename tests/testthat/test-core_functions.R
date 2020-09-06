##################################################
## Project: driveR
## Script purpose: Testthat testing script for
## core functions
## Date: Sep 6, 2020
## Author: Ege Ulgen
##################################################

# create_features_df ------------------------------------------------------
test_that("`create_features_df` works", {
    path2annovar_csv <- system.file("extdata/example.hg19_multianno.csv",
                                    package = "driveR")
    path2phenolyzer_out <- system.file("extdata/example.annotated_gene_list",
                                       package = "driveR")

    # prep_phenolyzer_input = FALSE
    expect_is(features_df <- create_features_df(annovar_csv_path = path2annovar_csv,
                                                scna_df = example_scna_table,
                                                phenolyzer_annotated_gene_list_path = path2phenolyzer_out),
              "data.frame")
    expect_equal(ncol(features_df), 27)

    # prep_phenolyzer_input = TRUE
    expect_is(create_features_df(annovar_csv_path = path2annovar_csv,
                                 scna_df = example_scna_table,
                                 prep_phenolyzer_input = TRUE),
              "character")
})
