##################################################
## Project: driveR
## Script purpose: Testthat testing script for
## core functions
## Date: Sep 13, 2020
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

test_that("`create_features_df` argument check works", {
    expect_error(create_features_df(prep_phenolyzer_input = "INVALID"),
                 "`prep_phenolyzer_input` should be logical")
})

# prioritize_driver_genes ----------------------------------------------
test_that("`prioritize_driver_genes` works", {
    expect_is(driverness_df <- prioritize_driver_genes(example_features_table, "LUAD"),
              "data.frame")
})

test_that("`prioritize_driver_genes` argument check works", {
    expect_error(prioritize_driver_genes(example_features_table, "INVALID"),
                 "`cancer_type` should be one of the short names in `MTL_submodel_descriptions`")

    expect_error(prioritize_driver_genes(matrix()),
                 "`features_df` should be a data frame")

    expect_error(prioritize_driver_genes(example_features_table[, -1], "LUAD"),
                 "`features_df` should contain exactly 27 columns")

    req_names <- c("gene_symbol", "metaprediction_score", "noncoding_score",
                   "scna_score", "hotspot_double_hit", "phenolyzer_score",
                   "hsa03320", "hsa04010", "hsa04020", "hsa04024", "hsa04060",
                   "hsa04066", "hsa04110", "hsa04115", "hsa04150", "hsa04151",
                   "hsa04210", "hsa04310", "hsa04330", "hsa04340", "hsa04350",
                   "hsa04370", "hsa04510", "hsa04512", "hsa04520", "hsa04630",
                   "hsa04915")
    tmp <- example_features_table
    colnames(tmp)[1] <- "INVALID"
    expect_error(prioritize_driver_genes(tmp, "LUAD"),
                 paste0("`features_df` should contain the following columns: ",
                        paste(dQuote(req_names), collapse = ", ")))
})
