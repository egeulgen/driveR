##################################################
## Project: driveR
## Script purpose: Testthat testing script for
## core functions
## Date: Oct 1, 2020
## Author: Ege Ulgen
##################################################

# predict_coding_impact ---------------------------------------------------s
test_that("`predict_coding_impact` works", {
    path2annovar_csv <- system.file("extdata/example.hg19_multianno.csv",
                                    package = "driveR")
    expect_is(metapred_df <- predict_coding_impact(path2annovar_csv), "data.frame")
    expect_equal(ncol(metapred_df), 2)

    # corner case
    tmp <- read.csv(path2annovar_csv)
    tmp <- tmp[tmp$CADD_phred == ".", ]
    path2corner <- tempfile()
    write.csv(tmp, path2corner, row.names = FALSE)
    expect_is(metapred_df <- predict_coding_impact(path2corner), "data.frame")
    expect_equal(nrow(metapred_df), 0)
    expect_equal(ncol(metapred_df), 2)
})

test_that("`predict_coding_impact` argument checks work", {
    expect_error(predict_coding_impact("invalid/path/to/csv"),
                 "The file used for `annovar_csv_path` does not exist")

    path2annovar_csv <- system.file("extdata/example.hg19_multianno.csv",
                                    package = "driveR")
    annovar_df <- read.csv(path2annovar_csv)
    annovar_df$Gene.refGene <- NULL
    path2annovar_csv <- tempfile()
    write.csv(annovar_df, path2annovar_csv, row.names = FALSE)
    nec_cols <- c("Gene.refGene",
                  "SIFT_score", "Polyphen2_HDIV_score", "LRT_score",
                  "MutationTaster_score", "MutationAssessor_score",
                  "FATHMM_score", "GERP.._RS", "phyloP7way_vertebrate",
                  "CADD_phred", "VEST3_score", "SiPhy_29way_logOdds",
                  "DANN_score")
    expect_error(predict_coding_impact(path2annovar_csv),
                 paste0("The table in `annovar_csv_path` should contain all of the following columns: ",
                        paste(dQuote(nec_cols), collapse = ", ")))

    path2annovar_csv <- system.file("extdata/example.hg19_multianno.csv",
                                    package = "driveR")
    expect_error(predict_coding_impact(path2annovar_csv,
                                       keep_highest_score = "INVALID"),
                 "`keep_highest_score` should be logical")

    expect_error(predict_coding_impact(path2annovar_csv,
                                       keep_single_symbol = "INVALID"),
                 "`keep_single_symbol` should be logical")
})

# create_features_df ------------------------------------------------------
test_that("`create_features_df` works", {
    # personalized analysis
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


    # batch analysis
    path2annovar_csv <- system.file("extdata/example_cohort.hg19_multianno.csv",
                                    package = "driveR")
    expect_is(features_df <- create_features_df(annovar_csv_path = path2annovar_csv,
                                                scna_df = example_cohort_scna_table,
                                                phenolyzer_annotated_gene_list_path = path2phenolyzer_out,
                                                batch_analysis = TRUE),
              "data.frame")
    expect_equal(ncol(features_df), 27)

    # corner case 1
    tmp <- read.csv(path2annovar_csv)
    tmp <- tmp[tmp$CADD_phred == ".", ]
    path2corner <- tempfile()
    write.csv(tmp, path2corner, row.names = FALSE)
    expect_is(features_df <- create_features_df(annovar_csv_path = path2corner,
                                                scna_df = example_scna_table,
                                                phenolyzer_annotated_gene_list_path = path2phenolyzer_out),
              "data.frame")
    expect_equal(ncol(features_df), 27)

    # corner case 2
    tmp <- tmp[1:2, ]
    tmp$cosmic91_coding <- tmp$cosmic91_noncoding <- "."
    path2corner <- tempfile()
    write.csv(tmp, path2corner, row.names = FALSE)
    tmp_scna_table <- example_scna_table[1:2, ]
    tmp_scna_table$log2ratio <- .1
    expect_warning(features_df <- create_features_df(annovar_csv_path = path2corner,
                                                     scna_df = tmp_scna_table,
                                                     phenolyzer_annotated_gene_list_path = path2phenolyzer_out))

})

test_that("`create_features_df` argument check works", {
    expect_error(create_features_df(prep_phenolyzer_input = "INVALID"),
                 "`prep_phenolyzer_input` should be logical")
})

# prioritize_driver_genes ----------------------------------------------
test_that("`prioritize_driver_genes` works", {
    expect_is(drivers_df <- prioritize_driver_genes(example_features_table, "LUAD"),
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
