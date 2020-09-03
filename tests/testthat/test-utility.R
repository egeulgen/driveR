test_that("metapredictor utility function works", {
    path2annovar_csv <- system.file("extdata/imielinski.hg19_multianno.csv",
                                    package = "driveR")
    expect_is(metapred_df <- driveR:::create_metaprediction_score_df(path2annovar_csv), "data.frame")
    expect_equal(ncol(metapred_df), 2)
})
