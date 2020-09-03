#' Create Coding Impact Meta-prediction Score Data Frame
#'
#' @param annovar_csv_path path to ANNOVAR csv output file
#'
#' @return data frame of meta-prediction scores containing 2 cclumns: \describe{
#'   \item{gene_symbol}{ID of the enriched term}
#'   \item{metaprediction_score}{the metapredictor impact score}
#' }
#'
#' @examples
#' path2annovar_csv <- system.file("extdata/imielinski.hg19_multianno.csv",
#'                                 package = "driveR")
#' metapred_df <- driveR:::create_metaprediction_score_df(path2annovar_csv)
create_metaprediction_score_df <- function(annovar_csv_path) {
    anno_bechmark_df <- utils::read.csv(annovar_csv_path)

    # filter out missing scores
    anno_bechmark_df <- anno_bechmark_df[anno_bechmark_df$SIFT_score != ".", ]
    anno_bechmark_df$SIFT_score <- as.numeric(anno_bechmark_df$SIFT_score)

    anno_bechmark_df <- anno_bechmark_df[anno_bechmark_df$Polyphen2_HDIV_score != ".", ]
    anno_bechmark_df$Polyphen2_HDIV_score <- as.numeric(anno_bechmark_df$Polyphen2_HDIV_score)

    anno_bechmark_df <- anno_bechmark_df[anno_bechmark_df$LRT_score != ".", ]
    anno_bechmark_df$LRT_score <- as.numeric(anno_bechmark_df$LRT_score)

    anno_bechmark_df <- anno_bechmark_df[anno_bechmark_df$MutationTaster_score != ".", ]
    anno_bechmark_df$MutationTaster_score <- as.numeric(anno_bechmark_df$MutationTaster_score)

    anno_bechmark_df <- anno_bechmark_df[anno_bechmark_df$MutationAssessor_score != ".", ]
    anno_bechmark_df$MutationAssessor_score <- as.numeric(anno_bechmark_df$MutationAssessor_score)

    anno_bechmark_df <- anno_bechmark_df[anno_bechmark_df$FATHMM_score != ".", ]
    anno_bechmark_df$FATHMM_score <- as.numeric(anno_bechmark_df$FATHMM_score)

    anno_bechmark_df <- anno_bechmark_df[anno_bechmark_df$GERP.._RS != ".", ]
    anno_bechmark_df$GERP.._RS <- as.numeric(anno_bechmark_df$GERP.._RS)

    anno_bechmark_df <- anno_bechmark_df[anno_bechmark_df$phyloP7way_vertebrate != ".", ]
    anno_bechmark_df$phyloP7way_vertebrate <- as.numeric(anno_bechmark_df$phyloP7way_vertebrate)

    anno_bechmark_df <- anno_bechmark_df[anno_bechmark_df$CADD_phred != ".", ]
    anno_bechmark_df$CADD_phred <- as.numeric(anno_bechmark_df$CADD_phred)

    anno_bechmark_df <- anno_bechmark_df[anno_bechmark_df$VEST3_score != ".", ]
    anno_bechmark_df$VEST3_score <- as.numeric(anno_bechmark_df$VEST3_score)

    anno_bechmark_df <- anno_bechmark_df[anno_bechmark_df$SiPhy_29way_logOdds != ".", ]
    anno_bechmark_df$SiPhy_29way_logOdds <- as.numeric(anno_bechmark_df$SiPhy_29way_logOdds)

    anno_bechmark_df <- anno_bechmark_df[anno_bechmark_df$DANN_score != ".", ]
    anno_bechmark_df$DANN_score <- as.numeric(anno_bechmark_df$DANN_score)

    # Predict metapredictor probabilities
    pred_df <- stats::predict(metapredictor_model, newdata = anno_bechmark_df, type = "prob")
    anno_bechmark_df$metaprediction_score <- pred_df$non.neutral

    metapred_scores_df <- anno_bechmark_df[, c("Gene.refGene", "metaprediction_score")]
    colnames(metapred_scores_df) <- c("gene_symbol", "metaprediction_score")

    # keep first symbol if multiple symbols exist
    metapred_scores_df$gene_symbol[grepl(";", metapred_scores_df$gene_symbol)] <- sapply(metapred_scores_df$gene_symbol[grepl(";", metapred_scores_df$gene_symbol)], function(x) unlist(strsplit(x, ";"))[1])

    # keep highest score
    metapred_scores_df <- metapred_scores_df[order(metapred_scores_df$metaprediction_score, decreasing = TRUE), ]
    metapred_scores_df <- metapred_scores_df[!duplicated(metapred_scores_df$gene_symbol), ]
    return(metapred_scores_df)
}
