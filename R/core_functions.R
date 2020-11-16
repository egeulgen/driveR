#' Create Coding Impact Meta-prediction Score Data Frame
#'
#' @param annovar_csv_path path to 'ANNOVAR' csv output file
#' @param keep_highest_score boolean to indicate whether to keep only the maximal
#' impact score per gene (default = \code{TRUE}). If \code{FALSE}, all scores
#' per each gene are returned
#' @param keep_single_symbol in ANNOVAR outputs, a variant may be annotated as
#' exonic in multiple genes. This boolean argument controls whether or not to
#' keep only the first encountered symbol for a variant (default = \code{TRUE})
#' @param na.string string that was used to indicate when a score is not available
#' during annotation with ANNOVAR (default = ".")
#'
#' @return data frame of meta-prediction scores containing 2 columns: \describe{
#'   \item{gene_symbol}{HGNC gene symbol}
#'   \item{metaprediction_score}{metapredictor impact score}
#' }
#'
#' @export
#'
#' @examples
#' path2annovar_csv <- system.file("extdata/example.hg19_multianno.csv",
#'                                 package = "driveR")
#' metapred_df <- predict_coding_impact(path2annovar_csv)
predict_coding_impact <- function(annovar_csv_path,
                                  keep_highest_score = TRUE,
                                  keep_single_symbol = TRUE,
                                  na.string = ".") {
    # argument checks
    if (!file.exists(annovar_csv_path))
        stop("The file used for `annovar_csv_path` does not exist")

    annovar_df <- utils::read.csv(annovar_csv_path)

    nec_cols <- c("Gene.refGene",
                  "SIFT_score", "Polyphen2_HDIV_score", "LRT_score",
                  "MutationTaster_score", "MutationAssessor_score",
                  "FATHMM_score", "GERP.._RS", "phyloP7way_vertebrate",
                  "CADD_phred", "VEST3_score", "SiPhy_29way_logOdds",
                  "DANN_score")
    if (!all(nec_cols %in% colnames(annovar_df)))
        stop("The table in `annovar_csv_path` should contain all of the following columns: ",
             paste(dQuote(nec_cols), collapse = ", "))

    if (!is.logical(keep_highest_score))
        stop("`keep_highest_score` should be logical")

    if (!is.logical(keep_single_symbol))
        stop("`keep_single_symbol` should be logical")

    # filter out missing scores
    annovar_df <- annovar_df[annovar_df$SIFT_score != na.string, ]
    annovar_df$SIFT_score <- as.numeric(annovar_df$SIFT_score)

    annovar_df <- annovar_df[annovar_df$Polyphen2_HDIV_score != na.string, ]
    annovar_df$Polyphen2_HDIV_score <- as.numeric(annovar_df$Polyphen2_HDIV_score)

    annovar_df <- annovar_df[annovar_df$LRT_score != na.string, ]
    annovar_df$LRT_score <- as.numeric(annovar_df$LRT_score)

    annovar_df <- annovar_df[annovar_df$MutationTaster_score != na.string, ]
    annovar_df$MutationTaster_score <- as.numeric(annovar_df$MutationTaster_score)

    annovar_df <- annovar_df[annovar_df$MutationAssessor_score != na.string, ]
    annovar_df$MutationAssessor_score <- as.numeric(annovar_df$MutationAssessor_score)

    annovar_df <- annovar_df[annovar_df$FATHMM_score != na.string, ]
    annovar_df$FATHMM_score <- as.numeric(annovar_df$FATHMM_score)

    annovar_df <- annovar_df[annovar_df$GERP.._RS != na.string, ]
    annovar_df$GERP.._RS <- as.numeric(annovar_df$GERP.._RS)

    annovar_df <- annovar_df[annovar_df$phyloP7way_vertebrate != na.string, ]
    annovar_df$phyloP7way_vertebrate <- as.numeric(annovar_df$phyloP7way_vertebrate)

    annovar_df <- annovar_df[annovar_df$CADD_phred != na.string, ]
    annovar_df$CADD_phred <- as.numeric(annovar_df$CADD_phred)

    annovar_df <- annovar_df[annovar_df$VEST3_score != na.string, ]
    annovar_df$VEST3_score <- as.numeric(annovar_df$VEST3_score)

    annovar_df <- annovar_df[annovar_df$SiPhy_29way_logOdds != na.string, ]
    annovar_df$SiPhy_29way_logOdds <- as.numeric(annovar_df$SiPhy_29way_logOdds)

    annovar_df <- annovar_df[annovar_df$DANN_score != na.string, ]
    annovar_df$DANN_score <- as.numeric(annovar_df$DANN_score)

    if (nrow(annovar_df) == 0) {
        empty_df <- as.data.frame(matrix("", nrow = 0, ncol = 2))
        colnames(empty_df) <- c("gene_symbol", "metaprediction_score")
        warning("There were no mutations with prediction scores for all tools")
        return(empty_df)
    }


    # Predict metapredictor probabilities
    pred_df <- stats::predict(metapredictor_model, newdata = annovar_df, type = "prob")
    annovar_df$metaprediction_score <- pred_df$non.neutral

    metapred_scores_df <- annovar_df[, c("Gene.refGene", "metaprediction_score")]
    colnames(metapred_scores_df) <- c("gene_symbol", "metaprediction_score")

    # keep first symbol if multiple symbols exist
    if (keep_single_symbol)
        metapred_scores_df$gene_symbol[grepl(";", metapred_scores_df$gene_symbol)] <- vapply(metapred_scores_df$gene_symbol[grepl(";", metapred_scores_df$gene_symbol)], function(x) unlist(strsplit(x, ";"))[1], "gene sym")

    # keep highest score
    metapred_scores_df <- metapred_scores_df[order(metapred_scores_df$metaprediction_score, decreasing = TRUE), ]
    if (keep_highest_score)
        metapred_scores_df <- metapred_scores_df[!duplicated(metapred_scores_df$gene_symbol), ]

    return(metapred_scores_df)
}



#' Create Data Frame of Features for Driver Gene Prioritization
#'
#' @inheritParams predict_coding_impact
#' @inheritParams create_gene_level_scna_df
#' @param phenolyzer_annotated_gene_list_path path to 'phenolyzer'
#' "annotated_gene_list" file
#' @param prep_phenolyzer_input boolean to indicate whether or not to create
#' a vector of genes for use as the input of 'phenolyzer' (default = \code{FALSE}).
#' If \code{TRUE}, the features data frame is not created and instead the vector
#' of gene symbols (union of all genes for which scores are available) is
#' returned.
#' @inheritParams create_SCNA_score_df
#' @inheritParams determine_hotspot_genes
#' @inheritParams determine_double_hit_genes
#' @param verbose boolean controlling verbosity (default = \code{TRUE})
#'
#' @return If \code{prep_phenolyzer_input=FALSE} (default), a data frame of
#' features for prioritizing cancer driver genes (\code{gene_symbol} as
#' the first column and 26 other columns containing features). If
#' \code{prep_phenolyzer_input=TRUE}, the functions returns a vector gene symbols
#' (union of all gene symbols for which scores are available) to be used as the
#' input for performing 'phenolyzer' analysis.
#'
#' The features data frame contains the following columns:
#' \describe{
#'   \item{gene_symbol}{HGNC gene symbol}
#'   \item{metaprediction_score}{the maximum metapredictor (coding) impact score for the gene}
#'   \item{noncoding_score}{the maximum non-coding PHRED-scaled CADD score for the gene}
#'   \item{scna_score}{SCNA proxy score. SCNA density (SCNA/Mb) of the minimal common region (MCR) in which the gene is located}
#'   \item{hotspot_double_hit}{boolean indicating whether the gene is a hotspot gene (indication of oncogenes) or subject to double-hit (indication of tumor-suppressor genes)}
#'   \item{phenolyzer_score}{'phenolyzer' score for the gene}
#'   \item{hsa03320}{boolean indicating whether or not the gene takes part in this KEGG pathway}
#'   \item{hsa04010}{boolean indicating whether or not the gene takes part in this KEGG pathway}
#'   \item{hsa04020}{boolean indicating whether or not the gene takes part in this KEGG pathway}
#'   \item{hsa04024}{boolean indicating whether or not the gene takes part in this KEGG pathway}
#'   \item{hsa04060}{boolean indicating whether or not the gene takes part in this KEGG pathway}
#'   \item{hsa04066}{boolean indicating whether or not the gene takes part in this KEGG pathway}
#'   \item{hsa04110}{boolean indicating whether or not the gene takes part in this KEGG pathway}
#'   \item{hsa04115}{boolean indicating whether or not the gene takes part in this KEGG pathway}
#'   \item{hsa04150}{boolean indicating whether or not the gene takes part in this KEGG pathway}
#'   \item{hsa04151}{boolean indicating whether or not the gene takes part in this KEGG pathway}
#'   \item{hsa04210}{boolean indicating whether or not the gene takes part in this KEGG pathway}
#'   \item{hsa04310}{boolean indicating whether or not the gene takes part in this KEGG pathway}
#'   \item{hsa04330}{boolean indicating whether or not the gene takes part in this KEGG pathway}
#'   \item{hsa04340}{boolean indicating whether or not the gene takes part in this KEGG pathway}
#'   \item{hsa04350}{boolean indicating whether or not the gene takes part in this KEGG pathway}
#'   \item{hsa04370}{boolean indicating whether or not the gene takes part in this KEGG pathway}
#'   \item{hsa04510}{boolean indicating whether or not the gene takes part in this KEGG pathway}
#'   \item{hsa04512}{boolean indicating whether or not the gene takes part in this KEGG pathway}
#'   \item{hsa04520}{boolean indicating whether or not the gene takes part in this KEGG pathway}
#'   \item{hsa04630}{boolean indicating whether or not the gene takes part in this KEGG pathway}
#'   \item{hsa04915}{boolean indicating whether or not the gene takes part in this KEGG pathway}
#' }
#'
#' @export
#'
#' @examples
#' \donttest{
#' path2annovar_csv <- system.file("extdata/example.hg19_multianno.csv",
#'                                 package = "driveR")
#' path2phenolyzer_out <- system.file("extdata/example.annotated_gene_list",
#'                                    package = "driveR")
#' features_df <- create_features_df(annovar_csv_path = path2annovar_csv,
#'                                   scna_df = example_scna_table,
#'                                   phenolyzer_annotated_gene_list_path = path2phenolyzer_out)
#' }
#'
#' @seealso \code{\link{prioritize_driver_genes}} for prioritizing cancer driver genes
create_features_df <- function(annovar_csv_path,
                               scna_df,
                               phenolyzer_annotated_gene_list_path,
                               batch_analysis = FALSE,
                               prep_phenolyzer_input = FALSE,
                               log2_ratio_threshold = 0.25,
                               gene_overlap_threshold = 25,
                               MCR_overlap_threshold = 25,
                               hotspot_threshold = 5L,
                               log2_hom_loss_threshold = -1,
                               verbose = TRUE,
                               na.string = ".") {
    ### argument check
    if (!is.logical(prep_phenolyzer_input))
        stop("`prep_phenolyzer_input` should be logical")

    ### determine individual features
    # coding variant impact metaprediction scores
    if (verbose)
        message("Predicting impact of coding variants")
    metaprediction_scores_df <- predict_coding_impact(annovar_csv_path = annovar_csv_path,
                                                      na.string = na.string)

    # non-coding variant impact metaprediction scores
    if (verbose)
        message("Predicting impact of non-coding variants")
    noncoding_scores_df <- create_noncoding_impact_score_df(annovar_csv_path = annovar_csv_path,
                                                            na.string = na.string)

    # gene-level SCNA df
    if (verbose)
        message("Determining gene-level SCNAs (This may take a while)")
    gene_SCNA_df <- create_gene_level_scna_df(scna_df = scna_df,
                                              gene_overlap_threshold = gene_overlap_threshold)

    # SCNA scores
    if (verbose)
        message("Scoring SCNA events")
    scna_scores_df <- create_SCNA_score_df(gene_SCNA_df = gene_SCNA_df,
                                           log2_ratio_threshold = log2_ratio_threshold,
                                           MCR_overlap_threshold = MCR_overlap_threshold)

    # hotspot or double-hit genes
    if (verbose)
        message("Determining hotspot/double-hit genes")
    hotspot_genes <- determine_hotspot_genes(annovar_csv_path = annovar_csv_path,
                                             hotspot_threshold = hotspot_threshold)
    double_hit_genes <- determine_double_hit_genes(annovar_csv_path = annovar_csv_path,
                                                   gene_SCNA_df = gene_SCNA_df,
                                                   log2_hom_loss_threshold = log2_hom_loss_threshold,
                                                   batch_analysis = batch_analysis)
    hotspot_dhit_genes <- unique(c(hotspot_genes, double_hit_genes))

    # return `all_genes` if phenolyzer input is required
    all_genes <- unique(c(metaprediction_scores_df$gene_symbol,
                          noncoding_scores_df$gene_symbol,
                          scna_scores_df$gene_symbol,
                          hotspot_dhit_genes))
    if (length(all_genes) == 0) {
        warning("No genes with suitable alterations identified")
        return(all_genes)
    }
    if (prep_phenolyzer_input)
        return(all_genes)

    # phenolyzer scores
    if (verbose)
        message("Parsing 'phenolyzer' gene scores")
    phenolyzer_df <- utils::read.delim(phenolyzer_annotated_gene_list_path)
    phenolyzer_df <- phenolyzer_df[phenolyzer_df$Gene %in% all_genes, ]

    ### combine all into features data frame
    # create df
    all_genes <- unique(c(metaprediction_scores_df$gene_symbol,
                          noncoding_scores_df$gene_symbol,
                          scna_scores_df$gene_symbol,
                          hotspot_dhit_genes,
                          phenolyzer_df$Gene))
    features_df <- data.frame(gene_symbol = all_genes,
                              metaprediction_score = 0,
                              noncoding_score = 0,
                              scna_score = 0,
                              hotspot_double_hit = all_genes %in% hotspot_dhit_genes,
                              phenolyzer_score = 0)

    # populate with scores
    features_df$metaprediction_score <- metaprediction_scores_df$metaprediction_score[match(features_df$gene_symbol,
                                                                                            metaprediction_scores_df$gene_symbol)]
    features_df$noncoding_score <- noncoding_scores_df$CADD_score[match(features_df$gene_symbol,
                                                                        noncoding_scores_df$gene_symbol)]
    features_df$scna_score <- scna_scores_df$SCNA_density[match(features_df$gene_symbol,
                                                                scna_scores_df$gene_symbol)]
    features_df$phenolyzer_score <- phenolyzer_df$Score[match(features_df$gene_symbol,
                                                              phenolyzer_df$Gene)]
    # replace NAs w/ 0s
    features_df$metaprediction_score <- ifelse(is.na(features_df$metaprediction_score), 0, features_df$metaprediction_score)
    features_df$noncoding_score <- ifelse(is.na(features_df$noncoding_score), 0, features_df$noncoding_score)
    features_df$scna_score <- ifelse(is.na(features_df$scna_score), 0, features_df$scna_score)
    features_df$phenolyzer_score <- ifelse(is.na(features_df$phenolyzer_score), 0, features_df$phenolyzer_score)

    ### add KEGG cancer pathway memberships as features
    if (verbose)
        message("Assessing memberships to KEGG - cancer-related pathways")
    tmp <- lapply(KEGG_cancer_pathways, function(x) features_df$gene_symbol %in% x)
    tmp <- as.data.frame(tmp)
    features_df <- cbind(features_df, tmp)

    return(features_df)
}

#' Prioritize Cancer Driver Genes
#'
#' @param features_df the features data frame for all genes, containing the following columns:
#' \describe{
#'   \item{gene_symbol}{HGNC gene symbol}
#'   \item{metaprediction_score}{the maximum metapredictor (coding) impact score for the gene}
#'   \item{noncoding_score}{the maximum non-coding PHRED-scaled CADD score for the gene}
#'   \item{scna_score}{SCNA proxy score. SCNA density (SCNA/Mb) of the minimal common region (MCR) in which the gene is located}
#'   \item{hotspot_double_hit}{boolean indicating whether the gene is a hotspot gene (indication of oncogenes) or subject to double-hit (indication of tumor-suppressor genes)}
#'   \item{phenolyzer_score}{'phenolyzer' score for the gene}
#'   \item{hsa03320}{boolean indicating whether or not the gene takes part in this KEGG pathway}
#'   \item{hsa04010}{boolean indicating whether or not the gene takes part in this KEGG pathway}
#'   \item{hsa04020}{boolean indicating whether or not the gene takes part in this KEGG pathway}
#'   \item{hsa04024}{boolean indicating whether or not the gene takes part in this KEGG pathway}
#'   \item{hsa04060}{boolean indicating whether or not the gene takes part in this KEGG pathway}
#'   \item{hsa04066}{boolean indicating whether or not the gene takes part in this KEGG pathway}
#'   \item{hsa04110}{boolean indicating whether or not the gene takes part in this KEGG pathway}
#'   \item{hsa04115}{boolean indicating whether or not the gene takes part in this KEGG pathway}
#'   \item{hsa04150}{boolean indicating whether or not the gene takes part in this KEGG pathway}
#'   \item{hsa04151}{boolean indicating whether or not the gene takes part in this KEGG pathway}
#'   \item{hsa04210}{boolean indicating whether or not the gene takes part in this KEGG pathway}
#'   \item{hsa04310}{boolean indicating whether or not the gene takes part in this KEGG pathway}
#'   \item{hsa04330}{boolean indicating whether or not the gene takes part in this KEGG pathway}
#'   \item{hsa04340}{boolean indicating whether or not the gene takes part in this KEGG pathway}
#'   \item{hsa04350}{boolean indicating whether or not the gene takes part in this KEGG pathway}
#'   \item{hsa04370}{boolean indicating whether or not the gene takes part in this KEGG pathway}
#'   \item{hsa04510}{boolean indicating whether or not the gene takes part in this KEGG pathway}
#'   \item{hsa04512}{boolean indicating whether or not the gene takes part in this KEGG pathway}
#'   \item{hsa04520}{boolean indicating whether or not the gene takes part in this KEGG pathway}
#'   \item{hsa04630}{boolean indicating whether or not the gene takes part in this KEGG pathway}
#'   \item{hsa04915}{boolean indicating whether or not the gene takes part in this KEGG pathway}
#' }
#' @param cancer_type short name of the cancer type. All available cancer types
#' are listed in \code{\link{MTL_submodel_descriptions}}
#'
#' @return data frame with 3 columns:
#' \describe{
#'   \item{gene_symbol}{HGNC gene symbol}
#'   \item{driverness_prob}{estimated probability for each gene in \code{features_df} of being a
#'    cancer driver. The probabilities are calculated using the selected (via
#'    \code{cancer_type}) cancer type's sub-model.}
#'   \item{prediction}{prediction based on the cancer-type-specific threshold (either "driver" or "non-driver")}
#' }
#'
#'
#' @export
#'
#' @examples
#' drivers_df <- prioritize_driver_genes(example_features_table, "LUAD")
#'
#' @seealso \code{\link{create_features_df}} for creating the features table.
#' \code{\link{TCGA_MTL_fit}} for details on the MTL model used for prediction.
prioritize_driver_genes <- function(features_df, cancer_type) {
    # argument checks
    if (!is.data.frame(features_df))
        stop("`features_df` should be a data frame")

    if (ncol(features_df) != 27)
        stop("`features_df` should contain exactly 27 columns")
    req_names <- c("gene_symbol", "metaprediction_score", "noncoding_score",
                   "scna_score", "hotspot_double_hit", "phenolyzer_score",
                   "hsa03320", "hsa04010", "hsa04020", "hsa04024", "hsa04060",
                   "hsa04066", "hsa04110", "hsa04115", "hsa04150", "hsa04151",
                   "hsa04210", "hsa04310", "hsa04330", "hsa04340", "hsa04350",
                   "hsa04370", "hsa04510", "hsa04512", "hsa04520", "hsa04630",
                   "hsa04915")
    if (!all(colnames(features_df) == req_names))
        stop("`features_df` should contain the following columns: ",
             paste(dQuote(req_names), collapse = ", "))

    if (!cancer_type %in% driveR::MTL_submodel_descriptions$short_name)
        stop("`cancer_type` should be one of the short names in `MTL_submodel_descriptions`")

    # determine which model to use
    idx <- which(driveR::MTL_submodel_descriptions$short_name == cancer_type)

    # calculate probabilities
    newX <- as.matrix(features_df[, -1])
    score <- newX %*% TCGA_MTL_fit$W[, idx] + TCGA_MTL_fit$C[idx]
    yhat <- exp(score)
    yhat <- yhat / (1 + yhat)

    # turn into data frame
    prob_df <- data.frame(gene_symbol = features_df$gene_symbol,
                          driverness_prob = yhat[, 1])
    prob_df <- prob_df[order(prob_df$driverness_prob, decreasing = TRUE), ]

    threshold <- driveR::specific_thresholds[cancer_type]
    prob_df$prediction <- ifelse(prob_df$driverness_prob >= threshold, "driver", "non-driver")

    return(prob_df)
}
