#' driveR: A package for Ranking Genes According to Cancer Driverness using Individual Sequencing Results
#'
#' Cancer genomes contain large numbers of somatic alterations but few genes
#' drive tumor development. Identifying molecular cancer driver genes is critical
#' for precision oncology. Most of current approaches either identify driver
#' genes based on mutational recurrence or using estimated scores predicting
#' the functional consequences of mutations.
#'
#' driveR is a tool for scoring and ranking genes from individual or batch somatic whole
#' exome/genome sequencing data according to cancer driverness. As features,
#' driveR uses coding impact metaprediction scores, non-coding impact scores,
#' somatic copy number alteration scores, hotspot gene/double-hit gene condition,
#' ‘phenolyzer’ gene scores and memberships to cancer-related KEGG pathways. It
#' uses these features to calculate cancer-type-specific driverness probabilities
#' for each gene using the related task of a multi-task learning model.
#'
#'
#' @seealso \code{\link{create_features_df}} for creating the features table to
#' be used to calculate cancer driverness probabilities.
#' See \code{\link{calculate_driverness_probs}} for calculating cancer driverness
#' probabilites for each gene.
#'
#' @docType package
#' @name driveR
NULL
