#' Random Forest Model for Coding Impact Metaprediction
#'
#' A Random Forest model object for metaprediction of coding variants' impact,
#' using 12 impact scores from different coding impact predictors. The model was
#' trained on 711 coding variants, with 10-folds repeated 3 times cross-validation.
#'
#' @format model object
"metapredictor_model"

#' Multi-Task Learning Model for Predicting Cancer Driver Genes
#'
#' A Multi-Task Learning (MTL) classification model object for determining
#' cancer driver genes based on 26 features. The model was trained using
#' TCGA data (obtained from ICGC release 28) with lasso regularization. It contains
#' 21 sub-models for different cancer types.
#'
#' @format MTL model object
#' @seealso \code{\link{MTL_submodel_descriptions}} for short names and descriptions
#' of all sub-models.
"TCGA_MTL_fit"

#' Tumor type specific probability thresholds
#'
#' Driver gene probability thresholds for all 21 cancer types (submodels).
#'
#' @format vector with 21 elements
#' @seealso \code{\link{TCGA_MTL_fit}} for the Multi-Task Learning model.
"specific_thresholds"

#' MTL Sub-model Descriptions
#'
#' A data frame containing descriptions for all sub-models of the MTL model.
#'
#' @format A data frame with 21 rows and 2 variables:
#' \describe{
#'   \item{short_name}{short name for the cancer type}
#'   \item{description}{description of the cancer type}
#' }
#' @seealso \code{\link{TCGA_MTL_fit}} for the MTL model.
"MTL_submodel_descriptions"

#' Table of Pan-Cancer Minimal Common Regions
#'
#' A data set containing the minimal common regions (MCRs) across all cancer
#' types studied in Kim TM, Xi R, Luquette LJ, Park RW, Johnson MD, Park PJ.
#' Functional genomic analysis of chromosomal aberrations in a compendium of
#' 8000 cancer genomes. Genome Res. 2013;23(2):217-27.
#'
#' @format A data frame with 165 rows and 5 variables:
#' \describe{
#'   \item{chr}{chromosome the MCR is located in}
#'   \item{start}{start position of the MCR}
#'   \item{end}{end position of the MCR}
#'   \item{MCR_type}{the type ("Amp" or "Del") of the MCR peak}
#'   \item{SCNA_density}{SCNA per Mb within the MCR}
#' }
#' @source \url{https://pubmed.ncbi.nlm.nih.gov/23132910/}
"MCR_table"

#' KEGG "Pathways in cancer"-related Pathways - Gene Sets
#'
#' A list containing the genes involved in each Homo sapiens KEGG "Pathways in
#' cancer" (hsa05200)-related Pathways. Each element is a vector of gene symbols
#' located in the given pathway. Names correspond to the KEGG ID of the pathway.
#' \emph{Generated on Nov 17, 2020.}
#'
#' @format list containing 21 vectors of gene symbols. Each vector corresponds
#'   to a pathway.
#' @seealso \code{\link{KEGG_cancer_pathways_descriptions}} for descriptions of
#'  KEGG "Pathways in cancer"-related pathways.
"KEGG_cancer_pathways"

#' KEGG "Pathways in cancer"-related Pathways - Descriptions
#'
#' A data frame containing descriptions for KEGG "Pathways in cancer"
#' (hsa05200)-related pathways.
#' \emph{Generated on Nov 17, 2020.}
#'
#' @format A data frame with 21 rows and 2 variables:
#' \describe{
#'   \item{id}{KEGG pathway ID}
#'   \item{description}{KEGG pathway description}
#' }
"KEGG_cancer_pathways_descriptions"

#' Example Somatic Copy Number Alteration Table
#'
#' A data set containing the somatic copy number alteration data for the lung
#' adenocarcinoma patient studied in Imielinski M, Greulich H, Kaplan B, et al.
#' Oncogenic and sorafenib-sensitive ARAF mutations in lung adenocarcinoma.
#' J Clin Invest. 2014;124(4):1582-6.
#'
#' @format A data frame with 3160 rows and 4 variables:
#' \describe{
#'   \item{chr}{chromosome the segment is located in}
#'   \item{start}{start position of the segment}
#'   \item{end}{end position of the segment}
#'   \item{log2ratio}{\ifelse{html}{\out{log<sub>2</sub>}}{\eqn{log_2}} ratio of
#'   the segment}
#' }
#' @source \url{https://pubmed.ncbi.nlm.nih.gov/24569458/}
"example_scna_table"

#' Example Features Table for Driver Prioritization
#'
#' The example dataset containing features for prioritizing cancer driver genes for
#' the lung adenocarcinoma patient studied in Imielinski M, Greulich H, Kaplan B, et al.
#' Oncogenic and sorafenib-sensitive ARAF mutations in lung adenocarcinoma.
#' J Clin Invest. 2014;124(4):1582-6.
#'
#' @format A data frame with 4901 rows and 27 variables:
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
#' @seealso \code{\link{KEGG_cancer_pathways_descriptions}} for descriptions of
#'  KEGG "Pathways in cancer"-related pathways.
"example_features_table"

#' Example Cohort-level Somatic Copy Number Alteration Table
#'
#' A data set containing the somatic copy number alteration data for 10 randomly
#' selected samples from TCGA's LAML (Acute Myeloid Leukemia) cohort
#'
#' @format A data frame with 126147 rows and 5 variables:
#' \describe{
#'   \item{chr}{chromosome the segment is located in}
#'   \item{start}{start position of the segment}
#'   \item{end}{end position of the segment}
#'   \item{log2ratio}{\ifelse{html}{\out{log<sub>2</sub>}}{\eqn{log_2}} ratio of
#'   the segment}
#'   \item{tumor_id}{ID for the tumor containing the SCNA segment}
#' }
#' @source \url{https://dcc.icgc.org/releases/release_28}
"example_cohort_scna_table"

#' Example Cohort-level Features Table for Driver Prioritization
#'
#' The example dataset containing features for prioritizing cancer driver genes for 10 randomly
#' selected samples from TCGA's LAML (Acute Myeloid Leukemia) cohort
#'
#' @format A data frame with 349 rows and 27 variables:
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
#' @seealso \code{\link{KEGG_cancer_pathways_descriptions}} for descriptions of
#'  KEGG "Pathways in cancer"-related pathways.
"example_cohort_features_table"
