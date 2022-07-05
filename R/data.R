#' Tumor type specific probability thresholds
#'
#' Driver gene probability thresholds for all 21 cancer types (submodels).
#'
#' @format vector with 21 elements
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
"MTL_submodel_descriptions"

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
