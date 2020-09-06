#' Random Forest Model for Coding Impact Metaprediction
#'
#' A Random Forest model object for metaprediction of coding variants' impact,
#' using 12 impact scores from different coding impact predictors. The model was
#' trained on 711 coding variants, with 10-folds repeated 3 times cross-validation.
#'
#' @format model object
"metapredictor_model"

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

#' Somatic Copy Number Alteration Table for Imielinski et al.
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
"imielinski_scna_table"
