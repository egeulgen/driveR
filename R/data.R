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
#'   \item{log2ratio}{\ifelse{html}{\out{log<sub>2</sub>}}{\eqn{log_2}} ratio of the segment}
#' }
#' @source \url{https://pubmed.ncbi.nlm.nih.gov/24569458/}
"imielinski_scna_table"