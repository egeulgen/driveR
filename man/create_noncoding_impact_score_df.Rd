% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utility.R
\name{create_noncoding_impact_score_df}
\alias{create_noncoding_impact_score_df}
\title{Create Non-coding Impact Score Data Frame}
\usage{
create_noncoding_impact_score_df(annovar_csv_path, na.string = ".")
}
\arguments{
\item{annovar_csv_path}{path to 'ANNOVAR' csv output file}

\item{na.string}{string that was used to indicate when a score is not available
during annotation with ANNOVAR (default = '.')}
}
\value{
data frame of meta-prediction scores containing 2 columns: \describe{
\item{gene_symbol}{HGNC gene symbol}
\item{CADD_phred}{PHRED-scaled CADD score}
}
}
\description{
Create Non-coding Impact Score Data Frame
}
