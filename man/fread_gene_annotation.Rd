% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/functions.R
\name{fread_gene_annotation}
\alias{fread_gene_annotation}
\title{Import gene annotations}
\usage{
fread_gene_annotation(
  file,
  select = 1:3,
  col.names = c("gene_id", "best human BLAST hit", "PFAM domain"),
  search.column = NULL
)
}
\arguments{
\item{file}{path to file to load}

\item{select}{vector of column names or numbers to keep, drop the rest (see select in
\code{?data.table::fread()})}

\item{col.names}{(see col.names in \code{?data.table::fread()})}
}
\value{
data.table with the selected columns and additional column containing
concatenated values from selected columns (used for table look-up)
}
\description{
Import gene annotations
}
