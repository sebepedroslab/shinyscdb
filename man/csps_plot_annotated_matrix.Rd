% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/functions.R
\name{csps_plot_annotated_matrix}
\alias{csps_plot_annotated_matrix}
\title{Plot annotated matrix
(https://github.com/sebepedroslab/metacell-downstream-functions/blob/98f1d7fdd36d982d28e8773f10bfe13e342f628c/Cross_species_functions.R#L1531)}
\usage{
csps_plot_annotated_matrix(
  mat,
  name = "data",
  heatmap_colors = c("white", "orange", "orangered2", "#520c52"),
  min_val = 0,
  max_val = quantile(mat, 0.99),
  use_raster = TRUE,
  row_title = NULL,
  col_title = NULL,
  row_labels = NULL,
  col_labels = NULL,
  max_length_labels = 40,
  fontsize = 5,
  row_annot = NULL,
  row_annot_cols = NULL,
  row_annot_legend = FALSE,
  row_cluster = FALSE,
  row_cluster_method = "pearson",
  col_annot = NULL,
  col_annot_cols = NULL,
  col_annot_legend = FALSE,
  col_cluster = FALSE,
  col_cluster_method = "pearson",
  cex_dotplot = 0.02,
  do_dotplot = FALSE
)
}
\arguments{
\item{mat}{any type of data matrix}

\item{name}{name of the type of data in the matrix (default: "data")}

\item{heatmap_colors}{vector of colors to map to the data}

\item{min_val, max_val}{min and max values of the colorscale}

\item{use_raster}{whether to rasterise}

\item{row_title, col_title}{titles for rows and columns}

\item{row_labels, col_labels}{ad-hoc labels for rows and columns (if set to NULL, they are taken from \code{mat} object)}

\item{max_length_labels}{truncate \code{row_labels} and \code{col_labels} to this maximum length, in characters (default 40)}

\item{fontsize}{size of labels (default 5 pts)}

\item{row_annot, col_annot}{either dataframes where the 1st column is a vector of categories for each row/column (same order is assumed) and 2nd is a vector of colors, or simply a vector of categories. Default is NULL, i.e. no annotations.}

\item{row_annot_cols, col_annot_cols}{named vector of colors, where names are categories that match the vector in \code{row_annot}/\code{col_annot} (not necessary if \code{row_annot}/\code{col_annot} are dataframes).}

\item{row_annot_legend, col_annot_legend}{whether to plot row/col annotation legends}

\item{row_cluster, col_cluster}{if TRUE/FALSE, whether to cluster rows/columns with default parameters. Other built-in options are "pearson", "euclidean", or a precomputed \code{hclust} object.}

\item{cex_dotplot}{transformation factor for dot size, if \code{do_dotplot=TRUE}}

\item{do_dotplot}{draw a dot plot instead of a heatmap}
}
\value{
a ComplexHeatmap object
}
\description{
Plot annotated matrix
(https://github.com/sebepedroslab/metacell-downstream-functions/blob/98f1d7fdd36d982d28e8773f10bfe13e342f628c/Cross_species_functions.R#L1531)
}
