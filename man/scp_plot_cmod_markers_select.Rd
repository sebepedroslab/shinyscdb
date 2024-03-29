% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/functions.R
\name{scp_plot_cmod_markers_select}
\alias{scp_plot_cmod_markers_select}
\title{Prepare heatmap of gene expression: select marker genes & create annotations
same as in Downstream_functions.R, just works on matrix and not on mc object}
\usage{
scp_plot_cmod_markers_select(
  mc_fp,
  black_list = c(),
  sub_list_mc = NULL,
  gene_list = NULL,
  order_genes = TRUE,
  gene_annot_file = NULL,
  annot_header = FALSE,
  gene_font_size = 4,
  clust_ord = NULL,
  per_clust_genes = 20,
  gene_min_fold = 2,
  transverality_N = ncol(mc_fp),
  transv_excluded_mc = NULL,
  output_file = NULL
)
}
\arguments{
\item{black_list}{character, blacklisted genes}

\item{gene_list}{character, list of genes to plot, if NULL (default) ...}

\item{order_genes}{logical, whether to cluster genes (default: TRUE)}

\item{gene_annot_file}{charcter, file path to the file containig gene annotations,
it should have three tab separated columns containing gene ID, pfam architecture
and any additional annotation in the last column}

\item{annot_header}{logical, gene annotation file has column names?}

\item{gene_font_size}{numeric, size of the gene names plotted as rownames
(default: 4)}

\item{clust_ord}{character, metacells in the order in which they should be
plotted; if cluster order is not specified (default: NULL), it is
determined by hierarchical clustering}

\item{per_clust_genes}{integer, how many genes per cluster to aim to show in the heatmap
(default: 20)}

\item{gene_min_fold}{numeric, minimum fold change for a gene to be considered for plotting
(default: 2)}

\item{transverality_N}{integer, number of metacells in which a gene can be highly expressed (>1.4)
to be considered for plotting, by default this is the total number of metacells}

\item{transv_excluded_mc}{character, metacells to be excluded in transversality calculation
(default: NULL)}

\item{output_file}{optionally, a path to RDS file to which the function output will be saved}

\item{mc_object}{loaded metacell object (\code{gMCCov} class)}
}
\description{
Prepares heatmap of gene expression fold change for metacells and single cells
(no plotting done, returns list with selected markers and prepared annotations),
}
