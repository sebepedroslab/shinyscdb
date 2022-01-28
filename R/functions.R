# Exported functions -----------------------------------------------------------

#' Import gene annotations
#'
#' @param file path to file to load
#' @param select vector of column names or numbers to keep, drop the rest (see select in
#'  `?data.table::fread()`)
#' @param col.names (see col.names in `?data.table::fread()`)
#'
#' @return data.table with the selected columns and additional column containing
#'  concatenated values from selected columns (used for table look-up)
#' @importFrom data.table :=
#' @export
fread_gene_annotation <- function(
  file,
  select = 1:3,
  col.names = c("gene_id","PFAM domain","best human BLAST hit"),
  search.column = NULL
) {
  GENE_ANNT <- data.table::fread(
    file = file,
    header = TRUE,
    sep = "\t",
    fill = TRUE,
    select = select,
    col.names = col.names
  )
  if (!is.null(search.column)) {
    GENE_ANNT[,search_id:=do.call(paste,.SD),.SDcols=search.column]
  }
  return(GENE_ANNT)
}

#' Import metacell annotations
#'
#' @param file path to file to load
#'
#' @return daata.table with the following columns: mc, cell_type, color
#'
#' @export
fread_cell_annotation <- function(file, select=NULL) {

  CELL_ANNT <- data.table::fread(file = file, sep = "\t", fill = TRUE, select=select)
  setnames(CELL_ANNT, colnames(CELL_ANNT)[1:3], c("mc","cell_type","color"))

  return(CELL_ANNT)
}

# Helper functions -------------------------------------------------------------

#' Reduce vector of metacells
red_mc_vector <- function(x,range_sep=":") {
  all_mcs <- sort(as.integer(x))
  all_ir <- IRanges::IRanges(start = all_mcs, end = all_mcs)
  red_ir <- IRanges::reduce(all_ir)
  starts <- IRanges::start(red_ir)
  ends <- IRanges::end(red_ir)

  endsout <- ends
  for(i in 1:length(starts)) {
    if (starts[i] == ends[i])
      endsout[[i]] <- ""
  }
  outb <- paste(starts,endsout,sep=range_sep)
  outv <- str_remove(outb,paste0(range_sep,"$"))
  paste(outv,collapse=",")
}

#' Summarize cell annotation
summarize_cell_annotation <- function(annt) {

  setDT(annt)
  setorder(annt,"mc")
  tanns <- tapply(annt$mc, annt$cell_type, red_mc_vector)

  dt <- data.table(
    'cell type' = unique(annt$cell_type),
    metacells = tanns[unique(annt$cell_type)],
    cols = annt$color[match(unique(annt$cell_type),annt$cell_type)]
  )
  dt[, colshex := col2hex(cols)]
  dt[, metacells := cell_spec(
      metacells, "html", color = "white", align = "c",
      background = factor(`cell type`, dt$`cell type`, dt$colshex)
  )]
  dtt <- dt[,c("cell type","metacells")]
  knitr::kable(dtt, escape = FALSE, align = "c")
}

# Plotting functions -----------------------------------------------------------

#' Plot an expression barplots of a gene
#' @param nmat mc_fp matrix
#' @param umat sc UMI matrix
sg_plot  <- function(
  nmat, umat, cttable, gid=NULL, sid=NULL,
  mdnorm, annt, ctpalette=NULL, mc_label_size=9
){
  # selected gene
  if (all(is.null(gid),is.null(sid)))
    stop("Need to specify either gid or sid!")
  if (is.null(gid))
    gid <- as.character(annt[search_id==sid,gene_id])
  # selected gene umi count
  gxp <- as.numeric(umat[gid,])*10
  if(mdnorm==TRUE)
    gxp <- log2(gxp/median(gxp))
  # selected gene lfc
  gxl <- as.numeric(nmat[gid,])

  gdata <- data.table(
    metacells=as.numeric(colnames(umat)),
    cell_type=cttable$cell_type,
    lfp=gxp, lfc=gxl
  )

  # cell type colours
  if (is.null(ctpalette)) {
    ctpalette <- unique(cttable$color)
    names(ctpalette) <- unique(cttable$cell_type)
  }

  # gene names
  gene_name <- as.character(annt[annt[[1]]==gid,1])
  gene_domain <- as.character(annt[annt[[1]]==gid,2])
  gene_hsap <- as.character(annt[annt[[1]]==gid,3])
  bptitle <- paste(gene_name, gene_domain)
  if (gene_name!=gene_hsap)
    bptitle <- paste(bptitle, gene_hsap)

  # umi frac barplot
  max_label <- sprintf(" %.2f", max(gxp))
  pst_x_max <- which( gxp == max(gxp) ) ; pst_x_max <- pst_x_max[1]

  gp_umi_frac <- ggplot2::ggplot(
      data=gdata,
      aes(x=factor(metacells), y=lfp, fill=factor(cell_type,levels=names(ctpalette)))
    ) +
    ggplot2::geom_bar(stat="identity", colour="grey") +
    ggplot2::geom_blank(aes(y=1.2*lfp)) +
    ggplot2::scale_fill_manual(values=ctpalette) +
    ggplot2::labs(
      x="metacells", y="UMI/10k",
      title=bptitle, fill="cell type"
    ) +
    ggplot2::scale_y_continuous(expand=c(0,0)) +
    ggplot2::geom_hline(yintercept=max(gxp), lty=2) +
    ggplot2::annotate("text", x=pst_x_max, y=1.08*max(gxp), label=max_label, size=5) +
    ggplot2::theme(
      legend.position="none",
      panel.background=element_blank(),
      axis.line=element_line(colour="black"),
      axis.text.x=element_text(angle=90, vjust=0.5, hjust=1, size=mc_label_size),
      text = element_text(size=18)
    )

  # log2fc barplot
  max_label <- sprintf(" %.2f", max(gxl))
  pst_x_max <- which( gxl == max(gxl) ) ; pst_x_max <- pst_x_max[1]

  gp_log_frac <- ggplot2::ggplot(
      data=gdata,
      aes(x=factor(metacells), y=lfc, fill=factor(cell_type,levels=names(ctpalette)))
    ) +
    ggplot2::geom_bar(stat="identity", colour="grey") +
    ggplot2::geom_blank(aes(y=1.2*lfc)) +
    ggplot2::scale_fill_manual(values=ctpalette) +
    ggplot2::labs(
      x="metacells", y="FC", fill="cell type",
      caption="Dashed line indicates the maximum observed value."
    ) +
    ggplot2::geom_hline(yintercept=max(gxl), lty=2) +
    ggplot2::scale_y_continuous(expand=c(0,0)) +
    ggplot2::annotate("text", x=pst_x_max, y=1.08*max(gxl), label=max_label, size=5) +
    ggplot2::theme(
      legend.position="bottom",
      panel.background=element_blank(),
      axis.line=element_line(colour="black"),
      axis.text.x=element_text(angle=90, vjust=0.5, hjust=1, size=mc_label_size),
      text = element_text(size=18)
    )
  ggp <- egg::ggarrange(gp_umi_frac, gp_log_frac, nrow=2, ncol=1)
  return(ggp)
}

#' Plot 2d expression of a gene
#' @param nmat mc_fp matrix
#' @param umat sc UMI matrix
#' @param mc2d mc2d object
scp_plot_sc_2d_gene_exp = function(
  mc2d,
  nmat,
  umat,
  sid, annt, gene_id=NULL,
  sc_vector = NULL,
  sc_scale = "unidir",
  sc_transform = NULL,
  sc_min = 0,
  sc_max = 8,
  sc_label = "cell UMI",
  sc_zero_color = "gray95",
  do_umifrac_sc = FALSE,
  mc_vector = NULL,
  mc_scale = "unidir",
  mc_transform = log2,
  mc_min = 0,
  mc_max = 2,
  mc_label = "log2(FC)",
  mc_zero_color = NULL,
  unidir_color_scale = c("gray90","orange","orangered2","#520c52"),
  bidir_color_scale = c("midnightblue","dodgerblue3","deepskyblue","#b8e0ed","gray90","#eccac0","#ff8d36","#f34312","#8e0631"),
  plot_edges=TRUE,
  plot_mcs=TRUE,
  plot_mc_name=TRUE,
  width=12,
  height=12,
  res=NA,
  cex_mc=3,
  cex_sc=0.75,
  do_axes=FALSE) {

  # get expression vectors (if not already provided)
  gene_id <- as.character(annt[search_id==sid,gene_id])

  if (is.null(sc_vector)) {
    gid <- match(gene_id,rownames(umat))
    sc_vector = umat[gid,]
  }
  if (is.null(mc_vector)) {
    gid <- match(gene_id,rownames(nmat))
    mc_vector = nmat[gid,]
  }

  # apply transformations
  if (do_umifrac_sc) {
    sc_vector = sc_vector / Matrix::colSums(umat) * 10000
    sc_label = "UMI/10k"
  }
  if (!is.null(sc_transform)) {
    sc_vector = sc_transform(sc_vector)
  }
  if (!is.null(mc_transform)) {
    mc_vector = mc_transform(mc_vector)
  }

  # if sc_max or mc_max values are NULL, get them from a high quantile in the mc_vector
  if (is.null(sc_max)) {
    sc_max = round(quantile(sc_vector[sc_vector>0], 0.75))
  }
  if (is.null(mc_max)) {
    mc_max = round(quantile(mc_vector[mc_vector>0], 0.75))
  }

  # apply min/max to sc vector
  sc_vector [ sc_vector < sc_min ] = sc_min
  sc_vector [ sc_vector > sc_max ] = sc_max
  # apply min/max to mc vector
  mc_vector [ mc_vector < mc_min ] = mc_min
  mc_vector [ mc_vector > mc_max ] = mc_max

  # create color palettes for single cells
  if (sc_scale == "unidir") {
    sc_color_fun = scales::col_numeric(palette=unidir_color_scale, domain=c(sc_min, sc_max))
  } else if (sc_scale == "bidir") {
    sc_color_fun = scales::col_numeric(palette=bidir_color_scale, domain=c(sc_min, sc_max))
  } else {
    message("`sc_scale` should be either `unidir` or `bidir`")
  }
  sc_color = sc_color_fun(sc_vector)
  sc_color_labels = seq(from=sc_min, to=sc_max, length.out=9)
  sc_color_legend = sc_color_fun(sc_color_labels)
  # zero color
  if (!is.null(sc_zero_color)) {
    sc_color [ sc_vector == 0 ] = sc_zero_color
    sc_color_legend [ sc_color_labels == 0 ] = sc_zero_color
  }

  # metacells
  if (mc_scale == "unidir") {
    mc_color_fun = scales::col_numeric(palette=unidir_color_scale, domain=c(mc_min, mc_max))
  } else if (mc_scale == "bidir") {
    mc_color_fun = scales::col_numeric(palette=bidir_color_scale, domain=c(mc_min, mc_max))
  } else {
    message("`sc_scale` should be either `unidir` or `bidir`")
  }
  mc_color = mc_color_fun(mc_vector)
  mc_color_labels = seq(from=mc_min, to=mc_max, length.out=9)
  mc_color_legend = mc_color_fun(mc_color_labels)
  # zero color
  if (!is.null(mc_zero_color)) {
    mc_color [ mc_vector == 0 ] = mc_zero_color
    mc_color_legend [ mc_color_labels == 0 ] = mc_zero_color
  }

  # determine plot max/min
  if (plot_mcs) {

    xlim=c(min(mc2d@sc_x, mc2d@mc_x, na.rm = TRUE), max(mc2d@sc_x, mc2d@mc_x, na.rm = TRUE))
    ylim=c(min(mc2d@sc_y, mc2d@mc_y, na.rm = TRUE), max(mc2d@sc_y, mc2d@mc_y, na.rm = TRUE))

  } else {

    xlim=c(min(mc2d@sc_x, na.rm = TRUE), max(mc2d@sc_x, na.rm = TRUE))
    ylim=c(min(mc2d@sc_y, na.rm = TRUE), max(mc2d@sc_y, na.rm = TRUE))

  }

  # plot individual cells
  plot(
    mc2d@sc_x,
    mc2d@sc_y,
    pch = 19,lwd = 0,
    cex = cex_sc,
    col = alpha(sc_color,0.7),
    xlim = xlim,
    ylim = ylim,
    xlab = NA, ylab = NA, axes=do_axes
  )
  legend("topleft",fill=sc_color_legend, legend=sprintf("%.2f",sc_color_labels), cex=0.8, title=sc_label)

  # plot edges between metacells?
  if (plot_edges) {
    fr = mc2d@graph$mc1
    to = mc2d@graph$mc2
    segments(
      x0=mc2d@mc_x[fr],
      y0=mc2d@mc_y[fr],
      x1=mc2d@mc_x[to],
      y1=mc2d@mc_y[to],
      col="gray70")
  }

  # plot metacells?
  if (plot_mcs) {
    points(
      mc2d@mc_x,
      mc2d@mc_y,
      pch=19,lwd=0.5,
      cex=cex_mc,
      col=alpha(mc_color,0.9))
    legend("topright",fill=mc_color_legend, legend=sprintf("%.2f",mc_color_labels), cex=0.8, title=mc_label)
  }

  # plot metacell ids?
  if (plot_mc_name) {
    text(
      mc2d@mc_x,
      mc2d@mc_y,
      #cex=cex_mc,
      labels=names(mc2d@mc_x))
  }
  title(
    main = gene_id,
    sub = sprintf(
      "2D projection\nn = %i cells | n = %i metacells",
      length(mc2d@sc_x),
      length(mc2d@mc_x)
    )
  )

}

#' Plot 2d proj colored by mc
mc_2d_plot = function(
  mc2d,
  mcsc,
  cttable,
  plot_edges=TRUE,
  plot_mcs=TRUE,
  plot_mc_name=TRUE,
  width=12,
  height=12,
  res=NA,
  cex_mc=3,
  cex_sc=0.75) {

  # get colors of metacells
  mc_colors = structure(cttable$color, names=cttable$mc)
  # get colors of individual cells
  cell_colors=mc_colors[as.character(mcsc[names(mc2d@sc_x)])]

  # determine plot max/min
  if (plot_mcs) {

    xlim=c(min(mc2d@sc_x, mc2d@mc_x, na.rm = TRUE), max(mc2d@sc_x, mc2d@mc_x, na.rm = TRUE))
    ylim=c(min(mc2d@sc_y, mc2d@mc_y, na.rm = TRUE), max(mc2d@sc_y, mc2d@mc_y, na.rm = TRUE))

  } else {

    xlim=c(min(mc2d@sc_x, na.rm = TRUE), max(mc2d@sc_x, na.rm = TRUE))
    ylim=c(min(mc2d@sc_y, na.rm = TRUE), max(mc2d@sc_y, na.rm = TRUE))

  }
  # plot individual cells?
  plot(
    mc2d@sc_x,
    mc2d@sc_y,
    pch=19,lwd=0,
    cex=cex_sc,
    col = alpha(cell_colors,0.4),
    xlim = xlim,
    ylim = ylim,
    xlab = NA, ylab = NA,
    xaxt ="n", yaxt ="n", bty = "n"
  )

  # plot edges between metacells?
  if (plot_edges) {
    fr = mc2d@graph$mc1
    to = mc2d@graph$mc2
    segments(
      x0=mc2d@mc_x[fr],
      y0=mc2d@mc_y[fr],
      x1=mc2d@mc_x[to],
      y1=mc2d@mc_y[to],
      col="gray70")
  }

  # plot metacells?
  if (plot_mcs) {
    points(
      mc2d@mc_x,
      mc2d@mc_y,
      pch=19,lwd=0.5,
      cex=cex_mc,
      col=alpha(mc_colors,0.8))
  }

  # plot metacell ids?
  if (plot_mc_name) {
    text(
      mc2d@mc_x,
      mc2d@mc_y,
      #cex=cex_mc,
      labels=names(mc2d@mc_x))
  }
}

#' Select genes for heatmap
genes_select_dt <- function(sterm, nmat, annt) {
  grep1 <- grep(sterm,annt[[1]],ignore.case=TRUE)
  grep2 <- grep(sterm,annt[[2]],ignore.case=TRUE)
  grep3 <- grep(sterm,annt[[3]],ignore.case=TRUE)
  rids <- sort(unique(c(grep1, grep2, grep3)))
  gs <- annt[rids][[1]]
  gsf <- gs[gs %in% rownames(nmat)]
  gsfid <- match(gsf,annt[[1]])
  annt[gsfid,1:3]
}
genes_upload_dt <- function(gs,annt) {
  gsfid <- match(gs,annt[[1]])
  annt[gsfid,1:3]
}
genes_select_names <- function(dt,rid) {
  dt[rid,][[1]]
}

#' Plot multigene heatmap to show expression of selected genes
mgenes_hmap_height <- function(nmat, gids, annt, min_expression_fc = NULL, max_expression_fc = NULL){
  tryCatch({
    gs <- intersect(gids,rownames(nmat))
    if (!is.null(max_expression_fc)) {
      flt <- apply(nmat[gs,], 1, function(x) !(sort(x,decreasing=TRUE,na.last=TRUE)[1]<min_expression_fc))
      gs <- gs[flt]
    }
    if (!is.null(min_expression_fc)) {
      flt <- apply(nmat[gs,], 1, function(x) !(sort(x,decreasing=TRUE,na.last=TRUE)[1]>max_expression_fc))
      gs <- gs[flt]
    }
    message("Length: ",length(gs))
    return(length(gs))
  }, error=function(e) {
    message(e)
    message("Using fixed length")
    return(2)
  })
}
mgenes_hmap <- function(
  nmat, annt, gids,
  min_expression_fc = NULL,  max_expression_fc = NULL, # gene filtering
  scale_expression_fc = 4, # trim expression values, i.e. set anything > scale_expression_fc to this value
  heatmap_colors = c("white","gray99","orange","orangered2","#520c52"),
  ct_table, cell_type_palette = NULL, cluster_genes = TRUE,
  mcid_font_size = 12, mc_annotaion_height = unit(2, "mm")
){

  # selected genes
  message("Selected genes: ",length(gids))
  gs <- intersect(gids,rownames(nmat))
  rid <- match(gs, annt[[1]])
  gns <- annt[rid][[2]]
  badgns <- gns=="" | is.na(gns)
  gns[badgns] <- gs[badgns]
  row_labels <- gns
  names(row_labels) <- gs

  message("Genes in heatmap: ",nrow(nmat[gs,]))
  hm <- as.matrix(nmat[gs,])

  # filter genes
  if (!is.null(min_expression_fc)) {
    message("Filtering genes by min FC")
    flt <- apply(hm, 1, function(x) !(sort(x,decreasing=TRUE,na.last=TRUE)[1]<min_expression_fc))
    hm <- hm[flt,]
  }
  if (!is.null(max_expression_fc)) {
    message("Filtering genes by max FC")
    flt <- apply(hm, 1, function(x) !(sort(x,decreasing=TRUE,na.last=TRUE)[1]>max_expression_fc))
    hm <- hm[flt,]
  }

  # expression matrix
  hm <- pmin(hm,scale_expression_fc)
  hm[is.na(hm)] <- 0
  hm[is.nan(hm)] <- 0
  hm[is.infinite(hm)] <- 0

  # order genes
  if (cluster_genes == TRUE) {
    message("Clustering genes")
    gord <- order(apply(hm, 1,function(x) which.max(rollmean(x,1))))
    hm <- hm[gord,]
  }

  # heatmap colors
  if (is.null(min_expression_fc)) min_expression_fc=0
  if (is.null(max_expression_fc)) max_expression_fc=5
  max_expression_fc <- pmin(scale_expression_fc, max_expression_fc)
  message("Scaling colors for heatmap between ", min_expression_fc , " and ", max_expression_fc)
  col_fun = circlize::colorRamp2(
    breaks = seq(from = min_expression_fc, to = max_expression_fc, length.out = length(heatmap_colors)),
    colors = heatmap_colors
  )

  # cell type colours
  if (is.null(cell_type_palette)) {
    message("Colors for cell types")
    cell_type_palette <- ct_table$color
    names(cell_type_palette) <- ct_table$cell_type
  }
  cell_types <- names(cell_type_palette)
  cell_colours <- unique(cell_type_palette)
  names(cell_colours) <- unique(names(cell_type_palette))

  # cell type annotiation bar
  message("Cell types annotations")
  ct_ann <- ComplexHeatmap::columnAnnotation(
    ct = cell_types, col = list(ct=cell_colours), height = mc_annotaion_height,
    gp = gpar(fontsize = mcid_font_size), border = FALSE,
    show_annotation_name = FALSE, show_legend = FALSE
  )

  # gene annotation
  gs <- rownames(hm)
  annid <- match(gs,annt[[1]])
  ganns <- annt[annid][[3]]
  ganns[is.na(ganns)] <- ""
  ganns_tr <- unlist(lapply(ganns, function(x){
    if(nchar(x)>30) {
      paste0(substr(x,1,27),"...")
    } else {
      x
    }
  }))
  names(ganns_tr) <- gs
  gs_ann <- ComplexHeatmap::HeatmapAnnotation(
    which = "row", gn = anno_text(ganns_tr,which="row"),
    show_annotation_name = FALSE, show_legend = FALSE
  )

  # expression heatmap
  message("Building heatmap")
  ht1 <- ComplexHeatmap::Heatmap(
    cbind(hm), name = "expression FC", show_heatmap_legend = TRUE,
    cluster_columns = FALSE, cluster_rows = FALSE,
    #row_labels = row_labels[rownames(hm)],
    show_column_dend = FALSE, show_row_dend = FALSE,
    show_column_names = TRUE, show_row_names = TRUE,
    row_names_side = "left", column_names_side = "bottom",
    column_names_gp = gpar(fontsize = mcid_font_size), border = TRUE,
    col = col_fun, rect_gp = gpar(col = "gray88", lwd = 0.1),
    bottom_annotation = ct_ann, top_annotation = ct_ann,
    right_annotation = gs_ann,
    heatmap_legend_param = list(
      legend_direction = "horizontal",
      legend_width = unit(5, "cm"),
      border = TRUE
    )
  )
  draw(ht1, heatmap_legend_side = "bottom")
}


#' Generate table for highly expressed genes in a group of metacells
#' @param method, method to use for retrieving gnees expressed in metacells,
#' either "absolute" or "median"; if "absolute", all metacells (- lky percentage)
#' must have gene fc above threshold; if "median", the median in all metacells
#' must be above threshold
mc_gene_summary <- function(
  mc_ids, fc=2, method=c("absolute","median"), methodbg=method, lky=0,
  usefcbg=FALSE, fcbg=fc, lkybg=lky,
  mc_fp, mc_counts, annt, tfannt
){
  mc_counts <- as.data.frame(mc_counts)
  mcs <- colnames(mc_fp)
  mc_ids <- as.character(mc_ids)
  mc_others <- mcs[!(mcs %in% mc_ids)]
  if (length(method)>1) method=method[1]
  if (!method %in% c("absolute","median")) stop("Method must be either 'absolute' or 'median'!")

  # horizontal UMI fraction
  #tot_umis <- rowSums(mc_counts)
  if (length(mc_ids)==1) {
    umisum_selected <- mc_counts[[mc_ids]]
    names(umisum_selected) <- rownames(mc_counts)
  } else {
    umisum_selected <- rowSums(cbind(mc_counts[,mc_ids]))
  }
  umisum_others <- rowSums(mc_counts[mc_others])
  umi_frac <- umisum_selected/(umisum_selected+umisum_others)
  umi_frac[is.nan(umi_frac)] <- 0

  # median gene FC
  median_fc <- apply(cbind(mc_fp[,mc_ids]),1,median,na.rm=TRUE)
  names(median_fc) <- rownames(mc_fp)

  # "gap genes" - specifically expressed in selected mcs
  if (method == "absolute") {

    ntarget <- length(mc_ids) - length(mc_ids)*lky
    abovefc <- rowSums(cbind(mc_fp[,mc_ids]) > fc)
    gap_genes_abovefc <- abovefc > ntarget | abovefc == ntarget
    if (usefcbg==TRUE) {
      if (methodbg == "absolute") {
        ntargetbg <- length(mc_others) - lkybg*length(mc_others)
        belowfc <- rowSums(cbind(mc_fp[,mc_others]) < fcbg)
        gap_genes_belowfc <- belowfc > ntargetbg | belowfc== ntargetbg
        gap_genes <- rownames(mc_fp)[which(gap_genes_abovefc & gap_genes_belowfc)]
      } else if ( methodbg == "median") {
        gap_genes_belowfc <-  apply(cbind(mc_fp[,mc_others]),1,function(x) median(x,na.rm=TRUE)<fcbg)
        gap_genes <- rownames(mc_fp)[which(gap_genes_abovefc & gap_genes_belowfc)]
      }
    } else {
      fcbg <- NULL
      gap_genes <- rownames(mc_fp)[which(gap_genes_abovefc==TRUE)]
    }
    gap_genes_annot <- intersect(gap_genes,annt$gene_id)
    anntids <- match(gap_genes_annot,annt$gene_id)

  } else if (method == "median") {

    gap_genes_abovefc <- apply(cbind(mc_fp[,mc_ids]),1,function(x) median(x,na.rm=TRUE)>fc)
    if (usefcbg==TRUE) {
      if (methodbg == "absolute") {
        ntargetbg <- length(mc_others) - lkybg*length(mc_others)
        belowfc <- rowSums(cbind(mc_fp[,mc_others]) < fcbg)
        gap_genes_belowfc <- belowfc > ntargetbg | belowfc== ntargetbg
        gap_genes <- rownames(mc_fp)[which(gap_genes_abovefc & gap_genes_belowfc)]
      } else if ( methodbg == "median") {
        gap_genes_belowfc <-  apply(cbind(mc_fp[,mc_others]),1,function(x) median(x,na.rm=TRUE)<fcbg)
        gap_genes <- rownames(mc_fp)[which(gap_genes_abovefc & gap_genes_belowfc)]
      }
    } else {
      fcbg <- NULL
      gap_genes <- rownames(mc_fp)[which(gap_genes_abovefc==TRUE)]
    }
    gap_genes_annot <- intersect(gap_genes,annt$gene_id)
    anntids <- match(gap_genes_annot,annt$gene_id)

  }

  # gene table
  gap_dt <- annt[anntids,1:3]
  gap_dt[,':='(
    tot_umi = umisum_selected[gap_genes_annot],
    umi_frac = round(umi_frac[gap_genes_annot]*100,2),
    median_fc = round(median_fc[gap_genes_annot],2),
    tf = ifelse(gene_id %in% tfannt$gene,"yes","no")
  )]
  gap_dt[,tf:=factor(tf,levels=c("yes","no"))]
  setorder(gap_dt,-median_fc)
  setnames(gap_dt,c("gene_id","tot_umi","umi_frac","median_fc"),c("gene ID","total UMIs","% UMIs","median fc"))

  # summary of the search and results
  search_dt <- data.table(
    `selected metacells` = paste( as.character(mc_ids), collapse=", "),
    `minimum fold change in selected metacells` = fc,
    `maximum fold change in non-selected metacells` = fcbg,
    `number of genes` = length(gap_genes_annot)
  )
  search_tdt <- data.table::transpose(search_dt, keep.names = "V0")
  setnames(search_tdt,c(""," "))
  # output
  list(gene_summary=gap_dt, summary=search_tdt)
}

#' Plot heatmap of gene expression for metacells
#'
scp_plot_cmod_markers_mc <- function(
  marker_data_list,
  output_file = NULL,
  height = 10,
  width = 5,
  res = NA,
  show_heatmap_legend = TRUE,
  highlight_genes = NULL,
  show_gene_names = FALSE,
  gene_font_size = 5,
  clust_col = NULL, clust_bars = NULL, clust_col_others = NULL, col_others = NULL,
  clust_anno_size = unit(1,"mm"),
  show_mc_names = TRUE, mc_font_size = 5,
  heatmap_colors = c("white","gray99","orange","orangered2","#520c52"),
  max_expression_fc = 5,
  gene_chr_limit = 70,
  verbose=FALSE,
  save_rds=TRUE,
  print_border=TRUE,
  show_clust_borders=TRUE
) {

  # PLOT METACELL PROFILE
  message("Plotting metacell expression")

  # get variables necessary for hm definition (from previous function call)
  genes = marker_data_list$genes
  gene_ord = marker_data_list$gene_ord
  clust_ord = marker_data_list$clust_ord
  niche_geomean_n = marker_data_list$niche_geomean_n
  gene_labels_1 = marker_data_list$gene_labels_1
  gene_labels_2 = marker_data_list$gene_labels_2
  gids = marker_data_list$gids
  gene_font_col = marker_data_list$gene_font_col


  # truncate left-side annotations
  gene_labels_1 = stringr::str_trunc(gene_labels_1, gene_chr_limit)
  gene_labels_1 = stringr::str_pad(gene_labels_1, gene_chr_limit, side="right")
  # truncate right annotations
  gene_labels_2 = stringr::str_trunc(gene_labels_2, gene_chr_limit)
  gene_labels_2 = stringr::str_pad(gene_labels_2, gene_chr_limit, side="left")

  # define matrix per metacell, based on geometric means
  mat1 = pmin( log2(niche_geomean_n[genes[gene_ord],as.character(clust_ord)] + 1), max_expression_fc )

  # create gene annotations
  if (show_gene_names) {
    if (length(highlight_genes) > 1) {
      message("Gene annots highlights")
      row_ha_right = ComplexHeatmap::HeatmapAnnotation(
        which = "row", simple_anno_size = unit(1,"mm"),
        gene = anno_mark(
          which="row", side="right", at=gids, labels=gene_labels_1[gids],
          labels_gp=gpar(fontsize = gene_font_size, col = gene_font_col[gids]),
          extend=unit(1, "mm")
        )
      )
      row_ha_left = ComplexHeatmap::HeatmapAnnotation(
        which = "row", simple_anno_size = unit(1,"mm"),
        gene = anno_mark(
          which="row", side="left", at=gids, labels=gene_labels_2[gids],
          labels_gp=gpar(fontsize = gene_font_size, col = gene_font_col[gids]),
          extend=unit(1, "mm")
        )
      )
    } else {
      message("Gene annots")
      if (verbose) message(paste(head(gene_labels_1),collapse=", "), ",...")
      row_ha_right = ComplexHeatmap::HeatmapAnnotation(
        which = "row",
        gene = anno_text(which = "row", gene_labels_1, location = 0, just = "left",
                         gp = gpar(fontsize = gene_font_size, col = gene_font_col))
      )
      if (verbose) message(paste(head(gene_labels_2),collapse=", "), ",...")
      row_ha_left = ComplexHeatmap::HeatmapAnnotation(
        which = "row",
        gene = anno_text(which = "row", gene_labels_2, location = 1, just = "right",
                         gp = gpar(fontsize = gene_font_size, col = gene_font_col))
      )
    }
  } else {
    row_ha_left = ComplexHeatmap::HeatmapAnnotation(
      which = "row", empty = anno_empty(which = "row", border = FALSE)
    )
    row_ha_right = row_ha_left
  }


  message("Expression colors...")
  col_fun = colorRampPalette(colors = heatmap_colors)
  shades = col_fun(40)

  # mc labels
  message("Metacell labels...")
  collabs <- colnames(mat1)
  if (!show_mc_names) collabs <- rep("",length(collabs))
  column_lab_ha = ComplexHeatmap::HeatmapAnnotation(
    which = "column",
    LAB = anno_text(which = "column", collabs, gp = gpar(fontsize = mc_font_size, rot=90)),
    height=unit(2,"mm")
  )
  top_column_ha = c(column_lab_ha)
  bottom_column_ha = c(column_lab_ha)

  # cluster colours
  if (!is.null(clust_col)){
    message("Columns...")
    if (is.null(names(clust_col))) {
      names(clust_col) <- clust_ord
    } else {
      if (!all(names(clust_col) %in% clust_ord))
        stop("Colour and cluster names do not match!")
      clust_col <- clust_col[clust_ord]
    }
    column_col_ha = ComplexHeatmap::HeatmapAnnotation(
      which = "column",
      'cell type' = colnames(mat1),
      col = list('cell type' = clust_col),
      border = TRUE,
      simple_anno_size = clust_anno_size,
      height = unit(1,"mm"),
      show_annotation_name = TRUE, show_legend = FALSE, gap = unit(5, "mm"),
      annotation_name_gp = gpar(fontsize = mc_font_size)
    )
    top_column_ha =    c(column_col_ha,column_lab_ha, gap = unit(1, "mm"))
    bottom_column_ha = c(column_lab_ha,column_col_ha, gap = unit(1, "mm"))

  }
  # barplots
  if (!is.null(clust_bars)) {

    if (is.null(names(clust_bars)))
      names(clust_bars) <- as.character(clust_ord)
    anno_bar <- clust_bars[as.character(clust_ord)]
    baxl <- range(clust_bars)

    column_bar_ha <- ComplexHeatmap::HeatmapAnnotation(
      which = "column",
      BAR = anno_barplot(
        anno_bar, height = 3 * clust_anno_size, bar_width = 0.9,
        gp = gpar(fill = clust_col, col = clust_col, fontsize = mc_font_size),
        axis_param = list(gp = gpar(fontsize = mc_font_size), at = baxl, labels = baxl)
      ),
      show_annotation_name = FALSE, show_legend = FALSE, gap = unit(5, "mm")
    )
    top_column_ha = c(column_bar_ha,top_column_ha, gap = unit(1, "mm"))

  }

  # additional color annotation
  if (!is.null(clust_col_others)) {
    message("Additional annots...")
    if (!is.null(col_others)) {
      column_col_ha <- HeatmapAnnotation(
        which = "column",
        df = clust_col_others, col = col_others,
        border = TRUE, simple_anno_size = clust_anno_size, title=NULL,
        show_annotation_name = TRUE, show_legend = TRUE, gap = unit(5, "mm"),
        annotation_name_gp = gpar(fontsize = mc_font_size)
      )
    } else {
      column_col_ha <- HeatmapAnnotation(
        which = "column",
        df = clust_col_others,
        border = TRUE, simple_anno_size = clust_anno_size, title=NULL,
        show_annotation_name = TRUE, show_legend = TRUE, gap = unit(5, "mm"),
        annotation_name_gp = gpar(fontsize = mc_font_size)
      )
    }
    top_column_ha <- c(column_col_ha,top_column_ha, gap = unit(1, "mm"))
    bottom_column_ha = c(bottom_column_ha,column_col_ha, gap = unit(1, "mm"))
  }

  # expression heatmap
  h1 = ComplexHeatmap::Heatmap(
    mat1, name = "expression", col = shades, #use_raster = TRUE,
    cluster_rows = FALSE, cluster_columns = FALSE,
    width = width,
    height = height,
    # column_split=clust_col,
    column_title = sprintf( "%i metacells", ncol(mat1) ),
    row_title = sprintf( "%i marker genes", nrow(mat1) ),
    show_column_names = FALSE,
    show_row_names = FALSE,
    right_annotation = row_ha_right,
    left_annotation = row_ha_left,
    top_annotation = top_column_ha,
    bottom_annotation = bottom_column_ha,
    column_names_gp = gpar(fontsize = mc_font_size),
    show_heatmap_legend = show_heatmap_legend,
    heatmap_legend_param = list(
      legend_height = unit(4, "cm"),
      border = print_border
    ),
    border = print_border
  )

  # save figure
  if (!is.null(output_file)) {
    extension <- stringr::str_extract(output_file,"(png|pdf)$")
    if (save_rds) {
      output_file_rds <- stringr::str_replace(output_file, "(png|pdg)$","rds")
      saveRDS(h1, output_file_rds)
    }

    # open graphics device
    if (extension == "png") {
      png(output_file, height = height, width = width, res=res)
    } else if (extension == "pdf") {
      pdf(output_file, height = height, width = width, useDingbats=TRUE)
    }
  }

  # draw heatmap
  draw(h1)

  # add cluster borders
  if (!is.null(clust_col) & show_clust_borders) {
    mat2 <- rbind(clust_col[match(clust_ord, names(clust_col))])
  } else {
    mat2 <- rbind(clust_ord)
  }
  if (show_clust_borders) {
    change_clust <- which(sapply(2:ncol(mat2), function(i) mat2[,i] != mat2[,i - 1]))
    decorate_heatmap_body("expression", {
      for (i in change_clust) {
        grid.lines(x = i / ncol(mat2), y = c(0,1), gp = gpar(lty = 1, lwd = 0.5))
      }
    })
  }

  if (!is.null(output_file)) dev.off()

  message("metacell heatmap done")

}

#' Plot heatmap of gene expression for single cells
#'
scp_plot_cmod_markers_sc <- function(
  marker_data_list,
  mcsc,
  umat,
  output_file = NULL,
  height = 10,
  width = 5,
  res = NA,
  show_heatmap_legend = TRUE,
  highlight_genes = NULL,
  show_gene_names = FALSE,
  gene_font_size = 5,
  clust_col = NULL, clust_bars = NULL, clust_col_others = NULL, col_others = NULL,
  clust_anno_size = unit(1,"mm"),
  show_mc_names = TRUE, mc_font_size = 5,
  heatmap_colors = c("white","gray99","orange","orangered2","#520c52"),
  gene_chr_limit = 70,
  verbose=FALSE,
  save_rds=FALSE,
  smoothen = 5,
  max_expression_fc = 5,
  max_expression_fc_sc = 5,
  print_border=TRUE,
  show_clust_borders = TRUE

) {

  # get variables necessary for hm definition (from previous function call)
  genes = marker_data_list$genes
  gene_ord = marker_data_list$gene_ord
  clust_ord = marker_data_list$clust_ord
  niche_geomean_n = marker_data_list$niche_geomean_n
  gene_labels_1 = marker_data_list$gene_labels_1
  gene_labels_2 = marker_data_list$gene_labels_2
  gids = marker_data_list$gids
  gene_font_col = marker_data_list$gene_font_col


  # truncate left-side annotations
  gene_labels_1 = stringr::str_trunc(gene_labels_1, gene_chr_limit)
  gene_labels_1 = stringr::str_pad(gene_labels_1, gene_chr_limit, side="right")
  # truncate right annotations
  gene_labels_2 = stringr::str_trunc(gene_labels_2, gene_chr_limit)
  gene_labels_2 = stringr::str_pad(gene_labels_2, gene_chr_limit, side="left")

  # directory to save output files to
  if (is.null(output_file)) {
    outdir <- getwd()
  } else {
    outdir <- dirname(output_file)
  }

  # PLOT SINGLE-CELL PROFILE
  message("Plotting single cell expression")
  cell_order=c()
  for (niche in clust_ord){
    cells=names(mcsc[which(mcsc == niche)])
    cell_order=c(cell_order,cells)
  }
  cluster_cell_count=as.matrix(table(mcsc))
  n_cells_cluster=cluster_cell_count[clust_ord,1]
  cells_clusts <- unlist(mapply(rep, clust_ord, n_cells_cluster, USE.NAMES=FALSE))

  umis=umat[,names(mcsc)]
  mat = umis[genes, ]
  totu = colSums(umis[, cell_order])
  mat = t(t(mat) / totu) * 800

  lus_1 = log2(1 + 7 * mat[genes[gene_ord], cell_order])
  lus = apply(lus_1 - apply(lus_1, 1, median),2, function(x) pmax(x,0))
  lus_smoo = t(apply(lus[genes[gene_ord],cell_order], 1, function(x) rollmean(x, smoothen, fill=0)))

  # heatmap per metacell
  mat1sc <- pmin(lus_smoo,max_expression_fc_sc)
  colnames(mat1sc) <- colnames(lus_smoo)

  # define matrix per metacell, based on geometric means
  mat1 <- pmin(log2(niche_geomean_n[genes[gene_ord],as.character(clust_ord)] + 1),max_expression_fc)

  # colors in heatmap
  col_fun = colorRampPalette(colors = heatmap_colors)
  shades = col_fun(40)

  # mc labels
  top_column_ha <- HeatmapAnnotation(
    mclabstop = anno_empty(which = "column", border = FALSE, height = unit(2,"mm"))
  )
  bottom_column_ha <- HeatmapAnnotation(
    mclabsbottom = anno_empty(which = "column", border = FALSE, height = unit(2,"mm"))
  )

  # color annotations
  if (!is.null(clust_col)){

    if (is.null(names(clust_col))) {
      names(clust_col) = clust_ord
    }
    column_col_ha <- HeatmapAnnotation(
      which = "column",
      'cluster' = cells_clusts, col = list('cluster' = clust_col),
      border = c(TRUE),
      simple_anno_size = clust_anno_size,
      show_annotation_name = TRUE, show_legend = FALSE, gap = unit(5, "mm"),
      annotation_name_gp = gpar(fontsize = mc_font_size)
    )
    top_column_ha <- c(column_col_ha,top_column_ha, gap = unit(1, "mm"))
    bottom_column_ha <- c(bottom_column_ha,column_col_ha, gap = unit(1, "mm"))

  }

  # barplot annotation
  if (!is.null(clust_bars)) {

    cell_cols <- unlist(mapply(rep, clust_col[clust_ord], n_cells_cluster, USE.NAMES=FALSE))
    cell_col <- structure(cell_cols, names=cell_order)

    if (is.null(names(clust_bars)))
      names(clust_bars) <- as.character(cell_order)

    anno_bar <- clust_bars[as.character(cell_order)]
    baxl <- range(clust_bars)

    column_bar_ha <- ComplexHeatmap::HeatmapAnnotation(
      which = "column",
      BAR = anno_barplot(
        anno_bar, height = 3 * clust_anno_size, bar_width = 1,
        gp = gpar(fill = cell_col, col = cell_col, fontsize = mc_font_size),
        axis_param = list(gp = gpar(fontsize = mc_font_size), at = baxl, labels = baxl)
      ),
      show_annotation_name = FALSE, show_legend = FALSE, gap = unit(5, "mm")
    )
    top_column_ha = c(column_bar_ha,top_column_ha)

  }
  # additional color annotation
  if (!is.null(clust_col_others)) {
    message("Additional annots...")

    clust_col_others <- clust_col_others[as.character(cell_order),colnames(clust_col_others),drop=FALSE]

    if (!is.null(col_others)) {
      column_col_ha <- HeatmapAnnotation(
        which = "column",
        df = clust_col_others, col = col_others,
        border = TRUE, simple_anno_size = clust_anno_size, title=NULL,
        show_annotation_name = TRUE, show_legend = TRUE, gap = unit(5, "mm"),
        annotation_name_gp = gpar(fontsize = mc_font_size)
      )
    } else {
      column_col_ha <- HeatmapAnnotation(
        which = "column",
        df = clust_col_others,
        border = TRUE, simple_anno_size = clust_anno_size, title=NULL,
        show_annotation_name = TRUE, show_legend = TRUE, gap = unit(5, "mm"),
        annotation_name_gp = gpar(fontsize = mc_font_size)
      )
    }
    top_column_ha <- c(column_col_ha,top_column_ha, gap = unit(1, "mm"))
    bottom_column_ha = c(bottom_column_ha,column_col_ha, gap = unit(1, "mm"))
  }

  # create gene annotations
  if (show_gene_names) {
    if (length(highlight_genes) > 1) {
      message("Gene annots highlights")
      row_ha_right <- ComplexHeatmap::HeatmapAnnotation(
        which = "row", simple_anno_size = unit(1,"mm"),
        gene = anno_mark(
          which="row", side="right", at=gids, labels=gene_labels_1[gids],
          labels_gp=gpar(fontsize = gene_font_size, col = gene_font_col[gids]),
          extend=unit(1, "mm")
        )
      )
      row_ha_left <- ComplexHeatmap::HeatmapAnnotation(
        which = "row", simple_anno_size = unit(1,"mm"),
        gene = anno_mark(
          which="row", side="left", at=gids, labels=gene_labels_2[gids],
          labels_gp=gpar(fontsize = gene_font_size, col = gene_font_col[gids]),
          extend=unit(1, "mm")
        )
      )
    } else {
      message("Gene annots")
      if (verbose) message(paste(head(gene_labels_1),collapse=", "), ",...")
      row_ha_right <- ComplexHeatmap::HeatmapAnnotation(
        which = "row", #simple_anno_size = unit(1,"mm"),
        gene = anno_text(which = "row", gene_labels_1, location = 0, just = "left", gp = gpar(
          fontsize = gene_font_size, col = gene_font_col
        ))
      )
      if (verbose) message(paste(head(gene_labels_2),collapse=", "), ",...")
      row_ha_left <- ComplexHeatmap::HeatmapAnnotation(
        which = "row", #simple_anno_size = unit(1,"mm"),
        gene = anno_text(which = "row", gene_labels_2, location = 1, just = "right", gp = gpar(
          fontsize = gene_font_size, col = gene_font_col
        ))
      )
    }
  } else {
    row_ha_left <- ComplexHeatmap::HeatmapAnnotation(
      which = "row", empty = anno_empty(which = "row", border = FALSE)
    )
    row_ha_right <- row_ha_left
  }


  # expression heatmap
  h1sc <- Heatmap(
    mat1sc, name = "sc_expression", col = shades, use_raster = TRUE,
    cluster_rows = FALSE,
    cluster_columns = FALSE,
    show_column_names = FALSE,
    show_row_names = FALSE,
    width = width,
    height = height,
    # column_split=clust_col,
    column_title = sprintf( "%i cells", ncol(mat1sc) ),
    row_title = sprintf( "%i marker genes", nrow(mat1sc) ),
    right_annotation = row_ha_right,
    left_annotation = row_ha_left,
    top_annotation = top_column_ha,
    bottom_annotation = bottom_column_ha,
    show_heatmap_legend = show_heatmap_legend,
    border = print_border
  )
  hlistsc <- h1sc

  # save figure
  if (!is.null(output_file)) {
    extension <- stringr::str_extract(output_file,"(png|pdf)$")
    if (save_rds) {
      output_file_rds <- stringr::str_replace(output_file, "(png|pdg)$","rds")
      saveRDS(h1, output_file_rds)
    }

    # open graphics device
    if (extension == "png") {
      png(output_file, height = height, width = width, res=res)
    } else if (extension == "pdf") {
      pdf(output_file, height = height, width = width, useDingbats=TRUE)
    }
  }

  # heatmap
  # draw(hlistsc, padding = unit(c(50, 50, 50, 50), "mm")) #bottom, left, top, right
  # drop all this artificial padding...
  draw(hlistsc)

  # add labels of metacells
  if (!is.null(clust_col)) {
    mat2sc <- rbind(clust_col[match(cells_clusts, names(clust_col))])
  } else {
    mat2sc <- rbind(cells_clusts)
    if (is.null(colnames(mat2sc)))
      colnames(mat2sc) <- cells_clusts
  }
  change_clust_sc <- which(sapply(2:ncol(mat2sc), function(i) mat2sc[,i] != mat2sc[,i - 1]))
  change_mc_sc <- c(
    which(sapply(2:ncol(mat2sc), function(i) colnames(mat2sc)[i] != colnames(mat2sc)[i - 1])),
    ncol(mat2sc)
  )
  if (show_clust_borders) {
    decorate_heatmap_body("sc_expression", {
      for (i in change_clust_sc) {
        grid.lines(x = i / ncol(mat2sc), y = c(0,1), gp = gpar(lty = 1, lwd = 0.5))
      }
      for (i in change_mc_sc) {
        grid.lines(x = i / ncol(mat2sc), y = c(0,1), gp = gpar(lty = 1, lwd = 0.5))
      }
    })
  }

  .add_mc_labels <- function(pos,labs) {
    for (j in 1:length(pos)) {
      i <- pos[j]
      iprev <- ifelse(j == 1,0,pos[j - 1])
      nt <- ncol(mat2sc)
      grid.text(label = labs[j], x = i / nt - (i / nt - iprev / nt) / 2, y = 0.5, gp = gpar(fontsize = mc_font_size), rot=90)
    }
  }
  if (show_mc_names) {
    decorate_annotation("mclabstop", .add_mc_labels(pos=change_mc_sc, labs=colnames(mat1)))
    decorate_annotation("mclabsbottom", .add_mc_labels(pos=change_mc_sc, labs=colnames(mat1)))
  }

  if (!is.null(output_file)) dev.off()
  ht_opt(RESET = TRUE)

  message("single-cell heatmap done")

}

# Cross species functions
