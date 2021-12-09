
loadUI <- function(id, label) {
  ns <- shiny::NS(id)
}

loadServer <- function(id, config_file="config.yaml", config_id) {
  shiny::moduleServer(
    id,
    function(input, output, session) {

      iconf <- yaml::yaml.load_file(config_file, eval.expr=TRUE)
      INPUT_DIR <- conf[['default']]$data_dir

      mc_fp <- conf[[config_id]]$mcfp_file
      umi_frac <- conf[[config_id]]$umifrac_file
      raw_counts <- conf[[config_id]]$umicount_file
      raw_counts_sc <- conf[[config_id]]$umicountsc_file
      mc2d_obj <-  conf[[config_id]]$mc2d_file
      cell_type_annotation <- conf[[config_id]]$ann_file


      MCFP <- readRDS(file = file.path(INPUT_DIR, mc_fp))
      UMIFRAC <- readRDS(file = file.path(INPUT_DIR, umi_frac))
      UMICOUNTS <- readRDS(file = file.path(INPUT_DIR, raw_counts))
      UMICOUNTSC <- readRDS(file = file.path(INPUT_DIR, raw_counts_sc))
      MC2D <- readRDS(file = file.path(INPUT_DIR, mc2d_obj))
      CELL_ANNT <- fread_cell_annotation(file = file.path(INPUT_DIR, cell_type_annotation))

      return(reactive({
        MCFP
        UMIFRAC
        UMICOUNTS
        UMICOUNTSSC
        CELL_ANNT
      }))
    }
  )
}

#' User interface for atlas introduction
#' @export
introUI <- function(id, label="Intro") {
  ns <- NS(id)

  shiny::tagList(
    shiny::fluidRow(
      shinydashboard::box(
        title="Metacell 2D projection", width = 8, height = 800, solidHeader = TRUE,
        shiny::plotOutput(ns("plot_2d_proj"), height = 700)
      ),
      shinydashboard::box(
        title = "Cell type annotation", width = 4, height = 800, solidHeader = TRUE,
        tableOutput(ns("cell_annot_table"))
      )
    ),
    shiny::fluidRow(
      shinydashboard::box(
        title="Gene expression heatmap", width = 12, height = 3000, solidHeader = TRUE,
        shiny::h5(
          "Heatmap: relative enrichment of UMIs for genes (rows)
              in metacells (column number IDs)
              grouped in cell types (column color bar annotation).
              Barplot: symbiodinium signal in metacells."
        ),
        #imageOutput(ns("plot_expression"))
        shiny::plotOutput(ns("plot_expression"))
      )
    )
  )
}

#' Server logic for atlas introduction
#' @export
introServer <- function(id, config_file="config.yaml", config_id) {
  shiny::moduleServer(
    id,

    function(input, output, session) {

      conf <- yaml::yaml.load_file(config_file, eval.expr=TRUE)
      INPUT_DIR <- conf[['default']]$data_dir

      raw_counts_sc <- conf[[config_id]]$umicountsc_file
      mc2d_obj <-  conf[[config_id]]$mc2d_file
      cellmc_obj <-  conf[[config_id]]$cellmc_file
      cell_type_annotation <- conf[[config_id]]$ann_file
      gene_expression_plot <- conf[[config_id]]$heatmap_file
      marker_data_list <- conf[[config_id]]$markers_file

      UMICOUNTSC <- readRDS(file = file.path(INPUT_DIR, raw_counts_sc))
      MC2D <- readRDS(file = file.path(INPUT_DIR, mc2d_obj))
      CELLMC <- readRDS(file = file.path(INPUT_DIR, cellmc_obj))
      CELL_ANNT <- fread_cell_annotation(file = file.path(INPUT_DIR, cell_type_annotation))
      MARKER_LIST <- readRDS(file = file.path(INPUT_DIR, marker_data_list))

      plot_expression_filepath <- file.path(INPUT_DIR, gene_expression_plot)

      # expression heatmap
      # output$plot_expression <- renderImage({
      #   dims_img <- dim(png::readPNG(plot_expression_filepath))
      #   height_img <- as.integer(dims_img[1])
      #   width_img  <- as.integer(dims_img[2])
      #   ratio <- height_img/width_img
      #   width  <- pmin(width_img,1200)
      #   height <- width*ratio
      #   list(
      #     src = plot_expression_filepath, alt="Expression heatmap",
      #     width = width, height = height
      #   )
      # }, deleteFile = FALSE)

      clust_col <- structure(CELL_ANNT[[3]], names=CELL_ANNT[[1]])
      output$plot_expression <- shiny::renderPlot(
        #scp_plot_cmod_markers_sc(marker_data_list = MARKER_LIST, mcsc = CELLMC, umat = UMICOUNTSC, clust_col=clust_col)
        scp_plot_cmod_markers_mc(
          marker_data_list = MARKER_LIST, clust_col=clust_col, clust_anno_size = unit(3,"mm"),
          show_mc_names=TRUE,  mc_font_size=6, show_gene_names=TRUE, gene_font_size=6
        ),
        height = pmin(length(MARKER_LIST$genes)*10, 2800)
      )

      # 2d mc projection
      output$plot_2d_proj <- shiny::renderPlot(
        mc_2d_plot(MC2D,CELLMC,CELL_ANNT),
        height = function() {
          height = session$clientData[[sprintf("output_%s-plot_2d_proj_width",id)]]
          pmin(height, 700)
        },
        width = function() {
          width = session$clientData[[sprintf("output_%s-plot_2d_proj_width",id)]]
          pmin(width, 700)
        }
      )
      output$cell_annot_table <- function() {
        summarize_cell_annotation(CELL_ANNT)
      }

      # Return the reactive
      # return()
    }
  )
}

#' User interface for visualizing single gene expression
#' @export
singleGeneUI <- function(id, label="Single gene expression") {
  ns <- NS(id)

  shiny::tagList(
    shiny::fluidRow(
      shinydashboard::box(
        title="Gene search", width=6, solidHeader = TRUE,
        shiny::selectInput(
          inputId=ns("search_id"),
          label="Search gene by name or id",
          choices=NULL,
          selectize = TRUE
        )
      ),
      shinydashboard::box(
        title="Your search:", width=6, solidHeader=TRUE,
        tableOutput(ns("single_gene"))
      )
    ),
    shiny::fluidRow(
      shinydashboard::box(
        title="2D projection", width=8, height = 1000, solidHeader=TRUE,
        tryCatch(shiny::plotOutput(ns("gene_2d"), height = 900), error=function(e) NULL)
      ),
      shinydashboard::box(
        title = "Cell type annotation", width = 4, height = 1000, solidHeader = TRUE,
        tableOutput(ns("cell_annot_table"))
      )
    ),
    shiny::fluidRow(
      shinydashboard::box(
        title="", width=12, height = 900, solidHeader=TRUE,
        shiny::plotOutput(ns("gene_barplot"), height = 700)
      )
    )
  )

}

#' Server logic for visualizing single gene expression
#' @export
singleGeneServer <- function(id,  config_file="config.yaml", config_id) {

  shiny::moduleServer(
    id,

    function(input, output, session) {

      conf <- yaml::yaml.load_file(config_file, eval.expr=TRUE)
      INPUT_DIR <- conf[['default']]$data_dir

      mc_fp <- conf[[config_id]]$mcfp_file
      raw_counts <- conf[[config_id]]$umicount_file
      umifrac <- conf[[config_id]]$umifrac_file
      raw_counts_sc <- conf[[config_id]]$umicountsc_file
      mc2d_obj <-  conf[[config_id]]$mc2d_file
      cell_type_annotation <- conf[[config_id]]$ann_file

      MCFP <- readRDS(file = file.path(INPUT_DIR, mc_fp))
      UMICOUNTS <- readRDS(file = file.path(INPUT_DIR, raw_counts))
      UMIFRAC <- readRDS(file = file.path(INPUT_DIR, umifrac))
      UMICOUNTSC <- readRDS(file = file.path(INPUT_DIR, raw_counts_sc))
      MC2D <- readRDS(file = file.path(INPUT_DIR, mc2d_obj))
      CELL_ANNT <- fread_cell_annotation(file = file.path(INPUT_DIR, cell_type_annotation))
      ALL_GENES <- rownames(MCFP)

      updateSelectizeInput(
        session, "search_id",
        choices = GENE_ANNT[GENE_ANNT[[1]] %in% ALL_GENES, search_id],
        server = TRUE
      )

      # mini table showing gene name and id
      output$single_gene <- renderTable(
        GENE_ANNT[search_id==input$search_id, 1:(ncol(GENE_ANNT)-1)])

      # barplot of gene umifrac
      output$gene_barplot <- shiny::renderPlot(
        sg_plot(
          nmat=MCFP, umat=UMIFRAC, cttable=CELL_ANNT,
          sid=input$search_id, mdnorm=FALSE, annt=GENE_ANNT
        ),
        height = 600
      )

      # 2d projection
      output$gene_2d <- shiny::renderPlot(
        tryCatch({
          scp_plot_sc_2d_gene_exp(
            mc2d=MC2D, nmat=MCFP, umat=UMICOUNTSC,
            sid=input$search_id, annt=GENE_ANNT,
            plot_mcs=TRUE, plot_edges=TRUE, plot_mc_name=TRUE,
            do_umifrac_sc=TRUE, sc_max=NULL, sc_zero_color = "aliceblue"
          )
        },  error = function(e) plot.new()),
        height = function() {
          height = session$clientData[[sprintf("output_%s-gene_2d_width",id)]]
          pmin(height, 900)
        },
        width = function() {
          width = session$clientData[[sprintf("output_%s-gene_2d_width",id)]]
          pmin(width, 900)
        }
      )
      output$cell_annot_table <- function() {
        summarize_cell_annotation(CELL_ANNT)
      }

    }
  )

}

#' User interface for visualizing multiple gene expression
#' @export
multiGeneUI <- function(id, label="Multi gene expression") {
  ns <- NS(id)

  shiny::tagList(
    shiny::fluidRow(
      shinydashboard::box(
        title="Gene selection", width=12, height=NULL, solidHeader = TRUE,
        shiny::h5("Choose genes for which to plot the expression across metacells.
              Either search for genes by typing (part of) the name or gene id in the search bar,
              or upload a text file with genes."),
        shiny::h5("Select genes to show on the heatmap either by clicking on the individual rows
              of the search results table followed by 'Add selected genes', or just by clicking
              'Add all genes' to include all the search results."
        ), br(),
        radioButtons(
          ns("geneselecttype"), label = "Choose genes",
          choices = list("search for genes" = "search",  "upload list of genes" = "upload"),
          selected = "search"
        ),
        conditionalPanel(
          condition = "input.geneselecttype == 'search'", ns = ns,
          searchInput(
            inputId = ns("free_genes"),
            label = "Search for gene:",
            placeholder = "Type gene name or ID",
            btnSearch = icon("search"),
            btnReset = icon("remove"),
            width = "50%"
          )
        ),
        conditionalPanel(
          condition = "input.geneselecttype == 'upload'", ns = ns,
          shiny::h5("Upload a text text file with genes. Each gene should be on the new line."),
          fileInput(
            ns("genefile"), "Choose gene file",
            multiple = FALSE,
            accept = c("text/tsv","text/csv","text/tab-separated-values,text/plain","text/comma-separated-values,text/plain")
          )
        ),
        DTOutput(ns("geneSelectDT")),
        hr(),
        actionButton(ns("add_genes"), "Add selected genes"),
        actionButton(ns("add_all_genes"), "Add all genes")
      )
    ),
    shiny::fluidRow(
      shinydashboard::box(
        title="Choosen genes:", width=12, solidHeader=TRUE,
        DTOutput(ns("selected_genes_dt")), hr(),
        actionButton(ns("clear_genes"), "Clear selection"),
        downloadButton(ns("download_genes_data"),"Download table")
      )
    ),
    shiny::fluidRow(
      shinydashboard::box(
        title="Heatmap", width=12, solidHeader=TRUE,
        shiny::fluidRow(
          column(
            width=4,
            sliderInput(
              inputId = ns("fcthrshmap"), label="Filter genes by min FC:",
              min=0, max=5, step=0.1, round=FALSE, value=0, width = "100%"
            )
          ),
          column(
            width=4,
            sliderInput(
              inputId = ns("maxvalhmap"), label="Scale color to maximum value:",
              min=1, max=5, step=0.1, round=FALSE, value=5, width = "100%"
            )
          ),
          column(width=2, switchInput(ns("clustergenes"), "Cluster genes", value=TRUE, inline=FALSE)),
          column(width=1, downloadButton(ns("download_genes_hmap"),"Download heatmap"))
        ),
        uiOutput(ns("ui_genes_heatmap")),
        imageOutput(ns("plot_ct_legend_horizontal_heatmap"))
      )
    )
  )

}

#' Server logic for visualizing multiple gene expression
#' @export
multiGeneServer <- function(id, config_file="config.yaml", config_id) {

  shiny::moduleServer(
    id,

    function(input, output, session) {

      conf <- yaml::yaml.load_file(config_file, eval.expr=TRUE)
      INPUT_DIR <- conf[['default']]$data_dir

      mc_fp <- conf[[config_id]]$mcfp_file
      cell_type_annotation <- conf[[config_id]]$ann_file

      MCFP <- readRDS(file = file.path(INPUT_DIR, mc_fp))
      CELL_ANNT <- fread_cell_annotation(file = file.path(INPUT_DIR, cell_type_annotation))
      ALL_GENES <- rownames(MCFP)

      # Select genes
      genes_dt <- reactive(
        if (input$geneselecttype=="search") {
          genes_select_dt(sterm=input$free_genes,nmat=MCFP,annt=GENE_ANNT)
        } else if (input$geneselecttype=="upload") {
          req(input$genefile)
          gs <- fread(input$genefile$datapath,header=FALSE)[[1]]
          genes_upload_dt(gs=gs, annt=GENE_ANNT)
        }
      )
      output$geneSelectDT <- renderDT(
        genes_dt(), rownames = FALSE,
        options = list(
          dom = 'tp', scrollX = TRUE, ordering = FALSE, pageLength = 10,
          columnDefs = list(list(className = 'dt-center', targets = 0:2))
        )
      )
      selected_genes <- reactiveValues()
      selected_genes$values <- c()
      observeEvent(input$add_genes, {
        g <- genes_select_names(dt=genes_dt(),rid=input$geneSelectDT_rows_selected)
        selected_genes$values <- c(selected_genes$values, g)
      })
      observeEvent(input$add_all_genes, {
        g <- genes_dt()$gene_id
        selected_genes$values <- c(selected_genes$values, g)
      })
      observeEvent(input$clear_genes, {
        selected_genes$values <- c()
      })
      genes_selected <- reactive({
        rid <- unique(match(selected_genes$values,GENE_ANNT$gene_id))
        tryCatch(GENE_ANNT[rid, 1:3], error=function(e) NULL)
      })
      output$selected_genes_dt <- renderDT(
        if (!is.null(selected_genes$values))
          genes_selected(),
        rownames = FALSE,
        options = list(
          dom = 'tp', scrollX = TRUE, ordering = FALSE, pageLength = 10,
          columnDefs = list(list(className = 'dt-center', targets = 0:2))
        ),
        selection = list(mode = 'none')
      )
      output$download_genes_data <- downloadHandler(
        filename = function() {
          "selected_genes.tsv"
        },
        content = function(file) {
          write.table(genes_selected(), file, sep = "\t", quote = FALSE)
        }
      )

      # Heatmap of selected genes
      plotting_f <- function() {
        tryCatch(mgenes_hmap(
          nmat=MCFP, annt=GENE_ANNT, gids=selected_genes$values, fcthrs=input$fcthrshmap,
          palette=heatmap_colors, maxval=pmax(input$maxvalhmap,input$fcthrshmap),
          ct_table=CELL_ANNT, cluster_genes=input$clustergenes, mcid_font_size=6
        ),  error = function(e) message(e))
      }
      output$genes_hmap <- shiny::renderPlot({
        if (!is.null(selected_genes$values)) {
          print(plotting_f())
        }
      })
      hmh <- reactiveValues()
      hmap_height <- function() {
        hmh$ng <- mgenes_hmap_height(
          nmat = MCFP, gids = selected_genes$values, annt=GENE_ANNT, fcthrs=input$fcthrshmap
        )
        if (hmh$ng<5) {
          sf <- 60
        } else if (hmh$ng<15) {
          sf <- 50
        } else if (hmh$ng<25) {
          sf <- 20
        } else {
          sf <- 15
        }
        return(hmh$ng*sf)
      }
      output$ui_genes_heatmap <- renderUI({
        ns <- session$ns
        if (!is.null(selected_genes$values))
          shiny::plotOutput(ns("genes_hmap"), height = hmap_height(), width = "100%")
      })
      output$download_genes_hmap <- downloadHandler(
        filename = "selected_genes_heatmap.png",
        content = function(file) {
          png(file,width=2400,height=hmap_height(),res=110)
          print(plotting_f())
          dev.off()
        }
      )


    })
}

#' User interface for summarizing expression in group of metacells
#' @export
summaryUI <- function(id, config_file="config.yaml", label="Metacell summary") {
  ns <- NS(id)

  shiny::tagList(
    shiny::fluidRow(
      shinydashboard::box(
        title="1. Metacells selection", width=4, solidHeader = TRUE,
        shiny::h5("Find genes that are specifically expressed in a group of metacells."),
        shiny::h5("It's possible to select metacells in two ways: either by
               inputing individual metacells IDs (select 'Input'), or by
               selecting all metacells annotated as a cell type (select 'Annotation')."),
        br(),
        radioButtons(
          ns("mcselecttype"), label = "Choose metacells from:",
          choices = list("Input" = "input",  "Annotation" = "file"),
          selected = "input"
        ),
        conditionalPanel(
          condition = "input.mcselecttype == 'input'", ns=ns,
          shiny::selectInput(
            inputId=ns("ids_mcs"),
            label="Select one or more individual metacells",
            choices=NULL, multiple=TRUE, selected=list(1,2),
            selectize = TRUE
          )
        ),
        conditionalPanel(
          condition = "input.mcselecttype == 'file'", ns=ns,
          shiny::h5("Select a cell type from existing annotations or upload your own annotation file."),
          materialSwitch(ns("uploadann"),"Upload annotation file",value=FALSE),
          conditionalPanel(
            condition = "input.uploadann", ns=ns,
            shiny::h5("Upload a tab-separated text file with the following columns: metacell, cell type, color"),
            fileInput(
              ns("annfile"), "Choose annotation file",
              multiple = FALSE,
              accept = c("text/tsv","text/tab-separated-values,text/plain",".tsv")
            )
          ),
          uiOutput(ns("ctselect"))
        )
      ),
      shinydashboard::box(
        title="2. Fold change selection.", width=4, solidHeader = TRUE,
        shiny::h5("Select minimum fold change threshold to be used for filtering genes. Genes that have fc
               above this value in all selected metacells (if using absolute method) or median fc
               above this value in all selected metacells (if using median method) will be shown in the
               summary table."),
        radioButtons(
          ns("fc_sum_method"),"Select method to use for summarizing gene expression in selected metacells:",
          choices = c("absolute","median"), selected = "median", inline = TRUE
        ),
        sliderInput(
          inputId = ns("fc_selection"), label="Choose minimum fold change for selected metacells",
          min=1.0, max=5, step=0.1, round=FALSE, value=2
        ),
        shiny::h5("It is also possible to set a maximum fc threshold for non-selected (background) metacells.
               In this case, the summary will only include genes that have fc above minimum threshold
               in selected metacells, and below maximum theshold in all other metacells."),
        checkboxInput(ns("fcbg"), "Choose maximum background fold change threshold", value = FALSE),
        conditionalPanel(
          condition = "input.fcbg", ns=ns,
          radioButtons(
            ns("fcbg_sum_method"),"Select method to use for summarizing gene expression in non-selected metacells:",
            choices = c("absolute","median"), selected = "absolute", inline = TRUE
          ),
          sliderInput(
            inputId = ns("fcbg_selection"), label="Choose maximum fold change for non-selected metacells",
            min=1.0, max=5, step=0.1, round=FALSE, value=2
          )
        )
      ),
      shinydashboard::box(
        title="3. Generate summary", width=4, solidHeader=TRUE,
        shiny::h5("Click to generate summary using specified parameters."),
        actionButton(ns("dosummary"), "Go!"),
        h4("Your search: "),
        tableOutput(ns("summary"))
      )
    ),
    shiny::fluidRow(
      shinydashboard::box(
        title="Gene table", width=12, solidHeader=TRUE,
        DTOutput(ns("genes_summary_table")), shiny::br(),
        shiny::h5("Total UMIs is the total UMI count for a gene in selected metacells."),
        shiny::h5("% UMIs is the percentage of UMIs for a gene that come from the selected metacells."),
        shiny::h5("Median fc is the median enrichment of UMIs for a gene in selected metacells."),
        shiny::br(),
        downloadButton(outputId=ns("download_table"), label="Download gene table")
      )
    )
  )
}

#' Server logic for summarizing expression in group of metacells
#' @export
summaryServer <- function(id, config_file="config.yaml", config_id) {
  shiny::moduleServer(
    id,

    function(input, output, session) {

      conf <- yaml::yaml.load_file(config_file, eval.expr=TRUE)
      INPUT_DIR <- conf[['default']]$data_dir

      mc_fp <- conf[[config_id]]$mcfp_file
      raw_counts <- conf[[config_id]]$umicount_file
      cell_type_annotation <- conf[[config_id]]$ann_file

      MCFP <- readRDS(file = file.path(INPUT_DIR, mc_fp))
      UMICOUNTS <- readRDS(file = file.path(INPUT_DIR, raw_counts))
      CELL_ANNT <- fread_cell_annotation(file = file.path(INPUT_DIR, cell_type_annotation))

      mcinputdt <- eventReactive(
        eventExpr = {
          input$annfile
          input$uploadann
        },
        valueExpr = {
          if(input$uploadann == TRUE) {
            req(input$annfile)
            fread_cell_annotation(input$annfile$datapath)
          } else {
            CELL_ANNT
          }
        }
      )
      output$anndt <- renderDT(mcinputdt())

      output$ctselect <- renderUI({
        ns <- session$ns
        shiny::selectInput(
          ns("ids_cts"), "Select cell type",
          choices=mcinputdt()[[2]],
          multiple=FALSE, selected=list(1)
        )
      })

      updateSelectizeInput(
        session, "ids_mcs",
        choices = CELL_ANNT[[1]],
        server = TRUE
      )

      selected_mcs <- reactive({
        if(input$mcselecttype == 'input'){
          input$ids_mcs
        } else if(input$mcselecttype == 'file'){
          mcinputdt()[cell_type==input$ids_cts,mc]
        }
      })

      gsmcs <- eventReactive(
        eventExpr = input$dosummary,
        valueExpr = {
          ns <- session$ns
          print(input$fc_sum_method)
          mc_gene_summary(
            mc_ids=selected_mcs(), method=input$fc_sum_method, fc=input$fc_selection,
            usefcbg=input$fcbg, fcbg=input$fcbg_selection, methodbg=input$fcbg_sum_method,
            mc_fp=MCFP, mc_counts=UMICOUNTS, annt=GENE_ANNT, tfannt=TF_ANNT
          )
        }
      )

      # table summarizing gap genes in queried metacells
      output$genes_summary_table <- DT::renderDataTable(
        datatable(gsmcs()$gene_summary) %>% formatStyle('tf', target='row', backgroundColor = styleEqual(c("yes","no"), c("AntiqueWhite","white"))),
        rownames = FALSE,
        options = list(
          dom = 'tp', scrollX = TRUE, ordering = TRUE, pageLength = 20,
          columnDefs = list(list(className = 'dt-center', targets = 0:2))
        ),
        selection = list(mode = 'none')
      )

      # table summarizing query and results
      output$summary  <- renderTable(
        gsmcs()$summary
      )

      # download multiple metacell query table
      output$download_table <- downloadHandler(
        filename <- function(){
          mc_ids_names <- red_mc_vector(selected_mcs(),range_sep="-")
          fcn <- paste0("fc",input$fc_selection)
          if (input$fcbg==TRUE)
            fcn <- paste0(fcn,"_fcbg",input$fcbg_selection)
          sprintf("mc_summary_%s_%s.tsv",mc_ids_names,fcn)
        },
        content <- function(file){
          write.table(
            gsmcs()$gene_summary,
            file, row.names=FALSE, col.names=TRUE, quote=FALSE, sep="\t"
          )
        }
      )
    }
  )

}
