#' User interface for atlas introduction
#' @export
introUI <- function(id, label="Intro") {
  ns <- NS(id)

  shiny::tagList(
    shiny::fluidRow(
      shinydashboard::box(
        title="Metacell 2D projection", width = 8, height = 1000, solidHeader = TRUE,
        withSpinner(
          shiny::plotOutput(ns("plot_2d_proj"), height = 800),
          type = 8, color = "lightgrey", size = 0.5, hide.ui = FALSE
        )
      ),
      shinydashboard::box(
        title = "Cell type annotation", width = 4, height = 1000, solidHeader = TRUE,
        withSpinner(
          tableOutput(ns("cell_annot_table")),
          type = 8, color = "lightgrey", size = 0.5, hide.ui = FALSE
        )
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
        withSpinner(
          shiny::plotOutput(ns("plot_expression")),
          type = 8, color = "lightgrey", size = 0.5, hide.ui = FALSE
        )
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
      INPUT_DIR <- file.path(
        conf[['default']]$data_dir,
        conf[[config_id]]$data_subdir
      )

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

      clust_col <- structure(CELL_ANNT[[3]], names=CELL_ANNT[[1]])
      output$plot_expression <- shiny::renderPlot(
        scp_plot_cmod_markers_mc(
          marker_data_list = MARKER_LIST, clust_col=clust_col, clust_anno_size = unit(3,"mm"),
          show_mc_names=FALSE,  mc_font_size=6, show_gene_names=TRUE, gene_font_size=6
        ),
        height = pmin(length(MARKER_LIST$genes)*10, 2800)
      )

      # 2d mc projection
      output$plot_2d_proj <- shiny::renderPlot(
        mc_2d_plot(MC2D,CELLMC,CELL_ANNT,plot_edges=FALSE),
        height = function() {
          height = session$clientData[[sprintf("output_%s-plot_2d_proj_width",id)]]
          pmin(height, 900)
        },
        width = function() {
          width = session$clientData[[sprintf("output_%s-plot_2d_proj_width",id)]]
          pmin(width, 900)
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
        #withSpinner(
          tableOutput(ns("single_gene"))
        #  type = 8, color = "lightgrey", size = 0.5, hide.ui = FALSE
        #)
      )
    ),
    shiny::fluidRow(
      shinydashboard::box(
        title="2D projection", width=8, height = 1000, solidHeader=TRUE,
        tryCatch(
          withSpinner(
            shiny::plotOutput(ns("gene_2d"), height = 900),
            type = 8, color = "lightgrey", size = 0.5, hide.ui = FALSE
          ), error=function(e) NULL
        )
      ),
      shinydashboard::box(
        title = "Cell type annotation", width = 4, height = 1000, solidHeader = TRUE,
        withSpinner(
          tableOutput(ns("cell_annot_table")),
          type = 8, color = "lightgrey", size = 0.5, hide.ui = FALSE
        )
      )
    ),
    shiny::fluidRow(
      shinydashboard::box(
        title="", width=12, height = 900, solidHeader=TRUE,
        shiny::fluidRow(
          shiny::column(
            width=12,
            withSpinner(
              shiny::plotOutput(ns("gene_barplot"), height = 700),
              type = 8, color = "lightgrey", size = 0.5, hide.ui = FALSE
            )
          )
        ),
        shiny::fluidRow(
          shiny::column(
            width = 3,
            shiny::selectInput(
              inputId=ns("order_bars_by"),
              label="Order bars by",
              choices=c(
                "cell type" = "cell_type",
                "metacell" = "metacell"
              ),
              selectize = FALSE
            )
          ),
          shiny::column(
            width = 3,
            shiny::sliderInput(
              inputId=ns("mc_label_size"),
              label = "Metacell lables size",
              min = 0, max = 10, step = 1, value = 5
            )
          )
        )
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
      INPUT_DIR <- file.path(
        conf[['default']]$data_dir,
        conf[[config_id]]$data_subdir
      )

      mc_fp <- conf[[config_id]]$mcfp_file
      raw_counts <- conf[[config_id]]$umicount_file
      umifrac <- conf[[config_id]]$umifrac_file
      raw_counts_sc <- conf[[config_id]]$umicountsc_file
      mc2d_obj <-  conf[[config_id]]$mc2d_file
      cell_type_annotation <- conf[[config_id]]$ann_file
      gene_annotation <-  conf[[config_id]]$gene_annot_file

      MCFP <- readRDS(file = file.path(INPUT_DIR, mc_fp))
      UMICOUNTS <- readRDS(file = file.path(INPUT_DIR, raw_counts))
      UMIFRAC <- readRDS(file = file.path(INPUT_DIR, umifrac))
      UMICOUNTSC <- readRDS(file = file.path(INPUT_DIR, raw_counts_sc))
      MC2D <- readRDS(file = file.path(INPUT_DIR, mc2d_obj))
      CELL_ANNT <- fread_cell_annotation(file = file.path(INPUT_DIR, cell_type_annotation))
      ALL_GENES <- rownames(MCFP)
      GENE_ANNT <- fread_gene_annotation(
        file = file.path(INPUT_DIR, gene_annotation),
        select = c(1:3),
        col.names = c("gene_id","gene name","PFAM domain"),
        search.column = c("gene_id","gene name","PFAM domain")
      )

      updateSelectizeInput(
        session, "search_id",
        choices = GENE_ANNT[GENE_ANNT[[1]] %in% ALL_GENES, search_id],
        server = TRUE
      )

      # mini table showing gene name and id
      output$single_gene <- renderTable(
        GENE_ANNT[search_id==input$search_id, 1:(ncol(GENE_ANNT)-1)])

      # barplot of gene umifrac
      cttable = as.data.frame(CELL_ANNT)
      print("input to sg_plot")
      print(head(cttable))
      print("...")
      output$gene_barplot <- shiny::renderPlot(
        sg_plot(
          nmat=MCFP, umat=UMIFRAC, cttable=cttable, order_by=input$order_bars_by,
          sid=input$search_id, mdnorm=FALSE, annt=GENE_ANNT,
          mc_label_size=input$mc_label_size
        ),
        height = 600
      )

      # 2d projection
      output$gene_2d <- shiny::renderPlot(
        tryCatch({
          scp_plot_sc_2d_gene_exp(
            mc2d=MC2D, nmat=MCFP, umat=UMICOUNTSC,
            sid=input$search_id, annt=GENE_ANNT,
            plot_mcs=TRUE, plot_edges=FALSE, plot_mc_name=FALSE,
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
          column(width=2),
          column(
            width=3,
            sliderInput(
              inputId = ns("min_expression_fc"), label="Filter genes by min FC:",
              min=0, max=5, step=0.1, round=FALSE, value=0, width = "100%"
            )
          ),
          column(
            width=3,
            sliderInput(
              inputId = ns("scale_expression_fc"), label="Scale to max FC value:",
              min=1, max=5, step=0.1, round=FALSE, value=5, width = "100%"
            )
          ),
          column(width=1, switchInput(ns("clustergenes"), "Cluster genes", value=TRUE, inline=FALSE)),
          column(width=1),
          column(width=1, downloadButton(ns("download_genes_hmap"),"Download heatmap"))
        ),
        withSpinner(
          uiOutput(ns("ui_genes_heatmap")),
          type = 8, color = "lightgrey", size = 0.5, hide.ui = FALSE
        ),
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
      INPUT_DIR <- file.path(
        conf[['default']]$data_dir,
        conf[[config_id]]$data_subdir
      )

      mc_fp <- conf[[config_id]]$mcfp_file
      cell_type_annotation <- conf[[config_id]]$ann_file
      gene_annotation <-  conf[[config_id]]$gene_annot_file

      MCFP <- readRDS(file = file.path(INPUT_DIR, mc_fp))
      CELL_ANNT <- fread_cell_annotation(file = file.path(INPUT_DIR, cell_type_annotation))
      ALL_GENES <- rownames(MCFP)
      GENE_ANNT <- fread_gene_annotation(
        file = file.path(INPUT_DIR, gene_annotation),
        select = c(1:3),
        col.names = c("gene_id","gene name","PFAM domain"),
        search.column = c("gene_id","gene name","PFAM domain")
      )

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
          nmat=MCFP, annt=GENE_ANNT, gids=selected_genes$values,
          min_expression_fc=input$min_expression_fc,
          scale_expression_fc=pmax(input$scale_expression_fc,input$min_expression_fc),
          cluster_genes=input$clustergenes,
          heatmap_colors=heatmap_colors,
          ct_table=CELL_ANNT,
          mcid_font_size=6
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
          nmat = MCFP, gids = selected_genes$values, annt=GENE_ANNT,
          min_expression_fc=input$min_expression_fc
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
          # shiny::selectInput(
          #   inputId=ns("ids_mcs"),
          #   label="Select one or more individual metacells",
          #   choices=NULL, multiple=TRUE, selected=list(1,2),
          #   selectize = TRUE
          # )
          shiny::textInput(
            inputId=ns("ids_mcs"),
            label="Select one or more individual metacells",
            value = "1,3-4"
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
        withSpinner(
          DTOutput(ns("genes_summary_table")),
          type = 8, color = "lightgrey", size = 0.5, hide.ui = FALSE
        ),
        shiny::br(),
        shiny::h5("Total UMIs is the total UMI count for a gene in selected metacells."),
        shiny::h5("% UMIs is the percentage of UMIs for a gene that come from the selected metacells."),
        shiny::h5("Median fc is the median enrichment of UMIs for a gene in selected metacells."),
        shiny::br(),
        downloadButton(outputId=ns("download_table"), label="Download gene table"),
        downloadButton(outputId=ns("download_subset_table"), label="Download gene table for selected genes")
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
      INPUT_DIR <- file.path(
        conf[['default']]$data_dir,
        conf[[config_id]]$data_subdir
      )

      mc_fp <- conf[[config_id]]$mcfp_file
      raw_counts <- conf[[config_id]]$umicount_file
      cell_type_annotation <- conf[[config_id]]$ann_file
      gene_annotation <-  conf[[config_id]]$gene_annot_file
      tfs_annotation <-  conf[[config_id]]$tf_annot_file

      MCFP <- readRDS(file = file.path(INPUT_DIR, mc_fp))
      UMICOUNTS <- readRDS(file = file.path(INPUT_DIR, raw_counts))
      CELL_ANNT <- fread_cell_annotation(file = file.path(INPUT_DIR, cell_type_annotation))
      GENE_ANNT <- fread_gene_annotation(
        file = file.path(INPUT_DIR, gene_annotation),
        select = c(1:3),
        col.names = c("gene_id","gene name","PFAM domain"),
        search.column = c("gene_id","gene name","PFAM domain")
      )
      TF_ANNT <- fread_gene_annotation(
        file = file.path(INPUT_DIR, tfs_annotation),
        select = c(1),
        col.names=c("gene_id")
      )

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
        if (input$mcselecttype == 'input'){
          input_text = strsplit(input$ids_mcs, ",")[[1]]
          smcs <- unlist(lapply(input_text, function(text_string) {
            if (grepl("-",text_string)) {
              input_range = strsplit(text_string, "-")[[1]]
              print(input_range)
              seq(
                from = as.integer(stringr::str_remove_all(input_range[1], " ")),
                to = as.integer(stringr::str_remove_all(input_range[2], " ")),
                by = 1
              )
            } else {
              as.integer(stringr::str_remove_all(text_string, " "))
            }
          }))
          # check that they are valid mcs
          smcs[smcs %in% CELL_ANNT$mc]
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


      # download subset table
      # selected_genes_table <- eventReactive(output$download_subset_table, {
      #   gsmcs()$gene_summary[input$genes_summary_table_rows_selected]
      # })
      output$download_subset_table <- downloadHandler(
        filename <- function(){
          mc_ids_names <- red_mc_vector(selected_mcs(),range_sep="-")
          fcn <- paste0("fc",input$fc_selection)
          if (input$fcbg==TRUE)
            fcn <- paste0(fcn,"_fcbg",input$fcbg_selection)
          sprintf("mc_summary_%s_%s_selected_genes.tsv",mc_ids_names,fcn)
        },
        content <- function(file){
          write.table(
            gsmcs()$gene_summary[input$genes_summary_table_rows_selected],
            file, row.names=FALSE, col.names=TRUE, quote=FALSE, sep="\t"
          )
        }
      )

      # download full table
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

#' Cross species comparison landing page
#' @export
comparaUI <- function(id, label="Cross species comparison") {
  ns <- NS(id)

  shiny::tagList(
    shiny::fluidRow(
      shinydashboard::box(

        title="Parameters selection", width=6, solidHeader = TRUE,

        # select level of grouping
        shiny::selectInput(
          ns("level"), "Level of single-cell grouping",
          choices = c(
            # "metacells" = "mcs",
            "cell types" = "cts"
            # "broad cell types" = "bct"
          ), multiple = FALSE
        ),

        # select metric
        shiny::selectInput(
          ns("metric"), "Similarity metric",
          choices = c(
            "weighted pearson correlation" = "wpearson",
            "weighted spearman correlation" = "wspearman",
            "Jaccard index" = "jaccard",
            "Shafer index" = "shaferindex",
            "Kullback–Leibler divergence" = "kld",
            "Jensen–Shannon divergence" = "jsd",
            "orthogonal least squares" = "onls"
          ), multiple = FALSE
        ),

        # select orthologs
        shiny::selectInput(
          ns("orthos"), "Orthologs",
          choices = c(
            "all orthologs" = "all",
            "ono-to-one orthologs" = "o2o",
            "transcription factors" = "tfs"
          ), multiple = FALSE
        ),

        # select fold-change
        shiny::selectInput(
          ns("fcthrs"), "Gene fold-change threshold",
          choices = c(1.1, 1.5, 2), multiple = FALSE
        )

      ),

      shinydashboard::box(
        "Your search:", width = 6, solidHeader = TRUE,
        h5(""),
        tableOutput(ns("compara_info")),
        h5("When you change parameters, regenerate results by clicking 'Compare species' button on the sidebar.")
      )
    ),

    # heatmap
    shiny::fluidRow(
      shinydashboard::box(
        width = 8, height = 800,
        withSpinner(
          plotOutput(ns("main_heatmap"), click = ns("ht_click")),
          type = 8, color = "lightgrey", size = 0.5, hide.ui = FALSE
        )
      ),

      shinydashboard::box(
        "Selected cell type pair:", width = 4, height = 800, solidHeader = FALSE,
        br(),
        h5("Click on the cell in heatmap to see details about cell type pair similarity."),
        # file
        # shiny::textOutput(ns("compara_file")),
        br(),
        # genes
        verbatimTextOutput(ns("ht_click_content")),
        br(),
        h5("See below the table with annotation of shared genes for selected cell type.")
      )
    ),

    shiny::fluidRow(
      shinydashboard::box(
        width = 12, solidHeader = FALSE,
        DT::dataTableOutput(ns("ht_click_table"))
      )
    )

  )
}

#' Server logic for cross-species coparison
#' @export
comparaServer <- function(id, config_file="config.yaml", config_id1, config_id2) {
  shiny::moduleServer(
    id,

    function(input, output, session) {

      # updateSelectizeInput(
      #   session, "comparasp1",
      #   choices = sps, # structure(setdiff(sps,config_id2), names=names(sps)[-match(config_id2,sps)]),
      #   server = TRUE
      # )
      # observeEvent(input$comparasp2, {
      #   spsf <- structure(setdiff(sps,input$comparasp2), names=names(sps)[-match(input$comparasp2,sps)])
      #   updateSelectizeInput(session, "comparasp1", choices = spsf, server = TRUE)
      # })
      # observeEvent(input$comparasp1, {
      #   spsf <- structure(setdiff(sps,input$comparasp1), names=names(sps)[-match(input$comparasp1,sps)])
      #   print(spsf)
      #   updateSelectizeInput(session, "comparasp2", choices = spsf, server = TRUE)
      # })

      conf <- yaml::yaml.load_file(config_file, eval.expr=TRUE)
      COMPARA_DIR <-file.path(
        conf[['default']]$data_dir,
        conf[['default']]$compara_dir
      )

      # construct name of the file to load form the input parameters
      csps_file <-  file.path(COMPARA_DIR, sprintf(
        "csps_icc.%s.%s.%s-%s.%s.fc%.2f.rds",
        input$level, input$orthos, config_id1, config_id2, input$metric, as.numeric(input$fcthrs)
      ))
      message("Reading ", csps_file)
      if (file.exists(csps_file)) {
        CSPS <- readRDS(csps_file)
      } else {
        warning("switching species order")
        config_id1_orig <- config_id1
        config_id1 <- config_id2
        config_id2 <- config_id1_orig
        csps_file <-  file.path(COMPARA_DIR, sprintf(
          "csps_icc.%s.%s.%s-%s.%s.fc%.2f.rds",
          input$level, input$orthos, config_id1, config_id2, input$metric, as.numeric(input$fcthrs)
        ))
        message("Reading ", csps_file)
        CSPS <- tryCatch(
          readRDS(csps_file),
          error = function(e) NULL
        )
      }

      if (!is.null(CSPS)) {
        output$compara_file <- shiny::renderText(sprintf("%s; class %s", csps_file, class(CSPS)))

        # annotation files for species
        INPUT_DIR1 <- file.path(
          conf[['default']]$data_dir,
          conf[[config_id1]]$data_subdir
        )
        cann1_file <- file.path(INPUT_DIR1, conf[[config_id1]]$ann_file)
        cann1 <- fread(cann1_file, header=TRUE)
        gann1_file <- file.path(INPUT_DIR1, conf[[config_id1]]$gene_annot_file)
        gann1 <- fread(gann1_file, header=FALSE)
        tfs1_file <- file.path(INPUT_DIR1, conf[[config_id1]]$tf_annot_file)
        tfs1 <- fread(tfs1_file, header=FALSE)

        INPUT_DIR2 <- file.path(
          conf[['default']]$data_dir,
          conf[[config_id2]]$data_subdir
        )
        cann2_file <- file.path(INPUT_DIR2, conf[[config_id2]]$ann_file)
        cann2 <- fread(cann2_file, header=TRUE)
        gann2_file <- file.path(INPUT_DIR2, conf[[config_id2]]$gene_annot_file)
        gann2 <- fread(gann2_file, header=FALSE)
        tfs2_file <- file.path(INPUT_DIR2, conf[[config_id2]]$tf_annot_file)
        tfs2 <- fread(tfs2_file, header=FALSE)

        ann_cols <- c("cell_type","color") # this should be interactively selected, but for now only cell types
        ann1 <- unique(cann1[,..ann_cols])[,cell_type:=paste(config_id1,cell_type,sep="|")]
        ann2 <- unique(cann2[,..ann_cols])[,cell_type:=paste(config_id2,cell_type,sep="|")]

        print(sprintf("species1: %s (rows)", config_id1))
        print(sprintf("species1 cell types: %s (rows)", nrow(ann1)))
        print(sprintf("species2: %s (columns)", config_id2))
        print(sprintf("species2 cell types: %s (rows)", nrow(ann2)))

        # table with summary of search params
        compara_info_dt <- data.table(
          "clustering level" = input$level,
          "orthologs" = input$orthos,
          "species1" = config_id1,
          "species2" = config_id2,
          "similarity metric" = input$metric,
          "fc_threshold" = as.numeric(input$fcthrs),
          "variable genes" = length(CSPS$var_genes)
        ) %>% t() %>% as.data.table(keep.rownames = "parameter")
        setnames(compara_info_dt, c("parameter","value"))
        output$compara_info <- renderTable(compara_info_dt)

        print(sprintf("dim: %s x %s", nrow(CSPS$cor_matrix), ncol(CSPS$cor_matrix)))

        # non-interactive heatmap
        cor_heatmap <- csps_plot_annotated_matrix(
          mat = CSPS$cor_matrix,
          name = CSPS$method,
          row_annot = ann1, col_annot = ann2,
          fontsize = 8
        )
        output$cor_hmap_simple <- shiny::renderPlot({
          ComplexHeatmap::draw(cor_heatmap)
        })

        # interactive heatmap
        output$main_heatmap <- renderPlot({
          shiny_env$ht = draw(cor_heatmap)
          shiny_env$ht_pos = ht_pos_on_device(shiny_env$ht)
        }, width = 750, height = 750)

        # clicked pair
        output$ht_click_content <- renderText({
          if (is.null(input$ht_click)) {
            "Not selected."
          } else {
            pos1 = ComplexHeatmap:::get_pos_from_click(input$ht_click)
            ht = shiny_env$ht
            pos = selectPosition(
              ht, mark = FALSE, pos = pos1,
              verbose = FALSE, ht_pos = shiny_env$ht_pos
            )
            if (!is.null(pos)) {
              row_index = pos[1, "row_index"]
              column_index = pos[1, "column_index"]
              m = ht@ht_list[[1]]@matrix
              v = m[row_index, column_index]
              rn = rownames(m)[row_index]
              cn = colnames(m)[column_index]
              olg = CSPS$overlap_genes[[1]][[rn]][[cn]]
              if (!is.null(olg)) {
                ng = length(olg)
              } else {
                ng = ""
              }
              print(sprintf("number of genes: %s", ng))
              glue::glue(
                "{rn} (row {row_index})",
                "{cn} (column {column_index})",
                "value: {v}",
                "gene pairs: {ng}",
                .sep = "\n"
              )
            } else { "Not selected." }
          }
        })

        # table for clicked pair
        output$ht_click_table <- DT::renderDataTable({

          req(input$ht_click)
          pos1 = ComplexHeatmap:::get_pos_from_click(input$ht_click)
          ht = shiny_env$ht
          pos = selectPosition(
            ht, mark = FALSE, pos = pos1,
            verbose = FALSE, ht_pos = shiny_env$ht_pos
          )
          if (!is.null(pos)) {
            row_index = pos[1, "row_index"]
            column_index = pos[1, "column_index"]
            m = ht@ht_list[[1]]@matrix
            v = m[row_index, column_index]
            rn = rownames(m)[row_index]
            cn = colnames(m)[column_index]
            genes1 = CSPS$overlap_genes[[1]][[rn]][[cn]]
            genes2 = CSPS$overlap_genes[[2]][[rn]][[cn]]

            print(length(genes1)); print(length(genes2))

            if(!is.null(genes1)) {

              #genes1 <- str_remove(genes1,"\\.[0-9]")
              #genes2 <- str_remove(genes2,"\\.[0-9]")

              ids1 <- match(genes1, gann1[[1]])
              ids2 <- match(genes2, gann2[[1]])

              gdt <- cbind.data.frame(gann1[ids1], gann2[ids2])
              setDT(gdt)
              add_col_names <- c("bbh","pfam")
              gdtcolnames <- c(config_id1,add_col_names,config_id2,add_col_names)
              setnames(gdt, gdtcolnames)

              gdt[,tf:="no"]
              gdt[get(config_id1) %in% tfs1[[1]], tf:="yes"]
              gdt[get(config_id2) %in% tfs2[[1]], tf:="yes"]
              gdt[,tf:=factor(tf,levels=c("yes","no"))]
              setorder(gdt,tf)
              datatable(
                gdt,
                extensions=c("Buttons",'Scroller'),
                options = list(
                  dom = 'tp', scrollX = TRUE, ordering = TRUE, pageLength = 15
                )
              ) %>% formatStyle(
                'tf',
                target = 'row',
                backgroundColor = styleEqual(c("yes","no"), c("AntiqueWhite","white"))
              )

            } else {
              NULL
            }
          }

        }, rownames = FALSE, selection = list(mode = 'none'))

    }

      # TO-DOs:
      # gene annotations tables

      # Return the reactive
      # return()
    }
  )
}

