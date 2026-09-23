# Server -----
server <- function(input, output, session) {
  
  # Download Button for Userguides -----
  output$download_userguide_analysis <- downloadHandler(
    filename = "userguide_analysis.pdf",
    content = function(file) {
      file.copy(from = "pdf_userguide_analysis.pdf", to = file, overwrite = TRUE)
    }
  )
  
  output$download_userguide_sourcedata <- downloadHandler(
    filename = "userguide_sourcedata.pdf",
    content = function(file) {
      file.copy(from = "pdf_userguide_sourcedata.pdf", to = file, overwrite = TRUE)
    }
  )
  
  # Dynamic Input Change: Navigation A (RNA Enrichment Analysis) -----
  observeEvent(input$rna_ea_algorithm, {
    
    # GSEA is only pre-computed in linear analysis
    updateRadioButtons(
      session,
      "rna_group_design",
      choices = choices_design[[input$rna_ea_algorithm]],
      selected = "linear"
    )
  },
  # run update even in initial state
  ignoreInit = FALSE
  )
  
  # "Old" age has to be older than "Young" age
  observeEvent(input$age1, {
    age1_num <- as.numeric(input$age1)
    age2_valid_choices <- ages_rna[ages_rna > age1_num]
    updateSelectInput(session, "age2",
                      choices = age2_valid_choices,
                      selected = max(18, min(age2_valid_choices)))
  })
  
  # Initial Volcano & Enrichment Dot Plot Panel Title ----
  output$rna_vp_tab_title <- renderUI({
    p("Rendering graphs may take several seconds. Please wait.")
  })
  
  output$rna_ea_tab_title <- renderUI({
    p("Rendering graphs may take several seconds. Please wait.")
  })
  
  output$prot_vp_tab_title <- renderUI({
    p("Rendering graphs may take several seconds. Please wait.")
  })
  
  output$prot_ea_tab_title <- renderUI({
    p("Rendering graphs may take several seconds. Please wait.")
  })
  
  output$met_vp_tab_title <- renderUI({
    p("Rendering graphs may take several seconds. Please wait.")
  })
  
  output$met_ea_tab_title <- renderUI({
    p("Rendering graphs may take several seconds. Please wait.")
  })
  
  # Navigation A: Show Volcano Plot label update button When volcano plot is shown -----
  output$rna_vp_update_button <- renderUI({
    req(!is.null(ea_rna_data$gg))
    actionButton(inputId = "rna_submit_label_vp", label = "Update Label")
  })
  
  output$prot_vp_update_button <- renderUI({
    req(!is.null(ea_prot_data$gg))
    actionButton(inputId = "prot_submit_label_vp", label = "Update Label")
  })
  
  output$vp_metab_update_button <- renderUI({
    req(!is.null(ea_metab_data$gg_metab))
    actionButton(inputId = "submit_label_vp_metab", label = "Update Label")
  })
  
  # Dynamic Input Change: Navigation B (Pathview) -----
  observeEvent(input$source_rect_pv, {
    # update tissue list in response to data source (RNA / Protein)
    updateSelectInput(
      session,
      "tissue_rect_pv",
      choices = choices_tissue[[input$source_rect_pv]],
      selected = choices_tissue[[input$source_rect_pv]][1]
    )
    
    # initialize age choice
    if (input$source_rect_pv == "rna") {
      updateSelectInput(session, "age1_pv", choices = ages_rna[1:8], selected = 3)
      updateSelectInput(session, "age2_pv",
                        choices = ages_rna[!(ages_rna %in% c(1, 3))],
                        selected = 18)
    }
  },
  # run update even in initial state
  ignoreInit = FALSE)
  
  # "Old" age has to be older than "Young" age
  observeEvent(input$age1_pv, {
    if (input$source_rect_pv != "rna") return()
    age1_num <- as.numeric(input$age1_pv)
    age2_valid_choices <- ages_rna[ages_rna > age1_num]
    updateSelectInput(session, "age2_pv",
                      choices = age2_valid_choices,
                      selected = max(18, min(age2_valid_choices)))
  })
  
  # when source of rectangle is RNA and use all groups, metabolite cannot be shown
  observeEvent(list(input$source_rect_pv, input$design_rect_pv), ignoreInit = TRUE, {
    
    if (input$source_rect_pv == "rna" & input$design_rect_pv == "all_group") {
      updateRadioButtons(session, "source_circ_pv",
                         choiceNames = list(HTML("Not Show")),
                         choiceValues = list("not_show"),
                         selected = "not_show"
      )
    } else {
      updateRadioButtons(session, "source_circ_pv",
                         choiceNames = list(HTML("Metabolites<br> (Jankowski et al 2025)"),
                                            HTML("Not Show")),
                         choiceValues = list("metab", "not_show"),
                         selected = "metab",
                         inline = FALSE
      )
    }
  })
  
  # Define all the reactive values -----
  # Navigation A: RNA/protein enrichment analysis
  ea_rna_data <- reactiveValues(
    dat_vp = NULL,
    gg = NULL,
    dp = NULL,
    dat_ea_show = NULL,
    metadata_table = NULL,
    # parameters used for exported file names
    design = NULL,
    tissue = NULL,
    ea_algorithm = NULL,
    database = NULL,
    age1 = NULL,
    age2 = NULL
  )
  
  ea_prot_data <- reactiveValues(
    dat_vp = NULL,
    gg = NULL,
    dp = NULL,
    dat_ea_show = NULL,
    metadata_table = NULL,
    # parameters used for exported file names
    tissue = NULL,
    ea_algorithm = NULL,
    database = NULL
  )
  
  # Navigation A'': metabolite enrichment analysis
  # define reactive value
  ea_metab_data <- reactiveValues(
    dat_metab_vp = NULL,
    gg_metab = NULL,
    dp_metab = NULL,
    dat_metab_ea_show = NULL,
    metadata_table = NULL,
    # parameters used for exported file names
    tissue = NULL,
    database = NULL
  )
  
  # Navigation B: Pathview
  file_pv <- reactiveVal(NULL)
  file_pv_legend <- reactiveVal(NULL)
  has_gene <- reactiveVal(TRUE)
  has_cpd <- reactiveVal(TRUE)
  img_bin <- reactiveVal(NULL)
  
  # Navigation A (RNA/Protein Enrichment Analysis) -----
  observeEvent(input$submit_rna_ea, {
    
    # show waiter
    w_rna_ea <- Waiter$new(html = spin_3(), color = transparent(0.5))
    w_rna_ea$show()
    
    dat_vp <- NULL
    # format fold change and p-value data
    dat_vp <- my_load_degs(tissue = input$rna_tissue_ea,
                           group_design = input$rna_group_design,
                           age_young = input$age1,
                           age_old = input$age2)
    
    # Enrichment Analysis -----
    # prepare null results
    res_ea <- NULL
    dp <- NULL
    dat_ea_show <- NULL
    
    if (input$rna_ea_algorithm == "ora"){
      # define DEGs/DEPs
      diff_exps <- dat_vp %>%
        filter(padj < input$rna_p_thres, abs(log2FoldChange) > input$rna_fc_thres)
      
      # define background genes
      bg_exp <- readRDS("data/Transcriptome/gene_symbol_entrezid_table_2026-09-12_exported.rds") %>%
        pull(ENTREZID)
      
      if (input$rna_de_usage == "both") {
        diff_exps <- diff_exps %>% pull(ENTREZID)
      } else if (input$rna_de_usage == "up-regulated") {
        diff_exps <- diff_exps %>% filter(log2FoldChange > 0) %>% pull(ENTREZID)
      } else if (input$rna_de_usage == "down-regulated") {
        diff_exps <- diff_exps %>% filter(log2FoldChange < 0) %>% pull(ENTREZID)
      }
      
      # execute ORA
      res_ea_list <- my_cluterprofiler_ora(differential_exps = diff_exps,
                                           background_exps = bg_exp,
                                           bioterm_database = input$rna_database,
                                           tissue_plot_title = input$rna_tissue_ea)
      
    } else if (input$rna_ea_algorithm == "gsea") {
      gsea_table <- readRDS(paste0("data/Transcriptome/GSEA_res/GSEA_",
                                   input$rna_tissue_ea,
                                   "_",
                                   sub("^GO\\((.*)\\)$", "GO_\\1", input$rna_database),
                                   "_2026_09_14_exported.rds"))
      gsea_dp_grob <- my_gsea_dp(gsea_res = gsea_table,
                                 bioterm_database = input$rna_database,
                                 tissue_plot_title = input$rna_tissue_ea)
    }
    
    # Navigation A Rendering -----
    # change reactive value
    ea_rna_data$design <- input$rna_group_design
    ea_rna_data$tissue <- input$rna_tissue_ea
    ea_rna_data$ea_algorithm <- input$rna_ea_algorithm
    ea_rna_data$database <- input$rna_database
    
    if (input$rna_group_design == "two_group") {
      ea_rna_data$age1 <- input$age1
      ea_rna_data$age2 <- input$age2
    }
    
    ea_rna_data$dat_vp <- dat_vp
    
    if (input$rna_ea_algorithm == "ora"){
      
      ea_rna_data$gg <- my_volcano(input_data = dat_vp,
                                   p_threshold = input$rna_p_thres,
                                   fc_threshold = input$rna_fc_thres)
      
      ea_rna_data$dp <- res_ea_list$grob_dot_plot
      ea_rna_data$dat_ea_show <- res_ea_list$df_ea_show
      
    } else if (input$rna_ea_algorithm == "gsea"){
      
      ea_rna_data$gg <- my_volcano(input_data = dat_vp,
                                   p_threshold = 0.05,
                                   fc_threshold = 0.5)
      
      ea_rna_data$dp <- gsea_dp_grob
      ea_rna_data$dat_ea_show <- gsea_table
    }
    
    # metadata table
    if (input$rna_group_design == "linear"){
      ea_rna_data$metadata_table <- metadata_table_rna %>%
        filter(tissue == input$rna_tissue_ea)
    } else if (input$rna_group_design == "two_group") {
      ea_rna_data$metadata_table <- metadata_table_rna %>%
        filter(tissue == input$rna_tissue_ea & month %in% paste0("month", c(input$age1, input$age2)))
    }
    
    
    # Output in Navigation A
    # tab title 
    output$rna_vp_tab_title <- renderUI({
      if (is.null(ea_rna_data$tissue)) {
        return(NULL)
      }
      title_text <- if (ea_rna_data$design == "linear") {
        paste0(ea_rna_data$tissue, " bulk RNA-seq across all ages")
      } else {
        paste0(ea_rna_data$tissue, " bulk RNA-seq ", ea_rna_data$age1, "-month vs ", ea_rna_data$age2, "-month")
      }
      p(strong(title_text))
    })
    
    # tab title 
    output$rna_ea_tab_title <- renderUI({
      if (is.null(ea_rna_data$tissue)) {
        return(NULL)
      }
      title_text <- if (ea_rna_data$design == "linear") {
        paste0(ea_rna_data$tissue, " bulk RNA-seq Enrichment Analysis across all ages")
      } else {
        paste0(ea_rna_data$tissue, " bulk RNA-seq Enrichment Analysis ", ea_rna_data$age1, "-month vs ", ea_rna_data$age2, "-month")
      }
      p(strong(title_text))
    })
    
    # interactive volcano plot
    output$rna_graph_vp <- ggiraph::renderGirafe({
      girafe(ggobj = ea_rna_data$gg,
             width_svg = 6,
             height_svg = 4.5,
             options = list(
               opts_sizing(rescale = FALSE),
               opts_toolbar(
                 saveaspng = FALSE,
                 hidden = c("selection", "zoom", "misc")
               )
             )
      )
    })
    
    # Update Selectize List for shown label in VP
    updateSelectizeInput(
      session,
      "rna_vp_show_labels",
      choices = ea_rna_data$gg[["data"]]$label_hover,
      selected = NULL,
      server = TRUE
    )
    
    # Download Button for VP
    output$rna_download_vp <- downloadHandler(
      filename = function() {
        if_else(ea_rna_data$design == "linear",
                paste0("bulkRNA_volcano_plot_", ea_rna_data$tissue, "_across_all_ages.png"),
                paste0("bulkRNA_volcano_plot_", ea_rna_data$tissue, "_", ea_rna_data$age1, "vs", ea_rna_data$age2, ".png")) 
      },
      content = function(file) {
        ggsave(file, plot = ea_rna_data$gg, width = 6, height = 4.5, dpi = 300)
      }
    )
    
    # Download Button for Result Table
    output$rna_download_de_table <- downloadHandler(
      filename = function() {
        if_else(ea_rna_data$design == "linear",
                paste0("bulkRNA_DE_result_", ea_rna_data$tissue, "_across_all_ages.csv"),
                paste0("bulkRNA_DE_result_", ea_rna_data$tissuea, "_", ea_rna_data$age1, "vs", ea_rna_data$age2, ".csv")) 
      },
      content = function(file) {
        write.csv(ea_rna_data$dat_vp, file, row.names = FALSE)
      }
    )
    
    # Table
    output$rna_table_vp <- renderReactable({
      if (is.null(ea_rna_data$dat_vp) || nrow(ea_rna_data$dat_vp) == 0) {
        reactable(data.frame(Message = "No results to display"))
      } else {
        reactable(
          ea_rna_data$dat_vp,
          searchable = TRUE,
          theme = reactableTheme(style = list(fontSize = "12px")),
          # reduce digit to look window tidy
          columns = list(
            log2FoldChange = colDef(name = "covariate-adjusted log2FC", format = colFormat(digits = 5)),
            pvalue = colDef(name = "p-val", format = colFormat(digits = 3)),
            padj = colDef(name = "FDR (BH-adjusted p)", format = colFormat(digits = 3)),
            ENTREZID = colDef(show = FALSE)
          )
        )
      }
    }) 
    
    # Enrichment Dot Plot
    output$rna_graph_dp <- renderPlot(
      {ea_rna_data$dp}
    )
    
    # Download Button for Dot Plot
    output$rna_download_dp <- downloadHandler(
      filename = function() {
        if_else(ea_rna_data$design == "linear",
                paste0("bulkRNA_", ea_rna_data$database, "_", ea_rna_data$ea_algorithm, "_enrichment_plot_", ea_rna_data$tissue, "_across_all_ages.png"),
                paste0("bulkRNA_", ea_rna_data$database, "_", ea_rna_data$ea_algorithm, "_enrichment_plot_", ea_rna_data$tissue, "_", ea_rna_data$age1, "vs", ea_rna_data$age2, ".png")) 
      },
      content = function(file) {
        ggsave(file, plot = ea_rna_data$dp, width = 9, height = 9, dpi = 300)
      }
    )
    
    # Download Button for Result Table
    output$rna_download_ea_table <- downloadHandler(
      filename = function() {
        if_else(ea_rna_data$design == "linear",
                paste0("bulkRNA_", ea_rna_data$database, "_", ea_rna_data$ea_algorithm, "_enrichment_result_", ea_rna_data$tissue, "_across_all_ages.csv"),
                paste0("bulkRNA_", ea_rna_data$database, "_", ea_rna_data$ea_algorithm, "_enrichment_result_", ea_rna_data$tissue, "_", ea_rna_data$age1, "vs", ea_rna_data$age2, ".csv")) 
      },
      content = function(file) {
        write.csv(ea_rna_data$dat_ea_show, file, row.names = FALSE)
      }
    )
    
    output$rna_table_ea <- renderReactable({
      if (is.null(ea_rna_data$dat_ea_show) || nrow(ea_rna_data$dat_ea_show) == 0) {
        reactable(data.frame(Message = "No enrichment results to display"))
      } else {
        reactable(
          ea_rna_data$dat_ea_show,
          searchable = TRUE,
          theme = reactableTheme(style = list(fontSize = "10px"))
        )
      }
    })
    
    # show metadata
    output$rna_metadata_table <- renderTable(ea_rna_data$metadata_table)
    
    # Download Button for Metadata Table
    output$rna_download_meta_table <- downloadHandler(
      filename = function() {
        if_else(ea_rna_data$design == "linear",
                paste0("bulkRNA_sample_size_", ea_rna_data$tissue, "_across_all_ages.csv"),
                paste0("bulkRNA_sample_size_", ea_rna_data$tissue, "_", ea_rna_data$age1, "vs", ea_rna_data$age2, ".csv")) 
      },
      content = function(file) {
        write.csv(ea_rna_data$metadata_table, file, row.names = FALSE)
      }
    )
    
    
    # hide waiter
    w_rna_ea$hide()
  })
  
  # update volcano plot if labels to show is updated -----
  observeEvent(input$rna_submit_label_vp, ignoreInit = TRUE, {
    ea_rna_data$gg <- my_volcano(input_data = ea_rna_data$dat_vp,
                                 p_threshold = ifelse(input$rna_ea_algorithm == "gsea", 0.05 ,input$rna_p_thres),
                                 fc_threshold = ifelse(input$rna_ea_algorithm == "gsea", 0.5 ,input$rna_fc_thres))
    
    if (!is.null(input$rna_vp_show_labels) && length(input$rna_vp_show_labels) > 0){
      ea_rna_data$gg <- ea_rna_data$gg +
        geom_point(data = ea_rna_data$dat_vp %>% filter(SYMBOL %in% input$rna_vp_show_labels),
                   aes(x = log2FoldChange, y = -log10(padj)),
                   size = 2, shape = 21, color = "black", fill = "magenta") +
        ggrepel::geom_label_repel(data = ea_rna_data$dat_vp %>% filter(SYMBOL %in% input$rna_vp_show_labels),
                                  aes(label = SYMBOL, x = log2FoldChange, y = -log10(padj)),
                                  alpha = 0.8,
                                  color = "magenta",
                                  size = 4,
                                  force = 4,
                                  max.overlaps = 10)
    }
  })
  
  
  # Navigation A' (Protein Enrichment Analysis) -----
  observeEvent(input$submit_prot_ea, {
    
    # show waiter
    w_prot_ea <- Waiter$new(html = spin_3(), color = transparent(0.5))
    w_prot_ea$show()
    
    dat_vp <- NULL
    # format fold change and p-value data
    dat_vp <- readRDS(paste0("data/Proteome/Proteome_DE_res/CellRep2023_DE_", input$prot_tissue_ea, "_2026-09-15.rds"))
    
    # Enrichment Analysis -----
    # prepare null results
    res_ea <- NULL
    dp <- NULL
    dat_ea_show <- NULL
    
    if (input$prot_ea_algorithm == "ora"){
      # define DEGs/DEPs
      diff_exps <- dat_vp %>%
        filter(padj < input$prot_p_thres, abs(log2FoldChange) > input$prot_fc_thres)
      
      # define background genes
      bg_exp <- dat_vp %>%
        pull(ENTREZID) %>%
        unique()
      
      if (input$prot_de_usage == "both") {
        diff_exps <- diff_exps %>% pull(ENTREZID) %>% unique()
      } else if (input$prot_de_usage == "up-regulated") {
        diff_exps <- diff_exps %>% filter(log2FoldChange > 0) %>% pull(ENTREZID) %>% unique()
      } else if (input$prot_de_usage == "down-regulated") {
        diff_exps <- diff_exps %>% filter(log2FoldChange < 0) %>% pull(ENTREZID) %>% unique()
      }
      
      # execute ORA
      res_ea_list <- my_cluterprofiler_ora(differential_exps = diff_exps,
                                           background_exps = bg_exp,
                                           bioterm_database = input$prot_database,
                                           tissue_plot_title = input$prot_tissue_ea)
      
    } else if (input$prot_ea_algorithm == "gsea") {
      gsea_table <- readRDS(paste0("data/Proteome/Proteome_GSEA_res/GSEA_proteome_",
                                   input$prot_tissue_ea,
                                   "_",
                                   sub("^GO\\((.*)\\)$", "GO_\\1", input$prot_database),
                                   "_2026_09_15_exported.rds"))
      gsea_dp_grob <- my_gsea_dp(gsea_res = gsea_table,
                                 bioterm_database = input$prot_database,
                                 tissue_plot_title = input$prot_tissue_ea)
    }
    
    # Navigation A Rendering -----
    # change reactive value
    ea_prot_data$tissue <- input$prot_tissue_ea
    ea_prot_data$ea_algorithm <- input$prot_ea_algorithm
    ea_prot_data$database <- input$prot_database
    
    ea_prot_data$dat_vp <- dat_vp
    
    if (input$prot_ea_algorithm == "ora"){
      
      ea_prot_data$gg <- my_volcano(input_data = dat_vp,
                                    p_threshold = input$prot_p_thres,
                                    fc_threshold = input$prot_fc_thres,
                                    tooltip_label = "symbol")
      
      ea_prot_data$dp <- res_ea_list$grob_dot_plot
      ea_prot_data$dat_ea_show <- res_ea_list$df_ea_show
      
    } else if (input$prot_ea_algorithm == "gsea"){
      
      ea_prot_data$gg <- my_volcano(input_data = dat_vp,
                                    p_threshold = 0.1,
                                    fc_threshold = 0,
                                    tooltip_label = "symbol")
      
      ea_prot_data$dp <- gsea_dp_grob
      ea_prot_data$dat_ea_show <- gsea_table
    }
    
    # metadata table
    ea_prot_data$metadata_table <- metadata_table_prot %>%
      filter(tissue == input$prot_tissue_ea)
    
    # Output in Navigation A
    # tab title
    output$prot_vp_tab_title <- renderUI({
      if (is.null(ea_prot_data$tissue)) {
        return(NULL)
      }
      title_text <- paste0(str_to_title(ea_prot_data$tissue), " Proteome 8-month vs 18-month")
      p(strong(title_text))
    })
    
    # tab title 
    output$prot_ea_tab_title <- renderUI({
      if (is.null(ea_prot_data$tissue)) {
        return(NULL)
      }
      title_text <- paste0(str_to_title(ea_prot_data$tissue), " Proteome Enrichment Analysis 8-month vs 18-month")
      p(strong(title_text))
    })
    
    # interactive volcano plot
    output$prot_graph_vp <- ggiraph::renderGirafe({
      girafe(ggobj = ea_prot_data$gg,
             width_svg = 6,
             height_svg = 4.5,
             options = list(
               opts_sizing(rescale = FALSE),
               opts_toolbar(
                 saveaspng = FALSE,
                 hidden = c("selection", "zoom", "misc")
               )
             )
      )
    })
    
    # Update Selectize List for shown label in VP
    updateSelectizeInput(
      session,
      "prot_vp_show_labels",
      choices = unique(ea_prot_data$gg[["data"]]$label_hover),
      selected = NULL,
      server = TRUE
    )
    
    # Download Button for VP
    output$prot_download_vp <- downloadHandler(
      filename = paste0("Proteome_volcano_plot_", ea_prot_data$tissue, ".png"),
      content = function(file) {
        ggsave(file, plot = ea_prot_data$gg, width = 6, height = 4.5, dpi = 300)
      }
    )
    
    # Download Button for Result Table
    output$prot_download_de_table <- downloadHandler(
      filename = paste0("Proteome_DE_result_", ea_prot_data$tissue, ".csv"),
      content = function(file) {
        write.csv(ea_prot_data$dat_vp, file, row.names = FALSE)
      }
    )
    
    # Table
    output$prot_table_vp <- renderReactable({
      if (is.null(ea_prot_data$dat_vp) || nrow(ea_prot_data$dat_vp) == 0) {
        reactable(data.frame(Message = "No results to display"))
      } else {
        reactable(
          ea_prot_data$dat_vp,
          searchable = TRUE,
          theme = reactableTheme(style = list(fontSize = "12px")),
          # reduce digit to look window tidy
          columns = list(
            log2FoldChange = colDef(name = "covariate-adjusted log2FC", format = colFormat(digits = 5)),
            pvalue = colDef(name = "p-val", format = colFormat(digits = 3)),
            padj = colDef(name = "FDR (BH-adjusted p)", format = colFormat(digits = 3)),
            ENTREZID = colDef(show = FALSE)
          )
        )
      }
    }) 
    
    # Enrichment Dot Plot
    output$prot_graph_dp <- renderPlot(
      {ea_prot_data$dp}
    )
    
    # Download Button for Dot Plot
    output$prot_download_dp <- downloadHandler(
      filename = paste0("Proteome_", ea_prot_data$database, "_", ea_prot_data$ea_algorithm, "_enrichment_plot_", ea_prot_data$tissue, ".png"),
      content = function(file) {
        ggsave(file, plot = ea_prot_data$dp, width = 9, height = 9, dpi = 300)
      }
    )
    
    # Download Button for Result Table
    output$prot_download_ea_table <- downloadHandler(
      filename = paste0("Proteome_", ea_prot_data$database, "_", ea_prot_data$ea_algorithm, "_enrichment_result_", ea_prot_data$tissue, ".csv"),
      content = function(file) {
        write.csv(ea_prot_data$dat_ea_show, file, row.names = FALSE)
      }
    )
    
    output$prot_table_ea <- renderReactable({
      if (is.null(ea_prot_data$dat_ea_show) || nrow(ea_prot_data$dat_ea_show) == 0) {
        reactable(data.frame(Message = "No enrichment results to display"))
      } else {
        reactable(
          ea_prot_data$dat_ea_show,
          searchable = TRUE,
          theme = reactableTheme(style = list(fontSize = "10px"))
        )
      }
    })
    
    # show metadata
    output$prot_metadata_table <- renderTable(ea_prot_data$metadata_table)
    
    # Download Button for Metadata Table
    output$prot_download_meta_table <- downloadHandler(
      filename = paste0("Proteome_sample_size_", ea_prot_data$tissue, ".csv"),
      content = function(file) {
        write.csv(ea_prot_data$metadata_table, file, row.names = FALSE)
      }
    )
    
    # hide waiter
    w_prot_ea$hide()
  })
  
  # update volcano plot if labels to show is updated -----
  observeEvent(input$prot_submit_label_vp, ignoreInit = TRUE, {
    ea_prot_data$gg <- my_volcano(input_data = ea_prot_data$dat_vp,
                                  p_threshold = ifelse(input$prot_ea_algorithm == "gsea", 0.1 ,input$prot_p_thres),
                                  fc_threshold = ifelse(input$prot_ea_algorithm == "gsea", 0 ,input$prot_fc_thres),
                                  tooltip_label = "symbol")
    
    if (!is.null(input$prot_vp_show_labels) && length(input$prot_vp_show_labels) > 0){
      ea_prot_data$gg <- ea_prot_data$gg +
        geom_point(data = ea_prot_data$dat_vp %>% filter(symbol %in% input$prot_vp_show_labels),
                   aes(x = log2FoldChange, y = -log10(padj)),
                   size = 2, shape = 21, color = "black", fill = "magenta") +
        ggrepel::geom_label_repel(data = ea_prot_data$dat_vp %>% filter(symbol %in% input$prot_vp_show_labels),
                                  aes(label = symbol, x = log2FoldChange, y = -log10(padj)),
                                  alpha = 0.8,
                                  color = "magenta",
                                  size = 4,
                                  force = 4,
                                  max.overlaps = 10)
    }
  })
  
  
  # Navigation A'' (Metabolite Enrichment) -----
  observeEvent(input$submit_metab_ea, {
    
    # show waiter
    w_metab_ea <- Waiter$new(html = spin_3(), color = transparent(0.5))
    w_metab_ea$show()
    
    # ORA -----
    dat_metab_vp <- dat_metab %>%
      filter(tissue == input$met_tissue_ea) %>%
      dplyr::select(-tissue) %>%
      left_join(dat_compound_id, by = "compound") %>%
      dplyr::select(compound, log2FC_Aged, pval_Aged, padj_Aged, HMDB, KEGG, PubChem, METLIN, kegg_name) %>%
      dplyr::rename(log2FoldChange = log2FC_Aged,
                    pvalue = pval_Aged,
                    padj = padj_Aged)
    
    # define differentially expressed metabolites and background metabolites and extract compound IDs
    if(input$metab_database == "KEGG"){
      metab_col_pick <- c("KEGG", "kegg_name")
    } else if (input$metab_database == "SMPDB"){
      metab_col_pick <- c("HMDB", "compound")
    }
    
    # extract background compound IDs
    bg_metabs <- dat_metab_vp %>%
      dplyr::select(all_of(metab_col_pick[1])) %>%
      setNames("id") %>%
      drop_na(id) %>%
      pull(id)
    
    # extract differentially regulated compound IDs
    diff_metabs <- dat_metab_vp %>%
      filter(padj < input$metab_p_thres, abs(log2FoldChange) > input$metab_fc_thres) %>%
      dplyr::select(all_of(metab_col_pick[1]), log2FoldChange, padj) %>%
      dplyr::rename(id = all_of(metab_col_pick[1])) %>%
      drop_na(id)
    
    # define background genes
    
    if (input$metab_de_usage == "both") {
      diff_metabs <- diff_metabs %>% pull(id)
    } else if (input$metab_de_usage == "up-regulated") {
      diff_metabs <- diff_metabs %>% filter(log2FoldChange > 0) %>% pull(id)
    } else if (input$metab_de_usage == "down-regulated") {
      diff_metabs <- diff_metabs %>% filter(log2FoldChange < 0) %>% pull(id)
    }
    
    # execute ORA
    # prepare null results
    res_metab_ea <- NULL
    dp_metab <- NULL
    dat_meatab_ea_show <- NULL
    
    # Perform ORA and Create dotplot object
    if(input$metab_universe == "measured"){
      res_metab_ea_list <- my_microbiomeprofiler_ora(differential_metabs = diff_metabs,
                                                     specify_background = TRUE,
                                                     background_metabs = bg_metabs,
                                                     bioterm_database = input$metab_database,
                                                     cid_table = dat_compound_id,
                                                     tissue_plot_title = input$met_tissue_ea)
    } else if (input$metab_universe == "all_known") {
      res_metab_ea_list <- my_microbiomeprofiler_ora(differential_metabs = diff_metabs,
                                                     specify_background = FALSE,
                                                     bioterm_database = input$metab_database,
                                                     cid_table = dat_compound_id,
                                                     tissue_plot_title = input$met_tissue_ea)
    }
    
    # Navigation A' Rendering -----
    ea_metab_data$tissue <- input$met_tissue_ea
    ea_metab_data$database <- input$metab_database
    
    # change reactive value
    ea_metab_data$gg_metab <- my_volcano(input_data = dat_metab_vp,
                                         p_threshold = input$metab_p_thres,
                                         fc_threshold = input$metab_fc_thres,
                                         tooltip_label = "compound")
    ea_metab_data$dat_metab_vp <- dat_metab_vp
    ea_metab_data$dp_metab <- res_metab_ea_list$grob_metab_dot_plot
    ea_metab_data$dat_metab_ea_show <- res_metab_ea_list$df_metab_ea_show
    
    # metadata table
    ea_metab_data$metadata_table <- metadata_table_met %>%
      filter(tissue == input$met_tissue_ea)
    
    
    # Output in Navigation A'
    # tab title
    output$met_vp_tab_title <- renderUI({
      if (is.null(ea_metab_data$tissue)) {
        return(NULL)
      }
      title_text <- paste0(ea_metab_data$tissue, " Metabolome 3-6 months vs ≥21 months")
      p(strong(title_text))
    })
    
    # tab title 
    output$met_ea_tab_title <- renderUI({
      if (is.null(ea_metab_data$tissue)) {
        return(NULL)
      }
      title_text <- paste0(ea_metab_data$tissue, " Metabolome Enrichment Analysis 3-6 months vs ≥21 months")
      p(strong(title_text))
    })
    
    # interactive volcano plot
    output$graph_metab_vp <- ggiraph::renderGirafe({
      girafe(ggobj = ea_metab_data$gg_metab,
             width_svg = 6,
             height_svg = 4.5,
             options = list(
               opts_sizing(rescale = FALSE),
               opts_toolbar(
                 saveaspng = FALSE,
                 hidden = c("selection", "zoom", "misc")
               )
             )
      )
    })
    
    # Update Selectize List for shown label in VP
    updateSelectizeInput(
      session,
      "vp_metab_show_labels",
      choices = ea_metab_data$dat_metab_vp$compound,
      selected = NULL,
      server = TRUE
    )
    
    # Download Button for VP
    output$download_metab_vp <- downloadHandler(
      filename = paste0("Metabolome_volcano_plot_", ea_metab_data$tissue, ".png"),
      content = function(file) {
        ggsave(file, plot = ea_metab_data$gg_metab, width = 6, height = 4.5, dpi = 300)
      }
    )
    
    # Download Button for Result Table
    output$download_metab_de_table <- downloadHandler(
      filename = paste0("Metabolome_DE_result_", ea_metab_data$tissue, ".csv"),
      content = function(file) {
        write.csv(ea_metab_data$dat_metab_vp, file, row.names = FALSE)
      }
    )
    
    # Table
    output$table_metab_vp <- renderReactable({
      if (is.null(ea_metab_data$dat_metab_vp) || nrow(ea_metab_data$dat_metab_vp) == 0) {
        reactable(data.frame(Message = "No results to display"))
      } else {
        reactable(
          ea_metab_data$dat_metab_vp,
          searchable = TRUE,
          theme = reactableTheme(style = list(fontSize = "12px")),
          # reduce digit to look window tidy
          columns = list(
            log2FoldChange = colDef(name = "covariate-adjusted log2FC", format = colFormat(digits = 5)),
            pvalue = colDef(name = "p-val", format = colFormat(digits = 3)),
            padj = colDef(name = "FDR (BH-adjusted p)", format = colFormat(digits = 3))
          )
        )
      }
    }) 
    
    # Enrichment Dot Plot
    output$graph_metab_dp <- renderPlot(
      {ea_metab_data$dp_metab}
    )
    
    # Download Button for Dot Plot
    output$download_metab_dp <- downloadHandler(
      filename = paste0("Metabolome_", ea_metab_data$database, "_ora_enrichment_plot_", ea_metab_data$tissue, ".png"),
      content = function(file) {
        ggsave(file, plot = ea_metab_data$dp_metab, width = 9, height = 9, dpi = 300)
      }
    )
    
    # Download Button for Result Table
    output$download_metab_ea_table <- downloadHandler(
      filename = paste0("Metabolome_", ea_metab_data$database, "_ora_enrichment_result_", ea_metab_data$tissue, ".csv"),
      content = function(file) {
        write.csv(ea_metab_data$dat_metab_ea_show, file, row.names = FALSE)
      }
    )
    
    output$table_metab_ea <- renderReactable({
      if (is.null(ea_metab_data$dat_metab_ea_show) || nrow(ea_metab_data$dat_metab_ea_show) == 0) {
        reactable(data.frame(Message = "No enrichment results to display"))
      } else {
        reactable(
          ea_metab_data$dat_metab_ea_show,
          searchable = TRUE,
          theme = reactableTheme(style = list(fontSize = "10px"))
        )
      }
    })
    
    # show metadata
    output$met_metadata_table <- renderTable(ea_metab_data$metadata_table)
    
    # Download Button for Metadata Table
    output$met_download_meta_table <- downloadHandler(
      filename = paste0("Metabolome_sample_size_", ea_metab_data$tissue, ".csv"),
      content = function(file) {
        write.csv(ea_metab_data$metadata_table, file, row.names = FALSE)
      }
    )
    
    # hide waiter
    w_metab_ea$hide()
  })
  
  # update volcano plot if metabolite labels to show is updated -----
  observeEvent(input$submit_label_vp_metab, ignoreInit = TRUE, {
    ea_metab_data$gg_metab <-  my_volcano(input_data = ea_metab_data$dat_metab_vp,
                                          p_threshold = input$metab_p_thres,
                                          fc_threshold = input$metab_fc_thres,
                                          tooltip_label = "compound")
    
    if (!is.null(input$vp_metab_show_labels) && length(input$vp_metab_show_labels) > 0){
      ea_metab_data$gg_metab <- ea_metab_data$gg_metab +
        geom_point(data = ea_metab_data$dat_metab_vp %>% filter(compound %in% input$vp_metab_show_labels),
                   aes(x = log2FoldChange, y = -log10(padj)),
                   size = 2, shape = 21, color = "black", fill = "magenta") +
        ggrepel::geom_label_repel(data = ea_metab_data$dat_metab_vp %>% filter(compound %in% input$vp_metab_show_labels),
                                  aes(label = compound, x = log2FoldChange, y = -log10(padj)),
                                  alpha = 0.6,
                                  color = "magenta",
                                  size = 4,
                                  force = 4,
                                  max.overlaps = 10)
    }
  })
  
  # Navigation B (Pathview)-----
  
  observeEvent(input$submit_pv, {
    # show waiter
    w_pv <- Waiter$new(html = spin_3(), color = transparent(0.5))
    w_pv$show()
    
    # initialize binary image
    img_bin(NULL)
    
    # define gene/protein and compound fold changes to put in pathview()
    # gene/protein rectangle fold changes
    if (input$source_rect_pv == "not_show") {
      rect_fc_show <- NULL
    } else if (input$source_rect_pv == "protein") {
      rect_fc_show <- readRDS(paste0("data/Proteome/Proteome_DE_res/CellRep2023_DE_", input$tissue_rect_pv, "_2026-09-15.rds")) %>%
        filter(!is.na(ENTREZID)) %>%
        # protein data sometimes has duplicated ENRTEZID (bc single Ensembl protein ID sometimes has multiple gene IDs)
        arrange(desc(abs(log2FoldChange))) %>%
        distinct(ENTREZID, .keep_all = TRUE) %>%
        pull(log2FoldChange, name = ENTREZID)
    } else if (input$source_rect_pv == "rna" & input$design_rect_pv == "two_group") {
      rect_fc_show <- my_load_degs(tissue = input$tissue_rect_pv,
                                   group_design = "two_group",
                                   age_young = input$age1_pv,
                                   age_old = input$age2_pv) %>%
        pull(log2FoldChange, name = ENTREZID)
    } else if (input$source_rect_pv == "rna" & input$design_rect_pv == "all_group") {
      rect_fc_show <- readRDS(paste0("data/Transcriptome/RNA_log2FC/log2FC_mtx_", input$tissue_rect_pv, "_2026_09_12_exported.rds"))
    }
    
    # compound fold changes
    if (input$source_circ_pv == "not_show" || (input$source_rect_pv == "rna" & input$design_rect_pv == "all_group")){
      circ_fc_show <- NULL
    } else if (input$source_circ_pv == "metab"){
      circ_fc_show <- dat_metab %>%
        filter(tissue == input$tissue_circ_pv) %>%
        dplyr::select(-tissue) %>%
        left_join(dat_compound_id, by = "compound") %>%
        filter(!is.na(KEGG)) %>%
        pull(log2FC_Aged, name = KEGG)
    }
    
    # Run pathview
    gene_limit_choice <- list(not_show = 1,
                              rna = 1,
                              protein = 0.5)
    pathview_res <- my_pathview(gene_fc = rect_fc_show,
                                cpd_fc = circ_fc_show,
                                kegg_id = substring(input$kegg_id_pv, 1, 8),
                                # Change color-scale range based on the data source (RNAseq : range -1to1, Protein: range -0.5to0.5)
                                gene_limit = gene_limit_choice[[input$source_rect_pv]])
    
    # detect whether the pathway has rectangle/circle to show
    has_gene(my_has_gene_rect(substring(input$kegg_id_pv, 1, 8), species = "mmu"))
    has_cpd(my_has_cpd_circ(substring(input$kegg_id_pv, 1, 8), species = "mmu"))
    
    # Remove the extra files from server
    unlink(paste0(substring(input$kegg_id_pv, 1, 8), ".xml"))  # Remove XML file
    if (!all(is.null(rect_fc_show), is.null(circ_fc_show))){
      unlink(paste0(substring(input$kegg_id_pv, 1, 8), ".png"))  # Remove PNG file
    }
    
    # specify .png file name for legend
    if (!is.null(rect_fc_show) && is.matrix(rect_fc_show)){
      file_pv_legend("pathview_legend_multi.png")
    } else {
      file_pv_legend("pathview_legend.png")
    }
    
    # clean up unnecessary objects
    rm(rect_fc_show, circ_fc_show)
    gc();gc()
    
    # change filename to output
    file_pv(pathview_res$pathview_filename)
    # make binary file of pathview image (for downloading)
    img_bin(readBin(file_pv(), "raw", n = file.info(file_pv())$size))
    
    # legend
    output$pathview_legend <- renderImage(
      {
        list(src = file_pv_legend(),
             alt = "Pathview Legend",
             height = "80%")
      },
      deleteFile = FALSE
    )
    
    output$pathview_message <- renderUI({
      msgs <- c()
      
      if (!has_gene()) {
        msgs <- c(msgs, "No gene/protein nodes (rectangle) exist in this pathway.")
      }
      if (!has_cpd()) {
        msgs <- c(msgs, "No metabolite nodes (small circle) exist in this pathway.")
      }
      
      if (length(msgs) == 0) { return(NULL) }
      
      div(
        style = "color:#D9534F; font-style: italic; padding-top:5px;",
        paste(msgs, collapse = " | ")
      )
    })
    
    # Download Button for legend and pathview
    output$download_pv_legend <- downloadHandler(
      filename = "pathview_legend.png",
      content = function(file) {
        file.copy(from = file_pv_legend(), to = file, overwrite = TRUE)
      }
    )
    
    output$download_pv_plot <- downloadHandler(
      filename = "pathview_image.png",
      content = function(file) {
        # make binary object
        writeBin(img_bin(), file)
      },
      contentType = "image/png"
    )
    
    # Update image output
    output$pB <- renderImage({
      list(src = file_pv(),
           contentType = "image/png",
           alt = "Pathview Image",
           height = "80%")
    }, deleteFile = TRUE)
    
    # hide waiter
    w_pv$hide()
  })
  
}