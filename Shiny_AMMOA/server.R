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
  
  # Dynamic Input Change: Navigation A (RNA/Protein Enrichment Analysis) -----
  observeEvent(input$datasource, {
    # update tissue list in response to data source (RNA / Protein)
    updateSelectInput(
      session,
      "tissue_ea",
      choices = choices_tissue[[input$datasource]],
      selected = choices_tissue[[input$datasource]][1]
    )
    
    # Update available statistical design (only 8month vs 18 month is available for protein)
    updateRadioButtons(
      session,
      "design",
      choices = choices_design[[input$datasource]],
      selected = "two_group"
    )
    
    # initialize age choice
    if (input$datasource == "protein") {
      updateSelectInput(session, "age1", choices = 8, selected = 8)
      updateSelectInput(session, "age2", choices = 18, selected = 18)
    } else if (input$datasource == "rna") {
      updateSelectInput(session, "age1", choices = ages_rna[1:8], selected = 3)
      updateSelectInput(session, "age2",
                        choices = ages_rna[!(ages_rna %in% c(1, 3))],
                        selected = 18)
    }
    
    # run update even in initial state
    ignoreInit = FALSE
  })
  
  # "Old" age has to be older than "Young" age
  observeEvent(input$age1, {
    req(input$datasource == "rna")
    age1_num <- as.numeric(input$age1)
    age2_valid_choices <- ages_rna[ages_rna > age1_num]
    updateSelectInput(session, "age2",
                      choices = age2_valid_choices,
                      selected = max(18, min(age2_valid_choices)))
  })
  
  # Navigation A and A': Show Volcano Plot label update button When volcano plot is shown -----
  output$vp_update_button <- renderUI({
    req(!is.null(ea_data$gg))
    actionButton(inputId = "submit_label_vp", label = "Update Label")
  })
  
  # bulk RNAseq only: whether or not show uncharacterized genes
  # if I use the condition 'input$datasource == "rna"', checkbox (dis)appears when sidebar input changes
  # this I use condition equation using ea_data reactive value because it is only updated after "submit" button is clicked
  output$vp_uncharacterized_check <- renderUI({
    req(!is.null(ea_data$gg) && !any(colnames(ea_data$dat_vp) == "ensembl_id"))
    checkboxInput(inputId = "vp_uncharacterized",
                  label = "Show Uncharacterized Genes",
                  value = ea_data$vp_unchar_state)
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
    
    # run update even in initial state
    ignoreInit = FALSE
  })
  
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
  ea_data <- reactiveValues(
    dat_vp = NULL,
    gg = NULL,
    dp = NULL,
    dat_ea_show = NULL,
    # show uncharacterized gene or not
    vp_unchar_state = FALSE
  )
  
  # Navigation A': metabolite enrichment analysis
  # define reactive value
  ea_metab_data <- reactiveValues(
    dat_metab_vp = NULL,
    gg_metab = NULL,
    dp_metab = NULL,
    dat_metab_ea_show = NULL
  )
  
  # Navigation B: Pathview
  file_pv <- reactiveVal(NULL)
  file_pv_legend <- reactiveVal(NULL)
  has_gene <- reactiveVal(TRUE)
  has_cpd <- reactiveVal(TRUE)
  img_bin <- reactiveVal(NULL)
  
  # Navigation A (RNA/Protein Enrichment Analysis) -----
  observeEvent(input$submit_ea, {
    library(ggiraph)
    
    # show waiter
    w_ea <- Waiter$new(html = spin_3(), color = transparent(0.5))
    w_ea$show()
    
    
    # ORA -----
    dat_vp <- NULL
    # format fold change and p-value data
    if (input$datasource == "rna"){
      # for RNA, perforn DESeq2 processing
      dat_vp <- my_deseq2_est(meta_data = rna_metadata,
                              target_tissue = input$tissue_ea,
                              est_design = input$design,
                              rawcount_dir = "data/RNA_raw",
                              age_1 = input$age1,
                              age_2 = input$age2) %>%
        left_join(., gene_ids, by = c("symbol" = "SYMBOL"))
    } else if (input$datasource == "protein"){
      # for promote data, filter prepossessed data into selected tissue
      dat_vp <- dat_prot %>%
        filter(tissue == input$tissue_ea) %>%
        dplyr::select(-tissue)
    }
    
    # define DEGs/DEPs
    diff_exps <- dat_vp %>%
      drop_na(ENTREZID) %>%
      filter(padj < input$p_thres, abs(log2FoldChange) > input$fc_thres)
    
    # define background genes
    bg_exp <- dat_vp %>%
      drop_na(ENTREZID) %>%
      pull(ENTREZID)
    
    if (input$deg_usage == "both") {
      diff_exps <- diff_exps %>% pull(ENTREZID)
    } else if (input$deg_usage == "up-regulated") {
      diff_exps <- diff_exps %>% filter(log2FoldChange > 0) %>% pull(ENTREZID)
    } else if (input$deg_usage == "down-regulated") {
      diff_exps <- diff_exps %>% filter(log2FoldChange < 0) %>% pull(ENTREZID)
    }
    
    # execute ORA
    # prepare null results
    res_ea <- NULL
    dp <- NULL
    dat_ea_show <- NULL
    
    res_ea_list <- my_cluterprofiler_ora(differential_exps = diff_exps,
                                         background_exps = bg_exp,
                                         bioterm_database = input$database)
    
    # Navigation A Rendering -----
    # change reactive value
    if (input$datasource == "rna"){
      ea_data$gg <- my_volcano(input_data = dat_vp %>% dplyr::filter(symbol %in% gene_ids$SYMBOL[gene_ids$characterized]),
                               p_threshold = input$p_thres,
                               fc_threshold = input$fc_thres)
    } else if (input$datasource == "protein"){
      ea_data$gg <- my_volcano(input_data = dat_vp,
                               p_threshold = input$p_thres,
                               fc_threshold = input$fc_thres)
    }

    ea_data$dat_vp <- dat_vp
    ea_data$dp <- res_ea_list$grob_dot_plot
    ea_data$dat_ea_show <- res_ea_list$df_ea_show
    
    
    # Output in Navigation A
    # interactive volcano plot
    output$graph_vp <- ggiraph::renderGirafe({
      girafe(ggobj = ea_data$gg,
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
      "vp_show_labels",
      choices = ea_data$gg[["data"]]$label_hover,
      selected = NULL,
      server = TRUE
    )
    
    # Download Button for VP
    output$download_vp <- downloadHandler(
      filename = "volcano_plot.png",
      content = function(file) {
        ggsave(file, plot = ea_data$gg, width = 6, height = 4.5, dpi = 300)
      }
    )
    
    # Download Button for Result Table
    output$download_de_table <- downloadHandler(
      filename = function() {
        if_else(input$design == "linear",
                paste0("result_", input$tissue_ea, "_", input$datasource, "_", input$design, ".csv"),
                paste0("result_", input$tissue_ea, "_", input$datasource, "_", input$age1, "vs", input$age2, ".csv")) 
      },
      content = function(file) {
        write.csv(ea_data$dat_vp, file, row.names = FALSE)
      }
    )
    
    # Table
    output$table_vp <- renderReactable({
      if (is.null(ea_data$dat_vp) || nrow(ea_data$dat_vp) == 0) {
        reactable(data.frame(Message = "No results to display"))
      } else {
        reactable(
          ea_data$dat_vp,
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
    output$graph_dp <- renderPlot(
      {ea_data$dp}
    )
    
    # Download Button for Dot Plot
    output$download_dp <- downloadHandler(
      filename = "enrichment_dot_plot.png",
      content = function(file) {
        ggsave(file, plot = ea_data$dp, width = 9, height = 9, dpi = 300)
      }
    )
    
    # Download Button for Result Table
    output$download_ea_table <- downloadHandler(
      filename = "all_enriched_terms.csv",
      content = function(file) {
        write.csv(ea_data$dat_ea_show, file, row.names = FALSE)
      }
    )
    
    output$table_ea <- renderReactable({
      if (is.null(ea_data$dat_ea_show) || nrow(ea_data$dat_ea_show) == 0) {
        reactable(data.frame(Message = "No enrichment results to display"))
      } else {
        reactable(
          ea_data$dat_ea_show,
          searchable = TRUE,
          theme = reactableTheme(style = list(fontSize = "10px"))
        )
      }
    })
    
    # hide waiter
    w_ea$hide()
  })
  
  # update volcano plot if labels to show is updated -----
  observeEvent(input$submit_label_vp, ignoreInit = TRUE, {
    if (input$datasource == "rna" && !ea_data$vp_unchar_state){
      ea_data$gg <- my_volcano(input_data = ea_data$dat_vp %>% dplyr::filter(symbol %in% gene_ids$SYMBOL[gene_ids$characterized]),
                               p_threshold = input$p_thres,
                               fc_threshold = input$fc_thres)
    } else {
      ea_data$gg <- my_volcano(input_data = ea_data$dat_vp,
                               p_threshold = input$p_thres,
                               fc_threshold = input$fc_thres)
    }
    
    if (!is.null(input$vp_show_labels) && length(input$vp_show_labels) > 0){
      ea_data$gg <- ea_data$gg +
        geom_point(data = ea_data$dat_vp %>% filter(symbol %in% input$vp_show_labels),
                   aes(x = log2FoldChange, y = -log10(padj)),
                   size = 2, shape = 21, color = "black", fill = "magenta") +
        ggrepel::geom_label_repel(data = ea_data$dat_vp %>% filter(symbol %in% input$vp_show_labels),
                                  aes(label = symbol, x = log2FoldChange, y = -log10(padj)),
                                  alpha = 0.8,
                                  color = "magenta",
                                  size = 4,
                                  force = 4,
                                  max.overlaps = 10)
    }
  })
  
  # bulkRNA-seq only : switch show/ not show uncharacterized gene names
  observeEvent(input$vp_uncharacterized, {
    ea_data$vp_unchar_state <- input$vp_uncharacterized
    
    if (input$datasource == "rna" && !ea_data$vp_unchar_state){
      ea_data$gg <- my_volcano(input_data = ea_data$dat_vp %>% dplyr::filter(symbol %in% gene_ids$SYMBOL[gene_ids$characterized]),
                               p_threshold = input$p_thres,
                               fc_threshold = input$fc_thres)
    } else {
      ea_data$gg <- my_volcano(input_data = ea_data$dat_vp,
                               p_threshold = input$p_thres,
                               fc_threshold = input$fc_thres)
    }
    
    # Update Selectize List for shown label in VP
    updateSelectizeInput(
      session,
      "vp_show_labels",
      choices = ea_data$gg[["data"]]$label_hover,
      selected = NULL,
      server = TRUE
    )
  })
  
  # Navigation A' (Metabolite Enrichment) -----
  observeEvent(input$submit_metab_ea, {
    
    # show waiter
    w_metab_ea <- Waiter$new(html = spin_3(), color = transparent(0.5))
    w_metab_ea$show()
    
    # load library for volcano plot
    library(ggiraph)
    
    # ORA -----
    dat_metab_vp <- dat_metab %>%
      filter(tissue == input$tissue_metab_ea) %>%
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
                                                     cid_table = dat_compound_id)
    } else if (input$metab_universe == "all_known") {
      res_metab_ea_list <- my_microbiomeprofiler_ora(differential_metabs = diff_metabs,
                                                     specify_background = FALSE,
                                                     bioterm_database = input$metab_database,
                                                     cid_table = dat_compound_id)
    }
    
    # Navigation A' Rendering -----
    # change reactive value
    ea_metab_data$gg_metab <- my_volcano(input_data = dat_metab_vp,
                                         p_threshold = input$metab_p_thres,
                                         fc_threshold = input$metab_fc_thres,
                                         tooltip_label = "compound")
    ea_metab_data$dat_metab_vp <- dat_metab_vp
    ea_metab_data$dp_metab <- res_metab_ea_list$grob_metab_dot_plot
    ea_metab_data$dat_metab_ea_show <- res_metab_ea_list$df_metab_ea_show
    
    # Output in Navigation A'
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
      filename = "volcano_plot.png",
      content = function(file) {
        ggsave(file, plot = ea_metab_data$gg_metab, width = 6, height = 4.5, dpi = 300)
      }
    )
    
    # Download Button for Result Table
    output$download_metab_de_table <- downloadHandler(
      filename = function() {
        paste0("result_", input$tissue_metab_ea, "_", input$metab_datasource, ".csv")
      },
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
      filename = "enrichment_dot_plot.png",
      content = function(file) {
        ggsave(file, plot = ea_metab_data$dp_metab, width = 9, height = 9, dpi = 300)
      }
    )
    
    # Download Button for Result Table
    output$download_metab_ea_table <- downloadHandler(
      filename = "all_enriched_terms.csv",
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
      rect_fc_show <- dat_prot %>%
        filter(tissue == input$tissue_rect_pv) %>%
        filter(!is.na(ENTREZID)) %>%
        pull(log2FoldChange, name = ENTREZID)
    } else if (input$source_rect_pv == "rna" & input$design_rect_pv == "two_group") {
      rect_fc_show <- my_deseq2_est(meta_data = rna_metadata,
                                    target_tissue = input$tissue_rect_pv,
                                    est_design = "two_group",
                                    rawcount_dir = "data/RNA_raw",
                                    age_1 = input$age1_pv,
                                    age_2 = input$age2_pv) %>%
        left_join(., gene_ids, by = c("symbol" = "SYMBOL")) %>%
        filter(!is.na(ENTREZID)) %>%
        pull(log2FoldChange, name = ENTREZID)
    } else if (input$source_rect_pv == "rna" & input$design_rect_pv == "all_group") {
      rect_fc_show <- readRDS(paste0("data/RNA_log2FC/bulkTMS_logFC_", input$tissue_rect_pv, "_exported_2025_12_17.rds"))
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
