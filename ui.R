# UI -----

ui <- page_navbar(
  title = "Shiny AMMOA (Local Ver 1.0)",
  # Choose appearance theme
  theme = bs_theme(bootswatch = "minty"),
  
  # ---- Navigation A: Protein/RNA Enrichment ----
  nav_panel("RNA/Protein Enrichment",
            id = "panel_ea",
            useWaiter(),
            
            # sidebar
            layout_sidebar(
              sidebar = sidebar(
                
                # Appearance of Sidebar
                open = "always", # Sidebar cannot be closed
                width = 300,
                
                # Action Button
                card(
                  helpText("After you choose your input, click 'Submit' button."),
                  actionButton(inputId = "submit_ea",
                               label = strong("Submit")) # strong() -> bold style
                ),
                
                # Users' Choice -----
                card(
                  card_header(strong("Choose Your Input")),
                  
                  # Data Source Select
                  card(
                    radioButtons(
                      label = strong("Data Source"),
                      inputId = "datasource",
                      choiceNames = list(HTML("bulk RNA-seq<br> (Schaum et al 2020)"),
                                         HTML("Proteomics<br> (Keele et al 2023)")),
                      choiceValues = c("rna", "protein"),
                      selected = "rna",
                      inline = FALSE
                    )
                  ),
                  
                  # Tissue Select
                  # Because available tissue panel differs depends on dataset referring,
                  # interactive update is set in the server below
                  card(
                    selectInput(
                      label = strong("Tissue"),
                      inputId = "tissue_ea",
                      choices = NULL,
                      selected = NULL
                    )
                  ),
                  
                  
                  # DESeq2 Design
                  # the below is initial value (datasource == "rna")
                  # as only two group comparison (8 month vs 18 month) is available for proteomic data,
                  # interactive update is set in the server below
                  card(
                    radioButtons(
                      inputId = "design",
                      label = strong("Design for DEG/DEP Detection"),
                      choices = c("linear model using all age groups" = "linear",
                                  "two age groups comparison" = "two_group"),
                      selected = "two_group",
                      inline = TRUE
                    ),
                    
                    # age groups to be compared (only applicable when design == two_group)
                    conditionalPanel(
                      condition = "input.design == 'two_group'",
                      p(strong("Age Group (month)")),
                      fluidRow(
                        column(6, selectInput("age1", "Young", choices = NULL)),
                        column(6, selectInput("age2", "Old", choices = NULL))
                      )
                    )
                  ),
                  
                  # Threshold to Define DEGs/DEPs
                  card(
                    p(strong("DE Threshold")),
                    
                    numericInput(inputId = "p_thres",
                                 label = "FDR (BH-adjusted p)",
                                 value = 0.05,
                                 min = 0,
                                 max = 1),
                    numericInput(inputId = "fc_thres",
                                 label = "Log2 Fold Change",
                                 value = 0.5,
                                 min = 0)
                  ),
                  
                  # DEG usage for Enrichment Analysis
                  card(
                    radioButtons(
                      inputId = "deg_usage",
                      label = strong("DEs Used for Enrichment Analysis"),
                      choices = c("both", "up-regulated", "down-regulated"),
                      selected = "both"
                    ),
                    
                    # Database
                    radioButtons(
                      inputId = "database",
                      label = strong("Pathway Database"),
                      choices = c("KEGG", "GO(BP)", "GO(MF)", "GO(CC)"),
                      selected = "KEGG"
                    )
                  )
                )
              ),
              
              # Main Panel -----
              navset_card_pill(
                # Tab A1: Volcano Plot
                nav_panel(
                  title = "Volcano Plot",
                  
                  div(
                    style = "display:flex; flex-direction:column; height:90vh;",
                    # volcano plot
                    div(style = "flex:4; display:flex; gap:20px;",
                        # left：Volcano plot
                        div(
                          style = "flex:3;",
                          girafeOutput("graph_vp", height = "100%", width = "100%")),
                        # right：Selectize Input
                        div(style = "flex:1;",
                            selectizeInput(
                              inputId = "vp_show_labels",
                              label = "Labels to Show:",
                              choices = NULL,
                              multiple = TRUE
                            ),
                            ## show label action button
                            div(
                              style = "margin-top: 10px;",
                              uiOutput("vp_update_button")
                            ),
                            
                            # show uncharacterized gene checkbox
                            div(
                              style = "margin-top: 10px;",
                              uiOutput("vp_uncharacterized_check")
                            )
                        )
                    ),
                    # download button for plot and table
                    div(
                      style = "display: flex; justify-content: left; gap: 15px; margin-bottom: 5px;",
                      downloadButton("download_vp", "Download Plot"),
                      downloadButton("download_de_table", "Download Table")
                    ),
                    # table showing gene name, fold change, etc.
                    div(style = "flex:5; overflow-y:auto;",
                        p(strong("Results Table")),
                        reactableOutput("table_vp", height = "100%")),
                  )
                ),
                
                # Tab A2: Enrichment Analysis
                nav_panel(
                  title = "Pathway Enrichment",
                  div(
                    style = "display:flex; flex-direction:column; height:120vh;",
                    # dot plot showing entiched biological terms
                    div(style = "flex:5;", plotOutput("graph_dp", height = "100%")),
                    # download button for plot and table
                    div(
                      style = "display: flex; justify-content: left; gap: 15px; margin-bottom: 5px;",
                      downloadButton("download_dp", "Download Plot"),
                      downloadButton("download_ea_table", "Download Table")
                    ),
                    # table showing all enriched terms
                    div(style = "flex:3; overflow-y:auto;",
                        p(strong("All Enriched Biological Terms")),
                        reactableOutput("table_ea", height = "100%"))
                  )
                )
              )
            )
  ),
  
  # ---- Navigation A' Metabolome Enrichment ----
  nav_panel("Metabolite Enrichment",
            id = "panel_metab_ea",
            useWaiter(),
            
            # sidebar
            layout_sidebar(
              sidebar = sidebar(
                
                # Appearance of Sidebar
                open = "always", # Sidebar cannot be closed
                width = 300,
                
                # Action Button
                card(
                  helpText("After you choose your input, click 'Submit' button."),
                  actionButton(inputId = "submit_metab_ea",
                               label = strong("Submit")) # strong() -> bold style
                ),
                
                # Users' Choice -----
                card(
                  card_header(strong("Choose Your Input")),
                  
                  # Data Source Select
                  # Actually it doesn't have choice, but to make UI look similar to RNA/Protein Enrichment, use Radio Button
                  card(
                    radioButtons(
                      label = strong("Data Source"),
                      inputId = "metab_datasource",
                      choiceNames = list(HTML("Metabolites<br> (Jankowski et al 2025)")),
                      choiceValues = "metab",
                      selected = "metab",
                      inline = FALSE
                    )
                  ),
                  
                  # Tissue Select
                  card(
                    selectInput(
                      label = strong("Tissue"),
                      inputId = "tissue_metab_ea",
                      choices = c("Brain", "Colon", "Diaphragm", "Eye", "Heart", "Jejunum", "Kidney", "Liver",
                                  "Lung", "Pancreas", "Quadriceps", "Serum", "Skin(Ear)", "Soleus", "Spleen"),
                      selected = "Kidney"
                    )
                  ),
                  
                  # Comparison Design
                  # Like datasource input, this section doesn't have multiple choices
                  card(
                    p(strong("Age Group")),
                    p("Young: 15-19 weeks"),
                    p("Old: 90-99 weeks")
                  ),
                  
                  # Threshold to Define DEGs/DEPs
                  card(
                    p(strong("DE Threshold")),
                    
                    numericInput(inputId = "metab_p_thres",
                                 label = "FDR (BH-adjusted p)",
                                 value = 0.05,
                                 min = 0,
                                 max = 1),
                    numericInput(inputId = "metab_fc_thres",
                                 label = "Log2 Fold Change",
                                 value = 0,
                                 min = 0)
                  ),
                  
                  # DEG usage for Enrichment Analysis
                  card(
                    radioButtons(
                      inputId = "metab_de_usage",
                      label = strong("DEs Used for Enrichment Analysis"),
                      choices = c("both", "up-regulated", "down-regulated"),
                      selected = "both"
                    ),
                    
                    # Database
                    radioButtons(
                      inputId = "metab_database",
                      label = strong("Pathway Database"),
                      choices = c("KEGG", "SMPDB"),
                      selected = "KEGG"
                    ),
                    
                    # How to set universe (i.e., background metabolites)
                    radioButtons(
                      inputId = "metab_universe",
                      label = strong("Background Metabolites"),
                      choiceNames = c("Measured metabolites only",
                                      "All metabolites registered in DB"),
                      choiceValues = c("measured", "all_known"),
                      selected = "measured"
                    ),
                    
                    helpText(
                      "* 'All metabolites registered in DB' uses all KEGG/SMPDB metabolites; ",
                      "it may detect more pathways but can increase false positives."
                    )
                  )
                )
              ),
              
              # Main Panel -----
              navset_card_pill(
                # Tab A1: Volcano Plot
                nav_panel(
                  title = "Metabolite Volcano Plot",
                  
                  div(
                    style = "display:flex; flex-direction:column; height:90vh;",
                    p(strong("Volcano Plot")),
                    # volcano plot
                    div(style = "flex:4; display:flex; gap:20px;",
                        # left：Volcano plot
                        div(
                          style = "flex:3;",
                          girafeOutput("graph_metab_vp", height = "100%", width = "100%")),
                        # right：Selectize Input
                        div(style = "flex:1;",
                            selectizeInput(
                              inputId = "vp_metab_show_labels",
                              label = "Labels to Show:",
                              choices = NULL,
                              multiple = TRUE
                            ),
                            # show label action button placeholder
                            uiOutput("vp_metab_update_button")
                        )
                    ),
                    div(
                      style = "display: flex; justify-content: left; gap: 15px; margin-bottom: 5px;",
                      downloadButton("download_metab_vp", "Download Plot"),
                      downloadButton("download_metab_de_table", "Download Table")
                    ),
                    # table showing gene name, fold change, etc.
                    div(style = "flex:5; overflow-y:auto;",
                        p(strong("Results Table")),
                        reactableOutput("table_metab_vp", height = "100%")),
                  )
                ),
                
                # Tab A2: Enrichment Analysis
                nav_panel(
                  title = "Metabolite Pathway Enrichment",
                  div(
                    style = "display:flex; flex-direction:column; height:120vh;",
                    # dot plot showing entiched biological terms
                    div(style = "flex:5;", plotOutput("graph_metab_dp", height = "100%")),
                    # download button for plot and table
                    div(
                      style = "display: flex; justify-content: left; gap: 15px; margin-bottom: 5px;",
                      downloadButton("download_metab_dp", "Download Plot"),
                      downloadButton("download_metab_ea_table", "Download Table")
                    ),
                    # table showing all enriched terms
                    div(style = "flex:3; overflow-y:auto;",
                        p(strong("All Biological Terms")),
                        reactableOutput("table_metab_ea", height = "100%"))
                  )
                )
              )
            )
  ),
  
  # ---- Navigation B ----
  nav_panel("Pathway Mapping",
            id = "panel_pv",
            useWaiter(),
            
            # Sidebar ----
            layout_sidebar(
              sidebar = sidebar(
                
                # Appearance of Sidebar
                open = "always", # Sidebar cannot be closed
                width = 300,
                
                # Action Button
                card(
                  helpText("After you choose your input, click 'Submit' button."),
                  
                  actionButton(inputId = "submit_pv",
                               label = strong("Submit"))
                ),
                
                # parameter
                # tissue
                card(
                  p(strong("Sourec of Gene/Protein Nodes")),
                  
                  radioButtons(
                    label = strong("Data Type"),
                    inputId = "source_rect_pv",
                    choiceNames = list(HTML("bulk RNA-seq<br> (Schaum et al 2020)"),
                                       HTML("Proteomics<br> (Keele et al 2023)"),
                                       HTML("Not Show")),
                    choiceValues = list("rna", "protein", "not_show"),
                    selected = "rna",
                    inline = FALSE
                  ),
                  
                  conditionalPanel(
                    condition = "input.source_rect_pv != 'not_show'",
                    selectInput(
                      label = strong("Tissue"),
                      inputId = "tissue_rect_pv",
                      choices = NULL,
                      selected = NULL
                    )
                  ),
                  
                  conditionalPanel(
                    condition = "input.source_rect_pv == 'rna'",
                    radioButtons(
                      inputId = "design_rect_pv",
                      label = strong("Age Groups Used"),
                      choices = c("all age groups" = "all_group",
                                  "two age groups"  = "two_group"),
                      selected = "two_group",
                      inline = TRUE
                    ),
                    helpText("When you choose 'all age groups', compound nodes cannot be shown.")
                  ),
                  
                  conditionalPanel(
                    condition = "input.source_rect_pv == 'rna' && input.design_rect_pv == 'two_group'",
                    p(strong("Age Group (month)")),
                    fluidRow(
                      column(6, selectInput("age1_pv", "Young(Ref)", choices = NULL)),
                      column(6, selectInput("age2_pv", "Old", choices = NULL))
                    ),
                  )
                ),
                
                card(
                  p(strong("Source of Metabolite Nodes")),
                  radioButtons(
                    label = strong("Data Type"),
                    inputId = "source_circ_pv",
                    choiceNames = list(HTML("Metabolites<br> (Jankowski et al 2025)"),
                                       HTML("Not Show")),
                    choiceValues = list("metab", "not_show"),
                    selected = "metab",
                    inline = FALSE
                  ),
                  
                  conditionalPanel(
                    condition = "input.source_circ_pv != 'not_show'",
                    selectInput(
                      label = strong("Tissue"),
                      inputId = "tissue_circ_pv",
                      choices = c("Brain", "Colon", "Diaphragm", "Eye", "Heart", "Jejunum", "Kidney", "Liver",
                                  "Lung", "Pancreas", "Quadriceps", "Serum", "Skin(Ear)", "Soleus", "Spleen"),
                      selected = "Kidney"
                    )
                  ),
                ),
                
                card(
                  # kegg pathway
                  selectInput(
                    label = strong("KEGG Pathway"),
                    inputId = "kegg_id_pv",
                    choices = kegg_ids$label,
                    selected = "mmu00020: Citrate cycle (TCA cycle)"
                  ),
                  # CSS
                  tags$style(HTML(paste(
                    "/* current selection */",
                    "#kegg_pathway + .selectize-control .selectize-input { font-size: 14px; }",
                    "/* drop-down selection list */",
                    "#kegg_pathway + .selectize-control .selectize-dropdown .option { font-size: 14px; }",
                    sep = "\n"
                  ))),
                  p(tags$span(
                    style = "font-size:12px;",
                    "List of KEGG pathways are available at ",
                    tags$a(href = "https://www.kegg.jp/kegg/pathway.html", 
                           "https://www.kegg.jp/kegg/pathway.html", 
                           target = "_blank")
                  )
                  )
                )
              ),
              
              # Navigation B Main Panel -----
              card(
                p(strong("Mapped Pathway")),
                # Legend
                div(style = "height: 150px;", imageOutput("pathview_legend", height="100%", width="auto")),
                # message
                div(style = "height: 25px;", uiOutput("pathview_message")),
                div(style = "height: 25px;font-size:12px;",
                    p("Pathway diagrams are based on KEGG (Kyoto Encyclopedia of Genes and Genomes). © Kanehisa Laboratories. All rights reserved.")),
                # download button
                div(
                  style = "display: flex; justify-content: left; gap: 15px; margin-bottom: 5px;",
                  downloadButton("download_pv_legend", "Download Legend"),
                  downloadButton("download_pv_plot", "Download Image")
                ),
                # pathview image
                div(style="flex:1; display:flex;", imageOutput("pB", height="100%", width="auto"))
              )
            )
  ),
  
  # ---- Navigation C (User Guide)----
  nav_panel("Userguide",
            id = "panel_userguide",
            navset_card_pill(
              nav_panel(
                title = "Pathway Analysis",
                div(
                  style = "width: 200px; margin-bottom: 10px;",
                  downloadButton("download_userguide_analysis", "Download as PDF")
                ),
                div(
                  style = "border: 1px solid #ddd;border-radius: 6px;padding: 8px;background-color: #fafafa;",
                  tags$iframe(
                    src = "userguide_analysis.html",
                    style = "width:100%; height: 85vh; border:none;"
                  )
                )
              ),
              nav_panel(
                title = "Source Data Description",
                div(
                  style = "width: 200px; margin-bottom: 10px;",
                  downloadButton("download_userguide_sourcedata", "Download as PDF")
                ),
                div(
                  style = "border: 1px solid #ddd;border-radius: 6px;padding: 8px;background-color: #fafafa;",
                  tags$iframe(
                    src = "userguide_sourcedata.html",
                    style = "width:100%; height: 85vh; border:none;"
                  )
                )
              )
            )
  ),
  
  # ---- Navigation D ----
  nav_panel("About",
            id = "panel_about",
            page_fillable(
              card(
                h2("About this App"),
                h3("Shiny AMMOA (Shiny Aging Murine Multi-Omic Analyzer)"),
                h4("Local Version (Full-Featured)"),
                
                h3("Version History"),
                tags$ul(
                  tags$li(
                    c("2026-04-24 Ver. 1.0 (First-Launched Build)")
                  )
                ),
                
                h3("Please Site"),
                tags$ul(
                  tags$li(
                    c("Coming Soon")
                  )
                ),
                
                h3("Source Code"),
                tags$ul(
                  tags$li(
                    a("GitHub Repository",
                      href = "https://github.com/M-Ninomiya-Kanda/Shiny_AMMOA",
                      target = "_blank")
                  )
                ),
                
                h3("Author"),
                p("Mayuka Ninomiya Kanda (QB3, University of California, Berkeley)"),
                p("Contact: ", a("m.k.ninomiya@berkeley.edu", href = "mailto:m.k.ninomiya@berkeley.edu")),
              )
            )
  )
  
)
