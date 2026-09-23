## Function to detect if KEGG pathway .xml file has gene rectangle / compound circle -----
# Define function to define whether indicated KEGG pathway has gene to visualize
# function to detect if there is gene in rectangle or not
my_has_gene_rect <- function(pathway_id, species = "mmu") {
  # path for KGML file
  kgml_file <- paste0(pathway_id, ".xml")
  
  # if the file does not exist, return FALSE
  if (!file.exists(kgml_file)) {
    warning(paste("KGML file not found:", kgml_file))
    return(FALSE)
  }
  
  # load XML
  kgml <- tryCatch({
    xmlParse(kgml_file)
  }, error = function(e) return(NULL))
  
  if (is.null(kgml)) return(FALSE)
  
  # search entry nodes
  entries <- getNodeSet(kgml, "//entry[@type='gene']")
  
  # detect rectangle nodes
  for (entry in entries) {
    g_nodes <- getNodeSet(entry, "graphics[@type='rectangle']")
    if (length(g_nodes) > 0) {
      return(TRUE)
    }
  }
  
  return(FALSE)
}

my_has_cpd_circ <- function(pathway_id, species = "mmu") {
  # path for KGML file
  kgml_file <- paste0(pathway_id, ".xml")
  
  # if the file does not exist, return FALSE
  if (!file.exists(kgml_file)) {
    warning(paste("KGML file not found:", kgml_file))
    return(FALSE)
  }
  
  # load XML
  kgml <- tryCatch({
    xmlParse(kgml_file)
  }, error = function(e) return(NULL))
  
  if (is.null(kgml)) return(FALSE)
  
  # search entry nodes
  entries <- getNodeSet(kgml, "//entry[@type='compound']")
  
  # detect rectangle nodes
  for (entry in entries) {
    g_nodes <- getNodeSet(entry, "graphics[@type='circle']")
    if (length(g_nodes) > 0) {
      return(TRUE)
    }
  }
  
  return(FALSE)
}

# Load DESeq2 Results ----
my_load_degs <- function(tissue,
                         group_design,
                         age_young = NULL,
                         age_old = NULL) {
  
  dat_ids <- readRDS("data/Transcriptome/gene_symbol_entrezid_table_2026-09-12_exported.rds")
  
  if (group_design == "linear"){
    de_filename <- paste0("data/Transcriptome/DE_res_", group_design,
                          "/DESeq_", tissue, "_linear_2026-09-12_exported.rds")
  } else if (group_design == "two_group"){
    de_filename <- paste0("data/Transcriptome/DE_res_", group_design,
                          "/", tissue,
                          "/DESeq_", tissue, "_m", age_young, "_vs_m", age_old, "_2026-09-12_exported.rds")
  }
  
  de_res <- readRDS(de_filename) %>%
    mutate(ENTREZID = as.character(ENTREZID)) %>%
    left_join(dat_ids, by = "ENTREZID") %>%
    dplyr::select(SYMBOL, everything())
  
  return(de_res)
}


# clusterProfiler Over Represenatation Analysis ----
my_cluterprofiler_ora <- function(differential_exps,
                                  background_exps,
                                  bioterm_database,
                                  tissue_plot_title = NULL){
  requireNamespace("clusterProfiler")
  requireNamespace("enrichplot")
  requireNamespace("org.Mm.eg.db")
  
  res_ora <- NULL
  df_ea_show <- NA
  dot_plot <- NULL
  
  # Perform ORA
  if (length(differential_exps) == 0) {
    # In case no DEs were extracted
    dot_plot <- ggplot() +
      annotate("text", x = 0.5, y = 0.5,
               label = paste0("No DEs were detected with the current thresholds."),
               size = 6, hjust = 0.5, color = "red3") +
      theme_void() +
      theme(text = element_text(family = "sans"))
    df_ea_show <- data.frame()
    
  } else if (length(differential_exps) > 0){
    if (bioterm_database == "KEGG"){
      res_ora <- clusterProfiler::enrichKEGG(gene = unique(differential_exps),
                                             universe = unique(background_exps),
                                             organism = "mmu",
                                             pvalueCutoff = 0.05,
                                             pAdjustMethod = "BH")
      if (!is.null(res_ora)){
        res_ora <- setReadable(res_ora, OrgDb = org.Mm.eg.db::org.Mm.eg.db, keyType = "ENTREZID")
      }
    } else {
      res_ora <- clusterProfiler::enrichGO(gene = unique(differential_exps),
                                           universe = unique(background_exps),
                                           OrgDb = org.Mm.eg.db::org.Mm.eg.db,
                                           ont = str_sub(bioterm_database, 4, -2),
                                           pAdjustMethod = "BH",
                                           pvalueCutoff  = 0.05,
                                           readable = TRUE)
    }
    
    
    # Plot
    if (is.null(res_ora) || (nrow(res_ora@result) > 0 && sum(res_ora@result$p.adjust < 0.05) == 0)) {
      # in case no enriched term was found 
      dot_plot <- ggplot() +
        annotate("text", x = 0.5, y = 0.5,
                 label = paste0("No significantly enriched biological terms were found.",
                                "\nConsider loosing the DE thresholds."),
                 size = 6, hjust = 0.5, color = "red3") +
        theme_void() +
        theme(text = element_text(family = "sans"))
      df_ea_show <- data.frame()
    } else if (!is.null(res_ora) & sum(res_ora@result$p.adjust < 0.05) > 0){
      # if <10 enrichment term found
      if (sum(res_ora@result$p.adjust < 0.05) <= 10){
        df_ea_show <- res_ora@result %>%
          filter(p.adjust < 0.05)
        
        dot_plot <- df_ea_show %>%
          arrange(p.adjust) %>%
          # wrap biological term because some of them are too long to display
          mutate(Description = str_wrap(Description, width = 60, indent = 0, exdent = 0)) %>%
          mutate(Description = factor(Description, levels = rev(.$Description))) %>%
          ggplot(aes(x = Description, y = -log10(p.adjust), size = Count)) +
          geom_point(color = brewer.pal(3, "Set2")[2]) +
          coord_flip() +
          scale_size_continuous(name = "DEs Count",
                                labels = scales::number_format(accuracy = 1),
                                limits = c(1, NA),
                                range = c(3, 8)) +
          xlab(NULL) +
          ylab("-log10(p.adjust)") +
          labs(title = paste0(bioterm_database, " ORA of ", tissue_plot_title)) +
          theme_light() +
          theme(text = element_text(family = "sans", size = 16),
                axis.text.y = element_text(color = "black", size = 11),
                axis.line = element_line(color = "black"))
      } else {
        # make the cluster of enriched terms
        # get distance matrix
        res_ora@result <- res_ora@result %>% drop_na(ID)
        dist_mat <- stats::as.dist(1 - enrichplot::pairwise_termsim(res_ora, showCategory = sum(res_ora@result$p.adjust < 0.05))@termsim)
        
        # hierarchical clustering
        dot_plot <- stats::cutree(hclust(dist_mat, method = "ward.D2"), k = 5) %>%
          # convert cluster information to data.frame
          data.frame() %>%
          setNames("cluster") %>%
          rownames_to_column("Description") %>%
          mutate(cluster = factor(cluster)) %>%
          left_join(res_ora@result, by = "Description")
        
        df_ea_show <- dot_plot
        
        # if more than 25 pathways were found, show top5 terms from each cluster
        
        if (nrow(dot_plot) > 25){
          dot_plot <- dot_plot %>%
            group_by(cluster) %>%
            arrange(p.adjust, .by_group = TRUE) %>%
            dplyr::slice_head(n = 5) %>%
            ungroup() 
        }
        
        dot_plot <- dot_plot %>%
          arrange(cluster, p.adjust) %>%
          # wrap biological term because some of them are too long to display
          mutate(Description = str_wrap(Description, width = 60, indent = 0, exdent = 0)) %>%
          mutate(Description = factor(Description, levels = rev(.$Description))) %>%
          # dot plot
          ggplot(aes(x = Description, y = -log10(pvalue), size = Count, color = cluster)) +
          geom_point() +
          coord_flip() +
          scale_color_brewer(palette = "Set2") +
          scale_size_continuous(name = "DE Count",
                                labels = scales::number_format(accuracy = 1),
                                limits = c(1, NA),
                                range = c(3, 8)) +
          xlab(NULL) +
          ylab("-log10(p.adjust)") +
          labs(title = paste0(bioterm_database, " ORA of ", tissue_plot_title)) +
          theme_light() +
          theme(text = element_text(family = "sans", size = 16),
                axis.text = element_text(color = "black"),
                axis.text.y = element_text(color = "black", size = 11),
                axis.line = element_line(color = "black"))
      }
    }
  }
  
  out_res <- list()
  out_res$df_ea_show <- df_ea_show
  out_res$grob_dot_plot <- dot_plot
  
  gc();gc()
  return(out_res)
}

# GSEA Visualization ----
# As GSEA results are pre-computed, only visualize on-demand
my_gsea_dp <- function(gsea_res,
                       bioterm_database,
                       tissue_plot_title = NULL){
 
  if (sum(gsea_res$p.adjust < 0.05) <= 10){
    grob <- gsea_res %>%
      filter(p.adjust < 0.05) %>%
      mutate(signed_log10_pavl = -sign(NES)*log10(p.adjust)) %>%
      arrange(desc(signed_log10_pavl)) %>%
      # wrap biological term because some of them are too long to display
      mutate(Description = str_wrap(Description, width = 60, indent = 0, exdent = 0)) %>%
      mutate(Description = factor(Description, levels = rev(.$Description))) %>%
      ggplot(aes(x = Description, y = signed_log10_pavl, size = abs(NES))) +
      geom_hline(yintercept = 0, color = "grey20", linetype = "dashed") +
      geom_point(color = brewer.pal(3, "Set2")[2]) +
      coord_flip() +
      scale_size_continuous(name = "|NES|",
                            range = c(3, 8)) +
      xlab(NULL) +
      ylab("signed log10(p.adjust)") +
      labs(title = paste0(bioterm_database, " GSEA of ", tissue_plot_title)) +
      theme_light() +
      theme(text = element_text(family = "sans", size = 16),
            axis.text.y = element_text(color = "black", size = 11),
            axis.line = element_line(color = "black"))
    
  } else {
    
    # if more than 25 pathways were found, show top5 terms from each cluster
    
    if (sum(gsea_res$p.adjust < 0.05) > 25){
      grob <- gsea_res %>%
        group_by(cluster) %>%
        arrange(p.adjust, .by_group = TRUE) %>%
        dplyr::slice_head(n = 5) %>%
        ungroup() 
    } else {
      grob <- gsea_res
    }
    
    grob <- grob %>%
      filter(p.adjust < 0.05) %>%
      mutate(signed_log10_pavl = -sign(NES)*log10(p.adjust)) %>%
      arrange(cluster, desc(signed_log10_pavl)) %>%
      # wrap biological term because some of them are too long to display
      mutate(Description = str_wrap(Description, width = 60, indent = 0, exdent = 0)) %>%
      mutate(Description = factor(Description, levels = rev(.$Description))) %>%
      ggplot(aes(x = Description, y = signed_log10_pavl, size = abs(NES))) +
      geom_hline(yintercept = 0, color = "grey20", linetype = "dashed") +
      geom_point(aes(color = cluster)) +
      scale_color_brewer(palette = "Set2") +
      coord_flip() +
      scale_size_continuous(name = "|NES|",
                            range = c(3, 8)) +
      xlab(NULL) +
      ylab("signed log10(p.adjust)") +
      labs(title = paste0(bioterm_database, " GSEA of ", tissue_plot_title)) +
      theme_light() +
      theme(text = element_text(family = "sans", size = 16),
            axis.text.y = element_text(color = "black", size = 11),
            axis.line = element_line(color = "black"))
  }
  
  return(grob)
}

# volcano plot ----
my_volcano <- function(input_data, p_threshold, fc_threshold, tooltip_label = "SYMBOL"){
  grob <- input_data %>%
    # remove missing value
    drop_na(log2FoldChange, padj) %>%
    # color label
    mutate(change = case_when(padj < p_threshold & log2FoldChange > fc_threshold ~ "up-regulated",
                              padj < p_threshold & log2FoldChange < -fc_threshold ~ "down-regulated",
                              .default = "others")) %>%
    mutate(change = factor(change, levels = c("up-regulated", "down-regulated", "others"))) %>%
    # tooltip label to show
    dplyr::rename(label_hover = all_of(tooltip_label)) %>%
    mutate(tooltip_id = row_number()) %>%
    ggplot() +
    geom_vline(xintercept = fc_threshold, linetype = "dashed", color = "grey20") +
    geom_vline(xintercept = -fc_threshold, linetype = "dashed", color = "grey20") +
    geom_hline(yintercept = -log10(p_threshold), linetype = "dashed", color = "grey20") +
    ggiraph::geom_point_interactive(aes(x = log2FoldChange, y = -log10(padj), color = change,
                                        tooltip = label_hover, data_id = tooltip_id),
                                    hover_nearest = TRUE,
                                    hover_css = "r:2pt;stroke:black;stroke-width:1px;fill:magenta;",
                                    alpha = 0.8,
                                    size = 1) +
    scale_color_manual(name = "Change with Age",
                       values = c("up-regulated" = "orange",
                                  "down-regulated" = "royalblue",
                                  "others" = "grey30")) +
    xlim(c(-max(abs(input_data$log2FoldChange[!is.na(input_data$log2FoldChange) & !is.na(input_data$padj)])),
           max(abs(input_data$log2FoldChange[!is.na(input_data$log2FoldChange) & !is.na(input_data$padj)])))) +
    xlab("covariate-adjusted log2FC") +
    ylab("-Log10 FDR") +
    theme_classic() +
    theme(text = element_text(family = "sans", size = 14),
          axis.text = element_text(color = "black"))
  
  return(grob)
}


# Annotate Compound Name to enrichResult object -----
# As enrichKEGG() and enrichHMDB() cannot convert KDGG/HMDB compound ID to common name, I make in-house function
my_metab_id_convert <- function(id_string, keytype, cid_table){
  # split "/" separated ID list
  tmp1 <- str_split(id_string, "/", simplify = TRUE) %>%
    as.vector() %>%
    data.frame(ID = .)
  
  if(keytype == "KEGG"){
    col_pick <- c("KEGG", "kegg_name")
  } else if (keytype == "HMDB"){
    col_pick <- c("HMDB", "compound")
  }
  
  tmp2 <- cid_table %>%
    dplyr::select(all_of(col_pick)) %>%
    setNames(c("ID", "name")) %>%
    left_join(tmp1, ., by = "ID") %>%
    pull(name) %>%
    str_flatten(collapse = "/")
  
  return(tmp2)
}

my_metab_setReadable <- function(res_object, cid_table){
  # detect which compound ID is used in the enrichment analysis result object
  if(str_sub(res_object@result$ID[1], 1, 3) == "SMP"){
    keytype = "HMDB"
  } else if (str_sub(res_object@result$ID[1], 1, 3) == "map"){
    keytype = "KEGG"
  }
  
  out_res_object <- res_object
  
  if (is(res_object) == "enrichResult"){
    tmp <- out_res_object@result %>%
      mutate(geneID = map_chr(geneID, ~ my_metab_id_convert(.x, keytype = keytype, cid_table = cid_table))) %>%
      # change column name as it is the same name as gene enrichment result(e.g., "GeneRatio") and confusing
      dplyr::rename(MetabRatio = GeneRatio,
                    CompoundName = geneID)
  } else if (is(res_object) == "gseaResult"){
    tmp <- out_res_object@result %>%
      mutate(core_enrichment = map_chr(core_enrichment, ~ my_metab_id_convert(.x, keytype = keytype, cid_table = cid_table)))
  }
  
  out_res_object@result <- tmp
  
  return(out_res_object)
}


# Metabolite Over Representation Analysis with MicrobiomeProfiler----

my_microbiomeprofiler_ora <- function(differential_metabs,
                                      specify_background = TRUE,
                                      background_metabs = NULL,
                                      bioterm_database = c("KEGG", "SMPDB"),
                                      cid_table,
                                      tissue_plot_title = NULL){
  requireNamespace("clusterProfiler")
  requireNamespace("MicrobiomeProfiler")
  
  res_metab_ora <- NULL
  df_metab_ea_show <- NULL
  metab_dot_plot <- NULL
  
  # Perform ORA
  if (length(differential_metabs) == 0) {
    # In case no DEs were extracted
    metab_dot_plot <- ggplot() +
      annotate("text", x = 0.5, y = 0.5,
               label = paste0("No differential metabolites were detected with the current thresholds."),
               size = 6, hjust = 0.5, color = "red3") +
      theme_void()
    df_metab_ea_show <- data.frame()
    
  } else if (length(differential_metabs) > 0){
    if (bioterm_database == "KEGG" & specify_background){
      # KEGG Enrichment & Specify background metabolites
      res_metab_ora <- clusterProfiler::enricher(gene = unique(differential_metabs),
                                                 universe = unique(background_metabs),
                                                 TERM2GENE = dat_pathway2compound,
                                                 TERM2NAME = dat_pathway2name,
                                                 pvalueCutoff = 1,
                                                 pAdjustMethod = "BH",
                                                 minGSSize = 5)
    } else if (bioterm_database == "KEGG" & !specify_background) {
      # KEGG Enrichment & NOT Specify background metabolites
      res_metab_ora <- clusterProfiler::enricher(gene = unique(differential_metabs),
                                                 universe = unique(dat_pathway2compound$gene),
                                                 TERM2GENE = dat_pathway2compound,
                                                 TERM2NAME = dat_pathway2name,
                                                 pvalueCutoff = 1,
                                                 pAdjustMethod = "BH",
                                                 minGSSize = 5)
    } else if (bioterm_database == "SMPDB" & specify_background){
      # SMPDB Enrichment & Specify background metabolites
      res_metab_ora <- MicrobiomeProfiler::enrichHMDB(metabo_list = unique(differential_metabs),
                                                      universe = unique(background_metabs),
                                                      pvalueCutoff = 1,
                                                      pAdjustMethod = "BH",
                                                      minGSSize = 5)
    } else if (bioterm_database == "SMPDB" & !specify_background){
      # SMPDB Enrichment & NOT Specify background metabolites
      res_metab_ora <- MicrobiomeProfiler::enrichHMDB(metabo_list = unique(differential_metabs),
                                                      pvalueCutoff = 0.05,
                                                      pAdjustMethod = "BH",
                                                      minGSSize = 5)
    }
    
    # if at least one pathway is discovered, convert compound ID to common name
    if (!is.null(res_metab_ora)){
      res_metab_ora <- my_metab_setReadable(res_object = res_metab_ora,
                                            cid_table = cid_table)
    }
    
    # Plot
    if (is.null(res_metab_ora) ||
        nrow(res_metab_ora@result) == 0 ||
        !any(res_metab_ora@result$Count > 0, na.rm = TRUE))  {
      # in case no enriched term was found 
      metab_dot_plot <- ggplot() +
        annotate("text", x = 0.5, y = 0.5,
                 label = paste0("No biological pathways were found.",
                                "\nConsider loosing the DE thresholds."),
                 size = 6, hjust = 0.5, color = "red3") +
        theme_void() +
        theme(text = element_text(family = "sans"))
      df_metab_ea_show <- data.frame()
    } else {
      df_metab_ea_show <- res_metab_ora@result %>% filter(Count > 0)
      
      # if more than 25 pathways are discovered, show top 25 pathways
      if (nrow(df_metab_ea_show) > 25){
        metab_dot_plot <- df_metab_ea_show %>%
          arrange(p.adjust) %>%
          dplyr::slice_head(n = 25)
        
        g_title <- paste0("Top 25 ", bioterm_database, " Pathways in ", tissue_plot_title)
      } else {
        metab_dot_plot <- df_metab_ea_show
        g_title <- paste0(bioterm_database, " Pathways in", tissue_plot_title)
      }
      
      metab_dot_plot <- metab_dot_plot %>%
        arrange(p.adjust) %>%
        # wrap biological term because some of them are too long to display
        mutate(Description = str_wrap(Description, width = 60, indent = 0, exdent = 0)) %>%
        mutate(Description = factor(Description, levels = rev(.$Description))) %>%
        ggplot(aes(x = Description, y = -log10(p.adjust), size = Count)) +
        geom_point(aes(fill = p.adjust), shape = 21, color = "grey20") +
        coord_flip() +
        scale_size_continuous(name = "Metabolites Count",
                              labels = scales::number_format(accuracy = 1),
                              limits = c(1, NA),
                              range = c(3, 8)) +
        xlab(NULL) +
        ylab("-p.adjust") +
        labs(title = g_title) +
        theme_light(base_family = "sans") +
        theme(text = element_text(size = 16),
              axis.text = element_text(color = "black"),
              axis.text.y = element_text(color = "black", size = 11),
              axis.line = element_line(color = "black"))
      
      # if non-significant pathway(s) are picked-up to show, show padj=.05 vertical line
      if(max(metab_dot_plot$data$p.adjust > 0.05)){
        metab_dot_plot <- metab_dot_plot +
          scale_fill_viridis(direction = -1, limit = c(min(0.05, min(metab_dot_plot$data$p.adjust)), NA)) +
          geom_hline(yintercept = -log10(0.05), color = "red3", linetype = "dashed") +
          annotate(
            "text",
            x = 0.5,
            y = -log10(0.05),
            label = "p.adjust = 0.05",
            vjust = -1,
            hjust = 1,
            color = "red3",
            size = 4
          )
      } else {
        metab_dot_plot <- metab_dot_plot +
          scale_fill_viridis(direction = -1)
      }
    }
  }
  
  out_res <- list()
  out_res$res_metab_ora <- res_metab_ora
  out_res$df_metab_ea_show <- df_metab_ea_show
  out_res$grob_metab_dot_plot <- metab_dot_plot
  
  return(out_res)
}

# cutomized pathview ----
my_pathview <- function(gene_fc = NULL, cpd_fc = NULL, kegg_id, gene_limit = 1){
  requireNamespace("pathview")
  requireNamespace("org.Mm.eg.db")
  # data(bods, package = "pathview")
  
  # initialize the object to put results
  out_res <- list()
  out_res$pathview_filename <- NULL
  
  # pattern1: rectangle = not_show, circle = not_show
  if(all(is.null(gene_fc), is.null(cpd_fc))) {
    download.kegg(pathway.id = substring(kegg_id, 4, 8), species = "mmu", file.type=c("xml", "png"))
    out_res$pathview_filename <- paste0(kegg_id, ".png")
    
  } else {
    # pattern2: at least one of gene or cpd is colored, single color per node
    pv <- pathview::pathview(
      gene.data = gene_fc,
      cpd.data = cpd_fc,
      pathway.id = kegg_id,
      species = "mmu",
      kegg.native = TRUE,
      same.layer = FALSE,
      node.sum = "mean",
      # na.col = "#555555",
      limit = list(gene = gene_limit, cpd = 0.5),
      low  = list(gene = "#4477AA", cpd = "#542788"),
      mid  = list(gene = "#C4B6CC", cpd = "#D8CADA"),
      high = list(gene = "#EE6677", cpd = "#FFFF00")
    )
    
    # clean up temporary files
    rm(pv)
    gc();gc()
    
    if (is.matrix(gene_fc)){
      out_res$pathview_filename <- paste0(kegg_id, ".pathview.multi.png")
    } else {
      # pattern3: using all groups information 
      out_res$pathview_filename <- paste0(kegg_id, ".pathview.png") 
    }
  }
  return(out_res)
}


