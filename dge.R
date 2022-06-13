# This scripts investigates differential gene expression between CXCR4 high and CXCR4 low tumors

# tools -----

  library(microViz)
  library(limma)
  library(SPIA)
  
  insert_head()

  source_all('./tools/project_tools.R', 
             message = TRUE, 
             crash = TRUE)
  
# data container -----

  dge <- list()

# globals: analysis table, stratification ------

  insert_msg('Globals setup')
  
  dge$analysis_tbl <- tcga$expression %>% 
    filter(tissue_type == 'Tumor') %>% 
    mutate(CXCR4_strata = cut(CXCR4, 
                              c(-Inf, cxcr_surv$cutpoint_objects$os$CXCR4$cutoff, Inf), 
                              c('low', 'high')))
  
  dge$n_numbers <- dge$analysis_tbl %>% 
    count(CXCR4_strata)
  
  dge$plot_tag <- paste0('\nlow: n = ', dge$n_numbers$n[1], 
                         ', high: n = ', dge$n_numbers$n[2])

# serial testing - two-tailed T test ------

  insert_msg('Serial testing')
  
  dge$test_results <- test_two_groups(data = dge$analysis_tbl, 
                                      variables = tcga$annotation$gene_symbol, 
                                      split_fct = 'CXCR4_strata', 
                                      .parallel = TRUE) %>% 
    mutate(n = nrow(dge$analysis_tbl), 
           regulation = ifelse(significant == 'no', 'ns', 
                               ifelse(estimate > log2(1.5), 
                                      'upregulated', 
                                      ifelse(estimate < -log2(1.5), 
                                             'downregulated', 
                                             'ns'))), 
           entrez_id = translate_gene(response), 
           gene_symbol = response)

# n numbers -----
  
  insert_msg('N numbers of the regulated genes')
  
  dge$gene_numbers <- dge$test_results %>% 
    filter(regulation != 'ns') %>% 
    dlply(.(regulation), nrow)

# Visualization of the dge results in a volcano plot ----
  
  insert_msg('Volcano plot')
  
  dge$volcano <- plot_volcano(data = dge$test_results, 
                              regulation_variable = 'estimate', 
                              p_variable = 'p_adjusted', 
                              signif_level = 0.05, 
                              regulation_level = log2(1.5), 
                              x_lab = expression('log'[2]*' regulation CXCR4 high vs low'), 
                              y_lab = expression('-log'[10]*' pFDR'), 
                              top_significant = 10, 
                              label_variable = 'gene_symbol', 
                              label_type = 'label', 
                              txt_size = 2.5, 
                              txt_face = 3, 
                              plot_title = 'Differential gene expression', 
                              plot_subtitle = 'Differential gene expression', 
                              fill_title = 'Regulation', 
                              cust_theme = globals$common_theme)
    
# Top up- and downregulated genes -----
  
  insert_msg('Plotting top regulated genes')
 
  dge$top_plot <- plot_top(data = dge$test_results, 
                           regulation_variable = 'estimate', 
                           label_variable = 'gene_symbol', 
                           p_variable = 'p_adjusted', 
                           signif_level = 0.05, 
                           regulation_level = 0, 
                           lower_ci_variable = 'lower_ci', 
                           upper_ci_variable = 'upper_ci', 
                           top_regulated = 20, 
                           plot_title = 'Differential gene expression: CXCR4 high vs low tumors', 
                           plot_subtitle = 'Top 20 differentially up- and downregulated genes', 
                           fill_title = 'Regulation', 
                           plot_tag = dge$plot_tag, 
                           x_lab = expression('log'[2]*' regulation CXCR4 high vs low'), 
                           cust_theme = globals$common_theme, 
                           show_txt = FALSE) + 
    theme(axis.text.y = element_text(face = 'italic'))
  
# additional analyses ----
  
  insert_msg('Additional analyses: GO, KEGG and pathway enrichment')
  
  source('./dge scripts/enrichment_analysis.R')

# END -----
  
  insert_tail()