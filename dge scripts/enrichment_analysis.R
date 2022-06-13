# This script performs enrichment analyses in the set of upregulated and downregulated genes

  insert_head()
  
# data containers -----
  
  dge_enrich <- list()
  
# globals: regulation vectors -----
  
  insert_msg('Globals setup')
  
  dge_enrich$reg_vctrs <- dge$test_results %>% 
    filter(regulation != 'ns') %>% 
    dlply(.(regulation), 
          function(x) set_names(x$estimate, 
                                x$entrez_id))
  
  dge_enrich$all_vctrs <- dge$test_results$entrez_id
  
  dge_enrich$gene_pathway <- read_tsv('./input data/pathway_links.csv', 
                                      col_names = FALSE) %>% 
    set_names(c("GeneID", "PathwayID")) %>% 
    mutate(GeneID = stri_replace(GeneID, fixed = 'hsa:', replacement = ''))

  dge_enrich$pathway_names <- read_tsv('./input data/pathway_names.csv', 
                                       col_names = FALSE) %>%
    set_names(c("PathwayID", "Description"))
  
# GO and KEGG enrichment analysis ----
  
  insert_msg('GO and KEGG enrichment analysis')
  
  plan('multisession')
  
  dge_enrich$go_results <- dge_enrich$reg_vctrs %>% 
    map(names) %>% 
    future_map(goana, 
               .options = furrr_options(seed = TRUE))
  
  dge_enrich$kegg_results <- dge_enrich$reg_vctrs %>% 
    map(names) %>% 
    future_map(kegga, 
               gene.pathway = dge_enrich$gene_pathway, 
               pathway.names = dge_enrich$pathway_names, 
               species = 'Hs', 
               .options = furrr_options(seed = TRUE))
  
  plan('sequential')
  
# Formatting the GO/KEGG enrichment results ------
  
  insert_msg('Formatting the GO/KEGG results')

  dge_enrich[c('go_results', 
               'kegg_results')] <- list(map(dge_enrich$go_results, 
                                            filter, 
                                            Ont == 'BP'), 
                                        map(dge_enrich$kegg_results, 
                                            mutate, 
                                            Term = Pathway)) %>% 
    map(~map(.x, rownames_to_column, 'ID') %>% 
          map(mutate, p_adj = p.adjust(P.DE, 'BH')) %>% 
          map(filter, p_adj < 0.05) %>% 
          map2_dfr(., names(.), ~mutate(.x, gene_group = .y)) %>% 
          as_tibble)

# Plotting the top 20 most significant GO/KEGG terms ----
  
  insert_msg('Plot with top significant factors')
  
  dge_enrich$enrich_plots <- dge_enrich[c('go_results', 
                                          'kegg_results')] %>% 
    map(~dlply(.x, 'gene_group')) %>% 
    unlist(recursive = FALSE)
  
  dge_enrich$enrich_plots <- list(data = dge_enrich$enrich_plots, 
                                  plot_title = c(rep('BP GO enrichment: top 20 terms', 2), 
                                                 rep('KEGG enrichment: top 20 terms', 2)), 
                                  plot_subtitle = rep(c('Genes downregulated in CXCR4 hi tumors', 
                                                        'Genes upregulated in CXCR4 hi tumors'), 2), 
                                  fill_scale = list(c(significant = 'steelblue', ns = 'gray60'), 
                                                    c(significant = 'coral3', ns = 'gray60'), 
                                                    c(significant = 'steelblue', ns = 'gray60'), 
                                                    c(significant = 'coral3', ns = 'gray60'))) %>% 
    pmap(plot_signifcant, 
         p_variable = 'p_adj', 
         label_variable = 'Term', 
         top_significant = 20, 
         signif_level = 0.05, 
         x_lab = expression('-log'[10]*' pFDR'),
         fill_title = '', 
         cust_theme = globals$common_theme)

# SPIA -----
  
  insert_msg('Pathway analysis with spia')
  
  dge_enrich$spia_results <- spia(de = reduce(dge_enrich$reg_vctrs, c), 
                                  all = dge_enrich$all_vctrs) %>% 
    as_tibble

# Plotting the SPIA results as a volcano plot ----
  
  insert_msg('Plotting the SPIA results as a volcano')
  
  dge_enrich$spia_volcano <- plot_volcano(data = dge_enrich$spia_results, 
                                          regulation_variable = 'tA', 
                                          p_variable = 'pGFdr', 
                                          signif_level = 0.05, 
                                          regulation_level = 0, 
                                          x_lab = 'Pathway activation CXCR4 high vs low, tA', 
                                          y_lab = expression('-log'[10]*' pFDR'), 
                                          top_significant = 10, 
                                          label_variable = 'Name', 
                                          label_type = 'text', 
                                          txt_size = 2.5, 
                                          plot_title = 'Pathway perturbation analysis', 
                                          plot_subtitle = 'CXCR4 high vs low tumors, SPIA', 
                                          cust_theme = globals$common_theme) + 
    geom_vline(xintercept = 0, linetype = 'dashed') + 
    scale_fill_manual(values = c('firebrick', 'steelblue', 'gray70'), 
                      labels = c('activated', 'inhibited', 'ns'), 
                      name = 'Pathway\nstatus')
  
  dge_enrich$spia_volcano <- dge_enrich$spia_volcano + 
    labs(tag = dge_enrich$spia_volcano$labels$tag %>% 
           stri_replace(fixed = 'downregulated', replacement = 'inhibited') %>% 
           stri_replace(fixed = 'upregulated', replacement = 'activated') %>% 
           paste0('\n', .))
  
# plotting the SPIA results as a scatter plot ----
  
  insert_msg('Plotting the SPIA results as a bar plot')

  dge_enrich$spia_barplot <- plot_top(data = dge_enrich$spia_results %>% 
                                        filter(pGFdr < 0.05), 
                                      regulation_variable = 'tA', 
                                      label_variable = 'Name', 
                                      p_variable = 'pGFdr', 
                                      regulation_level = 0, 
                                      signif_level = 0.05, 
                                      top_regulated = 20, 
                                      plot_title = 'Pathway perturbation analysis: top 20 pathways', 
                                      plot_subtitle = 'CXCR4 high vs low tumors, SPIA', 
                                      x_lab = 'Pathway activation CXCR4 high vs low, tA', 
                                      cust_theme = globals$common_theme, 
                                      show_txt = FALSE) + 
    scale_color_manual(values = c('firebrick', 'steelblue', 'gray70'), 
                       labels = c('activated', 'inhibited', 'ns'), 
                       name = 'Pathway\nstatus') + 
    guides(fill = 'none')

# END -----
  
  insert_tail()