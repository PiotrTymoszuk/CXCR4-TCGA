# This script generates report figures: main and the supplements

  insert_head()
  
# data containers -----
  
  report_figures <- list()
  suppl_figures <- list()
  
# Figure 1: participant clustering by immune features -----
  
  insert_msg('Figure 1: participants clusters')
  
  report_figures$clusters <- plot_grid(immune_clust$hm_plot + 
                                         theme(plot.tag = element_blank()), 
                                       plot_grid(ggdraw(), 
                                                 immune_clust$cxcr4_plot + 
                                                   theme(legend.position = 'none'), 
                                                 ggdraw(), 
                                                 ncol = 3, 
                                                 rel_widths = c(0.15, 0.6, 0.25)), 
                                       nrow = 2, 
                                       rel_heights = c(0.6, 0.4), 
                                       labels = LETTERS, 
                                       label_size = 10) %>% 
    as_figure(label = 'figure_1_clustering_cxcr4', 
              w = 180, 
              h = 180)
  
# Figure 2: correlation of CXCR4 expression with the QuanTIseq immune scores ----

  insert_msg('Figure 2: correlations of CXCR4 with the QuanTIseq immune scores')
  
  report_figures$correlation <- immuno_corr$plots[c('Macrophage M2_QUANTISEQ', 
                                                    'T cell regulatory (Tregs)_QUANTISEQ', 
                                                    'T cell CD8+_QUANTISEQ', 
                                                    'uncharacterized cell_QUANTISEQ')] %>% 
    map(~.x + 
          theme(plot.tag = element_blank()) + 
          labs(y = 'Quantiseq estimate')) %>% 
    plot_grid(plotlist = ., 
              ncol = 1) %>% 
    plot_grid(plot_grid(immuno_corr$forest_plot + 
                          theme(legend.position = 'none') +
                          labs(title = 'Correlation with CXCR4 expression', 
                               subtitle = 'Quantiseq estimates'), 
                        ggdraw(), 
                        nrow = 2, 
                        rel_heights = c(0.85, 0.15)) + 
                theme(plot.margin = ggplot2::margin(r = 5, unit = 'mm')), 
              ., 
              ncol = 2, 
              rel_widths = c(0.65, 0.35), 
              labels = LETTERS, 
              label_size = 10) %>% 
    as_figure(label = 'figure_2_infiltration_correlation', 
              w = 180, 
              h = 210)
  
# Figure 3: correlations of the pooled immune signatures and CXCR4 with survival -----
  
  insert_msg('Figure 3: immune features, CXCR4 and survival')
  
  report_figures$survival <- cxcr_surv$km_plots %>% 
    map(~.x[['CXCR4']] + 
          theme(legend.position = 'none')) %>%
    c(., list(legend = get_legend(cxcr_surv$km_plots$os[[1]] + 
                                    theme(legend.position = 'right') + 
                                    labs(color = NULL)))) %>% 
    plot_grid(plotlist = ., 
              ncol = 2, 
              align = 'hv', 
              axis = 'tblr') %>% 
    as_figure(label = 'figure_3_survival_km', 
              w = 180, 
              h = 160)

# Figure 4: expression of immune scores in the CXCR4high and CXCR4low tumors ------
  
  insert_msg('Figure 4: expression of immune signatures in the CXCR4 strata')
  
  report_figures$immune_signatures <- imm_function$plots[c('immune score_XCELL', 
                                                           'microenvironment score_XCELL', 
                                                           'cytotoxicity score_MCPCOUNTER', 
                                                           'stroma score_XCELL', 
                                                           'exp_sign_genes', 
                                                           'exhaustion_genes', 
                                                           'tis_genes', 
                                                           'ifn_sign_genes', 
                                                           'cytotoxic_genes')] %>% 
    map(~.x + 
          theme(plot.tag = element_blank(), 
                axis.text.x = element_blank(), 
                legend.position = 'none') + 
          labs(x = '', 
               y = 'GSVA estimate')) %>% 
    c(list(legend = plot_grid(get_legend(imm_function$plots[[1]]), 
                              ggdraw() + 
                                draw_text(imm_function$plots[[1]]$labels$tag, 
                                          size = 8, 
                                          hjust = 0, 
                                          x = 0.25), 
                              nrow = 2))) %>% 
    plot_grid(plotlist = ., 
              ncol = 4, 
              align = 'hv') %>% 
    as_figure(label = 'figure_4_immune_signatures', 
              w = 180, 
              h = 180)
  
# Figure 5: expression of T cell-related genes in the CXCRhigh and CXCR4low tumors -----
  
  insert_msg('Figure 5: expression of T cell-related genes in CXCR4high and CXCR4low tumors')
  
  report_figures$t_cell_genes <- plot_grid(imm_genes$forest_plot + 
                                             theme(legend.position = 'none', 
                                                   plot.tag.position = 'right')) %>% 
    as_figure(label = 'figure_5_t_cell_genes', 
              w = 180, 
              h = 210)
  
# Figure 6: differential gene expression in CXCR4 expression strata ----
  
  insert_msg('Figure 6: differential gene expression in the CXCR4 expression strata')
  
  report_figures$dge <- dge$volcano + 
    theme(legend.position = 'bottom')
  
  report_figures$dge <- report_figures$dge %>% 
    as_figure(label = 'figure_6_diff_gene_exprs', 
              w = 180, 
              h = 180)
  
# Figure 7: results of pathway modulation analysis ----
  
  insert_msg('Figure 7: pathway modulation analysis')
  
  report_figures$spia <- dge_enrich$spia_barplot + 
    theme(plot.tag = element_blank())
  
  report_figures$spia <- report_figures$spia %>% 
    as_figure(label = 'figure_7_signaling', 
              w = 180, 
              h = 210)

# Figure S1: cluster development -----
  
  insert_msg('Figure S1: cluster development')
  
  suppl_figures$cluster_dev <- plot_grid(plot(immune_clust$clust_results, 
                                              'training', 
                                              cust_theme = globals$common_theme)$observation + 
                                           theme(plot.tag = element_blank()) + 
                                           labs(title = 'Participant SOM clustering', 
                                                subtitle = 'Training process'), 
                                         ggdraw(), 
                                         immune_clust$clust_qc_plots$node$wss + 
                                           labs(title = 'Clustering of the SOM nodes'), 
                                         immune_clust$clust_qc_plots$node$dendrogram + 
                                           expand_limits(y = -0.2) + 
                                           labs(title = 'Clustering of the SOM nodes', 
                                                subtitle = immune_clust$clust_qc_plots$node$wss$labels$subtitle), 
                                         ncol = 2, 
                                         align = 'hv', 
                                         axis = 'tblr', 
                                         labels = c('A', '', 'B', 'C'), 
                                         label_size = 10) %>% 
    as_figure(label = 'figure_s1_cluster_development', 
              w = 180, 
              h = 180)

# Figure S2: levels of key immune cell population and bona-fide tumor cells in the clusters -----
  
  insert_msg('Figure S2: infiltration and clusters')
  
  suppl_figures$infil_clusters <- immune_clust$violin_plots[c('Macrophage M1_QUANTISEQ', 
                                                              'Macrophage M2_QUANTISEQ', 
                                                              'T cell CD4+ (non-regulatory)_QUANTISEQ', 
                                                              'T cell regulatory (Tregs)_QUANTISEQ', 
                                                              'T cell CD8+_QUANTISEQ', 
                                                              'uncharacterized cell_QUANTISEQ')] %>% 
    map(~.x + 
          theme(legend.position = 'none') + 
          expand_limits(y = 0)) %>% 
    plot_grid(plotlist = ., 
              ncol = 2, 
              align = 'hv', 
              axis = 'tblr') %>% 
    as_figure('figure_s2_infiltration_clusters', 
              w = 180, 
              h = 210)
  
# Figure S3: differential gene expression, top regulated genes -----
  
  insert_msg('Figure S2: Differential gene expression, top genes')
  
  suppl_figures$top_dge <- dge$top_plot %>% 
    as_figure(label = 'figure_s3_top_diff_regulated', 
              w = 180, 
              h = 210)
  
# Figure S4: GO and KEGG enrichment analysis ----
  
  insert_msg('Figure S4: GO and KEGG term enrichment analysis')

  suppl_figures$enrichment <- plot_grid(dge_enrich$enrich_plots$go_results.downregulated + 
                                          theme(legend.position = 'none'), 
                                        dge_enrich$enrich_plots$go_results.upregulated + 
                                          theme(legend.position = 'none'), 
                                        nrow = 2, 
                                        align = 'hv', 
                                        rel_heights = c(0.57, 0.43), 
                                        labels = LETTERS, 
                                        label_size = 10) %>% 
    plot_grid(., 
              get_legend(dge_enrich$enrich_plots$go_results), 
              ncol = 2, 
              rel_widths = c(0.8, 0.2)) %>% 
    as_figure(label = 'figure_s4_go_kegg_enrich', 
              w = 180, 
              h = 220)
  
# Saving the figures -----
  
  insert_msg('Saving the figures')
  
  report_figures %>% 
    walk(save_figure, 
         path = './report/figures', 
         device = cairo_pdf)
  
  suppl_figures %>% 
    walk(save_figure, 
         path = './report/supplementary figures', 
         device = cairo_pdf)
  
# END ----
  
  insert_tail()