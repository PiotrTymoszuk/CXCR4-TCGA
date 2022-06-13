# This script checks the expression of CXCR4 in diverse immune infiltrates identified by Quantiseq
# by Spearman's correlation

  insert_head()
  
# data container ----  
  
  immuno_corr <- list()
  
# globals: cell type variables -----
  
  insert_msg('Globals setup')
  
  immuno_corr$variables <- globals$infiltration %>% 
    filter(algorithm == 'QUANTISEQ', 
           score == 'no')
  
  immuno_corr$var_pairs <- immuno_corr$variables$variable %>% 
    map(~c('CXCR4', .x)) %>% 
    set_names(immuno_corr$variables$variable)
  
# normality of the infiltration variables -----
  
  insert_msg('Normality')
  
  immuno_corr$normality <- explore(tcga$expression %>% 
                                     filter(tissue_type == 'Tumor'), 
                                   variables = c(immuno_corr$variables$variable, 'CXCR4'), 
                                   what = 'normality', 
                                   pub_styled = TRUE)
  
# correlations -----

  insert_msg('Serial correlation with Spearman test')
  
  immuno_corr$correlations <- immuno_corr$var_pairs %>% 
    map_dfr(~correlate_variables(tcga$expression %>% 
                                   filter(tissue_type == 'Tumor'), 
                                 variables = .x, 
                                 what = 'correlation', 
                                 type = 'spearman', 
                                 adj_method = 'BH', 
                                 pub_styled = TRUE)) %>% 
    mutate(plot_cap = paste(eff_size, significance, sep = ', '), 
           plot_cap = stri_replace(plot_cap, fixed = 'rho', replacement = '\u03C1'))

# serial plotting: single plots -----
  
  insert_msg('Serial plotting of single correlation plots')
  
  immuno_corr$plots <- list(variables = immuno_corr$var_pairs, 
                            y_lab = paste(immuno_corr$variables$cell_type, 'Quantiseq estimate', sep = '\n'), 
                            plot_subtitle = immuno_corr$correlations$plot_cap, 
                            plot_title = immuno_corr$variables$cell_type) %>% 
    pmap(plot_correlation, 
         tcga$expression %>% 
           filter(tissue_type == 'Tumor'),  
         type = 'correlation', 
         point_alpha = 0.8, 
         x_lab = expression('log'[2]*italic(' CXCR4')*' expression'), 
         cust_theme = globals$common_theme)
  
# summary plotting: forest plots for the most significant correlations for each algorithm -----
  
  insert_msg('Summary forest plot')
  
  immuno_corr$summary_plot_tbl <- immuno_corr$var_pairs %>% 
    map_dfr(~correlate_variables(tcga$expression %>% 
                                   filter(tissue_type == 'Tumor'), 
                                 variables = .x, 
                                 what = 'correlation', 
                                 type = 'spearman', 
                                 adj_method = 'BH', 
                                 pub_styled = FALSE)) %>% 
    mutate(significant = ifelse(p_adjusted < 0.05, 'significant', 'ns'), 
           regulation = ifelse(significant == 'ns', 
                               'ns', 
                               ifelse(estimate > 0, 'positive', 'negative')))
  
  ## correlation estimate forest plot
  
  immuno_corr$forest_plot <- plot_top(data = immuno_corr$summary_plot_tbl %>% 
                                        mutate(variable2 = stri_replace(variable2, 
                                                                        fixed = '_QUANTISEQ', 
                                                                        replacement = ''), 
                                               variable2 = stri_replace(variable2, 
                                                                        fixed = 'non-regulatory', 
                                                                        replacement = 'non-Treg')), 
                                      regulation_variable = 'estimate', 
                                      label_variable = 'variable2', 
                                      p_variable = 'p_adjusted', 
                                      signif_level = 0.05, 
                                      regulation_level = 0, 
                                      lower_ci_variable = 'lower_ci', 
                                      upper_ci_variable = 'upper_ci', 
                                      top_regulated = 100, 
                                      cust_theme = globals$common_theme, 
                                      plot_title = 'Correlation of CXCR4 expression with immune infiltration', 
                                      plot_subtitle = 'Quantiseq infiltration estimates, Spearman correlaton', 
                                      plot_tag = paste('\nn =', immuno_corr$summary_plot_tbl$n[1]), 
                                      x_lab = expression(rho), 
                                      show_txt = TRUE, 
                                      show_ci_txt = TRUE) + 
    scale_x_continuous(limits = c(-1.1, 1.1), 
                       breaks = seq(-1, 1, by = 0.25))
    
# END -----
  
  insert_msg()