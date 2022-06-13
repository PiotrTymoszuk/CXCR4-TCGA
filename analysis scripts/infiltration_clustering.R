# This script clusters the study participants in respect to the Quantiseq infiltration estimates

  insert_head()

# container ----
  
  immune_clust <- list()
  
# globals: analysis table and SOM grids -----
  
  insert_msg('Globals setup')
  
  ## clustering variables
  
  immune_clust$variables <- globals$infiltration %>% 
    filter(algorithm == 'QUANTISEQ', 
           score == 'no') 
  
  ## analysis table, min/max normalization
  
  immune_clust$analysis_tbl <- tcga$expression %>% 
    filter(tissue_type == 'Tumor') %>%
    select(patient_id, 
           all_of(immune_clust$variables$variable)) %>% 
    column_to_rownames('patient_id') %>% 
    min_max(complete_cases = TRUE)
  
  ## cluster colors
  
  immune_clust$clust_colors <- c('#1' = 'gray60', 
                                 '#2' = 'cornflowerblue', 
                                 '#3' = 'coral3')

# clustering the participants and immune features with SOM -----  
  
  insert_msg('Participnant clustering')

  immune_clust$clust_results <- combi_cluster(data = immune_clust$analysis_tbl, 
                                              distance_som = 'euclidean', 
                                              xdim = 7, 
                                              ydim = 7, 
                                              topo = 'hexagonal', 
                                              neighbourhood.fct = 'gaussian', 
                                              toroidal = FALSE, 
                                              rlen = 1000, 
                                              node_clust_fun = hcluster, 
                                              distance_nodes = 'euclidean', 
                                              k = 3, 
                                              seed = 123)
  
  ## cluster name recoding
  
  immune_clust$clust_results$clust_assignment <- immune_clust$clust_results$clust_assignment %>% 
    mutate(clust_id = car::recode(clust_id, "'1' = '#1'; '2' = '#2'; '3' = '#3'"))
  
  ## cluster quality control: plots, variance and cross-validation
  
  immune_clust$clust_qc_plots <- plot(immune_clust$clust_results, 
                                      type = 'diagnostic', 
                                      cust_theme = globals$common_theme)
  
  immune_clust$clust_variance <- var(immune_clust$clust_results)
  
  immune_clust$clust_pca <- plot(immune_clust$clust_results, 
                                 type = 'components', 
                                 with = 'data', 
                                 k = 3, 
                                 red_fun = 'pca', 
                                 cust_theme = globals$common_theme)
  
  immune_clust$clust_pca[c('node', 'final')] <- 
    immune_clust$clust_pca[c('node', 'final')] %>% 
    map(~.x + scale_fill_manual(values = immune_clust$clust_colors))
  
  immune_clust$clust_umap <- plot(immune_clust$clust_results, 
                                  type = 'components', 
                                  with = 'data', 
                                  k = 2, 
                                  red_fun = 'umap', 
                                  cust_theme = globals$common_theme, 
                                  random_state = 123)
    
  
  immune_clust$clust_cv <- cv(immune_clust$clust_results, 
                              nfolds = 10, 
                              kNN = 5, 
                              seed = 123, 
                              .parallel = FALSE)

# appending the analysis table with the cluster assignment information and CXCR4 expression ----
  
  insert_msg('Appending the analsis table with the cluster assignment information')
  
  immune_clust$analysis_tbl <- left_join(immune_clust$analysis_tbl %>% 
                                           rownames_to_column('patient_id'), 
                                         immune_clust$clust_results$clust_assignment %>% 
                                           set_names(c('patient_id', 'node', 'clust_id')), 
                                         by = 'patient_id') %>% 
    left_join(tcga$expression[c('patient_id', 'CXCR4')], 
              by = 'patient_id') %>% 
    as_tibble
  
# checking the immune estimate  and CXCR4 level differences between the clusters via Kruskal-Wallis test ----
  
  insert_msg('Checking the differences in immune infiltration and CXCR4 expression between the clusters')
  
  ## feature distribution and eov
  
  immune_clust$infil_dist <- explore(immune_clust$analysis_tbl, 
                                     variables = c(immune_clust$variables$variable, 'CXCR4'), 
                                     what = 'normality', 
                                     pub_styled = TRUE)
  
  immune_clust$infil_eov <- compare_variables(immune_clust$analysis_tbl, 
                                              variables = c(immune_clust$variables$variable, 'CXCR4'), 
                                              split_factor = 'clust_id', 
                                              what = 'variance', 
                                              pub_styled = TRUE)
  
  ## descriptive stats
  
  immune_clust$infil_desc <- explore(immune_clust$analysis_tbl, 
                                     variables = c(immune_clust$variables$variable, 'CXCR4'), 
                                     split_factor = 'clust_id', 
                                     what = 'table', 
                                     pub_styled = TRUE) %>% 
    reduce(left_join, by = 'variable') %>% 
    set_names(c('variable', '#1', '#2', '#3'))
  
  ## testing
  
  immune_clust$infil_test <- compare_variables(immune_clust$analysis_tbl, 
                                               variables = c(immune_clust$variables$variable, 'CXCR4'), 
                                               split_factor = 'clust_id', 
                                               what = 'test', 
                                               types = c(rep('kruskal_test', length(immune_clust$variables$variable)), 
                                                           'anova'), 
                                               ci = FALSE, 
                                               adj_method = 'BH')
  
# heat map of the clustering features -----
  
  insert_msg('Heat map of the clustering features')
  
  immune_clust$hm_plot <- plot_clust_hm(x_object = immune_clust$clust_results, 
                                        plot_title = 'Sample clustering by immune estimates', 
                                        plot_subtitle = 'Quantiseq infiltration estimates, SOM/HCl', 
                                        x_lab = 'PDAC sample', 
                                        cust_theme = globals$common_theme) + 
    labs(fill = 'Min/max estimate') +
    scale_y_discrete(labels = set_names(immune_clust$variables$cell_type, 
                                        immune_clust$variables$variable))
  
  immune_clust$hm_plot <- immune_clust$hm_plot + 
    labs(tag = paste0('\n', immune_clust$hm_plot$labels$tag))

# plotting the CXCR4 expression in the clusters ----
  
  immune_clust$cxcr4_plot <- plot_variable(immune_clust$analysis_tbl, 
                                           variable = 'CXCR4', 
                                           split_factor = 'clust_id', 
                                           type = 'box', 
                                           point_alpha = 0.8, 
                                           cust_theme = globals$common_theme, 
                                           x_lab = 'Cluster', 
                                           y_lab = expression('log'[2]*italic(' CXCR4')*' expression'), 
                                           plot_subtitle = immune_clust$infil_test %>% 
                                             filter(variable == 'CXCR4') %>% 
                                             .$significance, 
                                           plot_title = 'CXCR4 in immune inflitration clusters') + 
    scale_fill_manual(values = immune_clust$clust_colors)
  
  immune_clust$cxcr4_plot <- immune_clust$cxcr4_plot + 
    labs(tag = paste0('\n', stri_replace_all(immune_clust$cxcr4_plot$labels$tag, 
                                             fixed = '\n', 
                                             replacement = ', ')))
  
# single violin plots for the immune infiltrates in the clusters -----
  
  insert_msg('Immune cell levels in violin plots')
  
  immune_clust$violin_plot_tbl <- tcga$expression %>% 
    filter(tissue_type == 'Tumor') %>%
    select(patient_id, 
           all_of(immune_clust$variables$variable)) %>% 
    left_join(immune_clust$clust_results$clust_assignment %>% 
                set_names(c('patient_id', 'node', 'clust_id')), 
              by = 'patient_id')
  
  immune_clust$violin_plots <- list(variable = immune_clust$variables$variable, 
                                    plot_title = immune_clust$variables$cell_type, 
                                    plot_subtitle = immune_clust$infil_test$significance[1:10]) %>% 
    pmap(plot_variable, 
         immune_clust$violin_plot_tbl, 
         split_factor = 'clust_id', 
         type = 'violin', 
         point_hjitter = 0, 
         cust_theme = globals$common_theme, 
         x_lab = 'Immune infiltration cluster') %>% 
    map(~.x + 
          labs(tag = stri_replace_all(.x$labels$tag, fixed = '\n', replacement = ', ') %>% 
                 paste0('\n', .)) + 
          scale_fill_manual(values = immune_clust$clust_colors)) %>% 
    set_names(immune_clust$variables$variable)
  
# END -----

  insert_tail()