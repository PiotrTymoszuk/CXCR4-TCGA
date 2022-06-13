# This script compares the immune signature gene members between the CXCR4 high and CXCR4 low tumors
# The cutoff is defined by the maximal difference in Tumor-Related Survival

  insert_head()
  
# data container -----
  
  imm_genes <- list()
  
# globals: variables to be analyzed, analysis table -----
  
  insert_msg('Globals setup')
  
  ## variables
  
  imm_genes$variables <- globals[globals$genesig_index$sign_var] %>% 
    reduce(union)

  ## analysis table, re-leveling the CXCR4 strata
  
  imm_genes$analysis_tbl <- tcga$expression %>% 
    filter(tissue_type == 'Tumor') %>% 
    mutate(CXCR4_strata = cut(CXCR4, 
                              c(-Inf, cxcr_surv$cutpoint_objects$os$CXCR4$cutoff, Inf), 
                              c('low', 'high'))) %>% 
    select(patient_id, 
           CXCR4_strata, 
           all_of(imm_genes$variables))
  
  ## n numbers
  
  imm_genes$n_nubers <- imm_genes$analysis_tbl %>% 
    count(CXCR4_strata)
  
  imm_genes$plot_tag <- paste0('low: n = ', imm_genes$n_nubers$n[1], 
                               '\nhigh: n = ', imm_genes$n_nubers$n[2])
  
# serial analysis ----
  
  insert_msg('Serial analysis, t test')
  
  ## regulation estimates
  
  imm_genes$summary <- test_two_groups(data = imm_genes$analysis_tbl, 
                                       split_fct = 'CXCR4_strata', 
                                       variables = imm_genes$variables, 
                                       type = 't', 
                                       adj_method = 'BH')

  ### effect size, Cohen's D
  
  imm_genes$eff_size <- compare_variables(imm_genes$analysis_tbl, 
                                             variables = imm_genes$variables, 
                                             split_factor = 'CXCR4_strata', 
                                             what = 'eff_size', 
                                             types = 'cohen_d', 
                                             ci = FALSE, 
                                             pub_styled = TRUE, 
                                             adj_method = 'BH')
  
# single plots, appending with the corrected p values -----
  
  insert_msg('Single plots')
  
  imm_genes$plots <- list(variable = imm_genes$variables, 
                          plot_title = imm_genes$variables, 
                          plot_subtitle = imm_genes$eff_size$significance) %>% 
    pmap(plot_variable, 
         imm_genes$analysis_tbl, 
         split_factor = 'CXCR4_strata', 
         type = 'box', 
         cust_theme = globals$common_theme, 
         x_lab = 'CXCR4 strata', 
         y_lab = expression('log'[2]*' expression')) %>% 
    map(~.x + 
          scale_fill_manual(values = globals$strata_colors, 
                            labels = globals$strata_labs)) %>% 
    set_names(imm_genes$variables)

# summary forest plot ----
  
  insert_msg('Summary forest plot')
  
  imm_genes$forest_plot <- plot_top(data = imm_genes$summary, 
                                    regulation_variable = 'estimate', 
                                    label_variable = 'response', 
                                    p_variable = 'p_adjusted', 
                                    signif_level = 0.05, 
                                    regulation_level = 0, 
                                    lower_ci_variable = 'lower_ci', 
                                    upper_ci_variable = 'upper_ci', 
                                    top_regulated = 100, 
                                    fill_title = 'Regulation', 
                                    plot_title = 'Immune gene expression in CXCR4 high vs low tumors', 
                                    plot_subtitle = 'T test', 
                                    plot_tag = imm_genes$plot_tag, 
                                    x_lab = expression('log'[2]*' regulation, CXCR4 high vs low'), 
                                    cust_theme = globals$common_theme, 
                                    show_txt = FALSE) + 
    theme(axis.text.y = element_text(size = 8, face = 'italic'))
    
# END -----
  
  insert_tail()