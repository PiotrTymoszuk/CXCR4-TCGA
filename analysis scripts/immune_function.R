# This script compares the immune signature scores between the CXCR4 high and CXCR4 low tumors
# The cutoff is defined by the maximal difference in Tumor-Related Survival

  insert_head()
  
# data container -----
  
  imm_function <- list()
  
# globals: variables to be analyzed, analysis table -----
  
  insert_msg('Globals setup')
  
  ## variables and their labels
  
  imm_function$variables <- c(globals$genesig_index$sign_var, 
                              globals$infiltration %>% 
                                filter(score == 'yes') %>% 
                                .$variable)
  
  imm_function$var_labels <- c(globals$genesig_index$sign_name,  
                               globals$infiltration %>% 
                                 filter(score == 'yes') %>% 
                                 .$plot_label) %>% 
    stri_replace(fixed = ': pooled', 
                 replacement = '') %>% 
    set_names(imm_function$variables)
  
  ## analysis table, re-leveling the CXCR4 strata
  
  imm_function$analysis_tbl <- tcga$expression %>% 
    filter(tissue_type == 'Tumor') %>% 
    mutate(CXCR4_strata = cut(CXCR4, 
                              c(-Inf, cxcr_surv$cutpoint_objects$os$CXCR4$cutoff, Inf), 
                              c('low', 'high'))) %>%
    select(patient_id, 
           CXCR4_strata, 
           all_of(imm_function$variables))
  
# Variable normality and EOV -----
  
  insert_msg('Normality and EOV')
  
  imm_function$normality <- explore(imm_function$analysis_tbl, 
                                    variables = imm_function$variables, 
                                    what = 'normality', 
                                    pub_styled = TRUE)
  
  imm_function$eov <- compare_variables(imm_function$analysis_tbl, 
                                        variables = imm_function$variables, 
                                        split_factor = 'CXCR4_strata', 
                                        what = 'variance', 
                                        pub_styled = TRUE)
  
# descriptive statistics -----
  
  insert_msg('Descriptive stats')
  
  imm_function$desc_stats <- explore(imm_function$analysis_tbl, 
                                     variables = imm_function$variables, 
                                     split_factor = 'CXCR4_strata', 
                                     what = 'table', 
                                     pub_styled = TRUE) %>% 
    reduce(left_join, by = 'variable') %>% 
    set_names(c('variable', 'cxcr4_low', 'cxcr4_hi'))
  
# serial analysis ----
  
  insert_msg('Serial analysis, t test')

  ## effect size, Cohen's D
  
  imm_function$eff_size <- compare_variables(imm_function$analysis_tbl, 
                                             variables = imm_function$variables, 
                                             split_factor = 'CXCR4_strata', 
                                             what = 'eff_size', 
                                             types = 'cohen_d', 
                                             ci = FALSE, 
                                             pub_styled = TRUE, 
                                             adj_method = 'BH')

# single plots, appending with the corrected p values -----
  
  insert_msg('Single plots')
  
  imm_function$plots <- list(variable = imm_function$variables, 
                             plot_title = imm_function$var_labels, 
                             y_lab = paste(imm_function$var_labels, 'GSVA estimate', sep = ', '), 
                             plot_subtitle = imm_function$eff_size$significance) %>% 
    pmap(plot_variable, 
         imm_function$analysis_tbl, 
         split_factor = 'CXCR4_strata', 
         type = 'box', 
         cust_theme = globals$common_theme, 
         x_lab = 'CXCR4 strata') %>% 
    map(~.x + 
          scale_fill_manual(values = globals$strata_colors, 
                            labels = globals$strata_labs)) %>% 
    set_names(imm_function$variables)
  
# END -----
  
  insert_tail()