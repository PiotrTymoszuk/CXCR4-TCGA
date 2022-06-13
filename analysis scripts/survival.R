# This script investigates survival differences in CXCR4 and immune signature low and high expressors
# the optimal expression cutoff determined iteratively.

  insert_head()
  
# data container ----
  
  cxcr_surv <- list()
  
# globals: analysis table and survival objects ------
  
  ## variables to model, labels
  
  cxcr_surv$variables <- c(globals$gene_interest, 
                           globals$genesig_index$sign_var, 
                           globals$infiltration %>% 
                             filter(algorithm == 'QUANTISEQ') %>% 
                             .$variable)
  
  cxcr_surv$var_labels <- c(globals$gene_interest, 
                            globals$genesig_index$sign_name, 
                            globals$infiltration %>% 
                              filter(algorithm == 'QUANTISEQ') %>% 
                              .$plot_label) %>% 
    stri_replace(fixed = ': QUANTISEQ', 
                 replacement = '')
  
  ## naming vectors
  
  cxcr_surv$response_names <- c('os', 'trs', 'rfs')
  
  cxcr_surv$response_labs <- c('Overall survival', 
                               'Tumor-related survival', 
                               'Relapse-free survival') %>% 
    set_names(cxcr_surv$response_names)
  
  cxcr_surv$time_labs <- c('OS, days',
                           'TRS, days', 
                           'RFS, days') %>% 
    set_names(cxcr_surv$response_names)

  ## COX modeling variables
  
  cxcr_surv$cox_vars <- paste(cxcr_surv$variables, 
                              'strata', sep = '_')
  
  cxcr_surv$cox_confounders <- c('age', 
                                 'sex', 
                                 'surgical_procedure', 
                                 'tumor_grade', 
                                 'tumor_size', 
                                 'residual_tumor')
  
  ## analysis table: cxcr4 expression and clinical features for multi-variate modeling
  
  cxcr_surv$analysis_table <- tcga$expression %>% 
    filter(tissue_type == 'Tumor') %>% 
    select(patient_id, 
           all_of(cxcr_surv$variables), 
           all_of(cxcr_surv$cox_confounders), 
           death, 
           death_panc, 
           os_days, 
           relapse, 
           rfs_days, 
           treatment_outcome, 
           treatment_outcome_score)
  
  names(cxcr_surv$analysis_table) <- make.names(names(cxcr_surv$analysis_table))
  
  cxcr_surv$variables <- make.names(cxcr_surv$variables)
  
  cxcr_surv$var_labels <- set_names(cxcr_surv$var_labels, 
                                    cxcr_surv$variables)

# Finding the optimal cutoffs for CXCR4 expression in OS, OS tumor and RFS ----
  
  insert_msg('Finding optimal CXCR4 expression cutoffs')
  
  plan('multisession')

  cxcr_surv$cutpoint_objects <- cxcr_surv$variables %>% 
    future_map(function(var) list(event = c('death', 'death_panc', 'relapse'), 
                                  time = c('os_days', 'os_days', 'rfs_days')) %>% 
                 pmap(find_cutoff, 
                      variable = var, 
                      data = cxcr_surv$analysis_table, 
                      min_n = floor(0.2*(nrow(cxcr_surv$analysis_table)))) %>% 
                 set_names(c('os', 'trs', 'rfs')) %>% 
                 compact) %>% 
    set_names(cxcr_surv$variables) %>% 
    compact %>% 
    transpose
  
  plan('sequential')
  
  cxcr_surv$optimal_cutpoints <- cxcr_surv$cutpoint_objects %>% 
    map(~map2_dfr(.x, names(.x), ~tibble(variable = .y, 
                                         cutoff = cutoff(.x)[[1]]))) %>% 
    reduce(full_join, by = 'variable') %>% 
    set_names(c('variable', names(cxcr_surv$cutpoint_objects)))

# KM modeling with the optimal expression cutoffs -----
  
  insert_msg('KM modeling with the optimal cutoffs')

  ## summaries, p value correction
  
  cxcr_surv$km_summaries <- cxcr_surv$cutpoint_objects %>% 
    map(~map2_dfr(.x, names(.x), 
                  ~mutate(.x$cutoff_stats[1, ], variable = .y))) %>% 
    map(mutate, 
        p_adjusted = p.adjust(p_value, 'BH'), 
        plot_lab = paste('p =', signif(p_adjusted, 2)))

  
  ## plots 
  
  cxcr_surv$km_plots <- list(obj = cxcr_surv$cutpoint_objects, 
                             resp = names(cxcr_surv$cutpoint_objects), 
                             stats = map(cxcr_surv$km_summaries, ~.x$plot_lab)) %>% 
    pmap(function(obj, resp, stats) list(x = obj, 
                                         title = cxcr_surv$var_labels[names(obj)], 
                                         xlab = cxcr_surv$time_labs[[resp]], 
                                         pval = stats) %>% 
           pmap(plot, 
                type = 'km', 
                palette = unname(globals$strata_colors), 
                pval.size = 2.75) %>% 
           map(~.x$plot + 
                 globals$common_theme + 
                 labs(subtitle = .x$plot$labels$subtitle %>% 
                        stri_replace(regex = 'p =.*', replacement = '') %>% 
                        stri_replace(fixed = 'ns(', replacement = '') %>% 
                        paste0(.x$plot$labels$tag)) + 
                 theme(plot.tag = element_blank())))

# END ----
  
  insert_tail()