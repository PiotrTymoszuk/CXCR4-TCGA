# This script creates report tables

  insert_head()
  
# data container -----
  
  suppl_tables <- list()
  
# table s1: correlations of CXCR4 with immune features -----
  
  insert_msg('Table S1: correlations')
  
  suppl_tables$correlations <- immuno_corr$correlations %>%
    mutate(variable2 = stri_split(variable2, fixed = '_', simplify = TRUE)[, 1], 
           eff_size = stri_replace(eff_size, fixed = 'rho', replacement = '\u03C1')) %>% 
    select(variable2, eff_size, significance) %>% 
    set_names(c('Feature', 'Corr. coefficient', 'pFDR'))

# table S2: differentially expressed genes ----
  
  insert_msg('Table S2: differentially expressed genes')
  
  suppl_tables$dge <- dge$test_results %>% 
    filter(regulation != 'ns') %>% 
    select(response, 
           entrez_id,
           n, 
           regulation, 
           estimate, 
           lower_ci, 
           upper_ci, 
           p_adjusted) %>% 
    arrange(desc(estimate)) %>% 
    mutate(estimate = signif(estimate, 3), 
           lower_ci = signif(lower_ci, 3), 
           upper_ci = signif(upper_ci, 3), 
           p_adjusted = signif(p_adjusted, 2)) %>% 
    set_names(c('Gene', 
                'Entrez_ID', 
                'N', 
                'Regulation', 
                'Log2 Regulation', 
                'Lower CI', 
                'Upper CI', 
                'pFDR'))
  
# table S3: significant results of SPIA ----
  
  insert_msg('Table S3: significant results of SPIA')
  
  suppl_tables$spia_results <- dge_enrich$spia_results %>% 
    select(Name, Status, tA, pGFWER) %>% 
    mutate(tA = signif(tA, 3), 
           pGFWER = signif(pGFWER, 2)) %>% 
    filter(pGFWER < 0.05) %>%
    set_names(c('Pathway', 
                'Status', 
                'Modulation tA', 
                'pGFDR'))
  
# saving the tables -----
  
  insert_msg('Saving the tables')
  
  suppl_tables %>% 
    set_names(c('table S1 correlations', 
                'table S2 diff gene expression', 
                'table S3 signaling pathways')) %>% 
    write_xlsx(path = './report/supplementary tables.xlsx')
  
# END ------
  
  insert_tail()