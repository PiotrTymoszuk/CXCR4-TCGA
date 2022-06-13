# this script calculates signature scores using gene set variation analysis tool (GSVA)

  insert_head()

# data container -----
  
  gene_sign <- list()
  
# gene signature calculation -----
  
  insert_msg('Calculating the gene signatures')
  
  gene_sign$signature_est <- calculate(x = globals[globals$genesig_index$sign_var], 
                                       data = tcga$expression[tcga$annotation$gene_symbol]) %>% 
    mutate(sample_id = tcga$expression$sample_id)
  
# signature calculation of the particular cell types, scaled variables -----
  
  insert_msg('Calculating the cell signatures')
  
  gene_sign$cell_tbl <- tcga$expression %>% 
    select(all_of(globals$infiltration$variable)) %>%
    map_dfc(~scale(.x)[, 1]) %>% 
    mutate(sample_id = tcga$expression$sample_id, 
           tissue_type = tcga$expression$tissue_type)
  
  gene_sign$cell_est <- globals$infiltration %>% 
    filter(score == 'no') %>% 
    dlply(.(cell_type)) %>% 
    map(~.x$variable) %>% 
    calculate(data = gene_sign$cell_tbl %>% 
                select( - sample_id, - tissue_type)) %>% 
    mutate(sample_id = gene_sign$cell_tbl$sample_id)
  
  names(gene_sign$cell_est)[names(gene_sign$cell_est) != 'sample_id'] <- 
    paste(names(gene_sign$cell_est)[names(gene_sign$cell_est) != 'sample_id'], 
          'pooled', 
          sep = '_')
  
# updating the globals ------
  
  insert_msg('Updating the cell infiltration table in globals')
  
  globals$infiltration <- tibble(variable = names(gene_sign$cell_est)[names(gene_sign$cell_est) != 'patient_id']) %>% 
    mutate(cell_type = stri_replace(variable, 
                                    fixed = '_pooled', 
                                    replacement = ''), 
           algorithm = 'pooled', 
           score = 'no', 
           plot_label = paste(cell_type, 
                              algorithm, 
                              sep = ': '), 
           plot_label_comp = paste(cell_type, 
                                   algorithm, 
                                   sep = '\n')) %>% 
    rbind(globals$infiltration, .)

# END -----
  
  insert_tail()