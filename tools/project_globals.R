# This script contains project globals

# libraries ---

  require(plyr)
  require(tidyverse)
  require(stringi)
  
# data container ------
  
  globals <- list()
  
# graphics -----
  
  globals$common_text <- element_text(size = 8, 
                                      face = 'plain', 
                                      color = 'black')
  
  globals$common_margin <- ggplot2::margin(t = 5, 
                                           l = 4, 
                                           r = 2, 
                                           unit = 'mm')
  
  globals$common_theme <- theme_classic() + theme(axis.text = globals$common_text, 
                                                  axis.title = globals$common_text, 
                                                  plot.title = element_text(size = 8, 
                                                                            face = 'bold'), 
                                                  plot.subtitle = globals$common_text, 
                                                  plot.tag = element_text(size = 8, 
                                                                          face = 'plain', 
                                                                          color = 'black', 
                                                                          hjust = 0, 
                                                                          vjust = 1), 
                                                  plot.tag.position = 'bottom', 
                                                  legend.text = globals$common_text, 
                                                  legend.title = globals$common_text, 
                                                  strip.text = globals$common_text,
                                                  strip.background = element_rect(fill = 'gray95', 
                                                                                  color = 'gray80'), 
                                                  plot.margin = globals$common_margin, 
                                                  panel.grid.major = element_line(color = 'gray90'))
# gene variables -----
  
  ## the main gene of interest
  
  globals$gene_interest <- 'CXCR4'
  
  ## T cell exhaustion genes, according to Woroniecka (DOI 10.1158/1078-0432.CCR-17-1846)
  
  globals$exhaustion_genes <- c('CD274', ## PD-L1
                                'PDCD1LG2', ## PD-L2
                                'PDCD1', ## PD-1
                                'HAVCR2', ## TIM3
                                'CTLA4', 
                                'TIGIT', 
                                'LAG3', 
                                'CD160', 
                                'CD244', 
                                'ENTPD1', 
                                'BTLA') 
  
  ## cytotoxic genes, according to DOI: 10.3389/fimmu.2021.622563/full
  
  globals$cytotoxic_genes <- c('CCL5', 
                               'GBP5', 
                               'GZMA', 
                               'GZMH', 
                               'IRF1', 
                               'LAG3', 
                               'NKG7', 
                               'PRF1', 
                               'PSMB10')
  
  ## IFNg signaling, DOI: 10.1172/JCI91190. 
  
  globals$ifn_sign_genes <- c('IDO1', 
                              'CXCL10', 
                              'CXCL9', 
                              'HLA-DRA', 
                              'STAT1', 
                              'IFNG')
  
  ## expanded immune signature, DOI: 10.1172/JCI91190.
  
  globals$exp_sign_genes <- c('CD3D', 
                              'IDO1', 
                              'CIITA', 
                              'CD3E', 
                              'CCL5', 
                              'GZMK', 
                              'CD2',
                              'HLA-DRA', 
                              'CXCL13', 
                              'IL2RG', 
                              'NKG7', 
                              'CXCR6', 
                              'LAG3', 
                              'TAGAP', 
                              'CXCL10', 
                              'STAT1', 
                              'GZMB')
  
  ## tumor inflammation signature (TIS), according to DOI: 10.1186/s40425-018-0367-1
  
  globals$tis_genes <- c('CCL5', 
                         'CD27', 
                         'CD274', 
                         'CD276', 
                         'CD8A', 
                         'CMKLR1', 
                         'CXCL9', 
                         'CXCR6', 
                         'HLA-DQA1', 
                         'HLA-DRB1', 
                         'HLA-E', 
                         'IDO1', 
                         'LAG3', 
                         'NKG7', 
                         'PDCD1LG2', 
                         'PSMB10', 
                         'STAT1', 
                         'TIGIT')
  
# gene signature names -----
  
  globals$genesig_index <- tibble(sign_var = c('exhaustion_genes', 
                                                'cytotoxic_genes', 
                                                'ifn_sign_genes', 
                                                'exp_sign_genes', 
                                                'tis_genes'), 
                                  sign_name = c('Exhaustion Sig', 
                                                'Cytotoxic Sig', 
                                                'IFN-\u03B3 Sig',
                                                'Exp. ImmSig', 
                                                'TISig'), 
                                  sig_color = c('olivegreen4', 
                                                'firebrick3', 
                                                'darkorange3', 
                                                'steelblue', 
                                                'indianred3'))
  
# immune infiltration estimates -----
  
  globals$infiltration <- timer_data$est_table %>% 
    select(- sample_id) %>% 
    names
  
  globals$infiltration <- tibble(variable = globals$infiltration, 
                                 cell_type = stri_split_fixed(globals$infiltration, 
                                                              pattern = '_', 
                                                              simplify = TRUE)[, 1], 
                                 algorithm = stri_split_fixed(globals$infiltration, 
                                                              pattern = '_', 
                                                              simplify = TRUE)[, 2]) %>% 
    mutate(score = ifelse(stri_detect(variable, fixed = 'score'), 'yes', 'no'), ## index, if the signature is a score
           plot_label = paste(cell_type, 
                              algorithm, 
                              sep = ': '), 
           plot_label = stri_replace(plot_label, fixed = ' score', replacement = ''), 
           plot_label = stri_replace(plot_label, fixed = 'microenvironment', replacement = 'TME'), 
           plot_label = stri_replace(plot_label, fixed = 'cytotoxicity', replacement = 'Cytotox.'), 
           plot_label = stri_replace(plot_label, fixed = 'immune', replacement = 'Imm.Score'), 
           plot_label = stri_replace(plot_label, fixed = 'stroma', replacement = 'Stroma'), 
           plot_label_comp = paste(cell_type, 
                                   algorithm, 
                                   sep = '\n'))
  
  ## skipping the progenitor cell types irrelevant for the tumor immunity, B cells
  ## mast cells and eosinophils
  
  globals$infiltration <- globals$infiltration %>% 
    filter(!cell_type %in% c('Common myeloid progenitor', 
                             'Common lymphoid progenitor', 
                             'Granulocyte-monocyte progenitor', 
                             'Hematopoietic stem cell', 
                             'B cell naive', 
                             'B cell memory', 
                             'B cell plasma', 
                             'Class-switched memory B cell',
                             'B cell',
                             'Mast cell', 
                             'Mast cell activated', 
                             'Mast cell resting', 
                             'Eosinophil'))
  
# expression strata graphics ----
  
  ## strata labels and colors
  
  globals$strata_labs <- c(low = 'CXCR4 low', 
                           high = 'CXCR4 high')
  
  globals$strata_colors <- c(low = 'cornflowerblue', 
                             high = 'coral3')
  
# END -----