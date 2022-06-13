# This script provides tools for project data analyses

# tools ----
  
  library(plyr)
  library(tidyverse)

# annotation functions -----

  translate_gene <- function(id, 
                             id_type = 'gene_symbol', 
                             output_type = 'entrez_id', 
                             dictionary = tcga$annotation) {
    
    ## gets gene identifier of interest

    naming_vec <- dictionary[[output_type]] %>% 
      set_names(dictionary[[id_type]])
    
    return(naming_vec[id])
        
  }
  
# Varia -----
  
  mm_inch <- function(input_mm) {
    
    0.0393700787 * input_mm
    
  }
  
# END -----