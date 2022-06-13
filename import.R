# This script imports TCGA RNA Seq data together with the immune inflitration
# estimates provided by TIMER

# toolbox ----

  library(plyr)
  library(tidyverse)
  library(stringi)
  library(soucer)
  library(gseaTools)

  source_all('./tools/project_tools.R', message = TRUE, crash = TRUE)
  
  insert_head()

# data containers -----

  timer_data <- list()
  tcga_expr <- list()
  tcga_clinics <- list()
  
  tcga <- list() ## for the final data container

# reading the Timer immune estimate table for all TCGA tumors -----
  
  insert_msg('Reading the timer data')
  
  timer_data$est_table <- read_csv('./input data/infiltration_estimation_for_tcga.csv.gz') %>% 
    mutate(sample_id = cell_type) %>% 
    select( - cell_type)

# downloading TCGA clinical and expression data if not present on the disc ------
  
  insert_msg('Reading the TCGA clinical and expression data')
  
  if(length(list.files('./input data/TCGA', pattern = 'PAAD__gene*')) == 0) {
    
    source('./data import scripts/data_download_TCGA.R')
    
    stop('Copy the expression and clinical files into ./input data/TCGA and run the pipeline again')
    
  }
  
# clearing the clinical and expression data data ------
  
  insert_msg('Clearning the clinical and expression data')
  
  c('./data import scripts/data_clearing_tcga_clinics.R', 
    './data import scripts/data_clearing_tcga_expression.R') %>% 
    source_all(source, message = TRUE, crash = TRUE)
  
# merging the expression data with the essential set of clinical information -----
  
  insert_msg('Appending the expression data with the essential clinical information')
  
  tcga[c('clinical', 
         'drugs', 
         'radiation')] <- tcga_clinics[c('clinical', 
                                         'drugs', 
                                         'radiation')]
  
  tcga[c('expression', 
         'annotation')] <- tcga_expr[c('exprs', 
                                       'annotation')]
  
  tcga$expression <- right_join(tcga$clinical, 
                                tcga$expression, 
                                by = 'patient_id')
  
# merging the timer estimates with the expression data set ----
  
  insert_msg('Mergining the Timer estimates with the expression data set')
  
  tcga$expression <- left_join(tcga$expression, 
                               timer_data$est_table, 
                               by = 'sample_id')
  
# globals setup -----
  
  insert_msg('Globals setup')
  
  source('./tools/project_globals.R')

# calculation of the signatures with GSVA, appending the expression table -----
  
  insert_msg('Calculation of the gene signatures with GSVA')
  
  source('./data import scripts/signature_calculation.R')
  
  tcga$expression <- left_join(tcga$expression, 
                               gene_sign$signature_est, 
                               by = 'sample_id') %>% 
    left_join(gene_sign$cell_est, 
              by = 'sample_id')
  
# restricting the sample set to the the ductal adenocarcinoma samples ------
  
  insert_msg('Ductal carcinoma samples only')
  
  tcga$pdac_ids <- tcga$clinical %>% 
    filter(histology == 'Adeno ductal') %>% 
    .$patient_id
  
  tcga[c('clinical', 
         'drugs', 
         'radiation', 
         'expression')] <- tcga[c('clinical', 
                                  'drugs', 
                                  'radiation', 
                                  'expression')] %>% 
    map(filter, 
        patient_id %in% tcga$pdac_id)
  
# END ----
  
  #rm(tcga_clinics, 
    # tcga_expr, 
    # timer_data)
  
  insert_tail()