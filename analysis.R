# This script executes data analysis scripts

  library(plyr)
  library(tidyverse)
  library(exda)
  library(clustTools)
  library(somKernels)
  library(microViz)
  library(kmOptimizer)
  library(survival)
  library(furrr)
  
  c('./tools/project_tools.R') %>%
    source_all(source, message = TRUE, crash = TRUE)

  insert_head()
  
# Executable scripts -----
  
  insert_msg('Deploying the analysis scripts')
  
  c('./analysis scripts/infiltration_clustering.R', ## participant clustering by immune infiltration/Quantiseq
    './analysis scripts/expression_localization.R', ## Spearman correlations
    './analysis scripts/survival.R', ## KM and Cox modeling
    './analysis scripts/immune_function.R', ## Immune function signatures
    './analysis scripts/immune_genes.R') %>% ## Immunity-related genes
    source_all(source, message = TRUE, crash = TRUE)
  
# END -----
  
  insert_tail()