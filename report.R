# This script executes scripts generating figures, tables and the report pdf

# tools ----

  library(cowplot)
  library(rmarkdown)
  library(knitr)
  library(kableExtra)
  library(bookdown)
  library(writexl)
  library(soucer)
  library(figur)

  insert_head()

# executable scripts ------

  insert_msg('Executing the report parts')
  
  c('./report scripts/report_figures.R', 
    './report scripts/report_tables.R', 
    './report scripts/deploy_report.R') %>% 
    source_all(source, message = TRUE, crash = TRUE)
  
# END ----
  
  insert_tail()