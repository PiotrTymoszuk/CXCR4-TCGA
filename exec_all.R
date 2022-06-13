# This program executes the scripts of the CXCR4 project

# libraries -----

  library(soucer)

  print(source_all(c('import.R',
                     'analysis.R', 
                     'dge.R', 
                     'report.R'), 
                   message = TRUE, crash = FALSE))

  save.image()
  
# END ----