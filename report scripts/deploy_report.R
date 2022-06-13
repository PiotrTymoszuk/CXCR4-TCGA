# This script renders the pdf report file

  insert_head()
  
# data container ----
  
  report_launch <- list()

# rendering the report -----
  
  insert_msg('Rendering the report')
  
  render('./report/report.Rmd', 
         output_format = pdf_document2(number_sections = F)) 
  
# END -----
  
  insert_tail()