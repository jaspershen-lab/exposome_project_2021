setwd(r4projects::get_project_wd())
library(tidyverse)
rm(list = ls())

load("2_data/genexpr.Rdata")

setwd("3_data_analysis/1_transcriptome_data_analysis/data_preparation")

expression_data = 
genexpr@assayData$exprs

sample_info =
  genexpr@phenoData@data %>%
  as.data.frame() %>%
  dplyr::rename(sample_id = ID)

variable_info =
  genexpr@featureData@data %>% 
  tibble::rownames_to_column(var = "variable_id")

rownames(variable_info) = NULL

sample_info = 
  sample_info %>% 
  dplyr::mutate(subject_id = sample_id) %>% 
  dplyr::mutate(sample_id = paste("child", sample_info$sample_id, sep = "_")) %>% 
  dplyr::select(sample_id, subject_id, everything())

colnames(expression_data) = sample_info$sample_id

save(expression_data, file = "expression_data")
save(variable_info, file = "variable_info")
save(sample_info, file = "sample_info")








