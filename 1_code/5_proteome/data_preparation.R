setwd(r4projects::get_project_wd())
library(tidyverse)
rm(list = ls())

load("2_data/proteome.Rdata")

setwd("3_data_analysis/8_proteome_data_analysis/data_preparation")

expression_data = 
proteome@assayData$exprs

sample_info =
  proteome@phenoData@data %>%
  as.data.frame() %>%
  dplyr::rename(sample_id = ID)

variable_info =
  proteome@featureData@data %>% 
  dplyr::rename(variable_id = Prot_ID) %>% 
  dplyr::mutate(Prot_ID = variable_id)

rownames(variable_info) = NULL

sample_info = 
sample_info %>% 
  dplyr::mutate(subject_id = sample_id) %>% 
  dplyr::mutate(sample_id = paste("child", sample_info$sample_id, sep = "_")) %>% 
  dplyr::select(sample_id, subject_id, everything())

colnames(expression_data) = 
sample_info$sample_id

save(expression_data, file = "expression_data")
save(variable_info, file = "variable_info")
save(sample_info, file = "sample_info")








