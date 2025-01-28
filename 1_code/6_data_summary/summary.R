setwd(r4projects::get_project_wd())
library(tidyverse)
rm(list = ls())

###load data
###exposome chemical
load("3_data_analysis/4_exposome_chemical_data_analysis/data_preparation/expression_data")
load("3_data_analysis/4_exposome_chemical_data_analysis/data_preparation/sample_info")
load("3_data_analysis/4_exposome_chemical_data_analysis/data_preparation/variable_info")

exposome_chemical_expression_data = expression_data
exposome_chemical_sample_info = sample_info %>% 
  dplyr::mutate(sample_id = as.character(sample_id))
exposome_chemical_variable_info = variable_info

###exposome air
load("3_data_analysis/5_exposome_air_data_analysis/data_preparation/expression_data")
load("3_data_analysis/5_exposome_air_data_analysis/data_preparation/sample_info")
load("3_data_analysis/5_exposome_air_data_analysis/data_preparation/variable_info")

exposome_air_expression_data = expression_data
exposome_air_sample_info = sample_info %>% 
  dplyr::mutate(sample_id = as.character(sample_id))
exposome_air_variable_info = variable_info

###exposome outdoor
load("3_data_analysis/6_exposome_outdoor_data_analysis/data_preparation/expression_data")
load("3_data_analysis/6_exposome_outdoor_data_analysis/data_preparation/sample_info")
load("3_data_analysis/6_exposome_outdoor_data_analysis/data_preparation/variable_info")

exposome_outdoor_expression_data = expression_data
exposome_outdoor_sample_info = sample_info %>% 
  dplyr::mutate(sample_id = as.character(sample_id))
exposome_outdoor_variable_info = variable_info

###exposome outdoor
load("3_data_analysis/7_lifecycle_outdoor_data_analysis/data_preparation/expression_data")
load("3_data_analysis/7_lifecycle_outdoor_data_analysis/data_preparation/sample_info")
load("3_data_analysis/7_lifecycle_outdoor_data_analysis/data_preparation/variable_info")

exposome_lifecycle_expression_data = expression_data
exposome_lifecycle_sample_info = sample_info %>% 
  dplyr::mutate(sample_id = as.character(sample_id))
exposome_lifecycle_variable_info = variable_info

####plasma proteome
load("3_data_analysis/8_proteome_data_analysis/data_preparation/expression_data")
load("3_data_analysis/8_proteome_data_analysis/data_preparation/sample_info")
load("3_data_analysis/8_proteome_data_analysis/data_preparation/variable_info")

proteome_expression_data = expression_data
proteome_sample_info = sample_info
proteome_variable_info = variable_info

###urine_metabolome
load("3_data_analysis/3_urine_metabolome_data_analysis/data_preparation/expression_data")
load("3_data_analysis/3_urine_metabolome_data_analysis/data_preparation/sample_info")
load("3_data_analysis/3_urine_metabolome_data_analysis/data_preparation/variable_info")

urine_metabolome_expression_data = expression_data
urine_metabolome_sample_info = sample_info
urine_metabolome_variable_info = variable_info

###serum_metabolome
load("3_data_analysis/2_serum_metabolome_data_analysis/data_preparation/expression_data")
load("3_data_analysis/2_serum_metabolome_data_analysis/data_preparation/sample_info")
load("3_data_analysis/2_serum_metabolome_data_analysis/data_preparation/variable_info")

serum_metabolome_expression_data = expression_data
serum_metabolome_sample_info = sample_info
serum_metabolome_variable_info = variable_info

###transcriptome
load("3_data_analysis/1_transcriptome_data_analysis/data_preparation/expression_data")
load("3_data_analysis/1_transcriptome_data_analysis/data_preparation/sample_info")
load("3_data_analysis/1_transcriptome_data_analysis/data_preparation/variable_info")

transcriptome_expression_data = expression_data
transcriptome_sample_info = sample_info
transcriptome_variable_info = variable_info

setwd("3_data_analysis/9_summary")


intersect(proteome_sample_info$sample_id,
          exposome_sample_info$subject_id)

setdiff(exposome_sample_info$subject_id, proteome_sample_info$sample_id)
setdiff(proteome_sample_info$sample_id, exposome_sample_info$subject_id)

temp = 
  proteome_sample_info %>% 
  dplyr::left_join(exposome_sample_info, by = "sample_id")




