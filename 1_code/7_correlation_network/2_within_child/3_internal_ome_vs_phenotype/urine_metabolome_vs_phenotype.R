##avoid source
# no_function()

setwd(r4projects::get_project_wd())
library(tidyverse)
rm(list = ls())
source("1_code/100_tools.R")

##load data
###child exposome outdoor
load(
  "3_data_analysis/6_exposome_outdoor_data_analysis/data_preparation/expression_data"
)
load("3_data_analysis/6_exposome_outdoor_data_analysis/data_preparation/sample_info")
load("3_data_analysis/6_exposome_outdoor_data_analysis/data_preparation/variable_info")

exposome_outdoor_variable_info <-
  variable_info

exposome_outdoor_sample_info =
  sample_info %>%
  dplyr::filter(stringr::str_detect(sample_id, pattern = "child"))

exposome_outdoor_expression_data =
  expression_data[, exposome_outdoor_sample_info$sample_id]

remain_idx = 
  exposome_outdoor_expression_data %>%
  apply(1, function(x) {
    sum(is.na(x))
  }) %>% 
  `==`(0) %>% 
  which()

exposome_outdoor_expression_data = exposome_outdoor_expression_data[remain_idx,]
exposome_outdoor_variable_info = exposome_outdoor_variable_info[remain_idx,]

##load data
###child urine_metabolome
load(
  "3_data_analysis/3_urine_metabolome_data_analysis/data_preparation/expression_data"
)
load("3_data_analysis/3_urine_metabolome_data_analysis/data_preparation/sample_info")
load("3_data_analysis/3_urine_metabolome_data_analysis/data_preparation/variable_info")

urine_metabolome_variable_info <-
  variable_info

urine_metabolome_sample_info =
  sample_info %>%
  dplyr::filter(stringr::str_detect(sample_id, pattern = "child"))

urine_metabolome_expression_data =
  expression_data[, urine_metabolome_sample_info$sample_id]

urine_metabolome_expression_data %>%
  apply(1, function(x) {
    sum(is.na(x))
  })

####only remain the overlapped samples
sample_id =
  intersect(exposome_outdoor_sample_info$sample_id,
            urine_metabolome_sample_info$sample_id)

exposome_outdoor_expression_data = 
  exposome_outdoor_expression_data[,sample_id]

exposome_outdoor_sample_info = 
  exposome_outdoor_sample_info[match(sample_id, exposome_outdoor_sample_info$sample_id),]

urine_metabolome_expression_data = 
  urine_metabolome_expression_data[,sample_id]

urine_metabolome_sample_info = 
  urine_metabolome_sample_info[match(sample_id, urine_metabolome_sample_info$sample_id),]

exposome_outdoor_sample_info$sample_id == urine_metabolome_sample_info$sample_id


setwd(r4projects::get_project_wd())
dir.create("3_data_analysis/correlation_network/within_child/urine_metabolome_vs_phenotype")
setwd("3_data_analysis/correlation_network/within_child/urine_metabolome_vs_phenotype")

#####glm
urine_metabolome_phenotype_glm = 
exposome_outdoor_sample_info[,c("Behavior", "Body.mass.index.z.score", "Intelligence.quotient")] %>% 
  purrr::map(function(x){
    cat(x[1], " ")
    temp_p = 
    urine_metabolome_expression_data %>% 
      t() %>% 
      as.data.frame() %>% 
      purrr::map(function(y){
        #x is the phenotype
        #y is the exposome
            temp_data =
              data.frame(x = x,
                         y = y,
                         exposome_outdoor_sample_info)
            glm_reg =
              glm(
                x ~ y + Child.sex + Year.of.birth + Maternal.BMI + Gestational.age.at.birth +
                  Maternal.age + Child.height + Child.weight + Birthweight,
                family = gaussian,
                temp_data
              )
            
            temp =
              summary(glm_reg)$coefficients %>%
              as.data.frame()
            temp$`Pr(>|t|)`[2]
      }) %>% 
      unlist()
  }) %>% 
  do.call(cbind, .) %>% 
  as.data.frame()

urine_metabolome_phenotype_glm

save(urine_metabolome_phenotype_glm, file = "urine_metabolome_phenotype_glm")
