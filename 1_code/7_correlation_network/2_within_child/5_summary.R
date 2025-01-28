##avoid source
no_function()

setwd(r4projects::get_project_wd())
library(tidyverse)
rm(list = ls())
source("1_code/tools.R")

##load data
###load exposome with phenotye
load("3_data_analysis/correlation_network/within_child/exposome_vs_phenotype/exposome_phenotype_glm")
exposome_phenotype_glm

###load exposome with internal-ome
load("3_data_analysis/correlation_network/within_child/exposome_vs_internal_omics/exposome_internal_omics_cor")

exposome_phenotype = 
exposome_phenotype_glm %>% 
  dplyr::filter(p.adjust < 0.05) %>% 
  dplyr::pull(variable_id) %>% 
  unique()

exposome_internal_ome = 
  exposome_internal_omics_cor %>% 
  dplyr::filter(p.adjust < 0.05) %>% 
  dplyr::pull(from) %>% 
  unique()

length(exposome_phenotype)
length(exposome_internal_ome)

dir.create("3_data_analysis/correlation_network/within_child/9_summary")
setwd("3_data_analysis/correlation_network/within_child/9_summary")

temp = list("Associate with phenotype" = exposome_phenotype,
            "Associate with internal-ome feature" = exposome_internal_ome)

plot = 
  ggvenn(
    data = temp,
    fill_color = c(ggsci::pal_aaas()(n=10)[1],ggsci::pal_aaas()(n=10)[2]),
    stroke_color = "white", set_name_size = 10, text_size = 6, fill_alpha = 0.9, 
    text_color = "white"
  )
plot

ggsave(plot, filename = "associate_with_phenotype_internal_ome.pdf", width = 7, height = 7)




