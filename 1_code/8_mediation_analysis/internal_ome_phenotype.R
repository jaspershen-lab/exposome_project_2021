##avoid source
# no_function()

setwd(r4projects::get_project_wd())
library(tidyverse)
rm(list = ls())

##load data
####exposome
load("3_data_analysis/mediation_analysis/data_preparation/exposome_expression_data")
load("3_data_analysis/mediation_analysis/data_preparation/exposome_sample_info")
load("3_data_analysis/mediation_analysis/data_preparation/exposome_variable_info")

####internal_ome
load("3_data_analysis/mediation_analysis/data_preparation/internal_ome_expression_data")
load("3_data_analysis/mediation_analysis/data_preparation/internal_ome_sample_info")
load("3_data_analysis/mediation_analysis/data_preparation/internal_ome_variable_info")

#####use the linear mixed model to find which exposome have association with phenotye
library(lme4)
library(lmerTest)

dir.create("3_data_analysis/mediation_analysis/internal_ome_phenotype", showWarnings = FALSE)
setwd("3_data_analysis/mediation_analysis/internal_ome_phenotype")

# Find intersected samples between exposome and internal omics data
intersect_name <- intersect(colnames(exposome_expression_data), colnames(internal_ome_expression_data))

# Subset expression data to intersected samples
exposome_expression_data <- exposome_expression_data[, intersect_name]
internal_ome_expression_data <- internal_ome_expression_data[, intersect_name]

# Subset sample info to intersected samples (matching order)
exposome_sample_info <- exposome_sample_info[match(intersect_name, exposome_sample_info$sample_id), ]
internal_ome_sample_info <- internal_ome_sample_info[match(intersect_name, internal_ome_sample_info$sample_id), ]

# Initialize parallel plan for faster computation
plan(multisession, workers = availableCores() - 1) 

# Define phenotypes of interest
phenotypes <- c("Behavior", "Body.mass.index.z.score", "Intelligence.quotient")

# Run GLM for each internal_ome variable against each phenotype with covariates
internal_phenotype_ome_glm_p <- internal_ome_expression_data %>%
  t() %>%  # transpose so variables are rows
  as.data.frame() %>%
  furrr::future_map(function(x) {
    p_vals <- exposome_sample_info[, phenotypes] %>%
      purrr::map_dbl(function(y) {
        temp_data <- data.frame(
          x = x,
          y = y,
          exposome_sample_info
        )
        glm_reg <- glm(
          y ~ x + Child.sex + Year.of.birth + Maternal.BMI + Gestational.age.at.birth +
            Maternal.age + Child.height + Child.weight + Birthweight,
          family = gaussian,
          data = temp_data
        )
        # Extract p-value for x (second coefficient)
        summary(glm_reg)$coefficients["x", "Pr(>|t|)"]
      })
    # Adjust p-values for multiple testing (BH) across the three phenotypes
    p.adjust(p_vals, method = "BH")
  }) %>%
  do.call(rbind, .) %>%
  as.data.frame()

# Assign column names for phenotypes
colnames(internal_phenotype_ome_glm_p) <- phenotypes

# Optional: inspect number of significant variables per phenotype (FDR < 0.05)
significant_counts <- colSums(internal_phenotype_ome_glm_p < 0.05)
print(significant_counts)

# Save the p-value matrix for downstream use
save(internal_phenotype_ome_glm_p, file = "internal_phenotype_ome_glm_p.RData")

exposome_phenotype_glm <- exposome_expression_data %>%
  t() %>%
  as.data.frame() %>%
  furrr::future_map(function(x) {
    p_vals <- exposome_sample_info[, phenotypes] %>%
      purrr::map_dbl(function(y) {
        temp_data <- data.frame(
          x = x,
          y = y,
          exposome_sample_info
        )
        glm_reg <- glm(
          y ~ x + Child.sex + Year.of.birth + Maternal.BMI + Gestational.age.at.birth +
            Maternal.age + Child.height + Child.weight + Birthweight,
          family = gaussian,
          data = temp_data
        )
        summary(glm_reg)$coefficients["x", "Pr(>|t|)"]
      })
    p.adjust(p_vals, method = "BH")
  }) %>%
  do.call(rbind, .) %>%
  as.data.frame()

colnames(exposome_phenotype_glm) <- phenotypes

rownames(exposome_phenotype_glm) <- exposome_variable_info$variable_id

rownames(internal_phenotype_ome_glm_p) <- internal_ome_variable_info$variable_id

# Extract significant variable IDs for each phenotype (FDR < 0.05)
sig_behavior_vars <- rownames(exposome_phenotype_glm)[which(exposome_phenotype_glm$Behavior < 0.05)]
sig_bmi_vars <- rownames(exposome_phenotype_glm)[which(exposome_phenotype_glm$`Body.mass.index.z.score` < 0.05)]
sig_iq_vars <- rownames(exposome_phenotype_glm)[which(exposome_phenotype_glm$Intelligence.quotient < 0.05)]

internal_sig_counts <- colSums(internal_phenotype_ome_glm_p < 0.05)
cat("Significant internal omics variables:\n")
cat("Behavior:", internal_sig_counts["Behavior"], "\n")
cat("BMI:", internal_sig_counts["Body.mass.index.z.score"], "\n")
cat("IQ:", internal_sig_counts["Intelligence.quotient"], "\n")

save(
  internal_phenotype_ome_glm_p,
  exposome_phenotype_glm,
  sig_behavior_vars, sig_bmi_vars, sig_iq_vars,
  file = "significant_vars_and_pvals.RData"
)




# ###pie to show the how many exposome variables are associated with outcome
data.frame(class = c("significant", "no"),
           number = c(
             sum(exposome_phenotype_glm$Behavior < 0.05),
             nrow(exposome_variable_info) - sum(exposome_phenotype_glm$Behavior < 0.05)
           )) %>%
  ggplot(aes(x = 2, y = number, fill = class)) +
  geom_bar(stat = "identity", color = "white") +
  coord_polar("y", start = 0) +
  theme_void() +
  scale_fill_manual(values = c(
    "significant" = "#1f77b4", 
    "no" = "#d3d3d3"
  )) +
  xlim(0.5, 2.5)




















