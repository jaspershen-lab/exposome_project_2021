##avoid source
no_function()

setwd(r4projects::get_project_wd())
library(tidyverse)
rm(list = ls())

## Load data
#### Exposome
load("3_data_analysis/mediation_analysis/data_preparation/exposome_expression_data")
load("3_data_analysis/mediation_analysis/data_preparation/exposome_sample_info")
load("3_data_analysis/mediation_analysis/data_preparation/exposome_variable_info")

#### Internal_ome
load("3_data_analysis/mediation_analysis/data_preparation/internal_ome_expression_data")
load("3_data_analysis/mediation_analysis/data_preparation/internal_ome_sample_info")
load("3_data_analysis/mediation_analysis/data_preparation/internal_ome_variable_info")

## Set up output directory
library(lme4)
library(lmerTest)
dir.create("3_data_analysis/mediation_analysis/exposome_internal_ome", showWarnings = FALSE)
setwd("3_data_analysis/mediation_analysis/exposome_internal_ome")

## Match sample columns across datasets
intersect_name <- intersect(colnames(exposome_expression_data),
                            colnames(internal_ome_expression_data))

exposome_expression_data <- exposome_expression_data[, intersect_name]
internal_ome_expression_data <- internal_ome_expression_data[, intersect_name]

exposome_sample_info <- exposome_sample_info[match(intersect_name, exposome_sample_info$sample_id), ]
internal_ome_sample_info <- internal_ome_sample_info[match(intersect_name, internal_ome_sample_info$sample_id), ]

## Parallel setup
library(future)
library(furrr)
plan("multicore") # On Windows, use plan("multisession") instead

## Analysis
phenotypes <- c("Behavior", "Body.mass.index.z.score", "Intelligence.quotient")
exposome_phenotype_glm <- list()

# Ensure rownames for exposome data
if (is.null(rownames(exposome_expression_data))) {
  rownames(exposome_expression_data) <- exposome_variable_info$variable_id
}

# Loop over phenotypes
for (pheno in phenotypes) {
  message("Processing ", pheno, "...")
  
  pvals <- future_map_dbl(1:nrow(exposome_expression_data), function(i) {
    y <- exposome_expression_data[i, ]
    temp_data <- data.frame(
      x = as.numeric(exposome_sample_info[[pheno]]),
      y = as.numeric(y),
      exposome_sample_info
    )
    
    fit <- try(glm(
      x ~ y + Child.sex + Year.of.birth + Maternal.BMI + Gestational.age.at.birth +
        Maternal.age + Child.height + Child.weight + Birthweight,
      data = temp_data,
      family = gaussian()
    ), silent = TRUE)
    
    if (inherits(fit, "try-error")) return(NA)
    summary(fit)$coefficients["y", "Pr(>|t|)"]
  })
  
  adj_pvals <- p.adjust(pvals, method = "BH")
  exposome_phenotype_glm[[paste0(pheno, "_glm_p.adjust")]] <- adj_pvals
  
  # Report significant variables
  sig_ids <- rownames(exposome_expression_data)[which(adj_pvals < 0.05)]
  message(length(sig_ids), " significant exposome variables for ", pheno)
}

# Save results
exposome_phenotype_glm <- as.data.frame(exposome_phenotype_glm)
exposome_phenotype_glm$variable_id <- rownames(exposome_expression_data)
save(exposome_phenotype_glm, file = "exposome_phenotype_glm.Rdata")

# Reload and check significant variables
load("exposome_phenotype_glm.Rdata")

exposome_phenotype_glm$variable_id[which(exposome_phenotype_glm$Body.mass.index.z.score_glm_p.adjust < 0.05)]
exposome_phenotype_glm$variable_id[which(exposome_phenotype_glm$Intelligence.quotient_glm_p.adjust < 0.05)]


###pie to show the how many exposome variables are associated with outcome
data.frame(class = c("significant", "no"),
           number = c(
             sum(exposome_phenotype_glm$behavior_glm_p.adjust < 0.05),
             nrow(exposome_variable_info) - sum(exposome_phenotype_glm$behavior_glm_p.adjust < 0.05)
           )) %>%
  ggplot(aes(x = 2, y = number, fill = class)) +
  geom_bar(stat = "identity") +
  coord_polar("y", start = 200)
  theme_void() +
  scale_fill_manual(values = c(
    significant
  )) +
  xlim(.2, 2.5)




















