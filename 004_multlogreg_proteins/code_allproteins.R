#########################################################################
#########################################################################
######## INTENSHIP PROJECT - proteins (all) code (FINAL VERSION) ########
#########################################################################
#########################################################################

## Ciloë Sans (M1 internship project, june-august 2025)

## 1. start:
rm(list = ls()) ## clear working env
cat("\014") ## clear console
set.seed(821)

# environment ====
output_dir <- 'myfilepath/'
dir.create(output_dir, showWarnings = FALSE)

## import necessary libraries
library(doParallel)
library(dplyr)
library(foreach)
library(ggplot2)
library(glmnet)
library(haven)
library(mgcv)
library(survival)
library(tidyr)
library(xtable)

source("/myfilepath")

# data ====
PRS_df <- data.table::fread("/myfilepath/mydata.txt")
PGS002013_df <- data.table::fread("/myfilepath/mydata2.txt")
data <- left_join(PGS002013_df, PRS_df, by = "Idepic")
# Identify columns
group_cols <- names(data)[grepl("gr_clrt", names(data))]

# parallel set-up ====
total_cores <- 10
group_parallel_cores <- min(length(group_cols), total_cores)
outer_cl <- parallel::makeCluster(group_parallel_cores)
doParallel::registerDoParallel(outer_cl)

# analysis: parallel ====
cat("\n--- START ---\n")
foreach(group_col = group_cols, .packages = c("dplyr", "glmnet", "caret", "parallel", "doParallel"), .export = "perform_case_cohort_anova") %dopar% {
  
  # data
  cat("\n--- data format ---\n")
  data <- data |> 
    dplyr::filter(!is.na(.data[[group_col]]))
  feature_cols <- names(data)[grepl("seq", names(data))]
  feature_data <- data[, ..feature_cols]
  risk_groups <- data[[group_col]]
  weights <- data$weights_somalogic_genetic
  
  # Step 2: Prepare data
  cat("\n--- data prep ---\n")
  X_top <- as.matrix(feature_data)
  y_top <- as.factor(risk_groups)
  
  complete_samples <- complete.cases(X_top, y_top)
  X_top <- X_top[complete_samples, ]
  y_top <- y_top[complete_samples]
  weights <- weights[complete_samples]
  
  # Step 3: Cross-validation folds
  cat("\n--- folds ---\n")
  NFOLDS <- 10
  folds <- caret::createFolds(y_top, k = NFOLDS, list = TRUE, returnTrain = FALSE)
  foldid_vector <- rep(NA, length(y_top))
  for (i in 1:NFOLDS) foldid_vector[folds[[i]]] <- i
  
  # Step 4: Inner parallel cluster for cv.glmnet
  cat("\n--- parallel ---\n")
  inner_cores <- max(1, floor(total_cores / group_parallel_cores))
  inner_cl <- parallel::makeCluster(inner_cores)
  doParallel::registerDoParallel(inner_cl)
  
  # Step 5: glmnet
  cat("\n--- glmnet ---\n")
  model_glm <- glmnet::cv.glmnet(x = X_top, y = y_top,
                                 family = "multinomial",
                                 alpha = 1,
                                 type.multinomial = "ungrouped",
                                 foldid = foldid_vector,
                                 weights = weights,
                                 parallel = TRUE)
  
  parallel::stopCluster(inner_cl)
  
  # Save
  model_file <- paste0(output_dir, group_col, ".rds")
  saveRDS(model_glm, file = model_file)
  return(NULL)
}

# Stop outer parallel cluster
parallel::stopCluster(outer_cl)
cat("\n--- END ---\n")



