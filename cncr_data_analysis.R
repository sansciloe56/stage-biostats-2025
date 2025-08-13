################################################
################################################
###### INTERNSHIP PROJECT - data analyses ######
################################################
################################################

## Ciloë Sans (M1 internship project, june-august 2025)

## 1. start:
rm(list = ls()) ## clear working env
cat("\014") ## clear console

output_dir <- 'mypath'
dir.create(output_dir, showWarnings = FALSE)

## import necessary libraries
library(dplyr)
library(ggplot2)
library(glmnet)
library(haven)
library(mgcv)
library(survival)
library(tidyr)
library(xtable)

setwd("/myfilepath/") ## set working directory

## PRS data
PRS_data <- data.table::fread("/myfilepath/mydata.txt")

## select only v few vars of the PRS data (PRS scores specifically & EPIC ID to be able to fuse somalogic & PRS databases)
PRS_data_selected <- PRS_data[, c("idepic", "PGS002013_hmPOS_GRCh37", "PGS003850_hmPOS_GRCh37")]
PRS_data_selected$Idepic <- PRS_data_selected$idepic


## 2. frequency table cancer cases
## for the data with PGS
cancer_data_prs <- PRS_data[, c("cncr_mal_anyc", "cncr_mal_blad", "cncr_mal_brea", "cncr_mal_clrt", "cncr_mal_clrt_colon", 
                                "cncr_mal_clrt_rectum", "cncr_mal_coru", "cncr_mal_glio", "cncr_mal_kidn", "cncr_mal_live", 
                                "cncr_mal_lung", "cncr_mal_lymp", "cncr_mal_mela", "cncr_mal_ovar", "cncr_mal_panc", "cncr_mal_pros", 
                                "cncr_mal_stom", "cncr_mal_thyr", "cncr_mal_uadt")]

cancer_names_vars_prs <- c("cncr_mal_anyc", "cncr_mal_blad", "cncr_mal_brea", "cncr_mal_clrt", "cncr_mal_clrt_colon", 
                           "cncr_mal_clrt_rectum", "cncr_mal_coru", "cncr_mal_glio", "cncr_mal_kidn", "cncr_mal_live", 
                           "cncr_mal_lung", "cncr_mal_lymp", "cncr_mal_mela", "cncr_mal_ovar", "cncr_mal_panc", "cncr_mal_pros", 
                           "cncr_mal_stom", "cncr_mal_thyr", "cncr_mal_uadt")


cncr_freq_prs <- matrix(0, nrow = length(cancer_names_vars_prs), ncol = 2)

df_cncr_freq_prs <- as.data.frame(cncr_freq_prs, row.names = cancer_names_vars_prs)

colnames(df_cncr_freq_prs) <- c("Non-case", "Incident")

for (k in cancer_names_vars_prs) {
  df_cncr_freq_prs[k, 1] <- sum(cancer_data_prs[, ..k] == 0)
  df_cncr_freq_prs[k, 2] <- sum(cancer_data_prs[, ..k] == 1)
}

xtable(df_cncr_freq_prs) ## shows the nb of cancer cases per individual

## sexe individu
table(PRS_data$sex)

## make confusion matrix btw genetic data and cancer type (colon-rectum specifically)
## 1st make a 0/1 var for whether SNPs are available
PRS_data$SNPs_binary <- ifelse(is.na(PRS_data$PGS002013_hmPOS_GRCh37), 0, 1) ## if NA = 0, if info present = 1

length(which(is.na(PRS_data$PGS002013_hmPOS_GRCh37)))

table(PRS_data$SNPs_binary, PRS_data$sex)
table(PRS_data$SNPs_binary, PRS_data$cncr_mal_anyc)
table(PRS_data$SNPs_binary, PRS_data$cncr_mal_clrt)


m <- length(cancer_names_vars_prs)
genetic_info_cncr_types <- matrix(0, nrow = m, ncol = 2)

df_genetic_info_cncr_types <- as.data.frame(genetic_info_cncr_types, row.names = cancer_names_vars_prs)
colnames(df_genetic_info_cncr_types) <- c("Pas de cancer", "Cancer")

for (i in 1:2) {
  for (j in 1:m) {
    df_genetic_info_cncr_types[j, i] <- table(cancer_data_prs[[cancer_names_vars_prs[j]]], PRS_data[, SNPs_binary])[2, i]
  }
}

xtable(df_genetic_info_cncr_types)


## 3. barplot cancer cases:
PRS_data_long <- pivot_longer(PRS_data, cols = c("cncr_mal_blad", "cncr_mal_brea", "cncr_mal_clrt", "cncr_mal_clrt_colon", "cncr_mal_clrt_rectum", "cncr_mal_coru", "cncr_mal_glio", "cncr_mal_kidn", "cncr_mal_live", "cncr_mal_lung", "cncr_mal_lymp", "cncr_mal_mela", "cncr_mal_ovar", "cncr_mal_panc", "cncr_mal_pros", "cncr_mal_stom", "cncr_mal_thyr", "cncr_mal_uadt"), names_to = "cancer_type", values_to = "frequency")
PRS_data_long_updated <- PRS_data_long %>% filter(frequency == 1)

cancer_names <- c("Bladder", "Breast", "Colon-Rectum", "Colon", "Rectum", "Corpus uteri", "Glioma", "Kidney", "Liver", "Lung", "Lymphoma", "Melanoma", "Ovary", "Pancreas", "Prostate", "Stomach", "Thyroid", "Upper aero-\ndigestive tract")

cncrcases_barplot <- ggplot(PRS_data_long_updated, aes(x = cancer_type, fill = cancer_type)) +
  geom_bar(show.legend = FALSE) +
  scale_x_discrete(name = "Cancer type", labels = cancer_names) + 
  scale_y_continuous(name = "Nb of cases") + 
  theme(axis.text.x = element_text(face = "plain", size = 8, angle = 30), title = element_text(size = 16)) + 
  labs(title = "Cancer types and cases", size = 20)
ggsave("cncrcases_barplot.tiff", plot = cncrcases_barplot, width = 12, height = 8, dpi = 300)

