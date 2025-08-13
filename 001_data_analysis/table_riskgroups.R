########################################################
########################################################
######## INTERNSHIP PROJECT - table risk groups ########
########################################################
########################################################

## Ciloë Sans (M1 internship project, june-august 2025)

## 1. start:
rm(list = ls()) ## clear working env
cat("\014") ## clear console

## environment to set
output_dir <- 'mypath'
dir.create(output_dir, showWarnings = FALSE)

## import necessary libraries
library(dplyr)
library(ggplot2)
library(tibble)
library(tidyr)
library(xtable)

setwd("/myfilepath/") ## set working directory


## load risk groups (gen info, PRS002013 ONLY)
PGS002013_df <- data.table::fread("/myfilepath/mytable.txt")


## get col names
group_cols <- names(PGS002013_df)[grepl("gr_clrt", names(PGS002013_df))]


## see risk grs repartition
n <- length(group_cols)

sink("table.txt")
for (i in 1:n) {
  cat(group_cols[i], ":", "\n")
  print(table(PGS002013_df[[group_cols[i]]]))
}
sink()

