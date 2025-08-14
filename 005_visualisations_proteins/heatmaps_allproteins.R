#########################################################
#########################################################
######## INTERNSHIP PROJECT - plots proteins ALL ########
#########################################################
#########################################################

## Ciloë Sans (M1 internship project, june-august 2025)

## 1. start:
rm(list = ls()) ## clear working env
cat("\014") ## clear console

## environment to set
output_dir <- 'myfilepath'
dir.create(output_dir, showWarnings = FALSE)

## import necessary libraries
library(dplyr)
library(ggplot2)
library(gplots)
library(glmnet)
library(haven)
library(mgcv)
library(nnet)
library(pheatmap)
library(reshape2)
library(survival)
library(tibble)
library(tidyr)
library(xtable)

setwd("/myfilepath/") ## set working directory


## load data
PRS_df <- data.table::fread("/myfilepath/mydata.txt")

proteins_info_file <- data.table::fread("/myfilepath/mydata2.txt")

## load risk groups (gen info, PRS002013 ONLY)
PGS002013_df <- data.table::fread("/myfilepath/mydata3.txt")

## fuse both databases together
df_ALL <- left_join(PGS002013_df, PRS_df, by = "Idepic")

## load model data
model_glm <- readRDS("/myfilepath/mymodel.rds")

plot(model_glm) ## check if model did a strong shrinkage or not => red curve should have a V/U-shape
ggsave("model_glm_fitplot.tiff", height = 16, width = 10, dpi = 300)

## see risk grs repartition
table(PGS002013_df$risk_gr_clrt_cncr_PGS_restr_90)
table(PGS002013_df$risk_gr_clrt_cncr_PGS_restr_10)

## save coefs GLM model
coefs_glm <- coef(model_glm, s = "lambda.1se")

## extract proteins names
proteins_varnames <- grep("seq", names(df_ALL), value = TRUE)

## create class names vector
crc_classes_names <- c("Low_risk//No_cancer", "Low_risk//Cancer", "High_risk//No_cancer", "High_risk//Cancer")  ## replace with lvl names

## select non-0 coefs from GLM model
proteins_nonzero <- lapply(crc_classes_names, function(cls) {
  M <- coefs_glm[[cls]]
  if(is.null(M)) return(NULL)
  
  mat <- as.matrix(M) ## convert to matrix
  prot_of_int <- intersect(rownames(mat), proteins_varnames)
  
  if (length(prot_of_int) == 0) return(NULL)
  
  sub_coefs <- mat[prot_of_int, , drop = FALSE]
  idx <- which(sub_coefs != 0)
  
  if (length(idx) == 0) return(NULL)
  
  data.frame(
    protein = rownames(sub_coefs)[idx],
    coef = as.numeric(sub_coefs[idx]),
    group = cls,
    stringsAsFactors = FALSE
  )
})

crc_df <- do.call(rbind, proteins_nonzero)


## -> top-20 proteins per cncr group (incl common ones across groups)
if(is.null(crc_df)){
  top_proteins <- NULL
} else {
  top_proteins <- crc_df %>%
    group_by(group) %>%
    arrange(desc(abs(coef)))
  
  ## common proteins
  common <- Reduce(intersect, lapply(proteins_nonzero, function(df){df$protein}))
  
  print(common)
  
  top_proteins$common_proteins <- top_proteins$protein %in% common
  
  
  ## -> plot proteins to see which ones are common across all groups
  ggplot(top_proteins, aes(x = reorder(protein, abs(coef)), y = coef, fill = common_proteins)) +
    geom_col(show.legend = TRUE) +
    coord_flip() +
    facet_wrap(~ group, scales = "free_y") +
    labs(
      title = "Top proteins cancer development",
      x = "Protein name",
      y = "Coefficient"
    ) +
    scale_fill_manual(values = c("TRUE" = "blue", "FALSE"= "darkgrey")) + 
    theme_minimal()
  ggsave("plot_top_10.tiff", width = 34, height = 26, dpi = 300)
  
  
  
  ## ->  heatmap not shared proteins across all groups
  unique_proteins <- setdiff(unique(top_proteins$protein), common)
  
  alone_proteins <- top_proteins %>%
    group_by(protein) %>%
    filter(n_distinct(group) == 1) %>%
    ungroup()
  
  ## reshape data to wide format
  heatmap_data <- alone_proteins %>%
    select(protein, group, coef) %>%
    pivot_wider(names_from = group, values_from = coef, values_fill = 0) %>%
    column_to_rownames("protein") %>%
    as.matrix()
  
  ## get names & replace proteins' IDs by their full names
  matched_names <- proteins_info_file$EntrezGeneSymbol[match(rownames(heatmap_data), proteins_info_file$AptName)]
  rownames(heatmap_data) <- ifelse(is.na(matched_names), rownames(heatmap_data), matched_names)
  
  ## reverse to get proteins' names on x-axis instead of y-axis
  heatmap_data_scaled <- t(heatmap_data)
  
  ## define max val to use for the limits in the color palette (& to then center it to 0)
  palette_col <- colorRampPalette(c("#006837", "#a6d96a", "white", "#f46d43", "#a50026"))(59)
  palette_lim <- seq(-0.5, 0.5, length.out = (length(palette_col) + 1))
  
  ## plot heatmap
  heatmap_plot <- pheatmap::pheatmap(heatmap_data_scaled,
                                     cluster_rows = FALSE,
                                     cluster_cols = FALSE,
                                     cellheight = 40,
                                     angle_col = "90",
                                     fontsize_col = 8,
                                     color = palette_col,
                                     breaks = palette_lim,
                                     border_color = "lightgrey",#"transparent",
                                     main = "Non-shared proteins (among top proteins) across risk groups")
  ggsave("heatmap_notshared_prots.tiff", plot = heatmap_plot, width = 16, height = 10, dpi = 300)
  
  
  ## heatmap common protein to see how it varies across the 4 groups
  shared_proteins <- crc_df[crc_df$protein %in% common, ]
  
  ## reshape data to wide format
  heatmap_data_common <- shared_proteins %>%
    select(protein, group, coef) %>%
    pivot_wider(names_from = group, values_from = coef, values_fill = 0) %>%
    column_to_rownames("protein") %>%
    as.matrix()
  
  ## get names & replace proteins' IDs by their full names
  matched_names_common <- proteins_info_file$EntrezGeneSymbol[match(rownames(heatmap_data_common), proteins_info_file$AptName)]
  rownames(heatmap_data_common) <- ifelse(is.na(matched_names_common), rownames(heatmap_data_common), matched_names_common)
  
  ## reverse to get proteins' names on x-axis instead of y-axis
  heatmap_data_common_scaled <- t(heatmap_data_common)
  
  ## plot heatmap
  heatmap_plot_shared <- pheatmap::pheatmap(heatmap_data_common_scaled,
                                            cluster_rows = FALSE,
                                            cluster_cols = FALSE,
                                            cellheight = 80,
                                            cellwidth = 100,
                                            angle_col = "90",
                                            color = palette_col,
                                            breaks = palette_lim,
                                            border_color = "lightgrey",#"transparent",
                                            main = "Common protein(s) across risk groups")
  ggsave("heatmap_shared_prots.tiff", plot = heatmap_plot_shared, width = 16, height = 10, dpi = 300)
  
  
  ## heatmap common protein to see how it varies across the 4 groups
  table_prots <- table(top_proteins$protein)
  common_2grs <- names(table_prots)[table_prots == 2]
  shared_proteins_2grs <- top_proteins[top_proteins$protein %in% common_2grs, ]
  
  ## reshape data to wide format
  heatmap_data_common_2grs <- shared_proteins_2grs %>%
    select(protein, group, coef) %>%
    pivot_wider(names_from = group, values_from = coef, values_fill = 0) %>%
    column_to_rownames("protein") %>%
    as.matrix()
  
  
  ## get names & replace proteins' IDs by their full names
  matched_names_common_2grs <- proteins_info_file$EntrezGeneSymbol[match(rownames(heatmap_data_common_2grs), proteins_info_file$AptName)]
  rownames(heatmap_data_common_2grs) <- ifelse(is.na(matched_names_common_2grs), rownames(heatmap_data_common_2grs), matched_names_common_2grs)
  
  ## reverse to get proteins' names on x-axis instead of y-axis
  heatmap_data_common_scaled_2grs <- t(heatmap_data_common_2grs)
  
  ## plot heatmap
  heatmap_plot_shared_2grs <- pheatmap::pheatmap(heatmap_data_common_scaled_2grs,
                                                 cluster_rows = FALSE,
                                                 cluster_cols = FALSE,
                                                 cellheight = 40,
                                                 angle_col = "90",
                                                 color = palette_col,
                                                 breaks = palette_lim,
                                                 border_color = "lightgrey",
                                                 main = "Common protein(s) across 2 groups")
  ggsave("heatmap_shared_2grs_prots.tiff", plot = heatmap_plot_shared_2grs, width = 16, height = 10, dpi = 300)
}



