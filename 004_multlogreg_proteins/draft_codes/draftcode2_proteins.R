###################################################################
###################################################################
###### INTERNSHIP PROJECT - draft 2 proteins (not necessary) ######
###################################################################
###################################################################

## Ciloë Sans (M1 internship project, june-august 2025)

## 1. start:
setwd("/myfilepath/") ## set working directory
rm(list = ls()) ## clear working env
cat("\014") ## clear console


## import necessary libraries
library(dplyr)
library(ggplot2)
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


## load data
PRS_df <- data.table::fread("/myfilepath/mydata.txt")
proteins_info_file <- data.table::fread("/myfilepath/mydata2.txt")


## load risk groups (gen info, PRS002013 ONLY)
PGS002013_df <- data.table::fread("/myfilepath/mydata3.txt")


## fuse both databases together
df_ALL <- left_join(PGS002013_df, PRS_df, by = "Idepic")



## 2. protein analysis:
proteins_varnames <- grep("seq", names(df_ALL), value = TRUE)

## check if any NAs
for (i in 1:length(proteins_varnames)) {
  any(which(is.na(df_ALL[, proteins_varnames[i]]) != 0))
}

## convert risks col to factor var:
df_ALL$risk_gr_clrt_cncr_PGS_relax <- as.factor(df_ALL$risk_gr_clrt_cncr_PGS_relax)
levels(df_ALL$risk_gr_clrt_cncr_PGS_relax)

df_ALL$risk_gr_clrt_cncr_PGS_relax <- factor(df_ALL$risk_gr_clrt_cncr_PGS_relax, levels = c("Low risk//No cancer", "Low risk//Cancer", "High risk//No cancer", "High risk//Cancer"))
levels(df_ALL$risk_gr_clrt_cncr_PGS_relax)



## 3. mult log reg model:
## -> run univariate ANOVA on proteins to select those with the lowest p-values before running glmnet()
pvals_proteins <- sapply(proteins_varnames, function(protein_name){
  model_anova <- aov(df_ALL[[protein_name]] ~ df_ALL$risk_gr_clrt_cncr_PGS_relax)
  summary(model_anova)[[1]][["Pr(>F)"]][1]
})

## FDR to perform on the p-vals of ANOVA to adj for multiple comparisons
pvals_proteins_adj <- p.adjust(pvals_proteins, method = "fdr")

## select only proteins with p-val < 0.5
proteins_sel <- names(pvals_proteins_adj)[pvals_proteins_adj < 0.65]


## -> using package glmnet
vars_sel <- c(proteins_sel, "age", "sex", "pa_index", "bmi_c")
model_glm <- cv.glmnet(as.matrix(df_ALL[, ..vars_sel]), df_ALL$risk_gr_clrt_cncr_PGS_relax, family = "multinomial", alpha = 0.5, type.multinomial = "ungrouped", nfolds = 10)
save(model_glm, file = "mult_log_model1.rda")

plot(model_glm) ## check if model did a strong shrinkage or not => red curve should have a V/U-shape
ggsave("model_glm_fitplot.tiff", height = 16, width = 10, dpi = 300)

## save coefs GLM model
coefs_glm <- coef(model_glm, s = "lambda.1se")

## create class names vector
crc_classes_names <- c("Low_risk//No_cancer", "Low_risk//Cancer", "High_risk//No_cancer", "High_risk//Cancer")  ## replace with lvl names

## select non-0 coefs from GLM model
proteins_nonzero <- lapply(crc_classes_names, function(cls) {
  M <- coefs_glm[[cls]]
  idx <- which(M != 0 & rownames(M) != "(Intercept)")
  
  if (length(idx) == 0) return(NULL)
  
  data.frame(
    protein = rownames(M)[idx],
    coef = as.numeric(M[idx]),
    group = cls,
    stringsAsFactors = FALSE
  )
})

crc_df <- do.call(rbind, proteins_nonzero)


## -> top-20 proteins per cncr group (incl common ones across groups)
top_proteins <- crc_df %>%
  group_by(group) %>%
  arrange(desc(abs(coef))) %>%
  slice_head(n = 20)

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
alone_proteins <- top_proteins[top_proteins$protein %in% unique_proteins, ]


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

## plot heatmap
pheatmap::pheatmap(heatmap_data_scaled,
                   cluster_rows = FALSE,
                   cluster_cols = FALSE,
                   cellheight = 40,
                   #cellwidth = 16,
                   angle_col = "90",
                   color = colorRampPalette(c("darkgreen", "white", "darkred"))(5),
                   border_color = "lightgrey",#"transparent",
                   main = "Non-shared proteins (among top proteins) across risk groups")
ggsave("heatmap_notshared_prots.tiff", width = 16, height = 10, dpi = 300)



