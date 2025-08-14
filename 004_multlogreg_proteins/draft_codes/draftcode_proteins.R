#####################################################################
#####################################################################
######## INTERNSHIP PROJECT - draft proteins (not necessary) ########
#####################################################################
#####################################################################

## Ciloë Sans (M1 internship project, june-august 2025)

## 1. start:
rm(list = ls()) ## clear working env
cat("\014") ## clear console

## import necessary libraries
library(dplyr)
library(ggplot2)
library(glmnet)
library(haven)
library(mgcv)
library(nnet)
library(reshape2)
library(survival)
library(tidyr)
library(xtable)

## load data
PRS_df <- data.table::fread("/myfilepath/mycode.txt")
proteins_info_file <- data.table::fread("/myfilepath/mycode2.txt")

## load risk groups (gen info, PRS002013 ONLY)
PGS002013_df <- data.table::fread("/myfilepath/risk_score.txt")

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


## 3. univariate models (for 5 =/= proteins, one at a time)
## -> protein name seq.10014.31:
#df_ALL$risk_gr_clrt_cncr_PGS_relax2 <- relevel(df_ALL$risk_gr_clrt_cncr_PGS_relax, ref = "Low_risk//No_cancer")
formula_1p_1 <- as.formula(paste("risk_gr_clrt_cncr_PGS_relax ~ ", paste(proteins_varnames[10], collapse = " + ")))
model_1p_1 <- multinom(formula_1p_1, data = df_ALL)
summary(model_1p_1) ## AIC \approx 2625.133

## compute p-vals:
z1 <- summary(model_1p_1)$coefficients/summary(model_1p_1)$standard.errors
p_vals1 <- 2*(1 - pnorm(abs(z1)))
p_vals1 ## only l. risk // no cncr is stat signif


## -> protein name seq.21887.2:
formula_1p_2 <- as.formula(paste("risk_gr_clrt_cncr_PGS_relax ~ ", paste(proteins_varnames[3452], collapse = " + ")))
model_1p_2 <- multinom(formula_1p_2, data = df_ALL)
summary(model_1p_2) ## AIC \approx 2633.301

## compute p-vals:
z2 <- summary(model_1p_2)$coefficients/summary(model_1p_2)$standard.errors
p_vals2 <- 2*(1 - pnorm(abs(z2)))
p_vals2 ## not stat signif for any group


## -> protein name seq.11993.227:
formula_1p_3 <- as.formula(paste("risk_gr_clrt_cncr_PGS_relax ~ ", paste(proteins_varnames[702], collapse = " + ")))
model_1p_3 <- multinom(formula_1p_3, data = df_ALL)
summary(model_1p_3) ## AIC \approx 2635.083

## compute p-vals:
z3 <- summary(model_1p_3)$coefficients/summary(model_1p_3)$standard.errors
p_vals3 <- 2*(1 - pnorm(abs(z3)))
p_vals3 ## not stat signif for any group


## -> protein name seq.8775.61:
formula_1p_4 <- as.formula(paste("risk_gr_clrt_cncr_PGS_relax ~ ", paste(proteins_varnames[6937], collapse = " + ")))
model_1p_4 <- multinom(formula_1p_4, data = df_ALL)
summary(model_1p_4) ## AIC \approx 2632.831

## compute p-vals:
z4 <- summary(model_1p_4)$coefficients/summary(model_1p_4)$standard.errors
p_vals4 <- 2*(1 - pnorm(abs(z4)))
p_vals4 ## not stat signif for any group


## -> protein name seq.24423.9:
formula_1p_5 <- as.formula(paste("risk_gr_clrt_cncr_PGS_relax ~ ", paste(proteins_varnames[4223], collapse = " + ")))
model_1p_5 <- multinom(formula_1p_5, data = df_ALL)
summary(model_1p_5) ## AIC \approx 2629.665

## compute p-vals:
z5 <- summary(model_1p_5)$coefficients/summary(model_1p_5)$standard.errors
p_vals5 <- 2*(1 - pnorm(abs(z5)))
p_vals5 ## stat signif at 10% for high risk//no cncr & low risk//no cncr groups



## 4. multivariate models
## -> with 100 proteins:
proteins_sel_100 <- proteins_varnames[2378:2477]

formula_100p <- as.formula(paste("risk_gr_clrt_cncr_PGS_relax ~ ", paste(proteins_sel_100, collapse = " + ")))
model_100p <- multinom(formula_100p, data = df_ALL)


## select proteins with p-values under 0K.1 (stat signif at 10%):
## compute p-values
z100 <- summary(model_100p)$coefficients/summary(model_100p)$standard.errors
p_vals100 <- 2*(1 - pnorm(abs(z100)))
p_vals100

## get coefs from model
coefs_multinom <- summary(model_100p)$coefficients

## convert coefs matrix to long format
coefs_df_long <- melt(coefs_multinom, varnames = c("risk_group", "protein_name"), value.name = "coef")
coefs_df_long <- coefs_df_long[coefs_df_long$protein_name != "(Intercept)", ] ## need to be removed as not a protein

## select top 10 proteins (have only 100 proteins here)
top_10_p <- coefs_df_long %>% group_by(risk_group) %>% arrange(desc(abs(coef))) %>% slice_head(n = 10)

## select shared proteins
shared_p <- top_10_p %>% group_by(protein_name) %>% summarise(n_grs = n_distinct(risk_group)) %>% filter(n_grs > 1)

top_10_p <- top_10_p %>% mutate(group_in_common = protein_name %in% shared_p$protein_name)

## plot top 10 proteins for each group
ggplot(top_10_p, aes(x = reorder(protein_name, abs(coef)), y = coef, fill = group_in_common)) + 
  geom_col(show.legend = TRUE) + 
  coord_flip() + 
  facet_wrap(~ risk_group) + 
  labs(title = "Top proteins for each risk group", x = "Protein name", y = "Coef") + 
  scale_fill_manual(values = c("TRUE" = "blue", "FALSE"= "darkgrey")) + 
  theme_minimal()


## run univariate ANOVA on proteins to select those with the lowest p-values before running glmnet()
pvals_proteins <- sapply(proteins_varnames, function(protein_name){
  model_anova <- aov(df_ALL[[protein_name]] ~ df_ALL$risk_gr_clrt_cncr_PGS_relax)
  summary(model_anova)[[1]][["Pr(>F)"]][1]
})

## FDR to perform on the p-vals of ANOVA to adj for multiple comparisons
pvals_proteins_adj <- p.adjust(pvals_proteins, method = "fdr")

## select only proteins with p-val < 0.2
proteins_sel <- names(pvals_proteins_adj)[pvals_proteins_adj <= 0.23]

## using package glmnet
vars_sel <- c(proteins_sel, "sex", "country")
model_glm <- cv.glmnet(as.matrix(df_ALL[, ..vars_sel]), df_ALL$risk_gr_clrt_cncr_PGS_relax, family = "multinomial", alpha = 0.5, type.multinomial = "ungrouped", nfolds = 10)
plot(model_glm) ## check if model did a strong shrinkage or not => red curve should have a V/U-shape

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


# top-20 proteins per cncr group (incl common ones across groups)
top_proteins <- crc_df %>%
  group_by(group) %>%
  arrange(desc(abs(coef))) %>%
  slice_head(n = 20)

# common proteins
common <- Reduce(intersect, lapply(proteins_nonzero, function(df){df$protein}))

print(common)

top_proteins$common_proteins <- top_proteins$protein %in% common

top_prot_plot <- ggplot(top_proteins, aes(x = reorder(protein, abs(coef)), y = coef, fill = common_proteins)) +
  geom_col(show.legend = TRUE) +
  coord_flip() +
  facet_wrap(~ group, scales = "free_y") +
  labs(
    title = "Top 20 proteins cancer development",
    x = "Protein name",
    y = "Coefficient"
  ) +
  scale_fill_manual(values = c("TRUE" = "blue", "FALSE"= "darkgrey")) + 
  theme_minimal()
ggsave("plotname.tiff", plot = top_prot_plot, width = 20, height = 10, dpi = 1000)



