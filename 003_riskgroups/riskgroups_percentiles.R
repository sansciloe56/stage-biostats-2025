######################################################################
######################################################################
###### INTERNSHIP PROJECT - risk groups (percentiles threshold) ######
######################################################################
######################################################################

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
library(survival)
library(tidyr)
library(xtable)

## load data
PRS_data_1 <- data.table::fread("/myfilepath/mydata.txt")


## 2. some data cleaning, updating, ...
## create censored variable (= 1 if individual has colon-rectum cancer; = 0 otherwise)
PRS_data_1$censored <- ifelse(PRS_data_1$cncr_mal_clrt == 1, 1, 0)


## convert sex to a factor variable & rename its categories (men = 1; women = 2)
PRS_data_1$sex <- factor(PRS_data_1$sex)
levels(PRS_data_1$sex) <- c("Men", "Women")


## make age groups (useful for when looking at risk of dvping cancer in terms of risk groups at the end of the code)
PRS_data_1$age_rounded <- round(PRS_data_1$age_blood)

PRS_data_1$age_cat <- case_when(PRS_data_1$age_rounded %in% 30:39 ~ "30-39",
                                PRS_data_1$age_rounded %in% 40:49 ~ "40-49",
                                PRS_data_1$age_rounded %in% 50:59 ~ "50-59",
                                PRS_data_1$age_rounded %in% 60:80 ~ "60+")


## convert covariates of interest as factors
PRS_data_1$country <- as.factor(PRS_data_1$country)
PRS_data_1$smoke_stat <- as.factor(PRS_data_1$smoke_stat)
PRS_data_1$hli_dietscore <- as.factor(PRS_data_1$hli_dietscore)
PRS_data_1$hli_dietscore_c <- as.factor(PRS_data_1$hli_dietscore_c)
PRS_data_1$pa_index <- as.factor(PRS_data_1$pa_index)


## create start & end dates
PRS_data <- PRS_data_1 %>% group_by(Idepic) %>% mutate(start_age = age_recr, stop_age = age_exit_cancer_1st)


## check nb inds from each country that have genetic info present
PRS_data$gen_info <- ifelse(!is.na(PRS_data$PGS002013_hmPOS_GRCh37), 1, 0) ## has gen info = 1; has no gen info = 0


## filter w participants lost to fup bef 65 to remove
PRS_data_u65 <- PRS_data %>% filter(age_exit_t2d >= 65)


## check nb inds from each country that have genetic info present
PRS_data_u65$gen_info <- ifelse(!is.na(PRS_data_u65$PGS002013_hmPOS_GRCh37), 1, 0) ## has gen info = 1; has no gen info = 0
table(PRS_data_u65$gen_info, PRS_data_u65$country) ## only individuals from italy and the uk have inds w/ gen info



## handle missing values (remove/replace NAs)
## a. check if any NAs:
which(is.na(PRS_data_u65$sex)) ## none
which(is.na(PRS_data_u65$country)) ## none

which(is.na(PRS_data_u65$alc_re))
which(is.na(PRS_data_u65$smoke_stat))
which(is.na(PRS_data_u65$bmi_c))
which(is.na(PRS_data_u65$hli_dietscore))
## all the 4 above have no NAs except hli_dietscore so nothing to do for these except the last var mentioned

PRS_data_u65$hli_dietscore <- as.numeric(PRS_data_u65$hli_dietscore)
PRS_data_u65$hli_dietscore[is.na(PRS_data_u65$hli_dietscore)] <- mean(PRS_data_u65$hli_dietscore, na.rm = TRUE)

PRS_data_u65$pa_index <- factor(PRS_data_u65$pa_index, levels = c("Inactive", "Moderately inactive", "Moderately active", "Active", "Missing"))
PRS_data_u65$smoke_stat <- factor(PRS_data_u65$smoke_stat, levels = c("Never", "Former", "Smoker"))

## d. get only df with no NAs for the variables used in the cox model
cox_model_data_geninfo <- PRS_data_u65[complete.cases(PRS_data_u65[, c("start_age", "stop_age", "censored", "cncr_mal_clrt", "sex", "country", "alc_re", "smoke_stat", "bmi_c", "hli_dietscore", "weights_somalogic_genetic", "PGS002013_hmPOS_GRCh37", "PGS003850_hmPOS_GRCh37")]), ]

## remove "Missing" cat from genetic info df (leads to issues and outliers in the residuals)
cox_model_data_geninfo <- cox_model_data_geninfo %>% filter(pa_index != "Missing")

## e. normalise PGS variables (leads to an extremely high HR & CI so not logical => rescaling needed)
cox_model_data_geninfo$PGS002013_hmPOS_GRCh37_norm <- scale(cox_model_data_geninfo$PGS002013_hmPOS_GRCh37)
cox_model_data_geninfo$PGS003850_hmPOS_GRCh37_norm <- scale(cox_model_data_geninfo$PGS003850_hmPOS_GRCh37)



## 3. implement the Cox PH models:
## -> all inds w/ weights (but no PGS)
cox_model_ALL_weighted <- coxph(Surv(start_age, stop_age, censored) ~ strata(sex, country) + pa_index + alc_re + smoke_stat + bmi_c, weights = weights_somalogic_genetic, data = cox_model_data_geninfo)

## -> only PGS + weights (2 dfs as 2 PGS)
cox_model_PGS_relax <- coxph(Surv(start_age, stop_age, censored) ~ strata(sex, country) + pa_index + alc_re + smoke_stat + bmi_c + PGS002013_hmPOS_GRCh37_norm, weights = weights_somalogic_genetic, data = cox_model_data_geninfo)
cox_model_PGS_restr <- coxph(Surv(start_age, stop_age, censored) ~ strata(sex, country) + pa_index + alc_re + smoke_stat + bmi_c + PGS003850_hmPOS_GRCh37_norm, weights = weights_somalogic_genetic, data = cox_model_data_geninfo)



## 4. group inds (per stratum combination) based on their risk score
## compute predicted survivals according to each stratum
cox_model_data_geninfo$pred_risk_ALL_weighted <- predict(cox_model_ALL_weighted, type = "lp")
cox_model_data_geninfo$pred_risk_PGS_relax <- predict(cox_model_PGS_relax, type = "lp")
cox_model_data_geninfo$pred_risk_PGS_restr <- predict(cox_model_PGS_restr, type = "lp")

## group inds based on their stratum category
cox_model_data_geninfo$strat_gr <- interaction(cox_model_data_geninfo$sex, cox_model_data_geninfo$country, drop = TRUE)


## make groupings by looping with =/= thresholds for the =/= %iles
percentile <- seq(0.1, 0.9, 0.1)

## -> ALL but weighted (i.e., no PGS):
for (i in seq_along(percentile)) {
  p <- percentile[i]
  name_col1 <- paste0("pred_risk_ALL_weighted_", p * 100)
  name_col_risk_gr1 <- paste0("risk_gr_clrt_cncr_PGS_ALL_", p * 100)
  
  vals <- cox_model_data_geninfo$pred_risk_ALL_weighted
  lower_cut <- quantile(vals, probs = p, na.rm = TRUE)
  upper_cut <- quantile(vals, probs = 1 - p, na.rm = TRUE)
  
  cox_model_data_geninfo <- cox_model_data_geninfo %>%
    group_by(strat_gr) %>%
    mutate(
      !!name_col1 := case_when(
        pred_risk_ALL_weighted <= lower_cut ~ "Low_risk",
        pred_risk_ALL_weighted >= upper_cut ~ "High_risk",
        TRUE ~ NA_character_
      )
    ) %>%
    ungroup() %>%
    mutate(
      !!name_col_risk_gr1 := ifelse(
        !is.na(.data[[name_col1]]),
        paste(.data[[name_col1]], ifelse(cncr_mal_clrt == 1, "Cancer", "No_cancer"), sep = "//"),
        NA_character_
      )
    )
}


## -> PGS relaxed version
for (i in seq_along(percentile)) {
  p <- percentile[i]
  name_col2 <- paste0("pred_risk_PGS_relax_", p * 100)
  name_col_risk_gr2 <- paste0("risk_gr_clrt_cncr_PGS_relax_", p * 100)
  
  vals <- cox_model_data_geninfo$pred_risk_PGS_relax
  lower_cut <- quantile(vals, probs = p, na.rm = TRUE)
  upper_cut <- quantile(vals, probs = 1 - p, na.rm = TRUE)
  
  cox_model_data_geninfo <- cox_model_data_geninfo %>%
    group_by(strat_gr) %>%
    mutate(
      !!name_col2 := case_when(
        pred_risk_PGS_relax <= lower_cut ~ "Low_risk",
        pred_risk_PGS_relax >= upper_cut ~ "High_risk",
        TRUE ~ NA_character_
      )
    ) %>%
    ungroup() %>%
    mutate(
      !!name_col_risk_gr2 := ifelse(
        !is.na(.data[[name_col2]]),
        paste(.data[[name_col2]], ifelse(cncr_mal_clrt == 1, "Cancer", "No_cancer"), sep = "//"),
        NA_character_
      )
    )
}


## -> PGS restricted version
for (i in seq_along(percentile)) {
  p <- percentile[i]
  name_col3 <- paste0("pred_risk_PGS_restr_", p * 100)
  name_col_risk_gr3 <- paste0("risk_gr_clrt_cncr_PGS_restr_", p * 100)
  
  vals <- cox_model_data_geninfo$pred_risk_PGS_restr
  lower_cut <- quantile(vals, probs = p, na.rm = TRUE)
  upper_cut <- quantile(vals, probs = 1 - p, na.rm = TRUE)
  
  cox_model_data_geninfo <- cox_model_data_geninfo %>%
    group_by(strat_gr) %>%
    mutate(
      !!name_col3 := case_when(
        pred_risk_PGS_restr <= lower_cut ~ "Low_risk",
        pred_risk_PGS_restr >= upper_cut ~ "High_risk",
        TRUE ~ NA_character_
      )
    ) %>%
    ungroup() %>%
    mutate(
      !!name_col_risk_gr3 := ifelse(
        !is.na(.data[[name_col3]]),
        paste(.data[[name_col3]], ifelse(cncr_mal_clrt == 1, "Cancer", "No_cancer"), sep = "//"),
        NA_character_
      )
    )
}



## make df & save it:
risk_grs_names <- grep("risk_gr_clrt_cncr_", names(cox_model_data_geninfo), value = TRUE)
df_risk_grs <- cox_model_data_geninfo[, c("Idepic", "pred_risk_ALL_weighted", "pred_risk_PGS_relax", "pred_risk_PGS_restr", risk_grs_names)]
write.table(df_risk_grs, "risk_score_ALL.txt", sep = "\t", row.names = FALSE)

#table(df_risk_grs$risk_gr_clrt_cncr_PGS_relax_70)



