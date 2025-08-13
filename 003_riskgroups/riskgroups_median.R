#############################################################
#############################################################
#### INTERNSHIP PROJECT - risk groups (median threshold) ####
#############################################################
#############################################################

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
cox_model_data_nogeninfo <- PRS_data_u65[complete.cases(PRS_data_u65[, c("start_age", "stop_age", "censored", "cncr_mal_clrt", "sex", "country", "alc_re", "smoke_stat", "bmi_c", "hli_dietscore")]), ]
cox_model_data_geninfo <- PRS_data_u65[complete.cases(PRS_data_u65[, c("start_age", "stop_age", "censored", "cncr_mal_clrt", "sex", "country", "alc_re", "smoke_stat", "bmi_c", "hli_dietscore", "weights_somalogic_genetic", "PGS002013_hmPOS_GRCh37", "PGS003850_hmPOS_GRCh37")]), ]

## remove "Missing" cat from genetic info df (leads to issues and outliers in the residuals)
cox_model_data_geninfo <- cox_model_data_geninfo %>% filter(pa_index != "Missing")

## e. normalise PGS variables (leads to an extremely high HR & CI so not logical => rescaling needed)
cox_model_data_geninfo$PGS002013_hmPOS_GRCh37_norm <- scale(cox_model_data_geninfo$PGS002013_hmPOS_GRCh37)
cox_model_data_geninfo$PGS003850_hmPOS_GRCh37_norm <- scale(cox_model_data_geninfo$PGS003850_hmPOS_GRCh37)



## 3. implement the Cox PH models:
## -> all individuals (no gen info, no weights)
cox_model_ALL <- coxph(Surv(start_age, stop_age, censored) ~ strata(sex, country) + pa_index + alc_re + smoke_stat + bmi_c, data = cox_model_data_nogeninfo)

## -> all inds w/ weights (but no PGS)
cox_model_ALL_weighted <- coxph(Surv(start_age, stop_age, censored) ~ strata(sex, country) + pa_index + alc_re + smoke_stat + bmi_c, weights = weights_somalogic_genetic, data = cox_model_data_geninfo)

## -> only PGS + weights (2 dfs as 2 PGS)
cox_model_PGS_relax <- coxph(Surv(start_age, stop_age, censored) ~ strata(sex, country) + pa_index + alc_re + smoke_stat + bmi_c + PGS002013_hmPOS_GRCh37_norm, weights = weights_somalogic_genetic, data = cox_model_data_geninfo)
cox_model_PGS_restr <- coxph(Surv(start_age, stop_age, censored) ~ strata(sex, country) + pa_index + alc_re + smoke_stat + bmi_c + PGS003850_hmPOS_GRCh37_norm, weights = weights_somalogic_genetic, data = cox_model_data_geninfo)



## 4. group inds (per stratum combination) based on their risk score
## compute predicted survivals according to each stratum
cox_model_data_nogeninfo$pred_risk_ALL <- predict(cox_model_ALL, type = "lp")
cox_model_data_geninfo$pred_risk_ALL_weighted <- predict(cox_model_ALL_weighted, type = "lp")
cox_model_data_geninfo$pred_risk_PGS_relax <- predict(cox_model_PGS_relax, type = "lp")
cox_model_data_geninfo$pred_risk_PGS_restr <- predict(cox_model_PGS_restr, type = "lp")

## group inds based on their stratum category
cox_model_data_nogeninfo$strat_gr <- interaction(cox_model_data_nogeninfo$sex, cox_model_data_nogeninfo$country, drop = TRUE)
cox_model_data_geninfo$strat_gr <- interaction(cox_model_data_geninfo$sex, cox_model_data_geninfo$country, drop = TRUE)


## make groupings
## -> ALL (i.e., no gen. info.):
cox_model_data_nogeninfo <- cox_model_data_nogeninfo %>% group_by(strat_gr) %>% mutate(risk_gr = ifelse(pred_risk_ALL >= median(pred_risk_ALL), "High_risk", "Low_risk"))
cox_model_data_nogeninfo <- cox_model_data_nogeninfo %>% mutate(risk_gr_clrt_cncr = paste(risk_gr, ifelse(cncr_mal_clrt == 1, "Cancer", "No_cancer"), sep = "//"))

risk_gr_l_noc <- cox_model_data_nogeninfo %>% filter(risk_gr_clrt_cncr == "Low risk//No cancer")
risk_gr_l_c <- cox_model_data_nogeninfo %>% filter(risk_gr_clrt_cncr == "Low risk//Cancer")
risk_gr_h_noc <- cox_model_data_nogeninfo %>% filter(risk_gr_clrt_cncr == "High risk//No cancer")
risk_gr_h_c <- cox_model_data_nogeninfo %>% filter(risk_gr_clrt_cncr == "High risk//Cancer")


## -> ALL but weighted (i.e., no PGS):
cox_model_data_geninfo <- cox_model_data_geninfo %>% group_by(strat_gr) %>% mutate(risk_gr_ALL_weighted = ifelse(pred_risk_ALL_weighted >= median(pred_risk_ALL_weighted), "High_risk", "Low_risk"))
cox_model_data_geninfo <- cox_model_data_geninfo %>% mutate(risk_gr_clrt_cncr_ALL_weighted = paste(risk_gr_ALL_weighted, ifelse(cncr_mal_clrt == 1, "Cancer", "No_cancer"), sep = "//"))


## -> PGS relaxed version
cox_model_data_geninfo <- cox_model_data_geninfo %>% group_by(strat_gr) %>% mutate(risk_gr_PGS_relax = ifelse(pred_risk_PGS_relax >= median(pred_risk_PGS_relax), "High_risk", "Low_risk"))
cox_model_data_geninfo <- cox_model_data_geninfo %>% mutate(risk_gr_clrt_cncr_PGS_relax = paste(risk_gr_PGS_relax, ifelse(cncr_mal_clrt == 1, "Cancer", "No_cancer"), sep = "//"))

risk_gr_l_noc <- cox_model_data_geninfo %>% filter(risk_gr_clrt_cncr_PGS_relax == "Low risk//No cancer")
risk_gr_l_c <- cox_model_data_geninfo %>% filter(risk_gr_clrt_cncr_PGS_relax == "Low risk//Cancer")
risk_gr_h_noc <- cox_model_data_geninfo %>% filter(risk_gr_clrt_cncr_PGS_relax == "High risk//No cancer")
risk_gr_h_c <- cox_model_data_geninfo %>% filter(risk_gr_clrt_cncr_PGS_relax == "High risk//Cancer")


## -> PGS restricted version
cox_model_data_geninfo <- cox_model_data_geninfo %>% group_by(strat_gr) %>% mutate(risk_gr_PGS_restr = ifelse(pred_risk_PGS_restr >= median(pred_risk_PGS_restr), "High_risk", "Low_risk"))
cox_model_data_geninfo <- cox_model_data_geninfo %>% mutate(risk_gr_clrt_cncr_PGS_restr = paste(risk_gr_PGS_restr, ifelse(cncr_mal_clrt == 1, "Cancer", "No_cancer"), sep = "//"))



## make dfs only with id, risk score & risk gr it corresponds to
df_risk_ALL <- cox_model_data_nogeninfo[, c("Idepic", "pred_risk_ALL", "risk_gr_clrt_cncr")]
df_risk_ALL_weighted <- cox_model_data_geninfo[, c("Idepic", "pred_risk_ALL_weighted", "risk_gr_clrt_cncr_ALL_weighted")]
df_risk_PGS_relax <- cox_model_data_geninfo[, c("Idepic", "pred_risk_PGS_relax", "risk_gr_clrt_cncr_PGS_relax")]
df_risk_PGS_restr <- cox_model_data_geninfo[, c("Idepic", "pred_risk_PGS_restr", "risk_gr_clrt_cncr_PGS_restr")]

## save dfs
write.table(df_risk_ALL, "risk_score_ALL.txt", sep = "\t", row.names = FALSE)
write.table(df_risk_ALL_weighted, "risk_score_ALL_weighted.txt", sep = "\t", row.names = FALSE)
write.table(df_risk_PGS_relax, "risk_score_PGS_relax.txt", sep = "\t", row.names = FALSE)
write.table(df_risk_PGS_restr, "risk_score_PGS_restr.txt", sep = "\t", row.names = FALSE)



## 5. analysis of the =/= risks
plot(df_risk_ALL$pred_risk_ALL, type = "l")
plot(df_risk_ALL_weighted$pred_risk_ALL_weighted, type = "l")
plot(df_risk_PGS_relax$pred_risk_PGS_relax, type = "l")
plot(df_risk_PGS_restr$pred_risk_PGS_restr, type = "l")

table(cox_model_data_geninfo$hli_smoke[cox_model_data_geninfo$cncr_mal_clrt == 1])

hist(cox_model_data_geninfo$PGS002013_hmPOS_GRCh37_norm, breaks = 50, probability = TRUE)
lines(density(cox_model_data_geninfo$PGS002013_hmPOS_GRCh37_norm), col = "red")

hist(cox_model_data_geninfo$PGS003850_hmPOS_GRCh37_norm, breaks = 50, probability = TRUE)
lines(density(cox_model_data_geninfo$PGS003850_hmPOS_GRCh37_norm), col = "red")

hist(df_risk_ALL$pred_risk_ALL, breaks = 50, probability = TRUE)
lines(density(df_risk_ALL$pred_risk_ALL), col = "red")

hist(df_risk_ALL_weighted$pred_risk_ALL_weighted, breaks = 50, probability = TRUE)
lines(density(df_risk_ALL_weighted$pred_risk_ALL_weighted), col = "red")

hist(df_risk_PGS_relax$pred_risk_PGS_relax, breaks = 50, probability = TRUE)
lines(density(df_risk_PGS_relax$pred_risk_PGS_relax), col = "red")

hist(df_risk_PGS_restr$pred_risk_PGS_restr, breaks = 50, probability = TRUE)
lines(density(df_risk_PGS_restr$pred_risk_PGS_restr), col = "red")



