##############################################
##############################################
###### INTERNSHIP PROJECT - cox model 1 ######
##############################################
##############################################

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

summary(PRS_data_1$sex) ## just to check

## create start & end dates
PRS_data <- PRS_data_1 %>% group_by(Idepic) %>% mutate(start_age = age_recr, stop_age = age_exit_cancer_1st)

summary(PRS_data[, c("start_age", "stop_age")]) ## just to check


## handle missing values (remove/replace NAs)
#sex + country + pa_index + hli_alcohol + hli_smoke + hli_bmic + hli_dietscore
## a. check if any NAs:
which(is.na(PRS_data$sex)) ## none
which(is.na(PRS_data$country)) ## none
which(is.na(PRS_data$pa_index)) ## none

which(is.na(PRS_data$hli_alcohol))
which(is.na(PRS_data$hli_smoke))
which(is.na(PRS_data$hli_bmic))
which(is.na(PRS_data$hli_dietscore))
## all the 4 above have approx 204 NAs (& the same individuals with missing vals on their HLI)

## b. replace NAs by + rep category
unique(PRS_data$hli_alcohol)


## make df without NAs for the specific variables and covariates that are selected for the cox model
cox_model_data_clean <- PRS_data[complete.cases(PRS_data[, c("start_age", "stop_age", "censored", "cncr_mal_clrt", "sex", "country", "pa_index", "hli_alcohol", "hli_smoke", "hli_bmic", "hli_dietscore")]), ]


## 3. implement the Cox PH model:
## -> unstratified model
## alc_pattern + smoke_stat
cox_model_unstrat <- coxph(Surv(start_age, stop_age, censored) ~ sex + country + pa_index + hli_alcohol + hli_smoke + hli_bmic + hli_dietscore, data = cox_model_data_clean)
summary(cox_model_unstrat)

## check PH assumption (all p-values > 0.05? => satisfied if so)
check_ph <- cox.zph(cox_model_unstrat)
print(check_ph)


## -> stratified model
cox_model_strat <- coxph(Surv(start_age, stop_age, censored) ~ strata(sex, country) + pa_index + hli_alcohol + hli_smoke + hli_bmic + hli_dietscore, data = cox_model_data_clean)
summary(cox_model_strat)

## check PH assumption (all p-values > 0.05? => satisfied if so)
check_ph <- cox.zph(cox_model_strat)
print(check_ph)



####################################################################################################################################################################################
## METHOD 1 - predict()
## predict absolute risk of dvping clrt cancer for each gender
#pred_men_strat <- predict(cox_model_strat, type = "survival", newdata = df_men, times = seq(min(df_men$start_age), max(df_men$stop_age), by = 1))
#pred_women_strat <- predict(cox_model_strat, type = "survival", newdata = df_women, times = seq(min(df_women$start_age), max(df_women$stop_age), by = 1))
#
## compute abs risk from predicted values
#abs_risk_men <- 1 - pred_men_strat
#abs_risk_women <- 1 - pred_women_strat
#
## plot absolute risk
#plot(abs_risk_men, type = "l", col = "cyan")
#plot(abs_risk_women, type = "l", col = "pink")
####################################################################################################################################################################################



## METHOD 2 - survfit()
surv_cncr_risk_unstrat <- survfit(cox_model_unstrat)
surv_cncr_risk <- survfit(cox_model_strat)


plot(surv_cncr_risk, fun = "event", col = 1:8, lty = 1, xlim = c(min(PRS_data$stop_age), max(PRS_data$stop_age)), ylim = c(0, 1), xlab = "Age", ylab = "Absolute Risk", main = "Absolute risk for colon-rectum cancer (by sex and country)")
grid()
legend("topleft", legend = names(surv_cncr_risk$strata), col = 1:8, cex = 0.6, lty = 1)



