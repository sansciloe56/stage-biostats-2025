##########################################################
##########################################################
###### INTERNSHIP PROJECT - cox model FINAL version ######
##########################################################
##########################################################

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

summary(PRS_data_1$sex) ## just to check

## make age groups (useful for when looking at risk of dvping cancer in terms of risk groups at the end of the code)
PRS_data_1$age_rounded <- round(PRS_data_1$age_blood)

PRS_data_1$age_cat <- case_when(PRS_data_1$age_rounded %in% 30:39 ~ "30-39",
                                PRS_data_1$age_rounded %in% 40:49 ~ "40-49",
                                PRS_data_1$age_rounded %in% 50:59 ~ "50-59",
                                PRS_data_1$age_rounded %in% 60:80 ~ "60+")


## convert covariates of interest as factors
PRS_data_1$country <- as.factor(PRS_data_1$country)
PRS_data_1$smoke_stat <- as.factor(PRS_data_1$hli_smoke)
PRS_data_1$hli_dietscore <- as.factor(PRS_data_1$hli_dietscore)

## create start & end dates
PRS_data <- PRS_data_1 %>% group_by(Idepic) %>% mutate(start_age = age_recr, stop_age = age_exit_cancer_1st)

summary(PRS_data[, c("start_age", "stop_age")]) ## just to check

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


## b. get only df with no NAs for the variables used in the cox model
cox_model_data <- PRS_data_u65[complete.cases(PRS_data_u65[, c("start_age", "stop_age", "censored", "cncr_mal_clrt", "sex", "country", "alc_re", "smoke_stat", "bmi_c", "weights_somalogic_genetic", "PGS002013_hmPOS_GRCh37")]), ]


## c. normalise PGS variables (leads to an extremely high HR & CI so not logical => rescaling needed)
cox_model_data$PGS002013_hmPOS_GRCh37_norm <- scale(cox_model_data$PGS002013_hmPOS_GRCh37)
cox_model_data$PGS003850_hmPOS_GRCh37_norm <- scale(cox_model_data$PGS003850_hmPOS_GRCh37)

cox_model_data_clean <- cox_model_data %>% filter(pa_index != "Missing")

cox_model_data_clean$pa_index <- factor(cox_model_data_clean$pa_index, levels = c("Inactive", "Moderately inactive", "Moderately active", "Active"))
cox_model_data_clean$smoke_stat <- factor(cox_model_data_clean$smoke_stat, levels = c("Never", "Former", "Smoker"))



## 3. implement the Cox PH model:
## -> unstratified model
cox_model_unstrat <- coxph(Surv(start_age, stop_age, censored) ~ sex + country + pa_index + alc_re + smoke_stat + bmi_c + PGS002013_hmPOS_GRCh37_norm, weights = weights_somalogic_genetic, data = cox_model_data_clean)
summary(cox_model_unstrat)

## check PH assumption (all p-values > 0.05? => satisfied if so)
check_ph_unstrat <- cox.zph(cox_model_unstrat)
print(check_ph_unstrat)


## -> stratified model
cox_model_strat <- coxph(Surv(start_age, stop_age, censored) ~ strata(sex, country) + pa_index + alc_re + smoke_stat + bmi_c + PGS002013_hmPOS_GRCh37_norm, weights = weights_somalogic_genetic, data = cox_model_data_clean)
summary(cox_model_strat)

## check PH assumption (all p-values > 0.05? => satisfied if so)
check_ph <- cox.zph(cox_model_strat)
print(check_ph)



## 4. compute absolute risk
## PRS analysis
## a. check distribution:
mean_PRS <- mean(cox_model_data_clean$PGS002013_hmPOS_GRCh37_norm)

ggplot(data = cox_model_data_clean, aes(x = PGS002013_hmPOS_GRCh37_norm)) + 
  geom_histogram(bins = 50, fill = "cyan", color = "black") + 
  geom_vline(aes(xintercept = mean_PRS), color = "red", linetype = "dashed")


## b. make new df argument => assign each covariate the + rep category
range_PRS <- quantile(cox_model_data_clean$PGS002013_hmPOS_GRCh37_norm, 0.8)

new_df_m_it <- data.frame(start_age = 40,
                          stop_age = 70,
                          censored = 0, 
                          sex = factor("Men", levels = levels(cox_model_data_clean$sex)), 
                          country = factor("Italy", levels = levels(cox_model_data_clean$country)), 
                          pa_index = factor("Moderately inactive", levels = levels(cox_model_data_clean$pa_index)), 
                          alc_re = round(mean(cox_model_data_clean$alc_re), 4), 
                          smoke_stat = factor("Never", levels = levels(cox_model_data_clean$smoke_stat)), 
                          bmi_c = round(mean(cox_model_data_clean$bmi_c), 2), 
                          PGS002013_hmPOS_GRCh37_norm = range_PRS)

new_df_w_it <- data.frame(start_age = 40,
                          stop_age = 70,
                          censored = 0, 
                          sex = factor("Women", levels = levels(cox_model_data_clean$sex)), 
                          country = factor("Italy", levels = levels(cox_model_data_clean$country)), 
                          pa_index = factor("Moderately inactive", levels = levels(cox_model_data_clean$pa_index)), 
                          alc_re = round(mean(cox_model_data_clean$alc_re), 4), 
                          smoke_stat = factor("Never", levels = levels(cox_model_data_clean$smoke_stat)), 
                          bmi_c = round(mean(cox_model_data_clean$bmi_c), 2), 
                          PGS002013_hmPOS_GRCh37_norm = range_PRS)

new_df_m_uk <- data.frame(start_age = 40,
                          stop_age = 70,
                          censored = 0, 
                          sex = factor("Men", levels = levels(cox_model_data_clean$sex)), 
                          country = factor("United Kingdom", levels = levels(cox_model_data_clean$country)), 
                          pa_index = factor("Moderately inactive", levels = levels(cox_model_data_clean$pa_index)), 
                          alc_re = round(mean(cox_model_data_clean$alc_re), 4), 
                          smoke_stat = factor("Never", levels = levels(cox_model_data_clean$smoke_stat)), 
                          bmi_c = round(mean(cox_model_data_clean$bmi_c), 2), 
                          PGS002013_hmPOS_GRCh37_norm = range_PRS)

new_df_w_uk <- data.frame(start_age = 40,
                          stop_age = 70,
                          censored = 0, 
                          sex = factor("Women", levels = levels(cox_model_data_clean$sex)), 
                          country = factor("United Kingdom", levels = levels(cox_model_data_clean$country)), 
                          pa_index = factor("Moderately inactive", levels = levels(cox_model_data_clean$pa_index)), 
                          alc_re = round(mean(cox_model_data_clean$alc_re), 4), 
                          smoke_stat = factor("Never", levels = levels(cox_model_data_clean$smoke_stat)), 
                          bmi_c = round(mean(cox_model_data_clean$bmi_c), 2), 
                          PGS002013_hmPOS_GRCh37_norm = range_PRS)

####################################################################################################################################################################################
## TO IGNORE FOR NOW
## METHOD 1 - predict()
## predict absolute risk of dvping clrt cancer for each gender
#pred_strat <- predict(cox_model_strat, type = "survival", newdata = new_df)
#
## compute abs risk from predicted values
#abs_risk <- 1 - pred_strat
#
## plot absolute risk
#plot(pred_strat, type = "l", col = "cyan")
####################################################################################################################################################################################


## c. check no NAs/extreme values/outliers
pred_surv <- predict(cox_model_strat, type = "lp")
pred_risk <- 1 - pred_surv

plot(pred_risk, type = "l")

## d. check norm of the resids w QQ-plot
qqnorm(pred_risk)
qqline(pred_risk, col = "red")

hist(pred_risk, breaks = 50, prob = TRUE)
lines(density(pred_risk), col ="red")
dev.off()

## METHOD 2 - survfit():
surv_cncr_risk_unstrat <- survfit(cox_model_unstrat, newdata = new_df_m_it)

surv_cncr_risk_m_it <- survfit(cox_model_strat, newdata = new_df_m_it)
surv_cncr_risk_w_it <- survfit(cox_model_strat, newdata = new_df_w_it)
surv_cncr_risk_m_uk <- survfit(cox_model_strat, newdata = new_df_m_uk)
surv_cncr_risk_w_uk <- survfit(cox_model_strat, newdata = new_df_w_uk)


## a. get abs risk
abs_cncr_risk_m_it <- 1 - surv_cncr_risk_m_it$surv
abs_cncr_risk_w_it <- 1 - surv_cncr_risk_w_it$surv
abs_cncr_risk_m_uk <- 1 - surv_cncr_risk_m_uk$surv
abs_cncr_risk_w_uk <- 1 - surv_cncr_risk_w_uk$surv

min_age <- min(cox_model_data_clean$start_age)
max_age <- max(cox_model_data_clean$stop_age)
max_pred <- max(abs_cncr_risk_m_it, abs_cncr_risk_w_it, abs_cncr_risk_m_uk, abs_cncr_risk_w_uk)

## b. plot abs risk
tiff("abs_risks.tiff", width = 2500, height = 2000, res = 400)

plot(surv_cncr_risk_m_it$time, abs_cncr_risk_m_it, type = "l", xlab = "Age", ylab = "Absolute Risk", main = "Absolute risk for colon-rectum cancer (by sex and country)", xlim = c(min_age, max_age), ylim = c(0, max_pred))
lines(surv_cncr_risk_w_it$time, abs_cncr_risk_w_it, col = "red")
lines(surv_cncr_risk_m_uk$time, abs_cncr_risk_m_uk, col = "cyan")
lines(surv_cncr_risk_w_uk$time, abs_cncr_risk_w_uk, col = "purple")
grid()
legend("topleft", legend = c("Men, IT", "Women, IT", "Men, UK", "Women, UK"), col = c("black", "red", "cyan", "purple"), cex = 0.6, lty = 1)



## 5. group inds (per stratum combination) based on their risk score
## compute predicted survivals according to each stratum
cox_model_data_clean$pred_risk <- predict(cox_model_strat, type = "lp", newdata = cox_model_data_clean)
cox_model_data_clean$strat_gr <- interaction(cox_model_data_clean$sex, cox_model_data_clean$country, drop = TRUE)

cox_model_data_clean <- cox_model_data_clean %>% group_by(strat_gr) %>% mutate(risk_gr = ifelse(pred_risk >= median(pred_risk, na.rm = TRUE), "High_risk", "Low_risk"))
cox_model_data_clean <- cox_model_data_clean %>% mutate(risk_gr_clrt_cncr = paste(risk_gr, ifelse(cncr_mal_clrt == 1, "Cancer", "No_cancer"), sep = "//"))

## plot barplots to see how each covariate affects each risk gr (i.e., for each stratum)
risk_gr_l_noc <- cox_model_data_clean %>% filter(risk_gr_clrt_cncr == "Low_risk//No_cancer")
risk_gr_l_c <- cox_model_data_clean %>% filter(risk_gr_clrt_cncr == "Low_risk//Cancer")
risk_gr_h_noc <- cox_model_data_clean %>% filter(risk_gr_clrt_cncr == "High_risk//No_cancer")
risk_gr_h_c <- cox_model_data_clean %>% filter(risk_gr_clrt_cncr == "High_risk//Cancer")

## plot density proportion of age for each risk gr (+ ref density of overall df to compare the proportions)
tiff("agecat_prop.tiff", width = 2000, height = 1500, dpi = 400)

plot(density(cox_model_data_clean$age_blood), xlab = "Age", main = "Proportion of age categories for each risk group", ylim = c(0, 0.08))
density(risk_gr_l_noc$age_blood) %>% lines(col = "red")
density(risk_gr_l_c$age_blood) %>% lines(col = "blue")
density(risk_gr_h_noc$age_blood) %>% lines(col = "orange")
density(risk_gr_h_c$age_blood) %>% lines(col = "purple")
grid()
legend("topleft", col = c("black", "red", "blue", "orange", "purple"), c("Reference", "Low risk, no cncr", "Low risk, cncr", "High risk, no cncr", "High risk, cncr"), lty = 1, cex = 0.6)
dev.off()
## conclusion: + ppl w clrt cncr fall in the 55-60 yrs old age range & those w/out clrt cncr follow a similar distr than the ref density (i.e., peak age around 57)

## if want to see exact age proportion per group
hist(cox_model_data_clean$age_blood)
hist(risk_gr_l_noc$age_blood, add = TRUE, col = "red")
hist(risk_gr_h_noc$age_blood, add = TRUE, col = "orange")
hist(risk_gr_l_c$age_blood, add = TRUE, col = "blue")
hist(risk_gr_h_c$age_blood, add = TRUE, col = "green")



