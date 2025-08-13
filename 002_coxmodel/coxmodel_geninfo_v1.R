##############################################
##############################################
###### INTERNSHIP PROJECT - cox model 2 ######
##############################################
##############################################

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
PRS_data_1$hli_alcohol <- as.factor(PRS_data_1$hli_alcohol)
PRS_data_1$hli_smoke <- as.factor(PRS_data_1$hli_smoke)
PRS_data_1$hli_bmic <- as.factor(PRS_data_1$hli_bmic)
PRS_data_1$hli_dietscore <- as.factor(PRS_data_1$hli_dietscore)
PRS_data_1$hli_dietscore_c <- as.factor(PRS_data_1$hli_dietscore_c)

## create start & end dates
PRS_data <- PRS_data_1 %>% group_by(Idepic) %>% mutate(start_age = age_recr, stop_age = age_exit_cancer_1st)

summary(PRS_data[, c("start_age", "stop_age")]) ## just to check



## check nb inds from each country that have genetic info present
PRS_data$gen_info <- ifelse(!is.na(PRS_data$PGS002013_hmPOS_GRCh37), 1, 0) ## has gen info = 1; has no gen info = 0
table(PRS_data$gen_info, PRS_data$country) ## only individuals from italy and the uk have inds w/ gen info



## handle missing values (remove/replace NAs)
## a. check if any NAs:
which(is.na(PRS_data$sex)) ## none
which(is.na(PRS_data$country)) ## none

which(is.na(PRS_data$hli_alcohol))
which(is.na(PRS_data$hli_smoke))
which(is.na(PRS_data$hli_bmic))
which(is.na(PRS_data$hli_dietscore))
## all the 4 above have approx 204 NAs (& the same individuals with missing vals on their HLI)


## b. replace NAs by + rep category
table(PRS_data$hli_alcohol)
table(PRS_data$hli_smoke)
table(PRS_data$hli_bmic)
table(PRS_data$hli_dietscore)
table(PRS_data$hli_dietscore_c)


## c. replace NAs by + rep cat
## first check + rep cat by looking at the proportions in each category:
table(PRS_data$hli_alcohol) ## change column with appropriate variable to check

PRS_data$hli_alcohol[is.na(PRS_data$hli_alcohol)] <- "0-6 (g/d)"
PRS_data$hli_smoke[is.na(PRS_data$hli_smoke)] <- "Never"
PRS_data$hli_bmic[is.na(PRS_data$hli_bmic)] <- "26-<30"
PRS_data$hli_dietscore[is.na(PRS_data$hli_dietscore)] <- 27
PRS_data$hli_dietscore_c[is.na(PRS_data$hli_dietscore_c)] <- 2


## d. get only df with no NAs for the variables used in the cox model
cox_model_data <- PRS_data[complete.cases(PRS_data[, c("start_age", "stop_age", "censored", "cncr_mal_clrt", "sex", "country", "hli_alcohol", "hli_smoke", "hli_bmic", "hli_dietscore_c", "weights_somalogic_genetic", "PGS002013_hmPOS_GRCh37")]), ]


## e. normalise PGS variables (leads to an extremely high HR & CI so not logical => rescaling needed)
cox_model_data$PGS002013_hmPOS_GRCh37_norm <- scale(cox_model_data$PGS002013_hmPOS_GRCh37)
cox_model_data$PGS003850_hmPOS_GRCh37_norm <- scale(cox_model_data$PGS003850_hmPOS_GRCh37)


## save df as it contains genetic info & can be useful to already have it for future work (has weights & PGS so a reduced df with only these present)
#write.table(cox_model_data, "df_only_genetic_info.txt", sep = "\t", row.names = FALSE)

cox_model_data_clean <- cox_model_data %>% filter(pa_index != "Missing")

cox_model_data_clean$pa_index <- factor(cox_model_data_clean$pa_index, levels = c("Inactive", "Moderately inactive", "Moderately active", "Active"))
cox_model_data_clean$hli_alcohol <- factor(cox_model_data_clean$hli_alcohol, levels = c("0-6 (g/d)", ">6-12 (g/d)", ">12-24 (g/d)", ">24-60 (g/d)", ">60 (g/d)"))
cox_model_data_clean$hli_smoke <- factor(cox_model_data_clean$hli_smoke, levels = c("Never", "Former, quit 10+ years", "Former, quit <= 10 years", "Current, 1-15 cig/day", "Current, 16+ cig/day"))
cox_model_data_clean$hli_bmic <- factor(cox_model_data_clean$hli_bmic, levels = c("<22", "22-<24", "24-<26", "26-<30", ">=30"))


## 3. implement the Cox PH model:
## -> unstratified model
## alc_pattern + smoke_stat
cox_model_unstrat <- coxph(Surv(start_age, stop_age, censored) ~ sex + country + pa_index + hli_alcohol + hli_smoke + hli_bmic + hli_dietscore_c + PGS002013_hmPOS_GRCh37_norm, weights = weights_somalogic_genetic, data = cox_model_data_clean)
summary(cox_model_unstrat)

## check PH assumption (all p-values > 0.05? => satisfied if so)
check_ph_unstrat <- cox.zph(cox_model_unstrat)
print(check_ph_unstrat)


## -> stratified model
cox_model_strat <- coxph(Surv(start_age, stop_age, censored) ~ strata(sex, country) + pa_index + hli_alcohol + hli_smoke + hli_bmic + hli_dietscore_c + PGS002013_hmPOS_GRCh37_norm, weights = weights_somalogic_genetic, data = cox_model_data_clean)
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


## make new df argument => assign each covariate the + rep category
range_PRS <- quantile(cox_model_data_clean$PGS002013_hmPOS_GRCh37_norm, 0.5)

new_df_m_it <- data.frame(start_age = 40,
                          stop_age = 80,
                          censored = 0, 
                          sex = factor("Men", levels = levels(cox_model_data_clean$sex)), 
                          country = factor("Italy", levels = levels(cox_model_data_clean$country)), 
                          pa_index = factor("Moderately inactive", levels = levels(cox_model_data_clean$pa_index)), 
                          hli_alcohol = factor("0-6 (g/d)", levels = levels(cox_model_data_clean$hli_alcohol)), 
                          hli_smoke = factor("Never", levels = levels(cox_model_data_clean$hli_smoke)), 
                          hli_bmic = factor("26-<30", levels = levels(cox_model_data_clean$hli_bmic)), 
                          hli_dietscore_c = factor("0", levels = levels(cox_model_data_clean$hli_dietscore_c)), 
                          PGS002013_hmPOS_GRCh37_norm = range_PRS)

new_df_w_it <- data.frame(start_age = 40,
                          stop_age = 80,
                          censored = 0, 
                          sex = factor("Women", levels = levels(cox_model_data_clean$sex)), 
                          country = factor("Italy", levels = levels(cox_model_data_clean$country)), 
                          pa_index = factor("Moderately inactive", levels = levels(cox_model_data_clean$pa_index)), 
                          hli_alcohol = factor("0-6 (g/d)", levels = levels(cox_model_data_clean$hli_alcohol)), 
                          hli_smoke = factor("Never", levels = levels(cox_model_data_clean$hli_smoke)), 
                          hli_bmic = factor("26-<30", levels = levels(cox_model_data_clean$hli_bmic)), 
                          hli_dietscore_c = factor("0", levels = levels(cox_model_data_clean$hli_dietscore_c)), 
                          PGS002013_hmPOS_GRCh37_norm = range_PRS)

new_df_m_uk <- data.frame(start_age = 40,
                          stop_age = 80,
                          censored = 0, 
                          sex = factor("Men", levels = levels(cox_model_data_clean$sex)), 
                          country = factor("United Kingdom", levels = levels(cox_model_data_clean$country)), 
                          pa_index = factor("Moderately inactive", levels = levels(cox_model_data_clean$pa_index)), 
                          hli_alcohol = factor("0-6 (g/d)", levels = levels(cox_model_data_clean$hli_alcohol)), 
                          hli_smoke = factor("Never", levels = levels(cox_model_data_clean$hli_smoke)), 
                          hli_bmic = factor("26-<30", levels = levels(cox_model_data_clean$hli_bmic)), 
                          hli_dietscore_c = factor("0", levels = levels(cox_model_data_clean$hli_dietscore_c)), 
                          PGS002013_hmPOS_GRCh37_norm = range_PRS)

new_df_w_uk <- data.frame(start_age = 40,
                          stop_age = 80,
                          censored = 0, 
                          sex = factor("Women", levels = levels(cox_model_data_clean$sex)), 
                          country = factor("United Kingdom", levels = levels(cox_model_data_clean$country)), 
                          pa_index = factor("Moderately inactive", levels = levels(cox_model_data_clean$pa_index)), 
                          hli_alcohol = factor("0-6 (g/d)", levels = levels(cox_model_data_clean$hli_alcohol)), 
                          hli_smoke = factor("Never", levels = levels(cox_model_data_clean$hli_smoke)), 
                          hli_bmic = factor("26-<30", levels = levels(cox_model_data_clean$hli_bmic)), 
                          hli_dietscore_c = factor("0", levels = levels(cox_model_data_clean$hli_dietscore_c)), 
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


## check no NAs/extreme values/outliers
pred_surv <- predict(cox_model_strat, type = "lp")
pred_risk <- 1 - pred_surv

plot(pred_risk, type = "l")

## check norm of the resids w QQ-plot
qqnorm(pred_risk)
qqline(pred_risk, col = "red")

hist(pred_risk, breaks = 50, prob = TRUE, ylim = c(0, 1))
lines(density(pred_risk), col ="red")



## METHOD 2 - survfit()
surv_cncr_risk_unstrat <- survfit(cox_model_unstrat, newdata = new_df_m_it)

surv_cncr_risk_m_it <- survfit(cox_model_strat, newdata = new_df_m_it)
surv_cncr_risk_w_it <- survfit(cox_model_strat, newdata = new_df_w_it)
surv_cncr_risk_m_uk <- survfit(cox_model_strat, newdata = new_df_m_uk)
surv_cncr_risk_w_uk <- survfit(cox_model_strat, newdata = new_df_w_uk)


## get abs risk
abs_cncr_risk_m_it <- 1 - surv_cncr_risk_m_it$surv
abs_cncr_risk_w_it <- 1 - surv_cncr_risk_w_it$surv
abs_cncr_risk_m_uk <- 1 - surv_cncr_risk_m_uk$surv
abs_cncr_risk_w_uk <- 1 - surv_cncr_risk_w_uk$surv

min_age <- min(cox_model_data_clean$start_age)
max_age <- max(cox_model_data_clean$stop_age)
max_pred <- max(abs_cncr_risk_m_it, abs_cncr_risk_w_it, abs_cncr_risk_m_uk, abs_cncr_risk_w_uk)

plot(surv_cncr_risk_m_it$time, abs_cncr_risk_m_it, type = "l", xlab = "Age", ylab = "Absolute Risk", main = "Absolute risk for colon-rectum cancer (by sex and country)", xlim = c(min_age, max_age), ylim = c(0, max_pred))
lines(surv_cncr_risk_w_it$time, abs_cncr_risk_w_it, col = "red")
lines(surv_cncr_risk_m_uk$time, abs_cncr_risk_m_uk, col = "cyan")
lines(surv_cncr_risk_w_uk$time, abs_cncr_risk_w_uk, col = "purple")
grid()
legend("topleft", legend = c("Men, IT", "Women, IT", "Men, UK", "Women, UK"), col = c("black", "red", "cyan", "purple"), cex = 0.6, lty = 1)



## 5. group inds (per stratum combination) based on their risk score
## compute predicted survivals according to each stratum
cox_model_data_clean$pred_risk <- predict(cox_model_strat, type = "lp")
cox_model_data_clean$strat_gr <- interaction(cox_model_data_clean$sex, cox_model_data_clean$country, drop = TRUE)

cox_model_data_clean <- cox_model_data_clean %>% group_by(strat_gr) %>% mutate(risk_gr = ifelse(pred_risk >= median(pred_risk), "High risk", "Low risk"))
cox_model_data_clean <- cox_model_data_clean %>% mutate(risk_gr_clrt_cncr = paste(risk_gr, ifelse(cncr_mal_clrt == 1, "Cancer", "No cancer"), sep = "//"))


## plot barplots to see how each covariate affects each risk gr (i.e., for each stratum)
risk_gr_l_noc <- cox_model_data_clean %>% filter(risk_gr_clrt_cncr == "Low risk//No cancer")
risk_gr_l_c <- cox_model_data_clean %>% filter(risk_gr_clrt_cncr == "Low risk//Cancer")
risk_gr_h_noc <- cox_model_data_clean %>% filter(risk_gr_clrt_cncr == "High risk//No cancer")
risk_gr_h_c <- cox_model_data_clean %>% filter(risk_gr_clrt_cncr == "High risk//Cancer")


## plot density proportion of age for each risk gr (+ ref density of overall df to compare the proportions)
plot(density(cox_model_data_clean$age_blood), xlab = "Age", main = "Proportion of age categories for each risk group", ylim = c(0, 0.06))
density(risk_gr_l_noc$age_blood) %>% lines(col = "red")
density(risk_gr_l_c$age_blood) %>% lines(col = "blue")
density(risk_gr_h_noc$age_blood) %>% lines(col = "orange")
density(risk_gr_h_c$age_blood) %>% lines(col = "purple")
grid()
legend("topleft", col = c("black", "red", "blue", "orange", "purple"), c("Reference", "Low risk, no cncr", "Low risk, cncr", "High risk, no cncr", "High risk, cncr"), lty = 1, cex = 0.6)
## conclusion: + ppl w clrt cncr fall in the 55-60 yrs old age range & those w/out clrt cncr follow a similar distr than the ref density (i.e., peak age around 57)



