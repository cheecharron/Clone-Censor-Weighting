library(dplyr)
library(sqldf)
library(emmeans)
library(lmerTest)
library(survival)
library(ggsurvfit)
library(LTRCtrees)
library(tictoc)
library(mgcv)
library(pec)
library(tidyr)

######### Censor-Clone Weighting Approach ##########

#create three dataset with censor date
early = long %>% filter(exp_group=='early' |
                          (exp_group!='early' & month<7)) %>%
  mutate(censor = ifelse(exp_group!='early' & month==6,1,0))

late = long %>% filter(exp_group=='late' |
                         (exp_group=='early' & month<=WLM_time) |
                         (exp_group=='never' & month<13)) %>%
  mutate(censor = ifelse(exp_group=='early' & month==WLM_time,1,
                         ifelse(exp_group=='never' & month==12,1,0)))

never = long %>% filter(exp_group=='never' |
                          (exp_group!='never' & month<=WLM_time)) %>%
  mutate(censor = ifelse(exp_group!='never',WLM,0))

#####################################
# model censor date and create IPCW #
#####################################

for (dataset in c('early','late','never')) {
#  dataset="early"
  
  df = get(dataset)

  mod = coxph(Surv(month_lag, month, censor)~first_weight+BMI+diabetes, data = df)
  
  # Obtain the baseline cumulative hazard function
  base_haz <- basehaz(mod, centered = F)

  # Merge the baseline cumulative hazard with the data
  # We'll use the `base_haz` and interpolate to get the cumulative hazard for 'start' and 'end' times
  # Ensure data is sorted by 'end' time before merging with baseline hazard

  # Interpolate the baseline cumulative hazard for each 'end' time
  base_haz_interp <- approx(base_haz$time, base_haz$haz, xout = df$month, rule = 2)$y
  
  # Calculate linear predictor (risk score) for each individual using the fitted Cox model
  df$linear_predictor <- predict(mod, newdata=df, type = "lp")
  
  # Calculate cumulative hazard for each individual
  df$cum_haz_individual <- base_haz_interp * exp(df$linear_predictor)

  # Calculate survival probability at each 'end' point
  df$value <- exp(-df$cum_haz_individual)
  df$IPW_denom = with(df, ifelse(censor==1, 1-value, 
                                 ifelse(censor==0,value,NA)))

  
  mod = coxph(Surv(month_lag, month, censor)~first_weight, data = df)
  
  # Obtain the baseline cumulative hazard function
  base_haz <- basehaz(mod, centered = F)

  # Merge the baseline cumulative hazard with the data
  # We'll use the `base_haz` and interpolate to get the cumulative hazard for 'start' and 'end' times
  # Ensure data is sorted by 'end' time before merging with baseline hazard

  # Interpolate the baseline cumulative hazard for each 'end' time
  base_haz_interp <- approx(base_haz$time, base_haz$haz, xout = df$month, rule = 2)$y
  
  # Calculate linear predictor (risk score) for each individual using the fitted Cox model
  df$linear_predictor <- predict(mod, newdata=df, type = "lp")
  
  # Calculate cumulative hazard for each individual
  df$cum_haz_individual <- base_haz_interp * exp(df$linear_predictor)

  # Calculate survival probability at each 'end' point
  df$value <- exp(-df$cum_haz_individual)
  df$IPW_num = with(df, ifelse(censor==1, 1-value, 
                                 ifelse(censor==0,value,NA)))
  
  #########################################################################
  
  df = df %>%
    mutate(IPCW = ifelse(IPW_num==IPW_denom,1,IPW_num/IPW_denom)) %>%
    filter(censor!=1)
  
  assign(dataset,df)
}

merged = rbind(early %>% dplyr::select(subject,month_lag,month,first_weight,wt_chg,weight_lead,IPCW) %>% mutate(exp="early"),
               late %>% dplyr::select(subject,month_lag,month,first_weight,wt_chg,weight_lead,IPCW) %>% mutate(exp="late"), 
               never %>% dplyr::select(subject,month_lag,month,first_weight,wt_chg,weight_lead,IPCW) %>% mutate(exp="never")) %>%
  mutate(exp = factor(exp))
