# The goal of this file is first to simulate longitudinal exposure data,
# whereby some exposure (x) can be predicted by covariates but in turn
# influences some outcome (y), which is also influenced by the 
# covariates. The second goal is to simulate two causal approaches to
# analyze the data: the weighted approach I came up with and the clone-
# censor weighting approach.

library(dplyr)
library(lme4)
library(sqldf)
library(emmeans)
library(lmerTest)
library(simr)
library(tidyr)

# The simulated dataset will have 10,000 patients who will each
# have 24 observations (corresponding to 1/month for 2 years). 
# Half of these patients will become exposed to treatment (x).
# Exposure to treatment x will be predicted by time-varying 
# covariates (a) and (b) via a monotonic function. 

set.seed(4321)

# Parameters
n <- 10000  # Number of participants
n_time <- 24  # Number of time points
BMI = c(n, mean=30, sd=2) #setting parameters for initial BMI

# Create a data frame
df <- data.frame(subject = rep(1:n))
df$BMI1 = rnorm(n,mean=28,sd=4)
df$height_in = rnorm(n,mean=71,sd=1.5)

## Diabetes Risk ##
# setting diabetes risk as very low (0.05%) for the lowest BMI (i.e., <=25).
# After that, the risk fo T2DM will increase by 130% for each BMI point.
# That is, BMI of 26 would have risk of (0.5%), and BMI of 45 would have risk 
# of 95.0%.

for (row in 1:n){
  print(row)
  df[row,'diabetes1'] = ifelse(df[row,'BMI1']<=25, rbinom(1, 1, .005), #setting risk of T2DM to .5% if BMI<=25
                      ifelse(df[row,'BMI1']>25, rbinom(1,1,.005*1.3^(df[row,'BMI1']-25)),NA))
  
  df[row,'WLM1'] = ifelse(df[row,'BMI1']>30, rbinom(1,1,.1),
                          ifelse(df[row,'BMI1']>30 & df[row,'diabetes1']==1, rbinom(1,1,.5), 0))
}

## WLM Exposure ##
# if BMI>30, likelihood is 10%.
# If diabetes=1 & BMI>30, likelihood is 50%.


# now need to loop over times to populate BMI, diabetes, WLM, and weight
# BMI will increase by an average of 1-point over the 2-year period (1/24 points each month)

for (col in c(2:24)){
  print(col)
  
  for (row in 1:n){
    #if previous WLM, BMI will drop by 1/12-point until it reaches approximately 24
    df[row,paste0('BMI',col)] = ifelse(df[row,paste0('WLM',col-1)]==1 & df[row,paste0('BMI',col-1)]>24,
                                       df[row,paste0('BMI',col-1)] - rnorm(1,mean=1/5,sd=1/10), 
                                       df[row,paste0('BMI',col-1)] + rnorm(1,mean=1/24,sd=1/100))

    df[row,paste0('diabetes',col)] = ifelse(df[row,paste0('diabetes',col-1)]==1,1,
                                            ifelse(df[row,paste0('BMI',col)]<=25, rbinom(1, 1, .005), #setting risk of T2DM to .5% if BMI<=25
                                                   ifelse(df[row,paste0('BMI',col)]>25, rbinom(1,1,.005*1.3^(df[row,paste0('BMI',col)]-25)),NA)))
    df[row,paste0('WLM',col)] = ifelse(df[row,paste0('WLM',col-1)]==1,1,
                                       ifelse(df[row,paste0('BMI',col)]>30, rbinom(1,1,.1),
                                              ifelse(df[row,paste0('BMI',col)]>30 & df[row,paste0('diabetes',col)]==1, rbinom(1,1,.5), 0)))

  }
}

## now need to transpose from wide to long ##
long = df %>% reshape(direction="long",
                      varying=colnames(df %>% dplyr::select(-c(subject,height_in))),
                      timevar="month",
                      times=c(1:24),
                      v.names=(c('BMI','diabetes','WLM')),
                      idvar=c('subject','height_in')) %>%
  mutate(weight = BMI*height_in^2/703) %>%
  arrange(subject,month) %>%
  group_by(subject) %>%
  mutate(weight_lead = lead(weight,1))

first = long %>% filter(month==1) %>% 
  dplyr::select(subject,weight) %>% 
  dplyr::rename(first_weight=weight)
long = merge(long,first,all.x=T) %>%
  mutate(wt_chg = weight_lead-first_weight)

#plot Kaplan-Meier curve for WLM exposure

temp1 = sqldf('select subject, min(month) as month, WLM
              from long
              where WLM==1
              group by subject')

temp2 = sqldf('select subject, 24 as month, WLM
              from long
              where subject not in
                (select subject
                from temp1)
              group by subject')

temp = rbind(temp1,temp2); rm(temp1); rm(temp2)

survfit2(Surv(month,WLM) ~ 1, data=temp) %>%
  ggsurvfit()

#############################################################
# establishing group membership based on timing of exposure #
#############################################################

temp$WLM_time = with(temp, ifelse(WLM==1,month,24))
temp$exp_group = with(temp, ifelse(WLM==1 & month<=6,'early',
                                   ifelse(WLM==1 & month<=12,'late','never')))

long = merge(long, temp %>% dplyr::select(c(subject,WLM_time,exp_group)), all.x=T) %>%
  group_by(subject) %>%
  mutate(month_lag = ifelse(month==1,0,lag(month,1)))
                           
rm(first); rm(temp); gc()
  