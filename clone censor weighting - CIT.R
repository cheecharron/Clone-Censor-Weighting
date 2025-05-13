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

#dataframe for storing integrated Brier skill scores
IBS = data.frame(matrix(ncol=5,nrow=0))
colnames(IBS) = c('group','denom_train','denom_test','num_train','num_test')
row=1

for (dataset in c('early','late','never')) {
#  dataset="early"
  IBS[row,1] = dataset
  
  df = get(dataset)
  
  ###############################################
  
  ### CREATE TRAINING AND VALIDATION SETS ###
  first = sqldf("select distinct subject from df")
  
  #randomly sample from list
  set.seed(1235)
  
  index = sample(1:nrow(first), .7*nrow(first))
  
  train_ids = first[index,]
  test_ids = first[-index,]
  
  train = df[df$subject %in% train_ids,]
  test = df[df$subject %in% test_ids,]
  
  #####################
  # Denominator Model #
  #####################
  
  tic('create model') #started 10/2 at 7:18pm
  mod = LTRCIT(Surv(month_lag, month, censor)~first_weight+BMI+diabetes, data = train)
  toc() #took 2.8 seconds 

#  df_nomiss = df[complete.cases(df),]
  
  tic('generate predictions')
  predobj = predict(mod, newdata=df, type='prob')
  toc() #took 0.17 seconds
  
  predobj = do.call(rbind, predobj) %>% data.frame
  
  #merging time predictions with observed time
  time = cbind(month=df$month, predobj) 
  
  #generating probability based on time index
  df$value = apply(time, 1, function(x) x$surv[which.min(abs(x$month-x$time))])
  df$IPW_denom = with(df, ifelse(censor==1, 1-value, 
                                 ifelse(censor==0,value,NA)))
  
  #######################################################################
  # need to generate Brier skill score for training and validation sets #
  #######################################################################
  
  ## FIRST TRAINING ##
  tic('Brier - training')
  
  trainb = sqldf("select * 
  			  from df
  			  where subject in 
  			  	(select subject 
  				from train)
                   order by subject, month")
  
  expprop = mean(trainb$censor)
  
  #calculate Brier scores
  trainb = trainb %>%
    mutate(brier = (month-month_lag)*(value-(1-censor))**2)
  
  #calculate Brier score of reference
  trainb$brier_ref = (trainb$month-trainb$month_lag)*(expprop-(1-trainb$censor))**2
  
  #tally Brier skill score
  IBS[row,2] = as.numeric(sqldf('select 1-avg(brier)/avg(brier_ref)
                      from
                        (select subject, sum(brier)/max(month) as brier, sum(brier_ref)/max(month) as brier_ref
                        from trainb
                        group by subject)'))
  
  toc()
  
  ## NOW THE REMAINDER OF THE DATA ##

  tic('Brier - validation')
  
  validb = sqldf("select * 
  			  from df
  			  where subject not in 
  			  	(select subject 
  				from train)")
  
  expprop = mean(validb$censor)
  
  #calculate Brier scores
  validb = validb %>%
    mutate(brier = (month-month_lag)*(value-(1-censor))**2)
  
  #calculate Brier score of reference
  validb$brier_ref = (validb$month-validb$month_lag)*(expprop-(1-validb$censor))**2
  
  #tally Brier skill score
  IBS[row,3] = as.numeric(sqldf('select 1-avg(brier)/avg(brier_ref)
                      from
                        (select subject, sum(brier)/max(month) as brier, sum(brier_ref)/max(month) as brier_ref
                        from validb
                        group by subject)'))
  
  toc()
  
  ###################
  # Numerator Model #
  ###################
  
  tic('create model') #started 10/2 at 7:18pm
  mod = LTRCIT(Surv(month_lag, month, censor)~first_weight, data = train)
  toc() #took 2.8 seconds 
  
  #  df_nomiss = df[complete.cases(df),]
  
  tic('generate predictions')
  predobj = predict(mod, newdata=df, type='prob')
  toc() #took 0.17 seconds
  
  predobj = do.call(rbind, predobj) %>% data.frame
  
  #merging time predictions with observed time
  time = cbind(month=df$month, predobj) 
  
  #generating probability based on time index
  df$value = apply(time, 1, function(x) x$surv[which.min(abs(x$month-x$time))])
  df$IPW_num = with(df, ifelse(censor==1, 1-value, 
                                 ifelse(censor==0,value,NA)))
  
  ############################################################################
  # need to generate integrated Brier score for training and validation sets #
  ############################################################################
  
  ## FIRST TRAINING ##
  tic('Brier - training')
  
  trainb = sqldf("select * 
  			  from df
  			  where subject in 
  			  	(select subject 
  				from train)
                   order by subject, month")
  
  expprop = mean(trainb$censor)
  
  #calculate Brier scores
  trainb = trainb %>%
    mutate(brier = (month-month_lag)*(value-(1-censor))**2)
  
  #calculate Brier score of reference
  trainb$brier_ref = (trainb$month-trainb$month_lag)*(expprop-(1-trainb$censor))**2
  
  #tally Brier skill score
  IBS[row,4] = as.numeric(sqldf('select 1-avg(brier)/avg(brier_ref)
                      from
                        (select subject, sum(brier)/max(month) as brier, sum(brier_ref)/max(month) as brier_ref
                        from trainb
                        group by subject)'))
  
  toc()
  
  ## NOW THE REMAINDER OF THE DATA ##
  
  tic('Brier - validation')
  
  validb = sqldf("select * 
  			  from df
  			  where subject not in 
  			  	(select subject 
  				from train)")
  
  expprop = mean(validb$censor)
  
  #calculate Brier scores
  validb = validb %>%
    mutate(brier = (month-month_lag)*(value-(1-censor))**2)
  
  #calculate Brier score of reference
  validb$brier_ref = (validb$month-validb$month_lag)*(expprop-(1-validb$censor))**2
  
  #tally Brier skill score
  IBS[row,5] = as.numeric(sqldf('select 1-avg(brier)/avg(brier_ref)
                      from
                        (select subject, sum(brier)/max(month) as brier, sum(brier_ref)/max(month) as brier_ref
                        from validb
                        group by subject)'))
  
  toc()
  
  #########################################################################
  
  df = df %>%
    mutate(IPCW = ifelse(IPW_num==IPW_denom,1,IPW_num/IPW_denom)) %>%
    filter(censor!=1)
  
  assign(dataset,df)
  
  row=row+1
}
  
write.csv(IBS, "IBS.csv")

merged = rbind(early %>% dplyr::select(subject,month_lag,month,first_weight,wt_chg,weight_lead,IPCW) %>% mutate(exp="early"),
               late %>% dplyr::select(subject,month_lag,month,first_weight,wt_chg,weight_lead,IPCW) %>% mutate(exp="late"), 
               never %>% dplyr::select(subject,month_lag,month,first_weight,wt_chg,weight_lead,IPCW) %>% mutate(exp="never")) %>%
  mutate(exp = factor(exp))
