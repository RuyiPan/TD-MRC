

temp <- readRDS("job_name=Gaussian_Predictionjob_num=1LPS.rds")
Gaussian<- c(temp$LPS,  temp$LPML, temp$WAIC$WAIC,temp$DIC)

temp <- readRDS("job_name=MixClayton_Predictionjob_num=1LPS_mix_new.rds")
MixC<- c(temp$LPS,  temp$LPML, temp$WAIC$WAIC,temp$DIC)
  

temp <- readRDS("job_name=MixClayton_Prediction_ct20_q7job_num=1LPS_mix_new.rds")
MixC_at20_q7<- c(temp$LPS,  temp$LPML, temp$WAIC$WAIC,temp$DIC)

temp <- readRDS("job_name=MixClayton_Prediction_ct30_q7job_num=1LPS_mix_new.rds")
MixC_at30_q7<- c(temp$LPS,  temp$LPML, temp$WAIC$WAIC,temp$DIC)

temp <- readRDS("job_name=SingleClayton_Predictionjob_num=1LPS_single.rds")
SingleC <- c(temp$LPS, temp$LPML, temp$WAIC$WAIC, temp$DIC)

round(rbind(MixC, MixC_at20_q7, MixC_at30_q7,Gaussian, SingleC),2)


