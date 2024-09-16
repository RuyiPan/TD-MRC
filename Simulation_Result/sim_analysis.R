

temp <- readRDS("job_name=Gaussian_Predictionjob_num=1LPS.rds")
Gaussian<- c(temp$LPS,  temp$LPML, temp$WAIC$WAIC,temp$DIC)

temp <- readRDS("job_name=MixClayton_Predictionjob_num=1LPS_mix_new.rds")
MixC<- c(temp$LPS,  temp$LPML, temp$WAIC$WAIC,temp$DIC)
  
temp <- readRDS("job_name=SingleClayton_Predictionjob_num=1LPS_single.rds")
SingleC <- c(temp$LPS, temp$LPML, temp$WAIC$WAIC, temp$DIC)

round(rbind(MixC, Gaussian, SingleC),2)


tmp <- readRDS("job_name=MixClayton_simulation_diff_parametersjob_num=1LPS_mix_new.rds")
tmp$LPML
