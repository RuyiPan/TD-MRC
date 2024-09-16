
Gaussian_LPS <- vector()
for (i in 1:6) {
  temp<- readRDS(paste0("Gaussian_Prediction/job_name=Gaussian_Predictionjob_num=",i, "LPS.rds"))
  Gaussian_LPS[i] <- temp$LPS
}
Gaussian_LPS

Mix_LPS <- vector()
for (i in 1:6){
  temp<- readRDS(paste0("MixC_ct=10_Prediction/job_name=MixClayton_Prediction_ak=10job_num=",i,"LPS_mix.rds"))
  Mix_LPS[i] <- temp$LPS
}
Mix_LPS

Single_LPS <- vector()
for (i in 1:6){
  temp<- readRDS(paste0("SingleClayton_Prediction/job_name=SingleClayton_Predictionjob_num=",i,"LPS_single.rds"))
  Single_LPS[i] <- temp$LPS
}
Single_LPS

rbind(Gaussian_LPS, Mix_LPS,Single_LPS )
round(Single_LPS , 2)


Guassian_CV_LPML <- vector()
Guassian_CV_DIC <- vector()
for (i in 0:9){
  temp <- readRDS(paste0("Gaussian_CV/job_name=Gaussian_CVjob_num=",i,"CV.rds"))
  Guassian_CV_LPML[i+1] <- temp$LPML
  Guassian_CV_DIC[i+1] <- temp$DIC
}
mean(Guassian_CV_LPML)
mean(Guassian_CV_DIC)
