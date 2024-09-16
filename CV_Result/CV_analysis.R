Gaussian <- NULL
for (i in 0:9) {
  temp <- readRDS(paste0("Gaussian/job_name=Gaussian_CVjob_num=",i,"CV.rds"))
  Gaussian <- rbind(Gaussian,
                    c(temp$MSE,temp$LPML,temp$WAIC$WAIC,temp$DIC))
}
apply(Gaussian, 2, mean)
 

SingleC <-NULL
for (i in 0:9) {
  temp <- readRDS(paste0("SingleC/job_name=SingleClayton_CVjob_num=",i,"LPS_single.rds"))
  SingleC  <- rbind(SingleC ,
                    c(temp$MSE,temp$LPML,temp$WAIC$WAIC,temp$DIC))
}
apply(SingleC , 2, mean)


MixC <-NULL
for (i in 0:9) {
  temp <- readRDS(paste0("MixC/job_name=MixClayton_CVjob_num=",i,"LPS_mix_new.rds"))
  MixC  <- rbind(MixC ,
                    c(temp$MSE,temp$LPML,temp$WAIC$WAIC,temp$DIC))
}
apply(MixC , 2, mean)


round(cbind(apply(MixC , 2, mean),
            apply(Gaussian, 2, mean),
            apply(SingleC , 2, mean)),
      2)
