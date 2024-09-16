
MixC_LPS <- vector()
MixC_LPML <- vector()
MixC_DIC <- vector()
MixC_WAIC <- vector()
for (i in 1:6) {
  temp <- readRDS(paste0("../Result_new/MixC_Prediction_MA0_S2/job_name=MixClayton_Prediction_MA0_S2job_num=",i,"LPS_mix_new.rds"))
  MixC_LPS[i] <- temp$LPS
  MixC_LPML[i] <- temp$LPML
  MixC_DIC[i] <- temp$DIC
  MixC_WAIC[i] <- temp$WAIC$WAIC
}
MixC_LPS
MixC_LPML
MixC_WAIC
MixC_DIC
#ct=1
#-3.388149  -4.522559 -12.325842  22.473366  25.834776  12.215031
#ct=2
#-4.217757  -4.214017 -13.425802  21.976766  27.099033  11.665025
#ct=5
#-9.033041  -4.557786 -11.197410  21.207068  26.406314  11.564243

Gaussian_DIC <- vector()
Gaussian_LPML <- vector()
Gaussian_WAIC <- vector()
Gaussian_LPS <- vector()
for (i in 1:6) {
  temp <- readRDS(paste0("../Result_new/Gaussian_Prediction/job_name=Gaussian_Predictionjob_num=",i,"LPS.rds"))
  Gaussian_LPS[i] <- temp$LPS
  Gaussian_DIC[i] <- temp$DIC
  Gaussian_LPML[i] <- temp$LPML
  Gaussian_WAIC[i] <- temp$WAIC$WAIC
}
Gaussian_LPS
Gaussian_LPML
Gaussian_WAIC
Gaussian_DIC

#5.193463 -1.645606 -1.278074 -7.122557 37.534893  6.122681
LPML_single <- vector()
WAIC_single <- vector()
DIC_single <- vector()
LPS_single <- vector()
for (i in 1:6) {
  temp <- readRDS(paste0("SingleC_Prediction/job_name=SingleClayton_Prediction_newjob_num=",i,"LPS_single.rds"))
  LPS_single[i] <- temp$LPS
  DIC_single[i] <- temp$DIC
  LPML_single[i] <- temp$LPML
  WAIC_single[i] <- temp$WAIC$WAIC
}
LPS_single
LPML_single
WAIC_single
DIC_single
#-79.53093 -170.49182 -153.62000 -258.18408 -246.04605 -199.73358
LPS <- rbind(MixC_LPS, Gaussian_LPS, LPS_single)
LPS <- round(LPS,2)
LPS
cbind(LPS,rowSums(LPS))

LPML <- rbind(MixC_LPML, Gaussian_LPML, LPML_single)
LPML <- round(LPML,2)
LPML
cbind(LPML,rowSums(LPML))

DIC <- rbind(MixC_DIC, Gaussian_DIC, DIC_single)
DIC <- round(DIC, 2)
DIC
cbind(DIC,rowSums(DIC))

WAIC <- rbind(MixC_WAIC, Gaussian_WAIC, WAIC_single)
WAIC <- round(WAIC, 2)
WAIC
cbind(WAIC,rowSums(WAIC))




MixC_MA0_S1 <- NULL
for (i in 1:6) {
  temp <- readRDS(paste0("../Result_new/MixC_Prediction_MA0_S1/job_name=MixC_Prediction_MA0_S1job_num=",i,"LPS_mix_new.rds"))
  MixC_MA0_S1 <- cbind(MixC_MA0_S1,
                       c(temp$LPS,temp$LPML, temp$WAIC$WAIC, temp$DIC))
}
round(MixC_MA0_S1,2)



MixC_MA0_S2 <- NULL
for (i in 1:6) {
  temp <- readRDS(paste0("../Result_new/MixC_Prediction_MA0_S2/job_name=MixC_Prediction_MA0_S2job_num=",i,"LPS_mix_new.rds"))
  MixC_MA0_S2 <- cbind(MixC_MA0_S2,
                       c(temp$LPS,temp$LPML, temp$WAIC$WAIC, temp$DIC))
}
round(MixC_MA0_S2,2)



MixC_MA1_S1 <- NULL
for (i in 1:6) {
  temp <- readRDS(paste0("../Result_new/MixC_Prediction_MA1_S1/job_name=MixC_Prediction_MA1_S1job_num=",i,"LPS_mix_new.rds"))
  MixC_MA1_S1 <- cbind(MixC_MA1_S1,
                       c(temp$LPS,temp$LPML, temp$WAIC$WAIC, temp$DIC))
}
round(MixC_MA1_S1,2)


MixC_MA1_S2 <- NULL
for (i in 1:6) {
  temp <- readRDS(paste0("../Result_new/MixC_Prediction_MA1_S2/job_name=MixC_Prediction_MA1_S2job_num=",i,"LPS_mix_new.rds"))
  MixC_MA1_S2 <- cbind(MixC_MA1_S2,
                       c(temp$LPS,temp$LPML, temp$WAIC$WAIC, temp$DIC))
}
round(MixC_MA1_S2,2)


MixC_MA2_S2 <- NULL
for (i in 1:6) {
  temp <- readRDS(paste0("../Result_new/MixC_Prediction_MA2_S2/job_name=MixC_Prediction_MA2_S2job_num=",i,"LPS_mix_new.rds"))
  MixC_MA2_S2 <- cbind(MixC_MA2_S2,
                       c(temp$LPS,temp$LPML, temp$WAIC$WAIC, temp$DIC))
}
round(MixC_MA2_S2,2)
LPS <- round(rbind(MixC_MA0_S1[1,],
      MixC_MA0_S2[1,],
      MixC_MA1_S1[1,],
      MixC_MA1_S2[1,],
      MixC_MA2_S2[1,]),2)
cbind(LPS,rowSums(LPS))

LPML <- round(rbind(MixC_MA0_S1[2,],
                   MixC_MA0_S2[2,],
                   MixC_MA1_S1[2,],
                   MixC_MA1_S2[2,],
                   MixC_MA2_S2[2,]),2)
LPML


WAIC <- round(rbind(MixC_MA0_S1[3,],
                    MixC_MA0_S2[3,],
                    MixC_MA1_S1[3,],
                    MixC_MA1_S2[3,],
                    MixC_MA2_S2[3,]),2)
WAIC




DIC <- round(rbind(MixC_MA0_S1[4,],
                    MixC_MA0_S2[4,],
                    MixC_MA1_S1[4,],
                    MixC_MA1_S2[4,],
                    MixC_MA2_S2[4,]),2)
DIC
