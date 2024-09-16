
library(tidyverse)
#
Parameters <- expand.grid(c(0,1,2,3,4,5),c(0,1,2,3,4), c(0,1,2))
gofs <- NULL
for (job_num in 1:nrow(Parameters)) {
  tmp <- readRDS(paste0("Result_new/job_name=MixClayton_fitting_realjob_num=",job_num,"LPS_mix_new.rds"))
  gofs <- rbind(gofs, 
                c(tmp$LPML,tmp$WAIC$WAIC,tmp$DIC))
}


res <- as.data.frame(cbind(Parameters,round(gofs,0)))
colnames(res) <- c("ct", "p", "q", "LPML", "WAIC", "DIC")
LPML_1 <- res %>% select(ct,p,q,LPML) %>% 
  filter(q==0) %>%
  select(ct,p,LPML) %>%
  pivot_wider(names_from =c("p"), values_from = "LPML")
LPML_1

LPML_2 <- res %>% select(ct,p,q,LPML) %>% 
  filter(q==1) %>%
  select(ct,p,LPML) %>%
  pivot_wider(names_from =c("p"), values_from = "LPML")
LPML_2

LPML_3 <-res %>% select(ct,p,q,LPML) %>% 
  filter(q==2) %>%
  select(ct,p,LPML) %>%
  pivot_wider(names_from =c("p"), values_from = "LPML")
LPML3
######################WAIC
res %>% select(ct,p,q,WAIC) %>% 
  filter(q==0) %>%
  select(ct,p,WAIC) %>%
  pivot_wider(names_from =c("p"), values_from = "WAIC")

res %>% select(ct,p,q,WAIC) %>% 
  filter(q==1) %>%
  select(ct,p,WAIC) %>%
  pivot_wider(names_from =c("p"), values_from = "WAIC")

res %>% select(ct,p,q,WAIC) %>% 
  filter(q==2) %>%
  select(ct,p,WAIC) %>%
  pivot_wider(names_from =c("p"), values_from = "WAIC")




######################DIC
res %>% select(ct,p,q,DIC) %>% 
  filter(q==0) %>%
  select(ct,p,DIC) %>%
  pivot_wider(names_from =c("p"), values_from = "DIC")

res %>% select(ct,p,q,DIC) %>% 
  filter(q==1) %>%
  select(ct,p,DIC) %>%
  pivot_wider(names_from =c("p"), values_from = "DIC")

res %>% select(ct,p,q,DIC) %>% 
  filter(q==2) %>%
  select(ct,p,DIC) %>%
  pivot_wider(names_from =c("p"), values_from = "DIC")


