Parameters <- rbind(expand.grid(c(1),c(0,1,3,5,10,20,30,40), c(0,1,2,3,4,5,6,7)),
                    expand.grid(c(1000), c(30),c(0,1,2,3,4,5,6,7)))


gofs <- NULL
for (job_num in 1:nrow(Parameters)) {
  tmp <- readRDS(paste0("MixC_diff_parameters_new/job_name=MixClayton_simulation_diff_parametersjob_num=",job_num,"LPS_mix_new.rds"))
  gofs <- rbind(gofs, 
                c(tmp$LPML,tmp$WAIC$WAIC,tmp$DIC))
}

res <- as.data.frame(cbind(Parameters,round(gofs,0)))
colnames(res) <- c("bk", "ct", "q", "LPML", "WAIC", "DIC")

res %>% select(bk, ct, q, LPML) %>%
  pivot_wider(names_from = c("q"), values_from = "LPML")


res %>% select(bk, ct, q, WAIC) %>%
  pivot_wider(names_from = c("q"), values_from = "WAIC")


res %>% select(bk, ct, q, DIC) %>%
  pivot_wider(names_from = c("q"), values_from = "DIC")
