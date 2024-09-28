burn.in <- 100
B <- 200
batch.size <- 50
range <- seq(burn.in*batch.size,B*batch.size,10)
est_tail_mean<- NULL
est_tail_lower <- NULL
est_tail_upper <- NULL
Parameters <- expand.grid(c(0,1,2,3,4,5),c(0,1,2,3,4), c(0,1,2))
job_num <- which(apply(Parameters, 1, function(row) all(row == c(1,1,2))))

temp_res <- readRDS(paste0("../Fitting_Result/Result_new/job_name=MixClayton_fitting_realjob_num=",job_num,"LPS_mix_new.rds"))
Theta <- temp_res$Theta
Pi <- temp_res$Pi
TT <- 36
for (t in 1:TT) {
  est_tail_mean <-  rbind(est_tail_mean,apply(do.call(rbind,lapply(range, function(i) {
     Pi[i,t,]*2^(-1/Theta[i,t,])
  })),2,mean))
  est_tail_lower <-  rbind(est_tail_lower,apply(do.call(rbind,lapply(range, function(i) {
    Pi[i,t,]*2^(-1/Theta[i,t,])
  })),2, quantile, probs=0.025))
  est_tail_upper <-  rbind(est_tail_upper,apply(do.call(rbind,lapply(range, function(i) {
    Pi[i,t,]*2^(-1/Theta[i,t,])
  })), 2, quantile, probs=0.975))
}
Date <- seq(as.Date("2017/01/01"), as.Date("2019/12/31"),by="months")
df_tail <- data.frame(
  Time = rep(1:TT, 4),
  Date = rep(Date[1:TT], 4),
  Mean = as.vector(est_tail_mean),
  Lower = as.vector(est_tail_lower),
  Upper = as.vector(est_tail_upper),
  Lambda = rep(paste0("Lambda", 1:4), each = TT)
)

df_tail$Lambda <- factor(df_tail$Lambda, levels = c("Lambda4", "Lambda3", "Lambda1", "Lambda2"))
library(ggplot2)
library(dplyr)
library(tidyr)
library(RColorBrewer)
custom_colors <- brewer.pal(n = 4, name = "Set1")
# Plot using ggplot2
fig_tail <- ggplot(df_tail, aes(x = Date)) +
  geom_ribbon(aes(ymin = Lower, ymax = Upper, fill = Lambda), alpha = 0.2) +
  geom_line(aes(y = Mean, color = Lambda), size = 1.5) +
  scale_y_continuous(limits = c(0, 0.75)) +
  scale_x_date(date_labels = "%Y-%m", date_breaks = "6 months") +
  scale_color_manual(values = custom_colors) +  # Set custom line colors
  scale_fill_manual(values = custom_colors) + 
  facet_wrap(~ Lambda, scales = "free_y") +
  labs(x = "", y = "") +  # Remove x and y labels
  # theme_minimal(base_size = 18) + 
  theme_bw(base_size = 18)+# Adjust text size
  theme(axis.text.x = element_text(angle = 45, size = 18, hjust = 1,face = "bold"),
        axis.text.y = element_text(size = 18, face = "bold"),
        strip.text = element_text(size = 18,face = "bold"),
        axis.title.y = element_text(size = 18, face = "bold",margin = margin(r = 10)),
        legend.position = "none")
fig_tail
ggsave("Figures/tail_coeff.pdf", plot = fig_tail, width = 10, height = 10)




fig_tail <- ggplot(df_tail, aes(x = Date)) +
  geom_ribbon(aes(ymin = Lower, ymax = Upper, fill = Lambda), alpha = 0.2) +
  geom_line(aes(y = Mean, color = Lambda), size = 1.5) +
  scale_y_continuous(limits = c(0, 1)) +
  scale_x_date(date_labels = "%Y-%m", date_breaks = "6 months") +
  scale_color_manual(values = custom_colors) +  # Set custom line colors
  scale_fill_manual(values = custom_colors) +  
  labs(x = "", y = "") +  # Remove x and y labels
  # theme_minimal(base_size = 18) + 
  theme_bw(base_size = 18)+# Adjust text size
  theme(axis.text.x = element_text(angle = 45, size = 18, hjust = 1,face = "bold"),
        axis.text.y = element_text(size = 18, face = "bold"),
        strip.text = element_text(size = 18,face = "bold"),
        axis.title.y = element_text(size = 18, face = "bold",margin = margin(r = 10)),
        legend.position = c(0.7, 0.8), 
        legend.title = element_blank(),  # Moves the legend to top center inside the plot
        legend.justification = c(0.4, 0.4))
fig_tail
ggsave("Figures/tail_full.pdf", plot = fig_tail, width = 10, height = 10)
