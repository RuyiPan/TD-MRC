Parameters <- rbind(expand.grid(c(1),c(0,1,3,5,10,20,30,40), c(0,1,2,3,4,5,6,7)),
                    expand.grid(c(1000), c(30),c(0,1,2,3,4,5,6,7)))
job_num <- which(apply(Parameters, 1, function(row) all(row == c(1,30,3))))

temp_res <- readRDS(paste0("MixC_diff_parameters_new/job_name=MixClayton_simulation_diff_parametersjob_num=",job_num,"LPS_mix_new.rds"))

burn.in <- 60
B <- 140
batch.size <- 50
range <- seq(burn.in*batch.size,B*batch.size,10)
Theta <- temp_res$Theta
Pi <- temp_res$Pi
Theta_mean <- NULL
Pi_mean <- NULL
for (j in c(1:4)) {
  t_mean <- rbind(colMeans(Theta[range,,j]), 
                  apply(Theta[range,,j], 2, quantile, probs=c(0.025, 0.975)))
  
  Theta_mean <- rbind(Theta_mean,
                      c(list(t_mean), j))
  
  p_mean <- rbind(colMeans(Pi[range,,j]), 
                  apply(Pi[range,,j], 2, quantile, probs=c(0.025, 0.975)))
  
  Pi_mean <- rbind(Pi_mean,
                   c(list(p_mean),j))
  
}
Theta_mean[1,]
Theta_mean[2,]
Theta_mean[3,]
Pi_mean[1,]
TT <- 20
a<-c(0.4,0.25,0.1,0.25)
Pi_true <- c(a[1]*(0.95)^{c(0:19)},
             a[2]*(1.05)^{c(0:19)},
             rep(a[3], 20),
             (1-a[3]-a[1]*(0.95)^{c(0:19)}-a[2]*(1.05)^{c(0:19)}))
th<-c(5,3,4,3) 
df_Pi <- data.frame(
  Time = rep(1:TT, 4),
  Date = rep(c(1:20), 4),
  Mean = unlist(lapply(1:4, function(k) Pi_mean[k][[1]][1, ])),
  Lower = unlist(lapply(1:4, function(k) Pi_mean[k][[1]][2, ])),
  Upper = unlist(lapply(1:4, function(k) Pi_mean[k][[1]][3, ])),
  Pi = rep(paste0("Pi", 1:4), each = TT),
  Pi_true = Pi_true 
)

df_Pi$Pi <- factor(df_Pi$Pi, levels = c("Pi4", "Pi3", "Pi1", "Pi2"))
library(ggplot2)
library(dplyr)
library(tidyr)
library(RColorBrewer)
custom_colors <- brewer.pal(n = 4, name = "Set1")
# Plot using ggplot2
fig_pi <- df_Pi %>% ggplot( aes(x = Time)) +
  geom_ribbon(aes(ymin = Lower, ymax = Upper, fill = Pi), alpha = 0.2) +
  geom_line(aes(y = Mean, color = Pi), size = 1.5) +
  scale_y_continuous(limits = c(0, 0.7)) +
  scale_x_continuous(breaks=seq(1,19,by=2)) +
  facet_wrap(~ Pi, scales = "free_y") +
  scale_color_manual(values = custom_colors) +  # Set custom line colors
  scale_fill_manual(values = custom_colors) +  
  labs(x = "", y = "") + 
  # Remove x and y labels
  # theme_minimal(base_size = 18) + 
  theme_bw(base_size = 18)+# Adjust text size
  theme(axis.text.x = element_text(angle = 0, size = 18,face = "bold"),
        axis.text.y = element_text(size = 18, face = "bold"),
        strip.text = element_text(size = 18,face = "bold"),
        axis.title.y = element_text(size = 18, face = "bold",margin = margin(r = 10)),
        legend.position = "none")+
  geom_line(aes(x = Time, y = Pi_true), color = "black", linetype = "dashed", size = 1.5)
fig_pi
ggsave("Figures/Pi.pdf", plot = fig_pi, width = 10, height = 10)


df_Theta <- data.frame(
  Time = rep(1:TT, 4),
  Date = rep(c(1:20), 4),
  Mean = unlist(lapply(1:4, function(k) Theta_mean[k][[1]][1, ])),
  Lower = unlist(lapply(1:4, function(k) Theta_mean[k][[1]][2, ])),
  Upper = unlist(lapply(1:4, function(k) Theta_mean[k][[1]][3, ])),
  Theta = rep(paste0("Theta", 1:4), each = TT),
  Theta_true=rep(th, each=TT)
)
df_Theta$Theta <- factor(df_Theta$Theta, levels = c("Theta4", "Theta3", "Theta1", "Theta2"))
library(ggplot2)
library(dplyr)
library(tidyr)
library(RColorBrewer)
custom_colors <- brewer.pal(n = 4, name = "Set1")
# Plot using ggplot2
fig_Theta <-df_Theta %>% ggplot( aes(x = Date)) +
  geom_ribbon(aes(ymin = Lower, ymax = Upper, fill = Theta), alpha = 0.2) +
  geom_line(aes(y = Mean, color = Theta), size = 1.5) +
  scale_y_continuous(limits = c(0, 17)) +
  scale_x_continuous(breaks=seq(1,19,by=2)) +
  facet_wrap(~ Theta, scales = "free_y") +
  scale_color_manual(values = custom_colors) +  # Set custom line colors
  scale_fill_manual(values = custom_colors) +  
  labs(x = "", y = "") +  # Remove x and y labels
  # theme_minimal(base_size = 18) + 
  theme_bw(base_size = 18)+# Adjust text size
  theme(axis.text.x = element_text(angle = 0, size = 18,face = "bold"),
        axis.text.y = element_text(size = 18, face = "bold"),
        strip.text = element_text(size = 18,face = "bold"),
        axis.title.y = element_text(size = 18, face = "bold",margin = margin(r = 10)),
        legend.position = "none")+
  geom_line(aes(x = Time, y = Theta_true), color = "black", linetype = "dashed", size = 1.5)
fig_Theta
ggsave("Figures/Theta.pdf", plot = fig_Theta, width = 10, height = 10)



#kappa
#acceptance rate
df_kappa <- data.frame(
  Date = 1:140,
  kappa=temp_res$kappa.all[, 1, 1]
)

fig_kappa <- df_kappa  %>% 
  ggplot(aes(x = Date, y = kappa)) +
  geom_line(linewidth = 1, alpha=1) +
  labs(y = "Kappa", x = "") +
  scale_x_continuous(breaks = c(0:7)*20)+
  theme_bw(base_size = 18)+# Adjust text size
  theme(axis.text.x = element_text(size = 18, face = "bold"),
        axis.text.y = element_text(size = 18, face = "bold"),
        strip.text = element_text(size = 18,face = "bold"),
        axis.title.y = element_text(size = 18, face = "bold",margin = margin(r = 10)),
        legend.position = "none")

fig_kappa

ggsave("Figures/kappa.pdf", plot = fig_kappa, width = 6, height = 6)
df_acc <- data.frame(
  Date = 1:140,
  acc=temp_res$acc.all[, 1, 1]
)

fig_acc<- df_acc  %>% 
  ggplot(aes(x = Date, y = acc)) +
  geom_line(linewidth = 1, alpha=1) +
  labs(y = "acceptance rate", x = "") +
  scale_x_continuous(breaks = c(0:7)*20)+
  theme_bw(base_size = 18)+# Adjust text size
  theme(axis.text.x = element_text(size = 18, face = "bold"),
        axis.text.y = element_text(size = 18, face = "bold"),
        strip.text = element_text(size = 18,face = "bold"),
        axis.title.y = element_text(size = 18, face = "bold",margin = margin(r = 10)),
        legend.position = "none")+
  geom_hline(yintercept = 0.3, linetype = "solid", color = "red", size = 1) +  # Line at y = 0.2
  geom_hline(yintercept = 0.4, linetype = "solid", color = "red", size = 1)  

fig_acc
ggsave("Figures/acc.pdf", plot = fig_acc, width = 6, height = 6)
