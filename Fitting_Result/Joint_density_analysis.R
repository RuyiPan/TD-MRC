source("../MixClayton/mcmc_helper.R")
Parameters <- expand.grid(c(0,1,2,3,4,5),c(0,1,2,3,4), c(0,1,2))
job_num <- which(apply(Parameters, 1, function(row) all(row == c(1,1,2))))

temp_res <- readRDS(paste0("Result_new/job_name=MixClayton_fitting_realjob_num=",job_num,"LPS_mix_new.rds"))

burn.in <- 100
B <- 200
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

df_Pi <- data.frame(
  Time = rep(1:TT, 4),
  Date = rep(Date[1:TT], 4),
  Mean = unlist(lapply(1:4, function(k) Pi_mean[k][[1]][1, ])),
  Lower = unlist(lapply(1:4, function(k) Pi_mean[k][[1]][2, ])),
  Upper = unlist(lapply(1:4, function(k) Pi_mean[k][[1]][3, ])),
  Pi = rep(paste0("Pi", 1:4), each = TT)
)
df_Pi$Pi <- factor(df_Pi$Pi, levels = c("Pi4", "Pi3", "Pi1", "Pi2"))
library(ggplot2)
library(dplyr)
library(tidyr)
library(RColorBrewer)
custom_colors <- brewer.pal(n = 4, name = "Set1")
# Plot using ggplot2
fig_pi <- ggplot(df_Pi, aes(x = Date)) +
  geom_ribbon(aes(ymin = Lower, ymax = Upper, fill = Pi), alpha = 0.2) +
  geom_line(aes(y = Mean, color = Pi), size = 1.5) +
  scale_y_continuous(limits = c(0, 1)) +
  facet_wrap(~ Pi, scales = "free_y") +
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
        legend.position = "none")
 
ggsave("Figures/Pi.pdf", plot = fig_pi, width = 10, height = 10)


df_Theta <- data.frame(
  Time = rep(1:TT, 4),
  Date = rep(Date[1:TT], 4),
  Mean = unlist(lapply(1:4, function(k) Theta_mean[k][[1]][1, ])),
  Lower = unlist(lapply(1:4, function(k) Theta_mean[k][[1]][2, ])),
  Upper = unlist(lapply(1:4, function(k) Theta_mean[k][[1]][3, ])),
  Theta = rep(paste0("Theta", 1:4), each = TT)
)
df_Theta$Theta <- factor(df_Theta$Theta, levels = c("Theta4", "Theta3", "Theta1", "Theta2"))
library(ggplot2)
library(dplyr)
library(tidyr)
library(RColorBrewer)
custom_colors <- brewer.pal(n = 4, name = "Set1")
# Plot using ggplot2
fig_Theta <- ggplot(df_Theta, aes(x = Date)) +
  geom_ribbon(aes(ymin = Lower, ymax = Upper, fill = Theta), alpha = 0.2) +
  geom_line(aes(y = Mean, color = Theta), size = 1.5) +
  scale_y_continuous(limits = c(0, 4)) +
  facet_wrap(~ Theta, scales = "free_y") +
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
        legend.position = "none")

ggsave("Figures/Theta.pdf", plot = fig_Theta, width = 10, height = 10)

Date <- seq(as.Date("2017/01/01"), as.Date("2019/12/31"),by="months")
format(Date, "%Y-%m")
u1 <- seq(0, 1, length.out=100)
u2 <- seq(0, 1, length.out=100)
Den <- list()
for (t in 1:TT) {
  temp <- NULL
  for (job_num in 0:19) {
    
    temp_res <- readRDS(paste0("Joint_density/job_name=Real_densityjob_num=",job_num,
                               "estimated_joint_dist.rds"))
    temp <- c(temp, temp_res$den_all[[t]])
  }
  Den[[t]] <- temp
}

for (t in 1:TT) {
  p <- expand.grid(PM2.5 = u1, O3 = u2) %>%
    mutate(Z = Den[[t]]) %>%
    ggplot(aes(PM2.5, O3, z = Z)) +
    geom_contour() +
    ggtitle(format(Date[t], "%Y-%m"))+
    geom_contour_filled(breaks = seq(0,3,by=0.5) ) +
    scale_fill_viridis_d(drop = FALSE) +
    theme_bw() +
    theme(plot.title = element_text(hjust = 0.5, size=40,face="bold"),aspect.ratio = 1) +
    theme(axis.text=element_text(size=22,face="bold"),
          axis.title=element_text(size=28,face="bold")) +
    theme(legend.position="none")
  
  fname <- paste0("Figures/den",t,".pdf")
  ggsave(fname, plot = p, width = 7.5, height = 7.5)
}

U_all <- readRDS("../Data/PM_O3_2017_2019.rds")
for (t in 1:TT) {
  temp <- as.data.frame(U_all[[t]])
  colnames(temp) <- c("PM2.5", "O3")
  p <- temp %>% ggplot(aes(x=PM2.5,y=O3))+
    geom_point(size=3)+
    ggtitle(format(Date[t], "%Y-%m"))+
    theme_bw() +
    theme(plot.title = element_text(hjust = 0.5, size=40,face="bold"),
          aspect.ratio = 1) +
    theme(axis.text=element_text(size=22,face="bold"),
          axis.title=element_text(size=28,face="bold"))
  fname <- paste0("Figures/real",t,".pdf")
  ggsave(fname, plot = p,width = 7.5, height = 7.5)
}
t <- 8 
# df_Pi %>% filter(Time==t)
as.vector(round(df_Pi %>% filter(Time==t) %>% select(Mean),2))

# df_Theta %>% filter(Time==t)
as.vector(round(df_Theta%>% filter(Time==t) %>% select(Mean),2))

