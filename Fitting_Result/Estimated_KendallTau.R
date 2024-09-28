real_data <- readRDS("../Data/PM_O3_2017_2019.rds")
kendall.tau <- vector()
spearman.rho <- vector()
for (t in 1:36) {
  temp <- real_data[[t]]
  kendall.tau[t] <- cor(temp[,1], temp[,2], method = "kendall")
  spearman.rho[t] <- cor(temp[,1], temp[,2], method = "spearman")
}

lam2rho <- function(lam) {
  (exp(2*lam)-1)/(exp(2*lam)+1)
}

burn.in <- 100
B <- 200
batch.size <- 50
range <- seq(burn.in*batch.size,B*batch.size,10)
est_tau_gaussian <- vector()
temp <- readRDS("../Fitting_Result/Result_new/job_name=Gaussian_fittingjob_num=1fitting_all.rds")
for (t in 1:36) {
  est_tau_gaussian[t] <- mean((2/pi)*asin(lam2rho(temp$lam_all[range,t])))
}

temp$lam_all[range]

Theta <- temp_res$Theta
Pi <- temp_res$Pi
est_tau_mixC <- vector()
Parameters <- expand.grid(c(0,1,2,3,4,5),c(0,1,2,3,4), c(0,1,2))
job_num <- which(apply(Parameters, 1, function(row) all(row == c(1,1,2))))

temp_res <- readRDS(paste0("../Fitting_Result/Result_new/job_name=MixClayton_fitting_realjob_num=",job_num,"LPS_mix_new.rds"))

for (t in 1:36) {
  est_tau_mixC[t] <- mean(do.call(rbind,lapply(range, function(i) {
    Theta[i,t,1]/(2+Theta[i,t,1])*Pi[i, t,1]-
      Theta[i,t,2]/(2+Theta[i,t,2])*Pi[i, t,2] +
      Theta[i,t,3]/(2+Theta[i,t,3])*Pi[i, t,3]-
      Theta[i,t,4]/(2+Theta[i,t,4])*Pi[i, t,4]
  })))
}

plot(time_series_data, est_tau_mixC, type="l")
lines(time_series_data, kendall.tau, col="red")
#empirical with estiatmed values
df <- data.frame(
  Date = time_series_data,
  est_tau_mixC = est_tau_mixC,
  kendall_tau = kendall.tau,
  est_tau_Gaussian= est_tau_gaussian
)
library(RColorBrewer)
custom_colors <- brewer.pal(n = 3, name = "Set1")
fig_tau <- df %>% 
  pivot_longer(c(est_tau_mixC, kendall_tau, est_tau_Gaussian), names_to = "Method", values_to = "Tau") %>%
  mutate(Method=ifelse(Method=="est_tau_mixC", "D.Clayton Mixture", ifelse(Method=="est_tau_Gaussian", "D.Gaussian", "Empirical Value")))%>%
  ggplot(aes(x = Date, y = Tau, color = Method)) +
  geom_line(linewidth = 1, alpha=1) +
  scale_y_continuous(
    limits = c(-0.65, 0.65),            # Set y-axis limits
    breaks = c(-0.6, -0.4, -0.2, 0, 0.2, 0.4, 0.6)  # Custom y-axis labels
  ) +
  labs(y = "Kendall's tau", x = "") +
  scale_x_date(date_labels = "%Y-%m", date_breaks = "6 months") +
  theme_bw(base_size = 18)+# Adjust text size
  theme(axis.text.x = element_text(angle = 45, size = 18, hjust = 1,face = "bold"),
        axis.text.y = element_text(size = 18, face = "bold"),
        strip.text = element_text(size = 18,face = "bold"),
        axis.title.y = element_text(size = 18, face = "bold",margin = margin(r = 10)),
        legend.position = c(0.7, 0.8), 
        legend.title = element_blank(),  # Moves the legend to top center inside the plot
        legend.justification = c(0.4, 0.4))+
  geom_hline(yintercept = 0, linetype = "dashed", color = "black")+
  scale_color_manual(values = c("D.Clayton Mixture" = "#4DAF4A", 
                                "D.Gaussian" = "#377EB8", 
                                "Empirical Value" = "#E41A1C"))
fig_tau
ggsave("Figures/Kendall_Tau.pdf", plot = fig_tau, width = 6, height = 6)


df <- data.frame(
  Date = time_series_data,
  kendall.tau=kendall.tau
)
custom_colors <- brewer.pal(n = 3, name = "Set1")
fig_tau <- df %>% pivot_longer(c(kendall.tau), names_to = "Method", values_to = "Tau") %>%
  ggplot(aes(x = Date, y = Tau, color = Method)) +
  geom_line(linewidth = 1, alpha=1) +
  scale_y_continuous(
    limits = c(-0.65, 0.65),            # Set y-axis limits
    breaks = c(-0.6, -0.4, -0.2, 0, 0.2, 0.4, 0.6)  # Custom y-axis labels
  ) +
  labs(y = "Kendall's tau", x = "") +
  scale_x_date(date_labels = "%Y-%m", date_breaks = "6 months") +
  theme_bw(base_size = 18)+# Adjust text size
  theme(axis.text.x = element_text(angle = 45, size = 18, hjust = 1,face = "bold"),
        axis.text.y = element_text(size = 18, face = "bold"),
        strip.text = element_text(size = 18,face = "bold"),
        axis.title.y = element_text(size = 18, face = "bold",margin = margin(r = 10)),
        legend.position = "none")+
  geom_hline(yintercept = 0, linetype = "dashed", color = "black")+
  scale_color_manual(values = custom_colors) +  # Set custom line colors
  scale_fill_manual(values = custom_colors) 
fig_tau 
ggsave("Figures/Kendall_Tau.pdf", plot = fig_tau, width = 6, height = 6)

df <- data.frame(
  Date = time_series_data,
  spearman.rho=spearman.rho
)
custom_colors <- brewer.pal(n = 3, name = "Set1")
fig_rho <- df %>% pivot_longer(c(spearman.rho), names_to = "Method", values_to = "Tau") %>%
  ggplot(aes(x = Date, y = Tau, color = Method)) +
  geom_line(linewidth = 1, alpha=1) +
  scale_y_continuous(
    limits = c(-0.65, 0.65),            # Set y-axis limits
    breaks = c(-0.6, -0.4, -0.2, 0, 0.2, 0.4, 0.6)  # Custom y-axis labels
  ) +
  labs(y = "Spearman's rho", x = "") +
  scale_x_date(date_labels = "%Y-%m", date_breaks = "6 months") +
  theme_bw(base_size = 18)+# Adjust text size
  theme(axis.text.x = element_text(angle = 45, size = 18, hjust = 1,face = "bold"),
        axis.text.y = element_text(size = 18, face = "bold"),
        strip.text = element_text(size = 18,face = "bold"),
        axis.title.y = element_text(size = 18, face = "bold",margin = margin(r = 10)),
        legend.position = "none")+
  geom_hline(yintercept = 0, linetype = "dashed", color = "black")+
  scale_color_manual(values = custom_colors) +  # Set custom line colors
  scale_fill_manual(values = custom_colors) 
fig_rho 
ggsave("Figures/Spearman_Rho.pdf", plot = fig_rho, width = 6, height = 6)
