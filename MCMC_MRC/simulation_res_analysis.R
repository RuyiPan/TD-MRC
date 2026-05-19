burn_in <- 100
B <- 160
range <- seq(burn_in*batch.size, B*batch.size, 10)

Theta <- temp_res$Theta
Pi <- temp_res$Pi
Theta_mean <- NULL
Pi_mean <- NULL
comp_num <- 2^ncol(temp_res$U_train[[1]])
for (j in c(1:comp_num)) {
  t_mean <- rbind(colMeans(Theta[range,,j]), 
                  apply(Theta[range,,j], 2, quantile, probs=c(0.025, 0.975)))
  
  Theta_mean <- rbind(Theta_mean,
                      c(list(t_mean), j))
  
  p_mean <- rbind(colMeans(Pi[range,,j]), 
                  apply(Pi[range,,j], 2, quantile, probs=c(0.025, 0.975)))
  
  Pi_mean <- rbind(Pi_mean,
                   c(list(p_mean),j))
  
}

Theta_mean[2,]
Theta_mean[3,]
TT <- length(U_train)
a<-c(0.4,0.25,0.1,0.25)
# Pi_true <- c(rep(1,TT),
#              rep(0,TT),
#              rep(0,TT),
#              rep(0,TT),
#              rep(0,TT),
#              rep(0,TT),
#              rep(0,TT),
#              rep(0,TT))
# th<-c(3,0,0,0,0,0,0,0) 
Pi_true <- c(a[1]*(0.95)^{c(0:19)},
             a[2]*(1.05)^{c(0:19)},
             rep(a[3], 20),
             (1-a[3]-a[1]*(0.95)^{c(0:19)}-a[2]*(1.05)^{c(0:19)}))
th<-c(5,3,4,3) 
df_Pi <- data.frame(
  Time = rep(1:TT, comp_num),
  Date = rep(c(1:TT), comp_num),
  Mean = unlist(lapply(1:comp_num, function(k) Pi_mean[k][[1]][1, ])),
  Lower = unlist(lapply(1:comp_num, function(k) Pi_mean[k][[1]][2, ])),
  Upper = unlist(lapply(1:comp_num, function(k) Pi_mean[k][[1]][3, ])),
  Pi = rep(paste0("Pi", 1:comp_num), each = TT),
  Pi_true = Pi_true 
)

df_Pi$Pi <- factor(df_Pi$Pi, levels =paste0("Pi",c(1:comp_num)))
library(ggplot2)
library(dplyr)
library(tidyr)
library(RColorBrewer)
custom_colors <- brewer.pal(n = comp_num, name = "Set1")
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


df_Theta <- data.frame(
  Time = rep(1:TT, comp_num),
  Date = rep(c(1:TT), comp_num),
  Mean = unlist(lapply(1:comp_num, function(k) Theta_mean[k][[1]][1, ])),
  Lower = unlist(lapply(1:comp_num, function(k) Theta_mean[k][[1]][2, ])),
  Upper = unlist(lapply(1:comp_num, function(k) Theta_mean[k][[1]][3, ])),
  Theta = rep(paste0("Theta", 1:comp_num), each = TT),
  Theta_true=rep(th, each=TT)
)
df_Theta$Theta <- factor(df_Theta$Theta, levels =paste0("Theta",c(1:comp_num)))
library(ggplot2)
library(dplyr)
library(tidyr)
library(RColorBrewer)
custom_colors <- brewer.pal(n = comp_num, name = "Set1")
# Plot using ggplot2
fig_Theta <-df_Theta %>% ggplot( aes(x = Date)) +
  geom_ribbon(aes(ymin = Lower, ymax = Upper, fill = Theta), alpha = 0.2) +
  geom_line(aes(y = Mean, color = Theta), size = 1.5) +
  scale_y_continuous(limits = c(-1, 15)) +
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



range <- seq(burn_in*batch.size,B*batch.size,10)
U_train <- temp_res$U_transformed_list
#WAIC
WAIC.value <- WAIC(U_train, Theta[range,,], Pi[range,,])
#LPML
LPML.value <- LPML(U_train, Theta[range,,], Pi[range,,])
#DIC
DIC.value <- DIC(U_train, Theta[range,,], Pi[range,,])

WAIC.value;DIC.value;LPML.value


apply(U_train[[t]][1:10,,], 1, function(row) sum(l.mg(row,Theta[1,1,])*Pi[1,1,]))
sum(l.mg(U_train[[t]][2,,],Theta[1,1,])*Pi[1,1,])
