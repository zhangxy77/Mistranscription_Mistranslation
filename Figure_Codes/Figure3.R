library(stringr)
library(tidyverse)
library(ggplot2)
library(cowplot)
library(ggpubr)
library(reshape2)
# ----------------------------------------------------------------------------------------- #
#                                             figure3. A                                   #
# ----------------------------------------------------------------------------------------- #
barcor_data_00 <- cbind(cor_00_data,0)
colnames(barcor_data_00)[[4]] <- "type" 
for(i in seq(1,nrow(barcor_data_00)))
{
  if(barcor_data_00[i,3]<=1e-05)
  {barcor_data_00[i,4] <- "1e-05"}
  if(barcor_data_00[i,3]>=1e-05 && barcor_data_00[i,3]<1.2e-05)
  {barcor_data_00[i,4] <- "1.2e-05"}
  if(barcor_data_00[i,3]>=1.2e-05 && barcor_data_00[i,3]<1.4e-05)
  {barcor_data_00[i,4] <- "1.4e-05"}
  if(barcor_data_00[i,3]>=1.4e-05 && barcor_data_00[i,3]<1.6e-05 )
  {barcor_data_00[i,4] <- "1.6e-05"}
  if(barcor_data_00[i,3]>=1.6e-05 && barcor_data_00[i,3]<3e-05 )
  {barcor_data_00[i,4] <- "3e-05"}
}
barcor_data_00$type <- factor(barcor_data_00$type, levels = c("1e-05","1.2e-05","1.4e-05","1.6e-05","3e-05"))
meannumber <- c("1e-05","1.2e-05","1.4e-05","1.6e-05","1.8e-05")
barcor_meandata_00 <- matrix(0,nrow = length(meannumber),ncol = 5)
colnames(barcor_meandata_00) <- c("gene","translation_error","transcription_error","type","se")
for(i in seq(1,length(meannumber)))
{
  t <- meannumber[i]
  tt <- barcor_data_00[barcor_data_00[,4]==t,]
  barcor_meandata_00[i,1] <- i
  barcor_meandata_00[i,2] <- mean(as.numeric(tt[,2]))
  barcor_meandata_00[i,3] <- mean(as.numeric(tt[,3]))
  barcor_meandata_00[i,4] <- t
  barcor_meandata_00[i,5] <- sd(tt[,2])/sqrt(length(tt[,2]))
}
barcor_meandata_00 <- as.data.frame(barcor_meandata_00)
barcor_meandata_00$type <- factor(barcor_meandata_00$type, levels = c("1e-05","1.2e-05","1.4e-05","1.6e-05","3e-05"))
barcor_meandata_00[,2] <- as.numeric(barcor_meandata_00[,2])
barcor_meandata_00[,3] <- as.numeric(barcor_meandata_00[,3])
barcor_meandata_00[,5] <- as.numeric(barcor_meandata_00[,5])
bar_100plot_00 <- ggplot(barcor_meandata_00, aes(x = as.numeric(type), y = translation_error)) +  
  geom_bar(stat = "identity", fill = "white", color = "grey35", width = 0.6, position = position_dodge(0.8)) + 
  geom_errorbar(aes(ymin = translation_error - se, ymax = translation_error + se), width = 0.2, color = "grey35", position = position_dodge(0.8), size = 1.5) +  
  geom_smooth(aes(x = as.numeric(type), y = translation_error), method = 'lm', se = FALSE, color = "black") +  
  scale_x_continuous(
    breaks = c(0.5, 1.5, 2.5, 3.5, 4.5,5.5), 
    labels = c("0.0","1.0", "1.2", "1.4", "1.6", "3.0"),
    limits = c(0.5, 5.7)
  ) +
  scale_y_continuous(labels = c(expression(italic(0)),   
                                expression(2),   
                                expression(4),   
                                expression(6)),   
                     expand = c(0, 0),   
                     breaks = c(0, 0.002, 0.004, 0.006),   
                     limits = c(0, 0.006)) +  
  labs(x = "", y = "", title = "") +  
  theme_bw(base_size = 16) +  
  theme_classic() +  
  theme(axis.text.x = element_text(size = 23), 
        axis.text.y = element_text(size = 23))
bar_100plot_00
# ----------------------------------------------------------------------------------------- #
#                                             figure3. B                                   #
# ----------------------------------------------------------------------------------------- #
barcor_data_00 <- cbind(cor_00_data,0)
colnames(barcor_data_00)[[4]] <- "type" 
for(i in seq(1,nrow(barcor_data_00)))
{
  if(barcor_data_00[i,3]<=5e-06)
  {barcor_data_00[i,4] <- "5e-06"}
  if(barcor_data_00[i,3]>=5e-06 && barcor_data_00[i,3]<8e-06)
  {barcor_data_00[i,4] <- "8e-06"}
  if(barcor_data_00[i,3]>=8e-06 && barcor_data_00[i,3]<1.1e-05)
  {barcor_data_00[i,4] <- "1.1e-05"}
  if(barcor_data_00[i,3]>=1.1e-05 && barcor_data_00[i,3]<1.4e-05)
  {barcor_data_00[i,4] <- "1.4e-05"}
  if(barcor_data_00[i,3]>=1.4e-05 && barcor_data_00[i,3]<1.7e-05)
  {barcor_data_00[i,4] <- "1.7e-05"}
}
barcor_data_00$type <- factor(barcor_data_00$type, levels = c("5e-06","8e-06","1.1e-05","1.4e-05","1.7e-05"))
meannumber <- c("5e-06","8e-06","1.1e-05","1.4e-05","1.7e-05")
barcor_meandata_00 <- matrix(0,nrow = length(meannumber),ncol = 5)
colnames(barcor_meandata_00) <- c("gene","translation_error","transcription_error","type","se")
for(i in seq(1,length(meannumber)))
{
  t <- meannumber[i]
  tt <- barcor_data_00[barcor_data_00[,4]==t,]
  barcor_meandata_00[i,1] <- i
  barcor_meandata_00[i,2] <- mean(as.numeric(tt[,2]))
  barcor_meandata_00[i,3] <- mean(as.numeric(tt[,3]))
  barcor_meandata_00[i,4] <- t
  barcor_meandata_00[i,5] <- sd(tt[,2])/sqrt(length(tt[,2]))
}
barcor_meandata_00 <- as.data.frame(barcor_meandata_00)
barcor_meandata_00$type <- factor(barcor_meandata_00$type, levels = c("5e-06","8e-06","1.1e-05","1.4e-05","1.7e-05"))
barcor_meandata_00[,2] <- as.numeric(barcor_meandata_00[,2])
barcor_meandata_00[,3] <- as.numeric(barcor_meandata_00[,3])
barcor_meandata_00[,5] <- as.numeric(barcor_meandata_00[,5])
bar_100plot_00 <- ggplot(barcor_meandata_00, aes(x = as.numeric(type), y = translation_error)) +  
  geom_bar(stat = "identity", fill = "white", color = "grey35", width = 0.6, position = position_dodge(0.8)) + 
  geom_errorbar(aes(ymin = translation_error - se, ymax = translation_error + se), width = 0.2, color = "grey35", position = position_dodge(0.8), size = 1.5) +  
  geom_smooth(aes(x = as.numeric(type), y = translation_error), method = 'lm', se = FALSE, color = "black") + 
  scale_x_continuous(
    breaks = c(0.5, 1.5, 2.5, 3.5, 4.5,5.5), 
    labels = c("0.0","0.5", "0.8", "1.1", "1.4","1.7"),
    limits = c(0.5, 5.7)
  ) +
  scale_y_continuous(labels = c(expression(italic(0)),   
                                expression(2),   
                                expression(4),   
                                expression(6)),   
                     expand = c(0, 0),   
                     breaks = c(0, 0.002, 0.004, 0.006),   
                     limits = c(0, 0.006)) +  
  labs(x = "", y = "", title = "") +  
  theme_bw(base_size = 16) +  
  theme_classic() +  
  theme(axis.text.x = element_text(size = 23), 
        axis.text.y = element_text(size = 23))

bar_100plot_00

# ----------------------------------------------------------------------------------------- #
#                                             figure3. C                                   #
# ----------------------------------------------------------------------------------------- #
barcor_data_00 <- cbind(cor_00_data,0)
colnames(barcor_data_00)[[4]] <- "type" 
for(i in seq(1,nrow(barcor_data_00)))
{
  if(barcor_data_00[i,3]<=5e-06)
  {barcor_data_00[i,4] <- "5e-06"}
  if(barcor_data_00[i,3]>=5e-06 && barcor_data_00[i,3]<1e-05)
  {barcor_data_00[i,4] <- "1e-05"}
  if(barcor_data_00[i,3]>=1e-05 && barcor_data_00[i,3]<1.5e-05)
  {barcor_data_00[i,4] <- "1.5e-05"}
  if(barcor_data_00[i,3]>=1.5e-05 && barcor_data_00[i,3]<2e-05)
  {barcor_data_00[i,4] <- "2e-05"}
  if(barcor_data_00[i,3]>=2e-05 && barcor_data_00[i,3]<3e-05)
  {barcor_data_00[i,4] <- "3e-05"}
}
barcor_data_00$type <- factor(barcor_data_00$type, levels = c("5e-06","1e-05","1.5e-05","2e-05","3e-05"))
meannumber <- c("5e-06","1e-05","1.5e-05","2e-05","3e-05")
barcor_meandata_00 <- matrix(0,nrow = length(meannumber),ncol = 5)
colnames(barcor_meandata_00) <- c("gene","translation_error","transcription_error","type","se")
for(i in seq(1,length(meannumber)))
{
  t <- meannumber[i]
  tt <- barcor_data_00[barcor_data_00[,4]==t,]
  barcor_meandata_00[i,1] <- i
  barcor_meandata_00[i,2] <- mean(as.numeric(tt[,2]))
  barcor_meandata_00[i,3] <- mean(as.numeric(tt[,3]))
  barcor_meandata_00[i,4] <- t
  barcor_meandata_00[i,5] <- sd(tt[,2])/sqrt(length(tt[,2]))
}
barcor_meandata_00 <- as.data.frame(barcor_meandata_00)
barcor_meandata_00$type <- factor(barcor_meandata_00$type, levels = c("5e-06","1e-05","1.5e-05","2e-05","3e-05"))
barcor_meandata_00[,2] <- as.numeric(barcor_meandata_00[,2])
barcor_meandata_00[,3] <- as.numeric(barcor_meandata_00[,3])
barcor_meandata_00[,5] <- as.numeric(barcor_meandata_00[,5])
bar_100plot_00 <- ggplot(barcor_meandata_00, aes(x = as.numeric(type), y = translation_error)) +  
  geom_bar(stat = "identity", fill = "white", color = "grey35", width = 0.6, position = position_dodge(0.8)) +
  geom_errorbar(aes(ymin = translation_error - se, ymax = translation_error + se), width = 0.2, color = "grey35", position = position_dodge(0.8), size = 1.5) +  
  geom_smooth(aes(x = as.numeric(type), y = translation_error), method = 'lm', se = FALSE, color = "black") +  
  scale_x_continuous(
    breaks = c(0.5, 1.5, 2.5, 3.5, 4.5,5.5), 
    labels = c("0.0","0.5", "1", "1.5", "2.0","3.0"),
    limits = c(0.5, 5.7)
  ) +
  scale_y_continuous(labels = c(expression(italic(0)), 
                                expression(1), 
                                expression(2), 
                                expression(3)), 
                     expand = c(0,0), 
                     breaks = c(0,0.001,0.002,0.003), 
                     limits = c(0,0.0035)) +  
  labs(x = "", y = "", title = "") +  
  theme_bw(base_size = 16) +  
  theme_classic() +  
  theme(axis.text.x = element_text(size = 23), 
        axis.text.y = element_text(size = 23))
bar_100plot_00

# ----------------------------------------------------------------------------------------- #
#                                             figure3. D                                   #
# ----------------------------------------------------------------------------------------- #
barcor_data_00 <- cbind(cor_00_data,0)
colnames(barcor_data_00)[[4]] <- "type" 
for(i in seq(1,nrow(barcor_data_00)))
{
  if(barcor_data_00[i,3]<=8e-06)
  {barcor_data_00[i,4] <- "8e-06"}
  if(barcor_data_00[i,3]>=8e-06 && barcor_data_00[i,3]<1.1e-05)
  {barcor_data_00[i,4] <- "1.1e-05"}
  if(barcor_data_00[i,3]>=1.1e-05 && barcor_data_00[i,3]<1.4e-05)
  {barcor_data_00[i,4] <- "1.4e-05"}
  if(barcor_data_00[i,3]>=1.4e-05 && barcor_data_00[i,3]<1.7e-05)
  {barcor_data_00[i,4] <- "1.7e-05"}
  if(barcor_data_00[i,3]>=1.7e-05 && barcor_data_00[i,3]<4e-05)
  {barcor_data_00[i,4] <- "4e-05"}
}
barcor_data_00$type <- factor(barcor_data_00$type, levels = c("8e-06","1.1e-05","1.4e-05","1.7e-05","4e-05"))
meannumber <- c("8e-06","1.1e-05","1.4e-05","1.7e-05","4e-05")
barcor_meandata_00 <- matrix(0,nrow = length(meannumber),ncol = 5)
colnames(barcor_meandata_00) <- c("gene","translation_error","transcription_error","type","se")
for(i in seq(1,length(meannumber)))
{
  t <- meannumber[i]
  tt <- barcor_data_00[barcor_data_00[,4]==t,]
  barcor_meandata_00[i,1] <- i
  barcor_meandata_00[i,2] <- mean(as.numeric(tt[,2]))
  barcor_meandata_00[i,3] <- mean(as.numeric(tt[,3]))
  barcor_meandata_00[i,4] <- t
  barcor_meandata_00[i,5] <- sd(tt[,2])/sqrt(length(tt[,2]))
}
barcor_meandata_00 <- as.data.frame(barcor_meandata_00)
barcor_meandata_00$type <- factor(barcor_meandata_00$type, levels = c("8e-06","1.1e-05","1.4e-05","1.7e-05","4e-05"))
barcor_meandata_00[,2] <- as.numeric(barcor_meandata_00[,2])
barcor_meandata_00[,3] <- as.numeric(barcor_meandata_00[,3])
barcor_meandata_00[,5] <- as.numeric(barcor_meandata_00[,5])
barcor_meandata_00[5,5] <- sd(barcor_meandata_00$translation_error)
bar_100plot_00 <- ggplot(barcor_meandata_00, aes(x = as.numeric(type), y = translation_error)) +  
  geom_bar(stat = "identity", fill = "white", color = "grey35", width = 0.6, position = position_dodge(0.8)) +
  geom_errorbar(aes(ymin = translation_error - se, ymax = translation_error + se), width = 0.2, color = "grey35", position = position_dodge(0.8), size = 1.5) +  
  geom_smooth(aes(x = as.numeric(type), y = translation_error), method = 'lm', se = FALSE, color = "black") + 
  scale_x_continuous(
    breaks = c(0.5, 1.5, 2.5, 3.5, 4.5,5.5),
    labels = c("0.0","0.8", "1.1", "1.4", "1.7","4.0"), 
    limits = c(0.5, 5.7)
  ) +
  scale_y_continuous(labels = c(expression(italic(0)),   
                                expression(2),   
                                expression(4),   
                                expression(6)),   
                     expand = c(0, 0),   
                     breaks = c(0, 0.002, 0.004, 0.006),   
                     limits = c(0, 0.006)) +  
  labs(x = "", y = "", title = "") +  
  theme_bw(base_size = 16) +  
  theme_classic() +  
  theme(axis.text.x = element_text(size = 23),
        axis.text.y = element_text(size = 23))
bar_100plot_00

# ----------------------------------------------------------------------------------------- #
#                                             figure3. E                                   #
# ----------------------------------------------------------------------------------------- #
barcor_data_00 <- cbind(cor_00_data,0)
colnames(barcor_data_00)[[4]] <- "type" 
for(i in seq(1,nrow(barcor_data_00)))
{
  if(barcor_data_00[i,3]<=2e-05)
  {barcor_data_00[i,4] <- "2e-05"}
  if(barcor_data_00[i,3]>=2e-05 && barcor_data_00[i,3]<2.2e-05)
  {barcor_data_00[i,4] <- "2.2e-05"}
  if(barcor_data_00[i,3]>=2.2e-05 && barcor_data_00[i,3]<2.4e-05)
  {barcor_data_00[i,4] <- "2.4e-05"}
  if(barcor_data_00[i,3]>=2.4e-05 && barcor_data_00[i,3]<2.6e-05)
  {barcor_data_00[i,4] <- "2.6e-05"}
  if(barcor_data_00[i,3]>=2.6e-05 && barcor_data_00[i,3]<2.8e-05)
  {barcor_data_00[i,4] <- "2.8e-05"}
  if(barcor_data_00[i,3]>=2.8e-05 && barcor_data_00[i,3]<3e-05)
  {barcor_data_00[i,4] <- "3e-05"}
}
barcor_data_00$type <- factor(barcor_data_00$type, levels = c("2e-05","2.2e-05","2.4e-05","2.6e-05","2.8e-05","3e-05"))
meannumber <- c("2e-05","2.2e-05","2.4e-05","2.6e-05","2.8e-05","3e-05")
barcor_meandata_00 <- matrix(0,nrow = length(meannumber),ncol = 5)
colnames(barcor_meandata_00) <- c("gene","translation_error","transcription_error","type","se")
for(i in seq(1,length(meannumber)))
{
  t <- meannumber[i]
  tt <- barcor_data_00[barcor_data_00[,4]==t,]
  barcor_meandata_00[i,1] <- i
  barcor_meandata_00[i,2] <- mean(as.numeric(tt[,2]))
  barcor_meandata_00[i,3] <- mean(as.numeric(tt[,3]))
  barcor_meandata_00[i,4] <- t
  barcor_meandata_00[i,5] <- sd(tt[,2])/sqrt(length(tt[,2]))
}
barcor_meandata_00 <- as.data.frame(barcor_meandata_00)
barcor_meandata_00$type <- factor(barcor_meandata_00$type, levels = c("2e-05","2.2e-05","2.4e-05","2.6e-05","2.8e-05","3e-05"))
barcor_meandata_00[,2] <- as.numeric(barcor_meandata_00[,2])
barcor_meandata_00[,3] <- as.numeric(barcor_meandata_00[,3])
barcor_meandata_00[,5] <- as.numeric(barcor_meandata_00[,5])
barcor_meandata_00[2,5] <- sd(barcor_meandata_00$translation_error)
bar_100plot_00 <- ggplot(barcor_meandata_00, aes(x = as.numeric(type), y = translation_error)) +  
  geom_bar(stat = "identity", fill = "white", color = "grey35", width = 0.6, position = position_dodge(0.8)) +  
  geom_errorbar(aes(ymin = translation_error - se, ymax = translation_error + se), width = 0.2, color = "grey35", position = position_dodge(0.8), size = 1.5) +  
  geom_smooth(aes(x = as.numeric(type), y = translation_error), method = 'lm', se = FALSE, color = "black") +   
  scale_x_continuous(
    breaks = c(0.5, 1.5, 2.5, 3.5, 4.5,5.5,6.5), 
    labels = c("0.0","2.0", "2.2", "2.4", "2.6","2.8", "3.0"), 
    limits = c(0.5, 6.7)
  ) +
  scale_y_continuous(labels = c(expression(italic(0)),   
                                expression(2),   
                                expression(4),   
                                expression(6)),   
                     expand = c(0, 0),   
                     breaks = c(0, 0.002, 0.004, 0.006),   
                     limits = c(0, 0.006)) +  
  labs(x = "", y = "", title = "") +  
  theme_bw(base_size = 16) +  
  theme_classic() +  
  theme(axis.text.x = element_text(size = 23), 
        axis.text.y = element_text(size = 23))

bar_100plot_00
# ----------------------------------------------------------------------------------------- #
#                                             figure3. F                                    #
# ----------------------------------------------------------------------------------------- #
##PLOT-bar-100
OR_m_100 <- list(OR_H_m_100,OR_Y_2_m_100,OR_M_2_m_100,OR_C_m_100,OR_D_m_100)
OR_m_100_plot <- data.frame(Species=c("Human","Yeast","Mouse","C.elegans","Drosophila"),
                            logOR=0,pvalue=0,CI_low=0,CI_high=0)
for(i in seq(1,length(OR_m_100)))
{
  OR_m_100_plot[i,2] <- log(OR_m_100[[i]]$measure[[2]])
  OR_m_100_plot[i,3] <- OR_m_100[[i]]$p.value[[4]]
  OR_m_100_plot[i,4] <- log(OR_m_100[[i]]$measure[[4]])
  OR_m_100_plot[i,5] <- log(OR_m_100[[i]]$measure[[6]])
}
OR_m_100_plot$pSig <- 0
for(i in seq(1,length(OR_m_100)))
{
  if(OR_m_100_plot[i,3] < 0.001){OR_m_100_plot[i,6] <- 0.001}
  if(OR_m_100_plot[i,3] > 0.001 && OR_m_100_plot[i,3] < 0.01){OR_m_100_plot[i,6] <- 0.01}
  if(OR_m_100_plot[i,3] > 0.01 &&OR_m_100_plot[i,3] < 0.05){OR_m_100_plot[i,6] <- 0.05}
}
OR_m_100_plot$Species <- factor(OR_m_100_plot$Species, 
                                levels = c("Human", "Mouse", "Drosophila", "C.elegans", "Yeast"),
                                labels = c("*H. sapiens*", "*M. musculus*", "*D. melanogaster*", 
                                           "*C. elegans*", "*S. cerevisiae*"))
OR_100_plot <- ggplot(OR_m_100_plot, aes(x = logOR, y = Species)) +
  geom_bar(stat = "identity", fill = "white", color = "grey35", width = 0.6, position = position_dodge(0.8)) +
  geom_errorbar(aes(xmin = CI_low, xmax = CI_high), width = 0.2, color = "black", position = position_dodge(0.8), linewidth = 1.5) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "black",linewidth = 1) +
  scale_y_discrete(position = "left") + # Move y-axis to the left
  scale_x_continuous(labels = c(expression(10),
                                expression(1),   
                                expression(0.01),
                                expression("0.0001"),
                                expression("0.000001")),   
                     expand = c(0, 0),   
                     breaks = c(1,0,-2,-4,-6),   
                     limits = c(-6.5, 1.5))+
  labs(title = "", x = "Odds Ratio", y = "Species") +
  theme_bw(base_size = 16) +
  theme_classic() +
  theme(
    axis.text.x = element_text(size = 28), 
    axis.title.y = element_text(size = 35, margin = margin(r = 20)), 
    axis.title.x = element_text(size = 35), 
    axis.text.y = element_markdown(size = 30),
    plot.title = element_text(size = 30, face = "bold", hjust = 0.5)
  )

OR_100_plot