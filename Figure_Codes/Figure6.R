library(stringr)
library(tidyverse)
library(ggplot2)
library(cowplot)
library(ggpubr)
library(reshape2)
# ----------------------------------------------------------------------------------------- #
#                                             figure6. A-B                                  #
# ----------------------------------------------------------------------------------------- #

human_TL <- human_TL[order(human_TL$Error_Rate,decreasing = T),]
human_high_gene <- human_TL[human_TL$Error_Rate>median(human_TL$Error_Rate),1]
human_low_gene <- human_TL[human_TL$Error_Rate<median(human_TL$Error_Rate),1]
human_sm_high <- human_sm[rownames(human_sm)%in% human_high_gene,]
human_sm_low <- human_sm[rownames(human_sm)%in% human_low_gene,]
human_sm_high$sim_real <- human_sm_high$`sim_N/S` /human_sm_high$`real_N/S`
human_sm_low$sim_real <- human_sm_low$`sim_N/S`/human_sm_low$`real_N/S`


yeast_high_gene <- yeast_TL[yeast_TL$Error_Rate>median(yeast_TL$Error_Rate),1]
yeast_low_gene <- yeast_TL[yeast_TL$Error_Rate<median(yeast_TL$Error_Rate),1]
yeast_sm_high <- yeast_sm[rownames(yeast_sm)%in% yeast_high_gene,]
yeast_sm_low <- yeast_sm[rownames(yeast_sm)%in% yeast_low_gene,]
yeast_sm_high$sim_real <- yeast_sm_high$`sim_N/S`/yeast_sm_high$`real_N/S`
yeast_sm_low$sim_real <- yeast_sm_low$`sim_N/S` / yeast_sm_low$`real_N/S` 


mouse_high_gene <- mouse_TL[mouse_TL$Error.rate>median(mouse_TL$Error.rate),1]
mouse_low_gene <- mouse_TL[mouse_TL$Error.rate<median(mouse_TL$Error.rate),1]
mouse_sm_high <- mouse_sm[rownames(mouse_sm)%in% mouse_high_gene,]
mouse_sm_low <- mouse_sm[rownames(mouse_sm)%in% mouse_low_gene,]
mouse_sm_high$sim_real <- mouse_sm_high$`sim_N/S`/mouse_sm_high$`real_N/S`
mouse_sm_low$sim_real <- mouse_sm_low$`sim_N/S`/mouse_sm_low$`real_N/S`


elegans_high_gene <- elegans_TL[elegans_TL$Error_Rate>median(elegans_TL$Error_Rate),1]
elegans_low_gene <- elegans_TL[elegans_TL$Error_Rate<median(elegans_TL$Error_Rate),1]
elegans_high_gene <- gsub("-", ".", elegans_high_gene)
elegans_low_gene <- gsub("-", ".", elegans_low_gene)
elegans_sm_high <- elegans_sm[rownames(elegans_sm)%in% elegans_high_gene,]
elegans_sm_low <- elegans_sm[rownames(elegans_sm)%in% elegans_low_gene,]
elegans_sm_high$sim_real <- elegans_sm_high$`sim_N/S` / elegans_sm_high$`real_N/S` 
elegans_sm_low$sim_real <- elegans_sm_low$`sim_N/S` / elegans_sm_low$`real_N/S`


drosophila_high_gene <- drosophila_TL[drosophila_TL$Error_Rate>median(drosophila_TL$Error_Rate),1]
drosophila_low_gene <- drosophila_TL[drosophila_TL$Error_Rate<median(drosophila_TL$Error_Rate),1]
drosophila_sm_high <- drosophila_sm[rownames(drosophila_sm)%in% drosophila_high_gene,]
drosophila_sm_low <- drosophila_sm[rownames(drosophila_sm)%in% drosophila_low_gene,]
drosophila_sm_high$sim_real <- drosophila_sm_high$`sim_N/S` / drosophila_sm_high$`real_N/S`
drosophila_sm_low$sim_real <- drosophila_sm_low$`sim_N/S` / drosophila_sm_low$`real_N/S`


#high
sm_gene_data_high <- data.frame(species = c("Human", "Yeast", "Mouse", "Drosophila", "C.elegans"),
                                Real = c(mean(human_sm_high$`real_N/S`),mean(yeast_sm_high$`real_N/S`),mean(mouse_sm_high$`real_N/S`),mean(drosophila_sm_high$`real_N/S`),mean(elegans_sm_high$`real_N/S`)),
                                Simulated = c(mean(human_sm_high$`sim_N/S`),mean(yeast_sm_high$`sim_N/S`),mean(mouse_sm_high$`sim_N/S`),mean(drosophila_sm_high$`sim_N/S`),mean(elegans_sm_high$`sim_N/S`)))
sm_gene_data_high_dt <- as.data.table(sm_gene_data_high)
sm_gene_plot_high <- melt(sm_gene_data_high_dt, id="species", variable.name="Attribute", value.name = "NonSyn/Syn")
mean_gene_high <- aggregate(sm_gene_plot_high$`NonSyn/Syn`, by=list(sm_gene_plot_high $species, sm_gene_plot_high $Attribute), FUN=mean)
sd_gene_high <- c(sd(elegans_sm_high$`real_N/S`),sd(drosophila_sm_high$`real_N/S`),sd(human_sm_high$`real_N/S`),sd(mouse_sm_high$`real_N/S`),sd(yeast_sm_high$`real_N/S`),sd(elegans_sm_high$`sim_N/S`),sd(drosophila_sm_high$`sim_N/S`),sd(human_sm_high$`sim_N/S`),sd(mouse_sm_high$`sim_N/S`),sd(yeast_sm_high$`sim_N/S`))
len_gene_high <- c(length(elegans_sm_high$`real_N/S`),length(drosophila_sm_high$`real_N/S`),length(human_sm_high$`real_N/S`),length(mouse_sm_high$`real_N/S`),length(yeast_sm_high$`real_N/S`),length(elegans_sm_high$`sim_N/S`),length(drosophila_sm_high$`sim_N/S`),length(human_sm_high$`sim_N/S`),length(mouse_sm_high$`sim_N/S`),length(yeast_sm_high$`sim_N/S`))
sm_gene_plot_res_high <- data.frame(mean_gene_high, sd=sd_gene_high, len=len_gene_high)
colnames(sm_gene_plot_res_high) = c("Species", "Attribute", "Mean", "Sd", "Count")
sm_gene_plot_res_high$Se <- sm_gene_plot_res_high$Sd/sqrt(sm_gene_plot_res_high$Count)
sm_gene_plot_res_high$Species <- factor(sm_gene_plot_res_high$Species, levels = c("Human", "Mouse", "Drosophila", "C.elegans", "Yeast"), 
                                        labels = c("*H. sapiens*", "*M. musculus*", "*D. melanogaster*", 
                                                   "*C. elegans*", "*S. cerevisiae*"))
high_plot <- ggplot(sm_gene_plot_res_high, aes(x=Species, y=Mean, fill=Attribute)) +
  geom_bar(stat="identity", position=position_dodge(), color="black", width=.8) +
  geom_errorbar(aes(ymin=Mean-Se, ymax=Mean +Se), position=position_dodge(.8), width=.2) +
  labs(title = "High Mistranslation", x = "", y = "") +
  theme_bw() +
  scale_fill_manual(values=c("grey35", "gray")) + 
  scale_y_continuous(expand=c(0,0)) +
  ylim(0, 7) +
  theme(axis.text.x = element_markdown(size = 23,face = "bold",angle = 45,vjust = 1, hjust = 1),
        axis.title.y = element_text(size = 25,face = "bold"), 
        axis.title.x = element_text(size = 25,face = "bold"),
        axis.text.y = element_text(size = 23),
        plot.title = element_text(size=25,hjust=0.5,face = "bold"),
        legend.title = element_text(size = 18,face = "bold") ,
        legend.text = element_text(size = 18,face = "bold"),
        legend.position = "none")
high_plot

#low
sm_gene_data_low <- data.frame(species = c("Human", "Yeast", "Mouse", "Drosophila", "C.elegans"),
                               Real = c(mean(human_sm_low$`real_N/S`),mean(yeast_sm_low$`real_N/S`),mean(mouse_sm_low$`real_N/S`),mean(drosophila_sm_low$`real_N/S`),mean(elegans_sm_low$`real_N/S`)),
                               Simulated = c(mean(human_sm_low$`sim_N/S`),mean(yeast_sm_low$`sim_N/S`),mean(mouse_sm_low$`sim_N/S`),mean(drosophila_sm_low$`sim_N/S`),mean(elegans_sm_low$`sim_N/S`)))
sm_gene_data_low_dt <- as.data.table(sm_gene_data_low)
sm_gene_plot_low <- melt(sm_gene_data_low_dt, id="species", variable.name="Attribute", value.name = "NonSyn/Syn")
mean_gene_low <- aggregate(sm_gene_plot_low$`NonSyn/Syn`, by=list(sm_gene_plot_low $species, sm_gene_plot_low $Attribute), FUN=mean)
sd_gene_low <- c(sd(elegans_sm_low$`real_N/S`),sd(drosophila_sm_low$`real_N/S`),sd(human_sm_low$`real_N/S`),sd(mouse_sm_low$`real_N/S`),sd(yeast_sm_low$`real_N/S`),sd(elegans_sm_low$`sim_N/S`),sd(drosophila_sm_low$`sim_N/S`),sd(human_sm_low$`sim_N/S`),sd(mouse_sm_low$`sim_N/S`),sd(yeast_sm_low$`sim_N/S`))
len_gene_low <- c(length(elegans_sm_low$`real_N/S`),length(drosophila_sm_low$`real_N/S`),length(human_sm_low$`real_N/S`),length(mouse_sm_low$`real_N/S`),length(yeast_sm_low$`real_N/S`),length(elegans_sm_low$`sim_N/S`),length(drosophila_sm_low$`sim_N/S`),length(human_sm_low$`sim_N/S`),length(mouse_sm_low$`sim_N/S`),length(yeast_sm_low$`sim_N/S`))
sm_gene_plot_res_low <- data.frame(mean_gene_low, sd=sd_gene_low, len=len_gene_low)
colnames(sm_gene_plot_res_low) = c("Species", "Attribute", "Mean", "Sd", "Count")
sm_gene_plot_res_low$Se <- sm_gene_plot_res_low$Sd/sqrt(sm_gene_plot_res_low$Count)
sm_gene_plot_res_low$Species <- factor(sm_gene_plot_res_low$Species, levels = c("Human", "Mouse", "Drosophila", "C.elegans", "Yeast"),  
                                       labels = c("*H. sapiens*", "*M. musculus*", "*D. melanogaster*", 
                                                  "*C. elegans*", "*S. cerevisiae*"))
low_plot <- ggplot(sm_gene_plot_res_low, aes(x=Species, y=Mean, fill=Attribute)) +
  geom_bar(stat="identity", position=position_dodge(),
           color="black", width=.8) +
  geom_errorbar(aes(ymin=Mean-Se, ymax=Mean +Se),
                position=position_dodge(.8), width=.2) +
  labs(title = "Low Mistranslation", x = "", y = "") +
  theme_bw()+
  scale_fill_manual(name = "Type of \nMistranscription",values=c("grey35", "gray")) + 
  scale_y_continuous(expand=c(0,0))+
  ylim(0, 7)+
  theme(axis.text.x = element_markdown(size = 23,face = "bold",angle = 45,vjust = 1, hjust = 1),
        axis.title.y = element_text(size = 25,face = "bold"), 
        axis.title.x = element_text(size = 25,face = "bold"),
        axis.text.y = element_text(size = 23),
        plot.title = element_text(size=25,hjust=0.5,face = "bold"),
        legend.title = element_text(size = 18,face = "bold") ,
        legend.text = element_text(size = 18,face = "bold")) 
low_plot

ggarrange(high_plot,low_plot ,ncol = 2, nrow =1,widths = c(4.5, 5.5))


# ----------------------------------------------------------------------------------------- #
#                                             figure6. C                                    #
# ----------------------------------------------------------------------------------------- #
#high-low
high_low_data <- data.frame(species = c("Human", "Yeast", "Mouse", "Drosophila", "C.elegans"),
                            high=c(median(human_sm_high$sim_real),median(yeast_sm_high$sim_real),median(mouse_sm_high$sim_real),median(drosophila_sm_high$sim_real),median(elegans_sm_high$sim_real)),
                            low=c(median(human_sm_low$sim_real),median(yeast_sm_low$sim_real),median(mouse_sm_low$sim_real),median(drosophila_sm_low$sim_real),median(elegans_sm_low$sim_real)))
pdf("highlow_circle.pdf", width = 8, height = 8)
ggplot(high_low_data, aes(x = low, y = high)) +  
  geom_point(  
    aes(  
      shape = ifelse(high < low, 1, 21),    
      fill = ifelse(high < low, NA, "gray30"),  
      color = ifelse(high < low, "gray60", "gray30")  
    ),  
    size = 10, stroke = 2, show.legend = FALSE    
  ) +  
  geom_abline(intercept = 0, slope = 1, color = "red") +  
  scale_shape_identity() +  
  scale_color_identity() +   
  scale_fill_identity() +    
  theme_bw() +
  labs(
    x = "Low Mistranslation",
    y = "High Mistranslation",
    title = "Relative efficiency of \npurging deleterious mistranscriptions"
  ) +
  theme(
    plot.title = element_text(size = 22,face = "bold",margin = margin(b = 20), hjust = 0.5),
    axis.text.x = element_text(size = 23),
    axis.text.y = element_text(size = 23),
    axis.title.x = element_text(size = 25, face = "bold", margin = margin(t = 20)),
    axis.title.y = element_text(size = 25, face = "bold", margin = margin(r = 20)),
    legend.position = "none",  
    plot.margin = unit(c(1, 3, 1, 1), "lines")  
  ) +
  scale_x_continuous(
    limits = c(1.0, 1.5),  
    expand = c(0, 0)        
  ) +
  scale_y_continuous(
    limits = c(1.0, 1.5),  
    expand = c(0, 0)
  )
dev.off()


# ----------------------------------------------------------------------------------------- #
#                                             figure6. D                                    #
# ----------------------------------------------------------------------------------------- #
#MHtest
TR_MH <-  array(c(H_OR_TR$data[[1]],H_OR_TR$data[[2]],H_OR_TR$data[[4]],H_OR_TR$data[[5]], 
                  Y_OR_TR$data[[1]],Y_OR_TR$data[[2]],Y_OR_TR$data[[4]],Y_OR_TR$data[[5]], 
                  M_OR_TL$data[[1]],M_OR_TL$data[[2]],M_OR_TL$data[[4]],M_OR_TL$data[[5]],
                  D_OR_TR$data[[1]],D_OR_TR$data[[2]],D_OR_TR$data[[4]],D_OR_TR$data[[5]],
                  E_OR_TR$data[[1]],E_OR_TR$data[[2]],E_OR_TR$data[[4]],E_OR_TR$data[[5]]), 
                dim = c(2,2,5),   
                dimnames = list(TR = c("TR<meann","TR>mean"),
                                TE = c("TE<mean","TE>mean"),
                                Speices = c("Human","Yeast","Mouse","Drosophila","Elegans")))
OR_TR_result <- mantelhaen.test(TR_MH)

#plot
OR_bootstap <- list(OR_TR_result,t_test_result_OR_TR_two_yeast,t_test_result_OR_TR_two_mouse,t_test_result_OR_TR_two_drosophila,t_test_result_OR_TR_two_elegans)
OR_bootstap_data <- data.frame(Species=c("Human","Yeast","Mouse","C.elegans","Drosophila","Combined"),
                               OR=0,pvalue=0,lower=0,upper=0)
for(i in seq(2,length(OR_bootstap)))
{
  OR_bootstap_data[i,2] <- OR_bootstap[[i]]$estimate[[1]]
  OR_bootstap_data[i,3] <- OR_bootstap[[i]]$p.value[[1]]
  OR_bootstap_data[i,4] <- OR_bootstap[[i]]$conf.int[[1]]
  OR_bootstap_data[i,5] <- OR_bootstap[[i]]$conf.int[[2]]
}
OR_bootstap_data$pSig <- 0
for(i in seq(1,6))
{
  if(OR_bootstap_data[i,3] < 0.001){OR_bootstap_data[i,6] <- 0.001}
  if(OR_bootstap_data[i,3] > 0.001 && OR_bootstap_data[i,3] < 0.01){OR_bootstap_data[i,6] <- 0.01}
  if(OR_bootstap_data[i,3] > 0.01 &&OR_bootstap_data[i,3] < 0.05){OR_bootstap_data[i,6] <- 0.05}
}
OR_bootstap_data$Species <- factor(OR_bootstap_data$Species, 
                                   levels = c("Human", "Mouse", "Drosophila", "C.elegans", "Yeast","Combined"), 
                                   labels = c("*H. sapiens*", "*M. musculus*", "*D. melanogaster*", 
                                              "*C. elegans*", "*S. cerevisiae*","Combined"))
pdf("TR_OR.pdf")
ggplot(OR_bootstap_data, aes(x = Species, y = OR)) +
  geom_point(size = 8,colour = "red") +  
  geom_errorbar(aes(ymin = lower, ymax = upper), width = 0.2,linewidth = 1.5) +  
  geom_hline(yintercept = 1, linetype = "dashed", color = "red",linewidth = 2) +  
  coord_flip() + 
  labs(x = "Species", y = "Odds ratio", title = "") +
  theme_bw(base_size = 16) +
  theme_classic() +
  theme(
    axis.text.x = element_text(size = 28,face = "bold"), 
    axis.title.y = element_text(size = 35, face = "bold",margin = margin(r = 20)), 
    axis.title.x = element_text(size = 35,face = "bold"), 
    axis.text.y = element_markdown(size = 30,face = "bold"),
    plot.title = element_text(size = 30, face = "bold", hjust = 0.5)
  )

dev.off()
