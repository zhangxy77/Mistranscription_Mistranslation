library(stringr)
library(tidyverse)
library(ggplot2)
library(cowplot)
library(ggpubr)
library(reshape2)
load('figure6.related.Rdata')
# ----------------------------------------------------------------------------------------- #
#                                             figure6. A-B                                  #
# ----------------------------------------------------------------------------------------- #
#利用TL数据计算high_low
human_TL <- human_TL[order(human_TL$Error_Rate,decreasing = T),]
human_high_gene <- human_TL[human_TL$Error_Rate>median(human_TL$Error_Rate),1]
human_low_gene <- human_TL[human_TL$Error_Rate<median(human_TL$Error_Rate),1]
human_sm_high <- human_sm[rownames(human_sm)%in% human_high_gene,]
human_sm_low <- human_sm[rownames(human_sm)%in% human_low_gene,]
human_sm_high$sim_real <- human_sm_high$`sim_N/S` /human_sm_high$`real_N/S`
human_sm_low$sim_real <- human_sm_low$`sim_N/S`/human_sm_low$`real_N/S`
wilcox.test(human_sm_high$`real_N/S`,human_sm_high$`sim_N/S`)
wilcox.test(human_sm_low$`real_N/S`,human_sm_low$`sim_N/S`)
wilcox.test(human_sm_high$sim_real,human_sm_low$sim_real)

yeast_high_gene <- yeast_TL[yeast_TL$Error_Rate>median(yeast_TL$Error_Rate),1]
yeast_low_gene <- yeast_TL[yeast_TL$Error_Rate<median(yeast_TL$Error_Rate),1]
yeast_sm_high <- yeast_sm[rownames(yeast_sm)%in% yeast_high_gene,]
yeast_sm_low <- yeast_sm[rownames(yeast_sm)%in% yeast_low_gene,]
yeast_sm_high$sim_real <- yeast_sm_high$`sim_N/S`/yeast_sm_high$`real_N/S`
yeast_sm_low$sim_real <- yeast_sm_low$`sim_N/S` / yeast_sm_low$`real_N/S` 
wilcox.test(yeast_sm_high$`real_N/S`,yeast_sm_high$`sim_N/S`)
wilcox.test(yeast_sm_low$`real_N/S`,yeast_sm_low$`sim_N/S`)
wilcox.test(yeast_sm_high$sim_real,yeast_sm_low$sim_real)

mouse_high_gene <- mouse_TL[mouse_TL$Error.rate>median(mouse_TL$Error.rate),1]
mouse_low_gene <- mouse_TL[mouse_TL$Error.rate<median(mouse_TL$Error.rate),1]
mouse_sm_high <- mouse_sm[rownames(mouse_sm)%in% mouse_high_gene,]
mouse_sm_low <- mouse_sm[rownames(mouse_sm)%in% mouse_low_gene,]
mouse_sm_high$sim_real <- mouse_sm_high$`sim_N/S`/mouse_sm_high$`real_N/S`
mouse_sm_low$sim_real <- mouse_sm_low$`sim_N/S`/mouse_sm_low$`real_N/S`
mouse_sm_low <- mouse_sm_low[mouse_sm_low$sim_real<3,]
wilcox.test(mouse_sm_high$`real_N/S`,mouse_sm_high$`sim_N/S`)
wilcox.test(mouse_sm_low$`real_N/S`,mouse_sm_low$`sim_N/S`)
wilcox.test(mouse_sm_high$sim_real,mouse_sm_low$sim_real)


elegans_high_gene <- elegans_TL[elegans_TL$Error_Rate>median(elegans_TL$Error_Rate),1]
elegans_low_gene <- elegans_TL[elegans_TL$Error_Rate<median(elegans_TL$Error_Rate),1]
elegans_high_gene <- gsub("-", ".", elegans_high_gene)
elegans_low_gene <- gsub("-", ".", elegans_low_gene)
elegans_sm_high <- elegans_sm[rownames(elegans_sm)%in% elegans_high_gene,]
elegans_sm_low <- elegans_sm[rownames(elegans_sm)%in% elegans_low_gene,]
elegans_sm_high$sim_real <- elegans_sm_high$`sim_N/S` / elegans_sm_high$`real_N/S` 
elegans_sm_low$sim_real <- elegans_sm_low$`sim_N/S` / elegans_sm_low$`real_N/S`
wilcox.test(elegans_sm_high$`real_N/S`,elegans_sm_high$`sim_N/S`)
wilcox.test(elegans_sm_low$`real_N/S`,elegans_sm_low$`sim_N/S`)
wilcox.test(elegans_sm_high$sim_real,elegans_sm_low$sim_real)

drosophila_high_gene <- drosophila_TL[drosophila_TL$Error_Rate>median(drosophila_TL$Error_Rate),1]
drosophila_low_gene <- drosophila_TL[drosophila_TL$Error_Rate<median(drosophila_TL$Error_Rate),1]
drosophila_sm_high <- drosophila_sm[rownames(drosophila_sm)%in% drosophila_high_gene,]
drosophila_sm_low <- drosophila_sm[rownames(drosophila_sm)%in% drosophila_low_gene,]
drosophila_sm_high$sim_real <- drosophila_sm_high$`sim_N/S` / drosophila_sm_high$`real_N/S`
drosophila_sm_low$sim_real <- drosophila_sm_low$`sim_N/S` / drosophila_sm_low$`real_N/S`
wilcox.test(drosophila_sm_high$`real_N/S`,drosophila_sm_high$`sim_N/S`)
wilcox.test(drosophila_sm_low$`real_N/S`,drosophila_sm_low$`sim_N/S`)
wilcox.test(drosophila_sm_high$sim_real,drosophila_sm_low$sim_real)

#high
sm_gene_data_high <- data.frame(species = c("Human", "Yeast", "Mouse", "Drosophila", "C.elegans"),
                                Real = c(mean(human_sm_high$`real_N/S`),mean(yeast_sm_high$`real_N/S`),mean(mouse_sm_high$`real_N/S`),mean(drosophila_sm_high$`real_N/S`),mean(elegans_sm_high$`real_N/S`)),
                                Simulated = c(mean(human_sm_high$`sim_N/S`),mean(yeast_sm_high$`sim_N/S`),mean(mouse_sm_high$`sim_N/S`),mean(drosophila_sm_high$`sim_N/S`),mean(elegans_sm_high$`sim_N/S`)))
#sm_gene_plot_high <- melt(sm_gene_data_high, id="species", variable.name="Attribute", value.name = "NonSyn/Syn")
# 将数据转换为 data.table
sm_gene_data_high_dt <- as.data.table(sm_gene_data_high)
# 使用 data.table 的 melt 函数
sm_gene_plot_high <- melt(sm_gene_data_high_dt, id="species", variable.name="Attribute", value.name = "NonSyn/Syn")
mean_gene_high <- aggregate(sm_gene_plot_high$`NonSyn/Syn`, by=list(sm_gene_plot_high $species, sm_gene_plot_high $Attribute), FUN=mean)
sd_gene_high <- c(sd(elegans_sm_high$`real_N/S`),sd(drosophila_sm_high$`real_N/S`),sd(human_sm_high$`real_N/S`),sd(mouse_sm_high$`real_N/S`),sd(yeast_sm_high$`real_N/S`),sd(elegans_sm_high$`sim_N/S`),sd(drosophila_sm_high$`sim_N/S`),sd(human_sm_high$`sim_N/S`),sd(mouse_sm_high$`sim_N/S`),sd(yeast_sm_high$`sim_N/S`))
len_gene_high <- c(length(elegans_sm_high$`real_N/S`),length(drosophila_sm_high$`real_N/S`),length(human_sm_high$`real_N/S`),length(mouse_sm_high$`real_N/S`),length(yeast_sm_high$`real_N/S`),length(elegans_sm_high$`sim_N/S`),length(drosophila_sm_high$`sim_N/S`),length(human_sm_high$`sim_N/S`),length(mouse_sm_high$`sim_N/S`),length(yeast_sm_high$`sim_N/S`))
sm_gene_plot_res_high <- data.frame(mean_gene_high, sd=sd_gene_high, len=len_gene_high)
colnames(sm_gene_plot_res_high) = c("Species", "Attribute", "Mean", "Sd", "Count")
sm_gene_plot_res_high$Se <- sm_gene_plot_res_high$Sd/sqrt(sm_gene_plot_res_high$Count)
sm_gene_plot_res_high$Species <- factor(sm_gene_plot_res_high$Species, levels = c("Human", "Mouse", "Drosophila", "C.elegans", "Yeast"),  # 按所需顺序排列
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
        legend.position = "none")  # 移除图例
high_plot

#low
sm_gene_data_low <- data.frame(species = c("Human", "Yeast", "Mouse", "Drosophila", "C.elegans"),
                               Real = c(mean(human_sm_low$`real_N/S`),mean(yeast_sm_low$`real_N/S`),mean(mouse_sm_low$`real_N/S`),mean(drosophila_sm_low$`real_N/S`),mean(elegans_sm_low$`real_N/S`)),
                               Simulated = c(mean(human_sm_low$`sim_N/S`),mean(yeast_sm_low$`sim_N/S`),mean(mouse_sm_low$`sim_N/S`),mean(drosophila_sm_low$`sim_N/S`),mean(elegans_sm_low$`sim_N/S`)))
#sm_gene_plot_low <- melt(sm_gene_data_low, id="species", variable.name="Attribute", value.name = "NonSyn/Syn")
# 将数据转换为 data.table
sm_gene_data_low_dt <- as.data.table(sm_gene_data_low)
# 使用 data.table 的 melt 函数
sm_gene_plot_low <- melt(sm_gene_data_low_dt, id="species", variable.name="Attribute", value.name = "NonSyn/Syn")
mean_gene_low <- aggregate(sm_gene_plot_low$`NonSyn/Syn`, by=list(sm_gene_plot_low $species, sm_gene_plot_low $Attribute), FUN=mean)
sd_gene_low <- c(sd(elegans_sm_low$`real_N/S`),sd(drosophila_sm_low$`real_N/S`),sd(human_sm_low$`real_N/S`),sd(mouse_sm_low$`real_N/S`),sd(yeast_sm_low$`real_N/S`),sd(elegans_sm_low$`sim_N/S`),sd(drosophila_sm_low$`sim_N/S`),sd(human_sm_low$`sim_N/S`),sd(mouse_sm_low$`sim_N/S`),sd(yeast_sm_low$`sim_N/S`))
len_gene_low <- c(length(elegans_sm_low$`real_N/S`),length(drosophila_sm_low$`real_N/S`),length(human_sm_low$`real_N/S`),length(mouse_sm_low$`real_N/S`),length(yeast_sm_low$`real_N/S`),length(elegans_sm_low$`sim_N/S`),length(drosophila_sm_low$`sim_N/S`),length(human_sm_low$`sim_N/S`),length(mouse_sm_low$`sim_N/S`),length(yeast_sm_low$`sim_N/S`))
sm_gene_plot_res_low <- data.frame(mean_gene_low, sd=sd_gene_low, len=len_gene_low)
colnames(sm_gene_plot_res_low) = c("Species", "Attribute", "Mean", "Sd", "Count")
sm_gene_plot_res_low$Se <- sm_gene_plot_res_low$Sd/sqrt(sm_gene_plot_res_low$Count)
sm_gene_plot_res_low$Species <- factor(sm_gene_plot_res_low$Species, levels = c("Human", "Mouse", "Drosophila", "C.elegans", "Yeast"),  # 按所需顺序排列
                                       labels = c("*H. sapiens*", "*M. musculus*", "*D. melanogaster*", 
                                                  "*C. elegans*", "*S. cerevisiae*"))
low_plot <- ggplot(sm_gene_plot_res_low, aes(x=Species, y=Mean, fill=Attribute)) +
  geom_bar(stat="identity", position=position_dodge(),
           color="black", width=.8) +
  geom_errorbar(aes(ymin=Mean-Se, ymax=Mean +Se),
                position=position_dodge(.8), width=.2) +
  labs(title = "Low Mistranslation", x = "", y = "") +
  theme_bw()+
  scale_fill_manual(name = "Type of \nMistranscription",values=c("grey35", "gray")) + # 更改柱状图颜色为黑色和灰色
  scale_y_continuous(expand=c(0,0))+
  ylim(0, 7)+
  theme(axis.text.x = element_markdown(size = 23,face = "bold",angle = 45,vjust = 1, hjust = 1),
        axis.title.y = element_text(size = 25,face = "bold"), 
        axis.title.x = element_text(size = 25,face = "bold"),
        axis.text.y = element_text(size = 23),
        plot.title = element_text(size=25,hjust=0.5,face = "bold"),
        legend.title = element_text(size = 18,face = "bold") ,
        legend.text = element_text(size = 18,face = "bold")) # 调整图例文本字号
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
      shape = ifelse(high < low, 1, 21),  # 使用 21 作为实心圆（填充），1 作为空白圆（没有填充）  
      fill = ifelse(high < low, NA, "gray30"),  
      color = ifelse(high < low, "gray60", "gray30")  
    ),  
    size = 10, stroke = 2, show.legend = FALSE  # 为虚线框设定线条宽度  
  ) +  
  geom_abline(intercept = 0, slope = 1, color = "red") +  # 添加对角线
  scale_shape_identity() +  # 使用原始形状编码
  scale_color_identity() +   # 使用原始颜色编码
  scale_fill_identity() +    # 使用原始填充编码
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
    legend.position = "none",  # 隐藏图例
    plot.margin = unit(c(1, 3, 1, 1), "lines")  # 为x轴标签预留空间
  ) +
  scale_x_continuous(
    limits = c(1.0, 1.5),   # 根据数据实际值调整x轴限制
    expand = c(0, 0)        # 移除边距扩展
  ) +
  scale_y_continuous(
    limits = c(1.0, 1.5),   # 根据数据实际值调整y轴限制
    expand = c(0, 0)
  )
dev.off()
# 新增high_low_ratio列
high_low_data$high_low_ratio <- high_low_data$high / high_low_data$low

# 统计high_low_ratio > 1的个数
num_greater_than_1 <- sum(high_low_data$high_low_ratio > 1)
total_count <- nrow(high_low_data)

# 执行二项检验，假设比例为0.5
binom_test_result <- binom.test(num_greater_than_1, total_count, p = 0.5, alternative = "greater", conf.level = 0.9)

# 查看结果
print(binom_test_result)

# 执行单样本t检验
t_test_result <- t.test(high_low_data$high_low_ratio, mu = 1, alternative = "greater")

# 查看结果
print(t_test_result)

# Wilcoxon符号秩检验
wilcox_test_result <- wilcox.test(high_low_data$high_low_ratio, mu = 1, alternative = "greater")

# 查看结果
print(wilcox_test_result)


# ----------------------------------------------------------------------------------------- #
#                                             figure6. D-E                                  #
# ----------------------------------------------------------------------------------------- #
OR_1 <- H_OR_TR$data[[1]]+ Y_OR_TR$data[[1]]+ M_OR_TL$data[[1]]+ D_OR_TR$data[[1]]+ E_OR_TR$data[[1]]
OR_2 <- H_OR_TR$data[[2]]+ Y_OR_TR$data[[2]]+ M_OR_TL$data[[2]]+ D_OR_TR$data[[2]]+ E_OR_TR$data[[2]]
OR_4 <- H_OR_TR$data[[4]]+ Y_OR_TR$data[[4]]+ M_OR_TL$data[[4]]+ D_OR_TR$data[[4]]+ E_OR_TR$data[[4]]
OR_5 <- H_OR_TR$data[[5]]+ Y_OR_TR$data[[5]]+ M_OR_TL$data[[5]]+ D_OR_TR$data[[5]]+ E_OR_TR$data[[5]]
OR_TR_matrix <- matrix(c(OR_1,OR_2,OR_4,OR_5),nrow = 2)
OR_TR_result <- oddsratio(OR_TR_matrix)
print(OR_TR_result)

OR_1_TL <- H_OR_TL$data[[1]]+ Y_OR_TL$data[[1]]+ M_OR_TR$data[[1]]+ D_OR_TL$data[[1]]+ E_OR_TL$data[[1]]
OR_2_TL <- H_OR_TL$data[[2]]+ Y_OR_TL$data[[2]]+ M_OR_TR$data[[2]]+ D_OR_TL$data[[2]]+ E_OR_TL$data[[2]]
OR_4_TL <- H_OR_TL$data[[4]]+ Y_OR_TL$data[[4]]+ M_OR_TR$data[[4]]+ D_OR_TL$data[[4]]+ E_OR_TL$data[[4]]
OR_5_TL <- H_OR_TL$data[[5]]+ Y_OR_TL$data[[5]]+ M_OR_TR$data[[5]]+ D_OR_TL$data[[5]]+ E_OR_TL$data[[5]]
OR_TL_matrix <- matrix(c(OR_1_TL,OR_2_TL,OR_4_TL,OR_5_TL),nrow = 2)
OR_TL_result <- oddsratio(OR_TL_matrix)
print(OR_TL_result)


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
TR_MH_result <- mantelhaen.test(TR_MH)
TR_MH_result 
TL_MH <-  array(c(H_OR_TL$data[[1]],H_OR_TL$data[[2]],H_OR_TL$data[[4]],H_OR_TL$data[[5]], 
                  Y_OR_TL$data[[1]],Y_OR_TL$data[[2]],Y_OR_TL$data[[4]],Y_OR_TL$data[[5]], 
                  M_OR_TR$data[[1]],M_OR_TR$data[[2]],M_OR_TR$data[[4]],M_OR_TR$data[[5]],
                  D_OR_TL$data[[1]],D_OR_TL$data[[2]],D_OR_TL$data[[4]],D_OR_TL$data[[5]],
                  E_OR_TL$data[[1]],E_OR_TL$data[[2]],E_OR_TL$data[[4]],E_OR_TL$data[[5]]), 
                dim = c(2,2,5),   
                dimnames = list(TL = c("TL<meann","TL>mean"),
                                TE = c("TE<mean","TE>mean"),
                                Speices = c("Human","Yeast","Mouse","Drosophila","Elegans")))
TL_MH_result <-mantelhaen.test(TL_MH)
TL_MH_result

#单样本t检测
TR_orvalue <- data.frame(species=c("Human","Yeast","Mouse","Drosophila","Elegans"),
                         TR_or=c(H_OR_TR$measure[[2]],Y_OR_TR$measure[[2]],M_OR_TL$measure[[2]],D_OR_TR$measure[[2]],E_OR_TR$measure[[2]]))
num_lower_than_1 <- sum(TR_orvalue$TR_or < 1)
total_count <- nrow(TR_orvalue)
# 执行二项检验，假设比例为0.5
binom_test_result_TR <- binom.test(num_lower_than_1, total_count, p = 0.5, alternative = "greater", conf.level = 0.9)
print(binom_test_result_TR)
# 执行单样本t检验
t_test_result_TR <- t.test(TR_orvalue$TR_or, mu = 1, alternative = "less")
print(t_test_result_TR)

#单样本t检测
TL_orvalue <- data.frame(species=c("Human","Yeast","Mouse","Drosophila","Elegans"),
                         TL_or=c(H_OR_TL$measure[[2]],Y_OR_TL$measure[[2]],M_OR_TR$measure[[2]],D_OR_TL$measure[[2]],E_OR_TL$measure[[2]]))
num_greater_than_1 <- sum(TL_orvalue$TL_or > 1)
total_count <- nrow(TL_orvalue)
# 执行二项检验，假设比例为0.5
binom_test_result_TL <- binom.test(num_greater_than_1, total_count, p = 0.5, alternative = "greater", conf.level = 0.9)
print(binom_test_result_TL)

# 执行单样本t检验
t_test_result_TL <- t.test(TL_orvalue$TL_or, mu = 1, alternative = "greater")
print(t_test_result_TL)

#plot
TR_OR_plot <- data.frame(
  Species= c("Human","Yeast","Mouse","Drosophila","C.elegans","Combined"),
  OR = c(TR_orvalue$TR_or,TR_MH_result$estimate),  # 每层的效应估计
  lower = c(H_OR_TR$measure[[4]], Y_OR_TR$measure[[4]], M_OR_TL$measure[[4]], D_OR_TR$measure[[4]], E_OR_TR$measure[[4]],TR_MH_result$conf.int[[1]]),  # 下界
  upper = c(H_OR_TR$measure[[6]], Y_OR_TR$measure[[6]], M_OR_TL$measure[[6]], D_OR_TR$measure[[6]], E_OR_TR$measure[[6]],TR_MH_result$conf.int[[2]])   # 上界
)
TR_OR_plot$Species <- factor(TR_OR_plot$Species, 
                             levels = c("Human", "Mouse", "Drosophila", "C.elegans", "Yeast","Combined"),  # 按所需顺序排列
                             labels = c("*H. sapiens*", "*M. musculus*", "*D. melanogaster*", 
                                        "*C. elegans*", "*S. cerevisiae*","Combined"))
TL_OR_plot <- data.frame(
  Species= c("Human","Yeast","Mouse","Drosophila","C.elegans","Combined"),
  OR = c(TL_orvalue$TL_or,TL_MH_result$estimate),  # 每层的效应估计
  lower = c(H_OR_TL$measure[[4]], Y_OR_TL$measure[[4]], M_OR_TR$measure[[4]], D_OR_TL$measure[[4]], E_OR_TL$measure[[4]],TL_MH_result$conf.int[[1]]),  # 下界
  upper = c(H_OR_TL$measure[[6]], Y_OR_TL$measure[[6]], M_OR_TR$measure[[6]], D_OR_TL$measure[[6]], E_OR_TL$measure[[6]],TL_MH_result$conf.int[[2]])   # 上界
)
TL_OR_plot$Species <- factor(TL_OR_plot$Species, 
                             levels = c("Human", "Mouse", "Drosophila", "C.elegans", "Yeast","Combined"),  # 按所需顺序排列
                             labels = c("*H. sapiens*", "*M. musculus*", "*D. melanogaster*", 
                                        "*C. elegans*", "*S. cerevisiae*","Combined"))
pdf("TR_TL_OR.pdf")
ggplot(TR_OR_plot, aes(x = Species, y = OR)) +
  geom_point(size = 8,colour = "red") +  # 效应估计点
  geom_errorbar(aes(ymin = lower, ymax = upper), width = 0.2,linewidth = 1.5) +  # 误差线
  geom_hline(yintercept = 1, linetype = "dashed", color = "red",linewidth = 2) +  # 中性值线
  coord_flip() +  # 横轴竖轴翻转
  labs(x = "Species", y = "Odds Ratio", title = "") +
  theme_bw(base_size = 16) +
  theme_classic() +
  theme(
    axis.text.x = element_text(size = 28,face = "bold"), 
    axis.title.y = element_text(size = 35, face = "bold",margin = margin(r = 20)), 
    axis.title.x = element_text(size = 35,face = "bold"), 
    axis.text.y = element_markdown(size = 30,face = "bold"),
    plot.title = element_text(size = 30, face = "bold", hjust = 0.5)
  )


ggplot(TL_OR_plot, aes(x = Species, y = OR)) +
  geom_point(size = 8,colour = "red") +  # 效应估计点
  geom_errorbar(aes(ymin = lower, ymax = upper), width = 0.2,linewidth = 1.5) +  # 误差线
  geom_hline(yintercept = 1, linetype = "dashed", color = "red",linewidth = 2) +  # 中性值线
  coord_flip() +  # 横轴竖轴翻转
  labs(x = "Species", y = "Odds Ratio", title = "") +
  scale_y_continuous(limits = c(0, 6), breaks = seq(0, 6, by = 0.5)) +
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
