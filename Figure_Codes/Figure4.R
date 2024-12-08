library(stringr)
library(tidyverse)
library(ggplot2)
library(cowplot)
library(ggpubr)
library(reshape2)
load('figure4.related.Rdata')
# ----------------------------------------------------------------------------------------- #
#                                             figure4. B                                    #
# ----------------------------------------------------------------------------------------- #
C0T7 <- data.frame(gene=fn_select$geno,COT7_fn=(fn_select$TUY_D7_rg),tag=0,gene1=0,gene2=0,fitness1=0,fitness2=0,fitness12=0)
for(i in seq(1,nrow(C0T7)))
{
  a <- unlist(strsplit(C0T7[i,1],split = ' '))
  if(length(a)==1){C0T7[i,3]=1}
  if(length(a)>1){C0T7[i,3]=2}
}
# 假设single和double mutation的fitness分别存储在以下向量中
single_fitness <- C0T7[C0T7[,3]==1,2]
double_fitness <- C0T7[C0T7[,3]==2,2]

for(i in seq(1,nrow(C0T7)))
{
  if(C0T7[i,3]==1)
  {C0T7[i,4] <- C0T7[i,1]}
  if(C0T7[i,3]==2)
  {
    a <- unlist(strsplit(C0T7[i,1],split = " "))
    C0T7[i,4] <- a[[1]]
    C0T7[i,5] <- a[[2]]
  }
}

for(i in seq(1,nrow(C0T7)))
{
  if(C0T7[i,3]==1)
  {C0T7[i,6] <- C0T7[i,2]}
  if(C0T7[i,3]==2)
  {
    C0T7[i,8] <- C0T7[i,2]
    a <- which(C0T7[,1]==C0T7[i,4])
    if(length(a)>0)
    { C0T7[i,6] <- C0T7[a,2]}
    b <- which(C0T7[,1]==C0T7[i,5])
    if(length(b)>0)
    {C0T7[i,7] <- C0T7[b,2]}
  }
}
fitness_double <- C0T7[C0T7[,3]==2,]
fitness_double <- fitness_double[fitness_double[,6]!=0&fitness_double[,7]!=0&fitness_double[,8]!=0,]
fitness_double$s1 <- (fitness_double$fitness1-1)
fitness_double$s2 <- (fitness_double$fitness2-1)
fitness_double$s12 <- (fitness_double$fitness12-1)
fitness_double$e <- fitness_double$s12-(fitness_double$s1+fitness_double$s2)
tr<- 0.3;tl<- 0.8
fitness_double$tre <- fitness_double$s1*tr*(1-tl)
fitness_double$tle <- fitness_double$s2*tl*(1-tr)
fitness_double$tr_tl <- fitness_double$tre+fitness_double$tle+fitness_double$s12*tr*tl

fitness_double$u1 <- (fitness_double$tre)+(fitness_double$tle)
fitness_double$u2 <- (fitness_double$tr_tl)
fitness_2 <- fitness_double[fitness_double$u1<0&fitness_double$u2<0,]

comparison_results <- fitness_2 %>%
  group_by(gene1) %>%
  summarise(U1_less_than_U2 = all(u1 > u2))

# 统计 U1 全部小于 U2 的基因组合数量
num_true <- sum(comparison_results$U1_less_than_U2)
num_total <- nrow(comparison_results)

cat("u1 全部小于 u2 的基因组合数量：", num_true, " / ", num_total, "\n")
# 计算每个基因组合的 U1 和 U2 的数量
sample_sizes <- fitness_2 %>%
  group_by(gene1) %>%
  summarise(count_U1 = n(), count_U2 = n())

# 过滤样本数不足的基因组合
valid_genes <- sample_sizes %>%
  filter(count_U1 > 1 & count_U2 > 1) %>%
  pull(gene1)

# 仅保留有效基因组合的数据
filtered_data <- fitness_2 %>%
  filter(gene1 %in% valid_genes)

# 对每个有效的基因组合进行配对 t 检验
test_results <- filtered_data %>%
  group_by(gene1) %>%
  summarise(p_value = t.test(u1, u2, paired = TRUE, alternative = "less")$p.value)
# 设置显著性水平
alpha <- 0.05

Siggene <- test_results[order(test_results$p_value,decreasing = F),] 
Siggenename <-  unlist(c(Siggene[1:100,1]))
# 统计显著小于的基因组合数量
num_significant <- sum(test_results$p_value < alpha)
num_total <- nrow(test_results)

cat("u1 显著小于 u2 的基因组合数量：", num_significant, " / ", num_total, "\n")
# 可视化显著性检验结果

fitness_plot <- fitness_2[fitness_2[,4]%in% Siggenename,]
colnames(fitness_plot)[[12]] <- "Single"
colnames(fitness_plot)[[13]] <- "Double"
# 将数据转换为长格式，以便绘制箱线图
long_data <- fitness_plot %>%
  pivot_longer(cols = c(Single, Double), names_to = "Variable", values_to = "Value")

count_fitness <- fitness_2 %>%
  group_by(gene1) %>%
  summarise(count_gene = n(),mean_single=mean(u1),mean_Double=mean(u2))
count_fitness <- merge(count_fitness,test_results,by = "gene1")
count_fitness_less <- count_fitness[count_fitness$gene1 %in% Siggenename,]

ggplot(count_fitness_less, aes(x = mean_single, y = mean_Double)) +
  # 绘制带红色边框的三角形
  geom_point(aes(
    shape = ifelse(mean_Double < mean_single, 25,17),
    size = count_gene + 1.5,
    color = ifelse(p_value < 0.05, "red", NA)
  ), stroke = 1.5, show.legend = FALSE) +  # 隐藏图例
  
  # 绘制填充颜色的三角形
  geom_point(aes(
    shape = ifelse(mean_Double < mean_single, 25,17),
    fill = ifelse(mean_Double < mean_single,  "gray30", "gray60"),
    size = count_gene,
    color = ifelse(mean_Double < mean_single,  "gray30", "gray60")
  ), stroke = 1, show.legend = FALSE) +  # 隐藏图例
  
  geom_abline(intercept = 0, slope = 1, color = "red") +
  scale_shape_identity() +
  scale_color_identity() +
  scale_fill_identity() +
  scale_size_continuous(name = "genes", range = c(2, 8)) +
  theme_bw() +
  labs(x = "Expectation of epistasis in Single Mutation", y = "Observation of epistasis in Double Mutation", title = "") +
  theme(
    plot.title = element_text(size = 22, hjust = 0.5),
    axis.text.x = element_text(size = 22),
    axis.text.y = element_text(size = 22),
    axis.title.x = element_text(size = 20),
    axis.title.y = element_text(size = 20, margin = margin(r = 20)),  # 向左移动y轴标题
    legend.position = "none",  # 隐藏图例
    plot.margin = unit(c(1, 3, 1, 1), "lines")  # 为x轴标签预留空间
  )  +
  scale_x_continuous(limits = c(-0.125,0))+
  scale_y_continuous(limits = c(-0.125,0))

special_gene <- fitness_double[fitness_double$gene1=="C284T",]
special_gene$label <- paste(special_gene$gene1,special_gene$gene2,sep = "+")
pdf("C284T.pdf", width = 10, height = 8)
ggplot(special_gene, aes(x = u1, y = u2)) +
  # 绘制填充颜色的三角形
  geom_point(aes(
    shape = ifelse(u2 < u1, 25,17),
    fill = ifelse(u2 < u1,  "gray30", "gray60"),
    color = ifelse(u2 < u1, "gray30", "gray60")
  ), size = 10, stroke = 1, show.legend = FALSE) +  # 直接在外部设置三角形大小
  # 添加基因名称标签，位置在三角形的正上方
  #geom_text(aes(label = label), hjust = -0.2, vjust = -0.8, size = 5) +  # 调整标签位置和大小
  geom_abline(intercept = 0, slope = 1, color = "red") +
  scale_shape_identity() +
  scale_color_identity() +
  scale_fill_identity() +
  theme_bw() +
  labs(x = "Fitness effect expected by additive model", y = "Observed fitness effect of\n Double mutants", title = "") +
  theme(
    plot.title = element_text(size = 22, hjust = 0.5),
    axis.text.x = element_text(size = 23),
    axis.text.y = element_text(size = 23),
    axis.title.x = element_text(size = 25,face = "bold",margin = margin(t = 20)),
    axis.title.y = element_text(size = 25, face = "bold",margin = margin(r = 20)),  # 向左移动y轴标题
    legend.position = "none",  # 隐藏图例
    plot.margin = unit(c(1, 3, 1, 1), "lines")  # 为x轴标签预留空间
  )  +
  scale_x_continuous(limits = c(-0.12,-0.01))+
  scale_y_continuous(limits = c(-0.12,-0.01))
dev.off()

# ----------------------------------------------------------------------------------------- #
#                                             figure4. C                                    #
# ----------------------------------------------------------------------------------------- #
genename <- c("AGS","AGY","AUS","AUU","AUY","TGS","TGY","TUS","TUU","TUY")
gene_fitness <- list()
for(i in seq(1,length(fitness_list)))
{
  tt <- load(fitness_list[[i]])
  gene_fitness[[i]] <- count_fitness[ sample(1:length(count_fitness[,1]),100),]
}
names(gene_fitness) <- genename
fitness <- matrix(0,nrow = 100,ncol = 10)
colnames(fitness) <- genename
fitness <- as.data.frame(fitness)
for(i in seq(1,10))
{
  tt <- gene_fitness[[i]]
  fitness[,i] <- tt[,4]-tt[,3]
}
fitness_long <- melt(fitness, variable.name = "Gene", value.name = "Fitness")
ggplot(fitness_long, aes(x = factor(Gene, levels = c("TGS","AGS","TUS","AUS","TGY","AGY","TUY","AUY","TUU","AUU")), y = Fitness)) +
  geom_boxplot() +   # 使用箱线图展示分布
  geom_jitter(width = 0.2, alpha = 0.5, color = "gray30") +  # 添加抖动点，展示每个fitness值
  geom_hline(yintercept = 0, linetype = "dashed", color = "red",linewidth = 0.8) +  # 添加红色虚线 
  theme_minimal() +
  labs(title = "", x = "", y = "") +
  theme(
    axis.title.x = element_text(size = 16),               # 调整x轴标题字体大小 
    axis.text.x = element_text(size = 14),                # 调整x轴标签字体大小
    axis.text.y = element_text(size = 14)                 # 调整y轴标签字体大小
  )


