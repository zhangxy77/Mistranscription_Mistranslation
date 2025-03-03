library(stringr)
library(tidyverse)
library(ggplot2)
library(cowplot)
library(ggpubr)
library(reshape2)

# ----------------------------------------------------------------------------------------- #
#                                             figure1. B                                    #
# ----------------------------------------------------------------------------------------- #
create_errortrend <- function(changetype_data, changecodon_data) {
  errortrend <- as.data.frame(matrix(0, nrow = length(changetype_data[, 1]), ncol = ncol(changecodon_data) - 1))
  colnames(errortrend) <- colnames(changecodon_data)[-1]
  return(data.frame(changetype = changetype_data, errortrend))
}

# Apply the function to each species
human_errortrend <- create_errortrend(changetype, Human_changecodon)
yeast_errortrend <- create_errortrend(changetype, yeast_changecodon)
mouse_errortrend <- create_errortrend(changetype, Mouse_changecodon)
drosophila_errortrend <- create_errortrend(changetype, drosophila_changecodon)
elegans_errortrend <- create_errortrend(changetype, elegans_changecodon)

update_errortrend <- function(tt, errortrend_data) {
  for(k in seq(1, nrow(tt))) {
    a <- unlist(strsplit(tt[k, 2], split = "-"))
    b <- unlist(strsplit(a[[1]], split = ""))
    c <- unlist(strsplit(a[[2]], split = ""))
    for(j in seq(1, 3)) {
      if(b[[j]] != c[[j]]) {
        tt[k, 1] <- paste(b[[j]], c[[j]], sep = "-")
      }
    }
  }
  
  for(k in seq(1, nrow(tt))) {
    for(j in seq(3, ncol(tt))) {
      if(tt[k, j] != 0) {
        loc_1 <- which(errortrend_data[, 2] == tt[k, 1])
        errortrend_data[loc_1, j] <- errortrend_data[loc_1, j] + tt[k, j]
      }
    }
  }
}

# Loop over species and apply the update function
for(i in seq(1, 5)) {
  tt <- species[[i]]
  tt <- cbind(0, tt)
  
  if(i == 1) {
    update_errortrend(tt, human_errortrend)
  } else if(i == 2) {
    update_errortrend(tt, yeast_errortrend)
  } else if(i == 3) {
    update_errortrend(tt, mouse_errortrend)
  } else if(i == 4) {
    update_errortrend(tt, drosophila_errortrend)
  } else if(i == 5) {
    update_errortrend(tt, elegans_errortrend)
  }
}

colnames(human_errortrend) <- gsub("\\.", "-", colnames(human_errortrend))
colnames(yeast_errortrend) <- gsub("\\.", "-", colnames(yeast_errortrend))
colnames(mouse_errortrend) <- gsub("\\.", "-", colnames(mouse_errortrend))
colnames(drosophila_errortrend) <- gsub("\\.", "-", colnames(drosophila_errortrend))
colnames(elegans_errortrend) <- gsub("\\.", "-", colnames(elegans_errortrend))
TR_H[,1] <- gsub("\\.", "-", TR_H[,1])
TR_Y[,1] <- gsub("\\.", "-", TR_Y[,1])
TR_M[,1] <- gsub("\\.", "-", TR_M[,1])
TR_D[,1] <- gsub("\\.", "-", TR_D[,1])
TR_E[,1] <- gsub("\\.", "-", TR_E[,1])
species_TR <- list(TR_H,TR_Y,TR_M,TR_D,TR_E)
species_errortrand <- list(human_errortrend,yeast_errortrend,mouse_errortrend,drosophila_errortrend,elegans_errortrend)
species_errortrand_bp <- list(human_errortrend,yeast_errortrend,mouse_errortrend,drosophila_errortrend,elegans_errortrend)
for(i in seq(1,5))
{
  tt <- species_errortrand[[i]]
  dd <- species_TR[[i]]
  for(j in seq(3,ncol(tt)))
  {
    loc_1 <- which(dd[,1]==colnames(tt)[[j]])
    b <- dd[loc_1,2]*10
    tt[,j] <- tt[,j]/b
  }
  species_errortrand_bp[[i]] <- tt
  
}
TR_H_mean <- species_errortrand_bp[[1]]
TR_Y_mean <- species_errortrand_bp[[2]]
TR_M_mean <- species_errortrand_bp[[3]]
TR_D_mean <- species_errortrand_bp[[4]]
TR_E_mean <- species_errortrand_bp[[5]]

calculate_stats <- function(df) {
  df$meanvalue <- apply(df[3:ncol(df)], 1, median, na.rm = TRUE)
  df$sd <- apply(df[3:ncol(df)], 1, sd, na.rm = TRUE)
  df$se <- apply(df[3:ncol(df)], 1, function(x) sd(x, na.rm = TRUE) / sqrt(sum(!is.na(x))))
  return(df)
}

# Calculate stats for all species
TR_H_mean <- calculate_stats(species_errortrand_bp[[1]])
TR_Y_mean <- calculate_stats(species_errortrand_bp[[2]])
TR_M_mean <- calculate_stats(species_errortrand_bp[[3]])
TR_D_mean <- calculate_stats(species_errortrand_bp[[4]])
TR_E_mean <- calculate_stats(species_errortrand_bp[[5]])

species_mean <- data.frame(type=TR_H_mean[,1],
                           Human=TR_H_mean$meanvalue,Yeast=TR_Y_mean$meanvalue,Mouse=TR_M_mean$meanvalue,
                           Drosophila=TR_D_mean$meanvalue,Elegans=TR_E_mean$meanvalue)
species_sd <- data.frame(type=TR_H_mean[,1],
                         Human_sd=TR_H_mean$sd,Yeast_sd=TR_Y_mean$sd,Mouse_sd=TR_M_mean$sd,
                         Drosophila_sd=TR_D_mean$sd,Elegans_sd=TR_E_mean$sd)
species_se <- data.frame(type=TR_H_mean[,1],
                         Human_se=TR_H_mean$se,Yeast_se=TR_Y_mean$se,Mouse_se=TR_M_mean$se,
                         Drosophila_se=TR_D_mean$se,Elegans_se=TR_E_mean$se)
species_mean_long <- melt(species_mean, id.vars = "type",
                          variable.name = "Species", value.name = "Value")
species_mean_se_long <- melt(species_se, id.vars = "type",
                             measure.vars = c("Human_se", "Yeast_se", "Mouse_se", "Drosophila_se", "Elegans_se"),
                             variable.name = "Species", value.name = "se")
species_mean_long$se <- species_mean_se_long$se
species_mean_long[,2] <- gsub("Elegans", "C.elegans", species_mean_long[,2])

# Plot
species_mean_long <- species_mean_long %>%
  mutate(Species = recode(Species,
                          "Human" = "H. sapiens",
                          "Mouse" = "M. musculus",
                          "Drosophila" = "D. melanogaster",
                          "C.elegans" = "C. elegans",
                          "Yeast" = "S. cerevisiae"))

species_mean_long$Species <- factor(species_mean_long$Species, 
                                    levels = c("H. sapiens", "M. musculus", "D. melanogaster", "C. elegans", "S. cerevisiae"))

ggplot(species_mean_long, aes(x = type, y = Value, fill = Species)) +  
  geom_bar(stat = "identity", position = "dodge") +  
  geom_errorbar(aes(ymin = Value - se, ymax = Value + se),   
                position = position_dodge(0.9), width = 0.25) +  
  theme_minimal() +  
  scale_y_continuous(labels = c(expression(italic(0)),   
                                expression(2%*%10^-5),   
                                expression(4%*%10^-5),   
                                expression(6%*%10^-5)),  
                     expand = c(0,0),   
                     breaks = c(0,0.00002,0.00004,0.00006),   
                     limits = c(0,0.00007)) +  
  labs(x = "Type of mistranscription", y = "Error rate / bp", title = "", fill = "Species") +  
  theme(  
    plot.title = element_text(size = 25, hjust = 0.5),  
    axis.text.x = element_text(size = 22),  
    axis.text.y = element_text(size = 22),  
    axis.title.x = element_text(size = 26, face = "bold"), 
    axis.title.y = element_text(size = 26, face = "bold", margin = margin(l = 20, r = 20)),
    legend.text = element_text(size = 20, face = "italic"),  
    legend.title = element_text(size = 20),  
    legend.position = c(0.05, 0.95),  
    legend.justification = c(0, 1)  
  ) +  
  scale_fill_manual(values = c("H. sapiens" = "gray20",  
                               "M. musculus" = "gray40",  
                               "D. melanogaster" = "gray60",  
                               "C. elegans" = "gray80",  
                               "S. cerevisiae" = "gray90"))

# ----------------------------------------------------------------------------------------- #
#                                             figure1. C                                    #
# ----------------------------------------------------------------------------------------- #
human_p1 <- ggplot(human_gene_transcription_error,aes(y = reorder(gene, Error.rate),x=Error.rate))+ 
  geom_point(size = 1)+  # Increase the size of the points 
  scale_x_continuous(labels = c(expression(italic(0)), 
                                expression(0.5), 
                                expression(1), 
                                expression(1.5), 
                                expression(2), 
                                expression(2.5)), 
                     expand = c(0,0), 
                     breaks = c(0,0.000005,0.00001,0.000015,0.00002,0.000025), 
                     limits = c(0,0.000028)) +
  labs(x = "", y = "", title = "")+
  theme(axis.title.x = element_text(size=30, margin = margin(r = 20)), 
        axis.text.x = element_text(size = 19), 
        axis.ticks.y = element_blank(),
        axis.title.y = element_text(size=30, margin = margin(r = 20)), # Increase the distance between y-axis title and labels
        axis.text.y = element_blank(),
        plot.title = element_blank())+
  coord_cartesian(clip = "off")
human_p1
ggMarginal(human_p1, type="histogram", fill = "grey48", xparams = list(bins=40), yparams = list(bins=40))

# ----------------------------------------------------------------------------------------- #
#                                             figure1.  E                                   #
# ----------------------------------------------------------------------------------------- #
human_p3 <- ggplot(human_gene_translation_error,aes(y = reorder(gene, Error_Rate),x=Error_Rate))+ 
  geom_point(size = 1)+  # Increase the size of the points
  scale_x_continuous(labels = c(expression(italic(0)), 
                                expression(2), 
                                expression(4), 
                                expression(6), 
                                expression(8), 
                                expression(10)), 
                     expand = c(0,0), 
                     breaks = c(0,0.002,0.004,0.006,0.008,0.01), 
                     limits = c(0,0.01062)) +
  labs(x = "", y = "", title = "")+
  theme(axis.title.x = element_text(size=20, margin = margin(r = 20)), 
        axis.text.x = element_text(size = 19), 
        axis.ticks.y = element_blank(),
        axis.title.y = element_text(size=30, margin = margin(r = 20)), # Increase the distance between y-axis title and labels
        axis.text.y = element_blank(),
        plot.title = element_blank())+
  coord_cartesian(clip = "off")
human_p3
ggMarginal(human_p3, type="histogram", fill = "grey48", xparams = list(bins=40), yparams = list(bins=40))