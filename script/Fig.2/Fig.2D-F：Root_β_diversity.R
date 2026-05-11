#library package
library(dplyr)
library(ggplot2)
library(RColorBrewer)
library(Rmisc)
library(vegan)
  

# Set Work Path
pwd <- dirname(rstudioapi::getSourceEditorContext()$path)
setwd(pwd)


#import data
metadata <- read.table("../../data/Metadata.txt",header = T, sep = "\t", fill = TRUE, na.strings = "")
row.names(metadata) = metadata$SampleID
bc_dist <- read.table("../../data/bray_curtis_distance_matrix.tsv",  header = T, row.names = 1, check.names=FALSE)


###Root Site ####
metadata_R <- subset(metadata, Area_Microbiome == "Root")
bc_dist_R <- bc_dist[row.names(bc_dist) %in% row.names(metadata_R), names(bc_dist) %in% row.names(metadata_R)]
bc_R_pcoa <- cmdscale(bc_dist_R, k=(10), eig=T)
bc_R_eig <- bc_R_pcoa$eig
bc_R_scores <- as.data.frame(bc_R_pcoa$points)
bc_R_pcoa_metadata <- cbind(metadata_R[,], bc_R_scores[,1:5])
names(bc_R_pcoa_metadata)
names(bc_R_pcoa_metadata)[which(names(bc_R_pcoa_metadata) %in% c("V1", "V2", "V3", "V4", "V5"))] <- c("PCoA1", "PCoA2", "PCoA3", "PCoA4", "PCoA5")


###Permutation multivariate non-parametric test(PERMANOVA)
set.seed(123)
bc_R_adonis_Site <- adonis2(bc_dist_R~Site, data=bc_R_pcoa_metadata, permutations = 999)
bc_R_adonis_Site
R_adonis_Site <- paste0(" Site adonis R2= ", round(bc_R_adonis_Site$R2, 3), "   P = ", bc_R_adonis_Site$`Pr(>F)`)


#visualization
bc_R_pcoa_metadata$Site <- factor(bc_R_pcoa_metadata$Site,levels = c("Village","Factory"))


#set colors
col_Site <- c("#8DD3C7","#FFFFB3")

bc_Root_pcoa_g12_ <- ggplot(bc_R_pcoa_metadata, aes(PCoA1, PCoA2, colour=Site, fill=Site, shape = Environment)) +
  geom_point(alpha=0.7,color = 'black',shape = 21, size=7) +  
  scale_shape_manual(values=c(21,24)) +  
  scale_fill_manual(values=col_Site) +
  scale_color_manual(values=col_Site) +
  labs(x=paste("PCoA 1 (", format(100 * bc_R_eig[1] / sum(bc_R_eig), digits=4), "%)", sep=""),
       y=paste("PCoA 2 (", format(100 * bc_R_eig[2] / sum(bc_R_eig), digits=4), "%)", sep=""),
       title = R_adonis_Site) +  
  theme_classic() +
  theme(
    plot.title = element_text(hjust = 0,vjust = 5, size = 10,face = "bold"), 
    legend.title = element_text(size=10,color="black", face = "bold"),
    legend.text= element_text(size=8,color="black"),
    axis.title = element_text(size = 10),  
    axis.text = element_text(size = 8),   
    plot.margin = margin(40, 20, 20, 20),
    legend.key.size = unit(6, "mm")   
  ) +
  guides(fill = guide_legend(
    override.aes = list(size = 5)  
  ))

bc_Root_pcoa_g12_
ggsave("../../figures/Fig.2/bc_dis_Root_Site.pdf", bc_Root_pcoa_g12_, width = 8, height =5.5)




############             Root Site Plants            ################
#Root Village Plants
metadata_R_V <- subset(metadata_R, Site == "Village")
table(metadata_R_V$Plant_name)

col_Plant <- c(
  "Ageratum conyzoides"        = "#FFFF99",
  "Bidens alba"                = "#999900",
  "Bidens pilosa"              = "#FFCC66",
  "Bidens pilosa var. radiata" = "#FF99FF",
  "Colocasia esculenta"        = "#9999FF",
  "Dimocarpus longan"          = "#FF6600",
  "Dryopteris crassirhizoma"   = "#FF00CC",
  "Eupatorium catarium"        = "#99FFFF",
  "Flacourtia indica"          = "#99FF66",
  "Jatropha curcas"            = "#99CCFF",
  "Murraya paniculata"         = "#3333FF",
  "Ophiopogon japonicus"       = "#33CC00",
  "Pteris multifida"           = "#330000",
  "Solanum nigrum"             = "#336600"
)

#PcoA
bc_dist_R_V <- bc_dist[row.names(bc_dist) %in% row.names(metadata_R_V), names(bc_dist) %in% row.names(metadata_R_V)]
bc_R_V_pcoa <- cmdscale(bc_dist_R_V, k=(10), eig=T)
bc_R_V_eig <- bc_R_V_pcoa$eig
bc_R_V_scores <- as.data.frame(bc_R_V_pcoa$points)
bc_R_V_pcoa_metadata <- cbind(metadata_R_V[,], bc_R_V_scores[,1:5])
names(bc_R_V_pcoa_metadata)
names(bc_R_V_pcoa_metadata)[which(names(bc_R_V_pcoa_metadata) %in% c("V1", "V2", "V3", "V4", "V5"))] <- c("PCoA1", "PCoA2", "PCoA3", "PCoA4", "PCoA5")


###Permutation multivariate non-parametric test(PERMANOVA)
set.seed(110)
bc_R_V_adonis_Plant <- adonis2(bc_dist_R_V~ScientificNames, data=metadata_R_V, permutations = 999)
bc_R_V_adonis_Plant
R_V_P_adonis <- paste0(" Plant adonis R2= ", round(bc_R_V_adonis_Plant$R2, 3), "   P = ", bc_R_V_adonis_Plant$`Pr(>F)`)

set.seed(123)
bc_R_V_adonis_Environment <- adonis2(bc_dist_R_V~Environment, data=metadata_R_V, permutations = 999)
bc_R_V_adonis_Environment
R_V_adonis_Environment <- paste0("Environment adonis R2= ", round(bc_R_V_adonis_Environment$R2, 3), "   P = ", bc_R_V_adonis_Environment$`Pr(>F)`, " ns")


#Visualization | ggplot 形状与颜色存在冲突，分两次组合图片
bc_R_V_pcoa_metadata$ScientificNames <- as.factor(bc_R_V_pcoa_metadata$ScientificNames)
Plant_color <- col_Plant[unique(bc_R_V_pcoa_metadata$ScientificNames)]


bc_R_V_pcoa_g12 <- ggplot(bc_R_V_pcoa_metadata, aes(PCoA1, PCoA2, fill=ScientificNames,  shape = Environment)) + 
  geom_point(alpha=0.7, color = "black", size=7) +  
  scale_fill_manual(values=Plant_color, name = 'Plant') +
  scale_shape_manual(values=c(21,24)) +  
  labs(x=paste("PCoA 1 (", format(100 * bc_R_V_eig[1] / sum(bc_R_V_eig), digits=4), "%)", sep=""),
       y=paste("PCoA 2 (", format(100 * bc_R_V_eig[2] / sum(bc_R_V_eig), digits=4), "%)", sep=""),
       title = paste(R_V_P_adonis, "\n",R_V_adonis_Environment)) +  
  theme_classic() +
  theme(
    plot.title = element_text(hjust = 0,vjust = 5, size = 10,face = "bold"), 
    legend.title = element_text(size=10,color="black", face = "bold"),
    legend.text= element_text(size=8,color="black"),
    axis.title = element_text(size = 10),  
    axis.text = element_text(size = 8),   
    plot.margin = margin(40, 20, 20, 20),
    legend.key.size = unit(6, "mm")
  ) +guides(
    shape = guide_legend(order = 2),  # 设置shape图例优先显示
    fill = guide_legend(order = 1)    # 设置fill图例其次显示
  )+
  guides(fill = guide_legend(order = 1,override.aes = list(size = 5)),
         shape =guide_legend(order = 2,override.aes = list(size = 5))
         )

bc_R_V_pcoa_g12 
ggsave("../../figures/Fig.2/bc_Root_Village_Plant_shape.pdf", bc_R_V_pcoa_g12, width = 9.5, height =6)


#again  This result has the correct color for shape
bc_R_V_pcoa_g12_ <- ggplot(bc_R_V_pcoa_metadata, aes(PCoA1, PCoA2, fill=ScientificNames, shape = Environment)) + 
  geom_point(alpha=0.7, color = "black",shape = 21, size=7) +  
  scale_fill_manual(values=Plant_color, name = 'Plant') +
  scale_shape_manual(values=c(21,24)) +  
  labs(x=paste("PCoA 1 (", format(100 * bc_R_V_eig[1] / sum(bc_R_V_eig), digits=4), "%)", sep=""),
       y=paste("PCoA 2 (", format(100 * bc_R_V_eig[2] / sum(bc_R_V_eig), digits=4), "%)", sep=""),
       title = paste(R_V_P_adonis, "\n",R_V_adonis_Environment)) +  
  theme_classic() +
  theme(
    plot.title = element_text(hjust = 0,vjust = 5, size = 10,face = "bold"), 
    legend.title = element_text(size=10,color="black", face = "bold"),
    legend.text= element_text(size=8,color="black"),
    axis.title = element_text(size = 10),  
    axis.text = element_text(size = 8),   
    plot.margin = margin(40, 20, 20, 20),
    legend.key.size = unit(6, "mm")   
  )+
  guides(fill = guide_legend(order = 1,override.aes = list(size = 5)))

bc_R_V_pcoa_g12_ 
ggsave("../../figures/Fig.2/bc_Root_Village_Plant_color.pdf", bc_R_V_pcoa_g12_, width = 9.5, height =6)



#Root Factory  Plants
metadata_R_H <- subset(metadata_R, Site == "Factory")
table(metadata_R_H$ScientificNames)
col_plant <- c("#FFFF99", "#FFCCFF", "#FFCC66", "#FF99FF", "#99CC33", "#FF6600", "#FF00CC", "#99FFFF", "#99FF66", "#99CCFF") 

#PcoA
bc_dist_R_H <- bc_dist[row.names(bc_dist) %in% row.names(metadata_R_H), names(bc_dist) %in% row.names(metadata_R_H)]
bc_R_H_pcoa <- cmdscale(bc_dist_R_H, k=(10), eig=T)
bc_R_H_eig <- bc_R_H_pcoa$eig
bc_R_H_scores <- as.data.frame(bc_R_H_pcoa$points)
bc_R_H_pcoa_metadata <- cbind(metadata_R_H[,], bc_R_H_scores[,1:5])
names(bc_R_H_pcoa_metadata)
names(bc_R_H_pcoa_metadata)[which(names(bc_R_H_pcoa_metadata) %in% c("V1", "V2", "V3", "V4", "V5"))] <- c("PCoA1", "PCoA2", "PCoA3", "PCoA4", "PCoA5")


###Permutation multivariate non-parametric test(PERMANOVA)
set.seed(123)
bc_R_H_adonis_Plant <- adonis2(bc_dist_R_H~ScientificNames, data=metadata_R_H, permutations = 999)
bc_R_H_adonis_Plant
R_H_P_adonis <- paste0(" Plant adonis R2= ", round(bc_R_H_adonis_Plant$R2, 3), "   P = ", bc_R_H_adonis_Plant$`Pr(>F)`)

set.seed(123)
bc_R_H_adonis_Environment <- adonis2(bc_dist_R_H~Environment, data=metadata_R_H, permutations = 999)
bc_R_H_adonis_Environment
R_H_adonis_Environment <- paste0("Environment adonis R2= ", round(bc_R_H_adonis_Environment$R2, 3), "   P = ", bc_R_H_adonis_Environment$`Pr(>F)`)


#Visualization
bc_R_H_pcoa_metadata$ScientificNames <- as.factor(bc_R_H_pcoa_metadata$ScientificNames)

bc_R_H_pcoa_g12 <- ggplot(bc_R_H_pcoa_metadata, aes(PCoA1, PCoA2, fill=ScientificNames, shape = Environment)) + 
  geom_point(alpha=0.7, color = "black", size=7) +  
  scale_fill_manual(values=col_plant, name = 'Plant') +
  scale_shape_manual(values=c(21,24)) +  
  labs(x=paste("PCoA 1 (", format(100 * bc_R_H_eig[1] / sum(bc_R_H_eig), digits=4), "%)", sep=""),
       y=paste("PCoA 2 (", format(100 * bc_R_H_eig[2] / sum(bc_R_H_eig), digits=4), "%)", sep=""),
       title = paste(R_H_P_adonis, "\n",R_H_adonis_Environment)) +  
  theme_classic() +
  theme(
    plot.title = element_text(hjust = 0,vjust = 5, size = 10,face = "bold"), 
    legend.title = element_text(size=10,color="black", face = "bold"),
    legend.text= element_text(size=8,color="black"),
    axis.title = element_text(size = 10),  
    axis.text = element_text(size = 8),   
    plot.margin = margin(40, 20, 20, 20),
    legend.key.size = unit(6, "mm")      
  )+ guides(
    fill = guide_legend(order = 1,override.aes = list(size = 5)),
    shape = guide_legend(order = 2,override.aes = list(size = 5))
  )

bc_R_H_pcoa_g12 
ggsave("../../figures/Fig.2/bc_Root_Factory_Plant_shape.pdf", bc_R_H_pcoa_g12, width = 9.5, height =6)


#again  This result has the correct color for shape
bc_R_H_pcoa_g12_ <- ggplot(bc_R_H_pcoa_metadata, aes(PCoA1, PCoA2, fill=ScientificNames, shape = Environment)) + 
  geom_point(alpha=0.7, color = "black",shape = 21, size=7) +  
  scale_fill_manual(values=col_plant, name = 'Plant') +
  scale_shape_manual(values=c(21,24)) +  
  labs(x=paste("PCoA 1 (", format(100 * bc_R_H_eig[1] / sum(bc_R_H_eig), digits=4), "%)", sep=""),
       y=paste("PCoA 2 (", format(100 * bc_R_H_eig[2] / sum(bc_R_H_eig), digits=4), "%)", sep=""),
       title = paste(R_H_P_adonis, "\n",R_H_adonis_Environment)) +  
  theme_classic() +
  theme(
    plot.title = element_text(hjust = 0,vjust = 5, size = 10,face = "bold"), 
    legend.title = element_text(size=10,color="black", face = "bold"),
    legend.text= element_text(size=8,color="black"),
    axis.title = element_text(size = 10),  
    axis.text = element_text(size = 8),   
    plot.margin = margin(40, 20, 20, 20),
    legend.key.size = unit(6, "mm")          
  )+ guides(
    fill = guide_legend(order = 1,override.aes = list(size = 5)))

bc_R_H_pcoa_g12_ 
ggsave("../../figures/Fig.2/bc_Root_Factory_Plant_color.pdf", bc_R_H_pcoa_g12_, width = 9.5, height =6)







