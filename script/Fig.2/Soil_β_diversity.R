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


###Soil Site ####
metadata_S <- subset(metadata, Area_Microbiome == "Soil")
bc_dist_S <- bc_dist[row.names(bc_dist) %in% row.names(metadata_S), names(bc_dist) %in% row.names(metadata_S)]
bc_S_pcoa <- cmdscale(bc_dist_S, k=(10), eig=T)
bc_S_eig <- bc_S_pcoa$eig
bc_S_scores <- as.data.frame(bc_S_pcoa$points)
bc_S_pcoa_metadata <- cbind(metadata_S[,], bc_S_scores[,1:5])
names(bc_S_pcoa_metadata)
names(bc_S_pcoa_metadata)[which(names(bc_S_pcoa_metadata) %in% c("V1", "V2", "V3", "V4", "V5"))] <- c("PCoA1", "PCoA2", "PCoA3", "PCoA4", "PCoA5")


###Permutation multivariate non-parametric test(PERMANOVA)
set.seed(123)
bc_S_adonis_Site <- adonis2(bc_dist_S~Site, data=bc_S_pcoa_metadata, permutations = 999)
bc_S_adonis_Site
S_adonis_Site <- paste0(" Site adonis R2= ", round(bc_S_adonis_Site$R2, 3), "   P = ", bc_S_adonis_Site$`Pr(>F)`)
#write.table(bc_S_adonis_Site,"Soil_Site_adonis_check.txt",sep = "\t",quote = FALSE)


#visualization
bc_S_pcoa_metadata$Site <- factor(bc_S_pcoa_metadata$Site,levels = c("Village","Factory"))
bc_S_pcoa_metadata$Environment <- as.factor(bc_S_pcoa_metadata$Environment)

#set colors
col_Site <- c("#8DD3C7","#FFFFB3")

bc_Soil_pcoa_g12_ <- ggplot(bc_S_pcoa_metadata, aes(PCoA1, PCoA2, colour=Site, fill=Site, shape = Environment)) +
  geom_point(alpha=0.7,color = 'black',shape = 21, size=7) +  
  scale_shape_manual(values=c(21,24)) +  
  scale_fill_manual(values=col_Site) +
  scale_color_manual(values=col_Site) +
  labs(x=paste("PCoA 1 (", format(100 * bc_S_eig[1] / sum(bc_S_eig), digits=4), "%)", sep=""),
       y=paste("PCoA 2 (", format(100 * bc_S_eig[2] / sum(bc_S_eig), digits=4), "%)", sep=""),
       title = S_adonis_Site) +  
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


bc_Soil_pcoa_g12_
ggsave("../../figures/Fig.2/bc_dis_Soil_Site.pdf", bc_Soil_pcoa_g12_, width = 8, height =5.5)



############             Soil Site Plants            ################
#Soil Village Plants
metadata_S_V <- subset(metadata_S, Site == "Village")
table(metadata_S_V$ScientificNames)
name <- unique(metadata_S_V$ScientificNames)


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
bc_dist_S_V <- bc_dist[row.names(bc_dist) %in% row.names(metadata_S_V), names(bc_dist) %in% row.names(metadata_S_V)]
bc_S_V_pcoa <- cmdscale(bc_dist_S_V, k=(10), eig=T)
bc_S_V_eig <- bc_S_V_pcoa$eig
bc_S_V_scores <- as.data.frame(bc_S_V_pcoa$points)
bc_S_V_pcoa_metadata <- cbind(metadata_S_V[,], bc_S_V_scores[,1:5])
names(bc_S_V_pcoa_metadata)
names(bc_S_V_pcoa_metadata)[which(names(bc_S_V_pcoa_metadata) %in% c("V1", "V2", "V3", "V4", "V5"))] <- c("PCoA1", "PCoA2", "PCoA3", "PCoA4", "PCoA5")


###Permutation multivariate non-parametric test(PERMANOVA)
set.seed(110)
bc_S_V_adonis_Plant <- adonis2(bc_dist_S_V~ScientificNames, data=metadata_S_V, permutations = 999)
bc_S_V_adonis_Plant
S_V_P_adonis <- paste0(" Plant adonis R2= ", round(bc_S_V_adonis_Plant$R2, 3), "   P = ", bc_S_V_adonis_Plant$`Pr(>F)`)
#write.table(bc_S_V_adonis_Plant,"Soil_Village_Plant_adonis_check.txt",sep = "\t",quote = FALSE)

set.seed(123)
bc_S_V_adonis_Environment <- adonis2(bc_dist_S_V~Environment, data=metadata_S_V, permutations = 999)
bc_S_V_adonis_Environment
S_V_adonis_Environment <- paste0("Environment adonis R2= ", round(bc_S_V_adonis_Environment$R2, 3), "   P = ", bc_S_V_adonis_Environment$`Pr(>F)`)
#write.table(bc_S_V_adonis_Environment,"Soil_Village_Environment_Plant_adonis_check.txt",sep = "\t",quote = FALSE)


#Visualization
bc_S_V_pcoa_metadata$ScientificNames <- as.factor(bc_S_V_pcoa_metadata$ScientificNames)
Plant_color <- col_Plant[name]


bc_S_V_pcoa_g12 <- ggplot(bc_S_V_pcoa_metadata, aes(PCoA1, PCoA2, fill=ScientificNames, shape = Environment)) + 
  geom_point(alpha=0.7, color = "black", size=7) +  
  scale_fill_manual(values=Plant_color, name = 'Plant') +
  scale_shape_manual(values=c(21,24)) +  
  labs(x=paste("PCoA 1 (", format(100 * bc_S_V_eig[1] / sum(bc_S_V_eig), digits=4), "%)", sep=""),
       y=paste("PCoA 2 (", format(100 * bc_S_V_eig[2] / sum(bc_S_V_eig), digits=4), "%)", sep=""),
       title = paste(S_V_P_adonis, "\n",S_V_adonis_Environment)) +  
  theme_classic() +
  theme(
    plot.title = element_text(hjust = 0,vjust = 5, size = 10,face = "bold"), 
    legend.title = element_text(size=10,color="black", face = "bold"),
    legend.text= element_text(size=8,color="black"),
    axis.title = element_text(size = 10),  
    axis.text = element_text(size = 8),   
    plot.margin = margin(40, 20, 20, 20),
    legend.key.size = unit(6, "mm")   
  ) + guides(fill = guide_legend(order = 1,override.aes = list(size = 5)),
             shape =guide_legend(order = 2,override.aes = list(size = 5))
  )

bc_S_V_pcoa_g12 

ggsave("../../figures/Fig.2/bc_Soil_Village_Plant_shape.pdf", bc_S_V_pcoa_g12, width = 9.5, height =6)


#again  This result has the correct color for shape
bc_S_V_pcoa_g12_ <- ggplot(bc_S_V_pcoa_metadata, aes(PCoA1, PCoA2, fill=ScientificNames, shape = Environment)) + 
  geom_point(alpha=0.7, color = "black",shape = 21, size=7) +  
  scale_fill_manual(values=Plant_color, name = 'Plant') +
  scale_shape_manual(values=c(21,24)) +  #  ,24 三角形
  labs(x=paste("PCoA 1 (", format(100 * bc_S_V_eig[1] / sum(bc_S_V_eig), digits=4), "%)", sep=""),
       y=paste("PCoA 2 (", format(100 * bc_S_V_eig[2] / sum(bc_S_V_eig), digits=4), "%)", sep=""),
       title = paste(S_V_P_adonis, "\n",S_V_adonis_Environment)) +  
  theme_classic() +
  theme(
    plot.title = element_text(hjust = 0,vjust = 5, size = 10,face = "bold"), 
    legend.title = element_text(size=10,color="black", face = "bold"),
    legend.text= element_text(size=8,color="black"),
    axis.title = element_text(size = 10),  
    axis.text = element_text(size = 8),   
    plot.margin = margin(40, 20, 20, 20),
    legend.key.size = unit(6, "mm")   
  ) + guides(fill = guide_legend(order = 1,override.aes = list(size = 5)))

bc_S_V_pcoa_g12_ 
ggsave("../../figures/Fig.2/bc_Soil_Village_Plant_color.pdf", bc_S_V_pcoa_g12_, width = 9.5, height =6)



#Soil Factory  Plants
metadata_S_H <- subset(metadata_S, Site == "Factory")
table(metadata_S_H$ScientificNames)
col_plant <- c("#FFFF99", "#FFCCFF", "#FFCC66", "#FF99FF", "#99CC33", "#FF6600", "#FF00CC", "#99FFFF", "#99FF66", "#99CCFF") 

#PcoA
bc_dist_S_H <- bc_dist[row.names(bc_dist) %in% row.names(metadata_S_H), names(bc_dist) %in% row.names(metadata_S_H)]
bc_S_H_pcoa <- cmdscale(bc_dist_S_H, k=(10), eig=T)
bc_S_H_eig <- bc_S_H_pcoa$eig
bc_S_H_scores <- as.data.frame(bc_S_H_pcoa$points)
bc_S_H_pcoa_metadata <- cbind(metadata_S_H[,], bc_S_H_scores[,1:5])
names(bc_S_H_pcoa_metadata)
names(bc_S_H_pcoa_metadata)[which(names(bc_S_H_pcoa_metadata) %in% c("V1", "V2", "V3", "V4", "V5"))] <- c("PCoA1", "PCoA2", "PCoA3", "PCoA4", "PCoA5")


###Permutation multivariate non-parametric test(PERMANOVA)
set.seed(123)
bc_S_H_adonis_Plant <- adonis2(bc_dist_S_H~ScientificNames, data=metadata_S_H, permutations = 999)
bc_S_H_adonis_Plant
S_H_P_adonis <- paste0(" Plant adonis R2= ", round(bc_S_H_adonis_Plant$R2, 3), "   P = ", bc_S_H_adonis_Plant$`Pr(>F)`)
#write.table(bc_S_H_adonis_Plant,"Soil_Factory_Plant_adonis_check.txt",sep = "\t",quote = FALSE)

set.seed(123)
bc_S_H_adonis_Environment <- adonis2(bc_dist_S_H~Environment, data=metadata_S_H, permutations = 999)
bc_S_H_adonis_Environment
S_H_adonis_Environment <- paste0("Environment adonis R2= ", round(bc_S_H_adonis_Environment$R2, 3), "   P = ", bc_S_H_adonis_Environment$`Pr(>F)`)
#write.table(bc_S_H_adonis_Environment,"Soil_Factory_Environment_Plant_adonis_check.txt",sep = "\t",quote = FALSE)


#Visualization
bc_S_H_pcoa_metadata$ScientificNames <- as.factor(bc_S_H_pcoa_metadata$ScientificNames)

bc_S_H_pcoa_g12 <- ggplot(bc_S_H_pcoa_metadata, aes(PCoA1, PCoA2, fill=ScientificNames, shape = Environment)) + 
  geom_point(alpha=0.7, color = "black", size=7) +  
  scale_fill_manual(values=col_plant, name = 'Plant') +
  scale_shape_manual(values=c(21,24)) +  
  labs(x=paste("PCoA 1 (", format(100 * bc_S_H_eig[1] / sum(bc_S_H_eig), digits=4), "%)", sep=""),
       y=paste("PCoA 2 (", format(100 * bc_S_H_eig[2] / sum(bc_S_H_eig), digits=4), "%)", sep=""),
       title = paste(S_H_P_adonis, "\n",S_H_adonis_Environment)) +  
  theme_classic() +
  theme(
    plot.title = element_text(hjust = 0,vjust = 5, size = 10,face = "bold"), 
    legend.title = element_text(size=10,color="black", face = "bold"),
    legend.text= element_text(size=8,color="black"),
    axis.title = element_text(size = 10),  
    axis.text = element_text(size = 8),   
    plot.margin = margin(40, 20, 20, 20),
    legend.key.size = unit(6, "mm")   
  )+ guides(fill = guide_legend(order = 1,override.aes = list(size = 5)),
            shape =guide_legend(order = 2,override.aes = list(size = 5))
  )

bc_S_H_pcoa_g12 

ggsave("../../figures/Fig.2/bc_Soil_Factory_Plant_shape.pdf", bc_S_H_pcoa_g12, width = 9.5, height =6)


#again  This result has the correct color for shape
bc_S_H_pcoa_g12_ <- ggplot(bc_S_H_pcoa_metadata, aes(PCoA1, PCoA2, fill=ScientificNames, shape = Environment)) + 
  geom_point(alpha=0.7, color = "black",shape = 21, size=7) +  
  scale_fill_manual(values=col_plant, name = 'Plant') +
  scale_shape_manual(values=c(21,24)) +  
  labs(x=paste("PCoA 1 (", format(100 * bc_S_H_eig[1] / sum(bc_S_H_eig), digits=4), "%)", sep=""),
       y=paste("PCoA 2 (", format(100 * bc_S_H_eig[2] / sum(bc_S_H_eig), digits=4), "%)", sep=""),
       title = paste(S_H_P_adonis, "\n",S_H_adonis_Environment)) +  
  theme_classic() +
  theme(
    plot.title = element_text(hjust = 0,vjust = 5, size = 10,face = "bold"), 
    legend.title = element_text(size=10,color="black", face = "bold"),
    legend.text= element_text(size=8,color="black"),
    axis.title = element_text(size = 10),  
    axis.text = element_text(size = 8),   
    plot.margin = margin(40, 20, 20, 20),
    legend.key.size = unit(6, "mm")      
  )+ guides(fill = guide_legend(order = 1,override.aes = list(size = 5)))

bc_S_H_pcoa_g12_ 

ggsave("../../figures/Fig.2/bc_Soil_Factory_Plant_color.pdf", bc_S_H_pcoa_g12_, width = 9.5, height =6)







