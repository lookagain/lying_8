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
metadata <- read.table("../../data/RNA-seq_metadata.txt",header = T, sep = "\t", fill = TRUE, na.strings = "")
row.names(metadata) = metadata$SampleID

tpm <- read.csv('metal_NIL_Syncom_DGE_sleuth_norm_tpm_genome.csv', row.names = 1, check.names = F)
tpm_t <- as.data.frame(t(tpm[, grep("CK_|Cd_|Mn_", colnames(tpm))]))


#euclidean distance
euc_dist <- vegdist(tpm_t, method = "euclidean") 

euc_pcoa <- cmdscale(euc_dist, k=(10), eig=T)
euc_eig <- euc_pcoa$eig
euc_scores <- as.data.frame(euc_pcoa$points)
euc_scores <- euc_scores[match(metadata$sample, rownames(euc_scores)), ]  #确保数据框行顺序一致

euc_pcoa_metadata <- cbind(metadata[,], euc_scores[,1:5])
names(euc_pcoa_metadata)
names(euc_pcoa_metadata)[which(names(euc_pcoa_metadata) %in% c("V1", "V2", "V3", "V4", "V5"))] <- c("PCoA1", "PCoA2", "PCoA3", "PCoA4", "PCoA5")



###Permutation multivariate non-parametric test(PERMANOVA)
set.seed(123)
euc_adonis_metal <- adonis2(euc_dist~metal, data=euc_pcoa_metadata, permutations = 9999)
euc_adonis_metal
adonis_metal <- paste0(" Metal adonis R2= ", round(euc_adonis_metal$R2, 3), "   P = ", euc_adonis_metal$`Pr(>F)`)

set.seed(123)
euc_adonis_condition <- adonis2(euc_dist~condition, data=euc_pcoa_metadata, permutations = 9999)
euc_adonis_condition
adonis_condition <- paste0("Inoculation adonis R2= ", round(euc_adonis_condition$R2, 3), "   P = ", euc_adonis_condition$`Pr(>F)`, " ns")


#visualization
euc_pcoa_metadata$metal <- factor(euc_pcoa_metadata$metal, levels = c("CK", "Cd", "Mn"))
euc_pcoa_metadata$condition <- factor(euc_pcoa_metadata$condition, levels =c("NIL", "SynCom") )


#set colors
col1 <- c("#C0E2FD","#FEC0C1","#CDC6FF")

RNAseq_pcoa <- ggplot(euc_pcoa_metadata, aes(PCoA1, PCoA2, colour= metal, fill= metal, shape = condition)) +
  scale_shape_manual(values=c(21,24), name = "Inoculation") +  
  scale_fill_manual(values=col1, name ="Treat") +
  scale_color_manual(values=col1) +
  labs(x=paste("PCoA 1 (", format(100 * euc_eig[1] / sum(euc_eig), digits=4), "%)", sep=""),
       y=paste("PCoA 2 (", format(100 * euc_eig[2] / sum(euc_eig), digits=4), "%)", sep=""),
       title = paste(adonis_metal, "\n",adonis_condition))+  
  theme_classic() +
  theme(
    plot.title = element_text(hjust = 0,vjust = 5, size = 10), 
    legend.title = element_text(size=10,color="black"),
    legend.text= element_text(size=8,color="black"),
    axis.title = element_text(size = 10),  
    axis.text = element_text(size = 8),   
    plot.margin = margin(40, 20, 20, 20),
    legend.key.size = unit(6, "mm")   
  ) +
  guides(fill = guide_legend(order = 1,override.aes = list(size = 6)) ,
         shape = guide_legend(order = 2,override.aes = list(size = 6))
  )

p1 <- RNAseq_pcoa + geom_point(alpha=0.7,color = 'black', size=7)  #,shape = 21
p1
p2 <- RNAseq_pcoa + geom_point(alpha=0.7,color = 'black', size=7, shape = 21)
p2

ggsave("../../figures/Fig.7/RNAseq_euc_dis_PCoA.pdf", p1, width = 8, height =5.5)
ggsave("../../figures/Fig.7/RNAseq_euc_dis_PCoA_colors.pdf", p2, width = 8, height =5.5)




