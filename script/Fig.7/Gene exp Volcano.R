##
library(ggplot2)
library(ggrepel)
library(openxlsx)

#
pwd <- dirname(rstudioapi::getSourceEditorContext()$path)
setwd(pwd)


#
CK <- read.csv("Metal_ CK_NIL vs CK_SynCom sleuth_gene_level.csv")
CK$log10qvalue <- -log10(CK$qval)

Cd <- read.csv("Metal_ Cd_NIL vs Cd_SynCom sleuth_gene_level.csv")
Cd$log10qvalue <- -log10(Cd$qval)

Mn <- read.csv("Metal_ Mn_NIL vs Mn_SynCom sleuth_gene_level.csv")
Mn$log10qvalue <- -log10(Mn$qval)

CK_Cd <- read.csv("Metal_ CK_NIL vs Cd_NIL sleuth_gene_level.csv")
CK_Cd$log10qvalue <- -log10(CK_Cd$qval)

CK_Mn <- read.csv("Metal_ CK_NIL vs Mn_NIL sleuth_gene_level.csv")
CK_Mn$log10qvalue <- -log10(CK_Mn$qval)

cut_off_fdr <- -log10(0.05)
cut_off_log2FC = 2


##  CK
CK$Sig = ifelse(CK$log10qvalue > cut_off_fdr &   
                    abs(CK$log2fc) >= cut_off_log2FC,  
                  ifelse(CK$log2fc > 0 ,'Up','Down'),'No change')
table(CK$Sig) 

CK_range <- range(CK$log2fc)
CK_yrange <- range(CK$log10qvalue)


###绘图——基础火山图###
p1 <- ggplot(CK, aes(x =log2fc, y=log10qvalue , colour=Sig)) + 
  geom_point(alpha=0.65, size=0.8) +  
  scale_color_manual(values=c("Down" = "#00CCFF", "No change" = "grey60","Up" = "#FF0000")) + xlim(c(-10, 10)) +ylim(c(0, 0.8)) +  
  labs(x="log2FC", y="-log10(qvalue)") +  
  ggtitle("CK_NIL vs CK_SynCom") +
  theme_bw() + 
  theme(legend.text= element_text(size=8,color="black"),
        axis.title = element_text(size = 12),
        axis.title.x = element_text(hjust = 0.5,vjust = 0),
        axis.text = element_text(size=8,color="black", hjust = 0.5),
        panel.grid.minor.y = element_blank(),
        panel.grid.major.y = element_blank(),
        panel.grid.minor.x = element_blank(), 
        panel.grid.major.x = element_blank(),
        plot.title = element_text(hjust = 0.5),
        legend.position="right", 
        plot.margin = margin(10, 10, 10, 10),
        legend.title = element_blank()) 
p1

ggsave("CK_Volcano_map.pdf", plot = p1, width = 6.5, height = 5.5)



##  Cd
Cd$Sig = ifelse(Cd$log10qvalue > cut_off_fdr &   
                  abs(Cd$log2fc) >= cut_off_log2FC, 
                ifelse(Cd$log2fc > 0 ,'Up','Down'),'No change')
table(Cd$Sig) 

Cd_range <- range(Cd$log2fc)
Cd_yrange <- range(Cd$log10qvalue)

###
p2 <- ggplot(Cd, aes(x =log2fc, y=log10qvalue , colour=Sig)) + 
  geom_point(alpha=0.65, size=0.8) +  
  scale_color_manual(values=c("Down" = "#00CCFF", "No change" = "grey60","Up" = "#FF0000")) + xlim(c(-10, 10)) + ylim(c(0, 0.8)) + 
  labs(x="log2FC", y="-log10(qvalue)") +  
  ggtitle("Cd_NIL vs Cd_SynCom") +
  theme_bw() + 
  theme(legend.text= element_text(size=8,color="black"),
        axis.title = element_text(size = 12),
        axis.title.x = element_text(hjust = 0.5,vjust = 0),
        axis.text = element_text(size=8,color="black", hjust = 0.5),
        panel.grid.minor.y = element_blank(),
        panel.grid.major.y = element_blank(),
        panel.grid.minor.x = element_blank(), 
        panel.grid.major.x = element_blank(),
        plot.title = element_text(hjust = 0.5),
        legend.position="right", 
        plot.margin = margin(10, 10, 10, 10),
        legend.title = element_blank()) 
p2

ggsave("Cd_Volcano_map.pdf", plot = p2, width = 6.5, height = 5.5)



##  Mn
Mn$Sig = ifelse(Mn$log10qvalue > cut_off_fdr &    
                  abs(Mn$log2fc) >= cut_off_log2FC,  
                ifelse(Mn$log2fc > 0 ,'Up','Down'),'No change')
table(Mn$Sig) 

Mn_range <- range(Mn$log2fc)
Mn_yrange <- range(Mn$log10qvalue)

#
p3 <- ggplot(Mn, aes(x =log2fc, y=log10qvalue , colour=Sig)) + 
  geom_point(alpha=0.65, size=0.8) + 
  scale_color_manual(values=c("Down" = "#00CCFF", "No change" = "grey60","Up" = "#FF0000")) +ylim(c(0, 3.5))+ 
  geom_vline(xintercept=c(-cut_off_log2FC,cut_off_log2FC),lty=5,col="black",lwd=0.1) + 
  geom_hline(yintercept = cut_off_fdr, lty=5,col="black",lwd=0.1) +  
  labs(x="log2FC", y="-log10(qvalue)") +  
  ggtitle("Mn_NIL vs Mn_SynCom") +
  theme_bw() + 
  theme(legend.text= element_text(size=8,color="black"),
        axis.title = element_text(size = 12),
        axis.title.x = element_text(hjust = 0.5,vjust = 0),
        axis.text = element_text(size=8,color="black", hjust = 0.5),
        panel.grid.minor.y = element_blank(),
        panel.grid.major.y = element_blank(),
        panel.grid.minor.x = element_blank(), 
        panel.grid.major.x = element_blank(),
        plot.title = element_text(hjust = 0.5),
        legend.position="right", 
        plot.margin = margin(10, 10, 10, 10),
        legend.title = element_blank()) 
p3

ggsave("Mn_Volcano_map.pdf", plot = p3, width = 6.5, height = 5.5)



###  CK vs Cd 
CK_Cd$Sig = ifelse(CK_Cd$log10qvalue > cut_off_fdr &    
                  abs(CK_Cd$log2fc) >= cut_off_log2FC,  
                ifelse(CK_Cd$log2fc > 0 ,'Up','Down'),'No change')
table(CK_Cd$Sig) 

CK_Cd_range <- range(CK_Cd$log2fc)
CK_Cd_yrange <- range(CK_Cd$log10qvalue)


#
p4 <- ggplot(CK_Cd, aes(x =log2fc, y=log10qvalue , colour=Sig)) + 
  geom_point(alpha=0.65, size=0.8) +  
  scale_color_manual(values=c("Down" = "#00CCFF", "No change" = "grey60","Up" = "#FF0000")) + xlim(c(-8, 8)) + ylim(c(0, 0.3))+  
  labs(x="log2FC", y="-log10(qvalue)") +  
  ggtitle("CK_NIL vs Cd_NIL") +
  theme_bw() + 
  theme(legend.text= element_text(size=8,color="black"),
        axis.title = element_text(size = 12),
        axis.title.x = element_text(hjust = 0.5,vjust = 0),
        axis.text = element_text(size=8,color="black", hjust = 0.5),
        panel.grid.minor.y = element_blank(),
        panel.grid.major.y = element_blank(),
        panel.grid.minor.x = element_blank(), 
        panel.grid.major.x = element_blank(),
        plot.title = element_text(hjust = 0.5),
        legend.position="right", 
        plot.margin = margin(10, 10, 10, 10),
        legend.title = element_blank()) 
p4

ggsave("CK_Cd_Volcano_map.pdf", plot = p4, width = 6.5, height = 5.5)



###  CK vs Mn
CK_Mn$Sig = ifelse(CK_Mn$log10qvalue > cut_off_fdr &   
                     abs(CK_Mn$log2fc) >= cut_off_log2FC,  
                   ifelse(CK_Mn$log2fc > 0 ,'Up','Down'),'No change')
table(CK_Mn$Sig) 

CK_Mn_range <- range(CK_Mn$log2fc)
CK_Mn_yrange <- range(CK_Mn$log10qvalue)


#
p5 <- ggplot(CK_Mn, aes(x =log2fc, y=log10qvalue , colour=Sig)) + 
  geom_point(alpha=0.65, size=0.8) +  
  scale_color_manual(values=c("Down" = "#00CCFF", "No change" = "grey60","Up" = "#FF0000")) +ylim(c(0, 3.5))+  
  geom_vline(xintercept=c(-cut_off_log2FC,cut_off_log2FC),lty=5,col="black",lwd=0.1) + 
  geom_hline(yintercept = cut_off_fdr, lty=5,col="black",lwd=0.1) + 
  labs(x="log2FC", y="-log10(qvalue)") + 
  ggtitle("CK_NIL vs Mn_NIL") +
  theme_bw() + 
  theme(legend.text= element_text(size=8,color="black"),
        axis.title = element_text(size = 12),
        axis.title.x = element_text(hjust = 0.5,vjust = 0),
        axis.text = element_text(size=8,color="black", hjust = 0.5),
        panel.grid.minor.y = element_blank(),
        panel.grid.major.y = element_blank(),
        panel.grid.minor.x = element_blank(), 
        panel.grid.major.x = element_blank(),
        plot.title = element_text(hjust = 0.5),
        legend.position="right", 
        plot.margin = margin(10, 10, 10, 10),
        legend.title = element_blank()) 
p5

ggsave("CK_Mn_Volcano_map.pdf", plot = p5, width = 6.5, height = 5.5)



##补充一个CK_NIL vs Mn_SynCom
CK_Mn_SynCom <- read.csv("Metal_ CK_NIL vs Mn_SynCom sleuth_gene_level.csv")
CK_Mn_SynCom$log10qvalue <- -log10(CK_Mn_SynCom$qval)



CK_Mn_SynCom$Sig = ifelse(CK_Mn_SynCom$log10qvalue > cut_off_fdr &    
                     abs(CK_Mn_SynCom$log2fc) >= cut_off_log2FC,  
                   ifelse(CK_Mn_SynCom$log2fc > 0 ,'Up','Down'),'No change')
table(CK_Mn_SynCom$Sig) 

CK_Mn_SynCom_range <- range(CK_Mn_SynCom$log2fc)
CK_Mn_SynCom_yrange <- range(CK_Mn_SynCom$log10qvalue)


#
p6 <- ggplot(CK_Mn_SynCom, aes(x =log2fc, y=log10qvalue , colour=Sig)) + 
  geom_point(alpha=0.65, size=0.8) +  
  scale_color_manual(values=c("Down" = "#00CCFF", "No change" = "grey60","Up" = "#FF0000")) +ylim(c(0, 2.5))+  
  geom_vline(xintercept=c(-cut_off_log2FC,cut_off_log2FC),lty=5,col="black",lwd=0.1) + 
  geom_hline(yintercept = cut_off_fdr, lty=5,col="black",lwd=0.1) + 
  labs(x="log2FC", y="-log10(qvalue)") + 
  ggtitle("CK_NIL vs Mn_SynCom") +
  theme_bw() + 
  theme(legend.text= element_text(size=8,color="black"),
        axis.title = element_text(size = 12),
        axis.title.x = element_text(hjust = 0.5,vjust = 0),
        axis.text = element_text(size=8,color="black", hjust = 0.5),
        panel.grid.minor.y = element_blank(),
        panel.grid.major.y = element_blank(),
        panel.grid.minor.x = element_blank(), 
        panel.grid.major.x = element_blank(),
        plot.title = element_text(hjust = 0.5),
        legend.position="right", 
        plot.margin = margin(10, 10, 10, 10),
        legend.title = element_blank()) 
p6

ggsave("CK_NIL_vs_Mn_SynCom_Volcano_map.pdf", plot = p6, width = 6.5, height = 5.5)



