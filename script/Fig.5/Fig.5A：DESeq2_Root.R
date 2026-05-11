#library package
library(microeco)
library(DESeq2)
library(ggplot2)
library(RColorBrewer)
library(stringr)

# Set Work Path
pwd <- dirname(rstudioapi::getSourceEditorContext()$path)
setwd(pwd)

#import
Table_ASV <- read.csv("../../data/Table_ASV.csv",check.names = F,row.names = 1)
Taxonomy_ASV <- read.csv("../../data/Taxonomy_ASV.csv",check.names = F,row.names = 1)
metadata <- read.table("../../data/Metadata.txt", sep = "\t",header = T, na.strings = "",row.names = 1)


#Root Site
metadata_Root <- metadata[metadata$Area_Microbiome == "Root",]
metadata_Root_H <- metadata_Root[metadata_Root$Site == "Factory",]
metadata_Root_V <- metadata_Root[metadata_Root$Site == "Village",] 
metadata_Root <- rbind(metadata_Root_H,metadata_Root_V)

ASV <- Table_ASV
ASV <- ASV[,colnames(ASV) %in% rownames(metadata_Root)]

tax <- Taxonomy_ASV
tax$ASV <- rownames(tax)


#Root Site DESeq2
asv_dataset <- microtable$new(sample_table = metadata_Root,
                              otu_table = ASV, 
                              tax_table = tax)   #15368 taxa with 0 abundance are removed from the otu_table
asv_dataset$tidy_dataset()
asv_dataset

t1 <- trans_diff$new(dataset = asv_dataset, method = "DESeq2", group = "Site",
                     filter_thres = 0.0001,
                     taxa_level = "ASV",
                     alpha = 0.01, p_adjust_method = "fdr",
                     remove_unknown = F)

t1$res_diff <- t1$res_diff[order(t1$res_diff$log2FoldChange,decreasing = F),]
t1_res_diff <- t1$res_diff
t1_res_diff$phylum <- gsub(".*?p__(.*?)(;|\\|).*", "\\1", t1_res_diff$Taxa)
t1_res_diff$ASV <-  gsub(".*(ASV.*)", "\\1", t1_res_diff$Taxa)
t1_res_diff$Family <- str_extract(t1_res_diff$Taxa, "(?<=f__)[^;|]+")



#Visualization
#biomakers
biomaker_F <- read.csv("../../data/rf_family_R_top30_res-abund-wilcox-sig.csv",check.names = F)
biomaker_F <- biomaker_F[biomaker_F$Group == 'Factory',]

sort(table(t1_res_diff$phylum), decreasing = TRUE)
t1_res_diff$hight_phylum <- t1_res_diff$phylum
hight_phylum <- c("Proteobacteria","Actinobacteriota", "Firmicutes","Bacteroidota","Myxococcota","Chloroflexi","Verrucomicrobiota")
t1_res_diff[!t1_res_diff$hight_phylum %in% hight_phylum,"hight_phylum"] <- "Others"


t1_res_diff$hight_phylum <-  factor(t1_res_diff$hight_phylum,
                                      levels=c(hight_phylum,"Others"), ordered=TRUE)

t1_res_diff <- t1_res_diff[order(t1_res_diff$hight_phylum),]

t1_res_diff[,"P.adj"][is.na(t1_res_diff[,"P.adj"])] <- 1

t1_res_diff$DCB <- as.factor(ifelse(t1_res_diff$P.adj >= 0.05,"NS",
                                                ifelse(t1_res_diff$log2FoldChange > 0,"Enriched","Depleted")))

t1_res_diff$DCB_lab <- as.factor(ifelse(t1_res_diff$P.adj >= 0.05,"",
                                    ifelse(t1_res_diff$log2FoldChange < 0, as.character(t1_res_diff$ASV), "")))



#sort
t1_res_diff$num <- as.numeric(gsub("ASV_", "", t1_res_diff$ASV))

t1_res_diff_Pro <- t1_res_diff[t1_res_diff$hight_phylum %in% "Proteobacteria",]
t1_res_diff_Pro <- t1_res_diff_Pro[order(t1_res_diff_Pro$num),]
t1_res_diff_Act <- t1_res_diff[t1_res_diff$hight_phylum %in% "Actinobacteriota",]
t1_res_diff_Act <- t1_res_diff_Act[order(t1_res_diff_Act$num),]
t1_res_diff_Fir <- t1_res_diff[t1_res_diff$hight_phylum %in% "Firmicutes",]
t1_res_diff_Fir <- t1_res_diff_Fir[order(t1_res_diff_Fir$num),]
t1_res_diff_Bac <- t1_res_diff[t1_res_diff$hight_phylum %in% "Bacteroidota",]
t1_res_diff_Bac <- t1_res_diff_Bac[order(t1_res_diff_Bac$num),]
t1_res_diff_Myx <- t1_res_diff[t1_res_diff$hight_phylum %in% "Myxococcota",]
t1_res_diff_Myx <- t1_res_diff_Myx[order(t1_res_diff_Myx$num),]
t1_res_diff_Chl <- t1_res_diff[t1_res_diff$hight_phylum %in% "Chloroflexi",]
t1_res_diff_Chl <- t1_res_diff_Chl[order(t1_res_diff_Chl$num),]
t1_res_diff_Ver <- t1_res_diff[t1_res_diff$hight_phylum %in% "Verrucomicrobiota",]
t1_res_diff_Ver <- t1_res_diff_Ver[order(t1_res_diff_Ver$num),]
t1_res_diff_Oth <- t1_res_diff[t1_res_diff$hight_phylum %in% "Others",]
t1_res_diff_Oth <- t1_res_diff_Oth[order(t1_res_diff_Oth$num),]

t1_res_diff_order <- rbind(t1_res_diff_Pro, t1_res_diff_Act, t1_res_diff_Fir, t1_res_diff_Bac, t1_res_diff_Myx, t1_res_diff_Chl, t1_res_diff_Ver, t1_res_diff_Oth)

num_order <- t1_res_diff_order$num

t1_res_diff_order$num <-factor(t1_res_diff_order$num,levels = num_order)


library(dplyr)
F_sig_ASV <- t1_res_diff[, c('Family', "DCB_lab")] %>% 
  filter(
    !is.na(Family) & Family != "" & Family %in% biomaker_F$Taxa,   
    !is.na(DCB_lab) & DCB_lab != "" 
  )

F_sig_ASV$ASV <- F_sig_ASV$DCB_lab
F_sig_ASV$DCB_lab <- gsub("_", "", F_sig_ASV$DCB_lab)
F_sig_ASV$label <- paste(F_sig_ASV$DCB_lab, F_sig_ASV$Family, sep = "|")
F_sig_ASV <- F_sig_ASV[,c("label","ASV")]


t1_res_diff_order <- merge(t1_res_diff_order, F_sig_ASV, "ASV", all=T)
t1_res_diff_order[is.na(t1_res_diff_order)] <- ''


#colors
col_tax <- c("Proteobacteria" = "#E41A1C",
             "Actinobacteriota" = "#596A98",
             "Firmicutes" = "#449B75",
             "Bacteroidota" = "#6B886D",
             "Acidobacteriota" = "#C66764",
             "Chloroflexi" = "#FF7F00",
             "Myxococcota" = "#AC5782",
             "Verrucomicrobiota" = "#FFE528",
             "Nitrospirota" = "#E485B7",
             "Planctomycetota" = "#00FF33",
             "Others" = "#C9992C",
             "Desulfobacterota" = "#3366FF",
             "Entotheonellaeota" = "#33FFFF",
             "Spirochaetota" = "#00CCCC",
             "Gemmatimonadota" = "#66CC00",
             "Dependentiae" = "#FF33FF",
             "Patescibacteria" = "#FFCCFF",
             "Bdellovibrionota" = "#CCFF99",
             "NB1-j" = "#CCFFFF"
)


P_manhattan_Site_Root <- 
  ggplot(t1_res_diff_order, aes(num, -log10(P.adj), color=hight_phylum, shape=DCB, fill=hight_phylum)) +
  geom_point(cex=log(t1_res_diff_order$baseMean + 1), alpha = 0.6) +
  scale_shape_manual(values=c(25,24,21)) +
  scale_fill_manual(values=col_tax) +
  scale_color_manual(values=col_tax) +
  geom_hline(aes(yintercept = 2), colour="red", linetype="dashed") +
  geom_text(label = t1_res_diff_order$label, size = log(t1_res_diff_order$baseMean + 1, base = 3), colour="grey50", vjust = 1.5,hjust = 0.3) + 
  theme_bw() + 
  labs(x=paste("Bacterial ASVs"),y=paste("-log10(P.adj)")) +
  ggtitle(paste("Root | Village vs Factory")) +
  theme_bw() + 
  theme(legend.position= "right") +
  theme(panel.grid = element_blank(), 
        panel.border = element_blank(),
        plot.title = element_text(size = 10, colour = "black"),
        legend.title = element_text(size = 12, colour = "black"),
        legend.text = element_text(size = 12, colour = "black"),
        axis.title = element_text(size = 12, colour = "black"),
        axis.ticks.x = element_blank(),
        axis.text.x = element_blank(),
        axis.text.y = element_text(size = 10, colour = "black"),
        axis.line.x = element_line(color = "black", linewidth = 1.5),
        axis.line.y = element_line(color = "black", linewidth = 0.8)) + 
  guides(colour = "none") + 
  guides(size = "none", fill = "none") + 
  guides(color = guide_legend(title = "Phylum"))

P_manhattan_Site_Root

ggsave("../../figures/Fig.5/P_manhattan_Village-vs-Factory_Root.pdf", P_manhattan_Site_Root, width = 8, height = 4)
write.csv(t1_res_diff_order,"Root_Site_DEseq2_sig_result_visual.csv",row.names = F)





