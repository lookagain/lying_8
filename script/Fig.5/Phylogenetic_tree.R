##
library(ape)
library(ggtree)
library(ggtreeExtra)
library(ggplot2)
library(ggnewscale)
library(dplyr)
library(showtext)
showtext_auto()


#
pwd <- dirname(rstudioapi::getSourceEditorContext()$path)
setwd(pwd)


#
candidate_ASV <- read.csv('../../data/DESeq2 Root Site important ASV.csv')
R_E_ASV <- read.csv('../../data/DESeq2 Factory Root Environment important 7 ASV overlap with 46 ASV.csv',check.names = F)
tree <- read.tree("../../data/Newick Export branch length.nwk")
abundance <- read.csv("../../data/Root Site 46ASV relative abundant.csv")


#
name1 <- c('Beijerinckiaceae', 'Comamonadaceae', 'Sphingomonadaceae', 'Hyphomicrobiaceae', 'Microscillaceae')
candidate_ASV$F_color <- ifelse(candidate_ASV$Family %in% name1, candidate_ASV$Family, "Others") 
candidate_ASV$F_color <- factor(candidate_ASV$F_color, levels = c(name1, "Others"))
candidate_ASV$phylum <- factor(candidate_ASV$phylum, levels = c('Proteobacteria', 'Actinobacteriota', 'Firmicutes', 'Bacteroidota'))

R_E_ASV <- R_E_ASV[,colnames(R_E_ASV) %in% c('ASV', 'Enriched')]
R_E_ASV
indoor <- c('ASV_224', 'ASV_615')
outdoor <- c('ASV_842', 'ASV_1267', 'ASV_1668', 'ASV_1751', 'ASV_1989')

candidate_ASV <- candidate_ASV %>%
  mutate(Environment = case_when(
    ASV %in% indoor  ~ 'indoor',
    ASV %in% outdoor ~ 'outdoor',
    TRUE             ~ 'ns'
  ))

data1 <- candidate_ASV[, colnames(candidate_ASV) %in% c('ASV', 'phylum', 'Family', 'F_color', 'log2FoldChange', 'Environment')]
data1 <- data1[match(tree$tip.label, data1$ASV), ]


#使用的7个菌
use_ASV <- c('ASV_575', 'ASV_1139', 'ASV_1200', 'ASV_2714', 'ASV_349', 'ASV_842', 'ASV_482')

data1$Inoculation <- ifelse(data1$ASV %in% use_ASV, "*", "")


##颜色
F_colors = c(
  "Beijerinckiaceae" = "#1f77b4",  
  "Sphingomonadaceae" = "#17becf",  
  "Microscillaceae" = "#7f7f7f",    
  "Comamonadaceae" = "#ff7f0e",    
  "Hyphomicrobiaceae" = "#8c564b",    
  "Others" = "#C9992C"
)

phylum_color <- c("Proteobacteria" = "#E41A1C",
                  "Actinobacteriota" = "#596A98",
                  "Firmicutes" = "#449B75",
                  "Bacteroidota" = "#6B886D")


##
p <- ggtree(tree, layout = "fan", branch.length="none",ladderize = F, linewidth = 0.6,open.angle = 20, color = "grey" )

p

p <- p %<+% data1[,]  #要保证data1的第一列与p中的label列相同,内容顺序倒是没有要求一致
head(p$data)
p$data$label <- gsub('_', '', p$data$label)

p1 <- p + geom_tree(aes(color = F_color), linewidth = 0.8) + 
  guides(color = guide_legend(title = "Family", override.aes = list(linewidth = 1.2)))+
  geom_tiplab(aes(label=label, color = F_color), show.legend = FALSE, size=4.5, hjust=-0.05, vjust = 0.5,offset=0, family="Arial", fontface=2)+
  scale_color_manual(values = F_colors)
p1

p2 <- p1 + geom_tiplab(aes(label=Inoculation), color = "red", show.legend = FALSE, size=7, hjust=-8, vjust = 0.8, offset=0, family="Arial", fontface=2)
p2


#
tipnode <- p$data[which(p$data$isTip==TRUE),]#仅保留外部节点
data1$nodeid <- tipnode$parent
unique_tem <- data1 %>% distinct(nodeid, .keep_all = TRUE)
df1 <- data.frame(Phylum =unique_tem$phylum, nodeid= unique_tem$nodeid)
df1$label1 <- ""


p3 <- p2 + new_scale_color() +
  geom_cladelab(data=df1,
                mapping =aes(node=nodeid, label= label, color=Phylum),
                extend = 0.5, barsize=4,offset = 8.5,
                show.legend = c(textcolor = FALSE)) +
  scale_color_manual(values = phylum_color,
                     guide = guide_legend(
                       override.aes = list(label = "")))
p3


ab_F <- abundance[abundance$Site == "Factory", c(7,1:6)]
ab_F$ASV_ID <- gsub('_', '', ab_F$ASV_ID)
test <- p$data


p4 <- p3 +
  geom_fruit(
    data = ab_F,
    geom = geom_bar,
    mapping = aes(y = ASV_ID, x = Mean),
    fill = "#FFFFB3",
    color= "black",
    pwidth = 0.7,
    width = 0.6,
    stat="identity",
    offset = 1.1
  ) 
p4


p5 <- p4+new_scale_fill()+             
  geom_fruit(geom = geom_tile,mapping = aes(fill = -log2FoldChange),pwidth = 1.5,width = 1.5,offset =0)+
  scale_fill_gradientn(colours = c( "#E0E0E0", "#CCD5DC","#B8CAD6","#A7BECF","#92B4CD", 
                                    "#7CAACC", "#649AC8","#478DC1","#0164AF"), 
                       breaks = seq(2,29,3), limits=c(2,29))+  
  guides(fill = guide_colorbar(reverse = T))+   
  theme( legend.background=element_rect(fill=NA), legend.text = element_text(color = "black",size = 10),
         legend.title =element_text(size = 11,colour = "black"),
         legend.key.width = unit(0.65, "cm"),
         legend.key.height = unit(0.62, "cm"),
         legend.key.size = unit(6, "mm"), 
         legend.spacing = unit(0.1, "mm"), 
         )+xlim(0,NA)
p5 

ggsave(filename="../../figures/Fig.5/candidate_46ASV_tree_legend.pdf",p5, height=24, width=32, units = "cm")


col_E <- c('indoor' = "#99CC00", 'outdoor' = "#FFCC00", ns = 'white')

p6 <- p5 + new_scale_fill()+ 
    geom_fruit(
      geom = geom_tile, mapping = aes(fill = Environment),
      color= "#009999", pwidth = 1.5, width = 1.5, offset = 0.15
    ) +
  scale_fill_manual(values = col_E)
p6

ggsave(filename="../../figures/Fig.5/candidate_46ASV_tree_figure.pdf",p6, height=24, width=32, units = "cm")


