#library package
library(ggplot2)
library(dplyr)
library(EasyStat)
library(Rmisc)
library(ggpubr)


# Set Work Path
pwd <- dirname(rstudioapi::getSourceEditorContext()$path)
setwd(pwd)


#import data
metadata <- read.table("../../data/Metadata.txt",header = T, sep = "\t", fill = TRUE, na.strings = "")
head(metadata)
row.names(metadata) = metadata$SampleID

shannon <- read.table("../../data/shannon_vector.tsv", header = T, sep ="\t", fill = T, na.strings = " ")
colnames(shannon)[1] <- c("SampleID")

alpha <- merge(metadata,shannon,"SampleID")


### Root       ####
metadata_R <- subset(alpha, Area_Microbiome == "Root")

#set colors and y-axis maximum
col_Site <- c("#8DD3C7","#FFFFB3","#FFFFB3")
y_max <- max(metadata_R$shannon_entropy) + 1


#Significance test
metadata_R_sig_test <- metadata_R[,c("SampleID","Site","shannon_entropy")]
colnames(metadata_R_sig_test)[2] <- "group" 
metadata_R_sig_test$group <- as.factor(metadata_R_sig_test$group)

res = KwWlx(data = metadata_R_sig_test, i= 3)
res[[1]]
res[[2]]


#sig label
R_sig <- data.frame(
  niche = "Root",
  group1 = "Village",
  group2 = "Factory",
  y.position = y_max,
  p.adj = res[[2]]$p.adj,
  p.adj.signif = paste0("P = ", res[[2]]$p.adj),
  xmin = c(1),
  xmax = c(2)
)

#visualization
metadata_R_SE <- summarySE(metadata_R, measurevar= "shannon_entropy", groupvars= "Site")
metadata_R_SE$Site <- factor(metadata_R_SE$Site, levels = c("Village","Factory"))
metadata_R_SE$niche <- "Root"



### Soil       ####
metadata_S <- subset(alpha, Area_Microbiome == "Soil")

#Significance test
metadata_S_sig_test <- metadata_S[,c("SampleID","Site","shannon_entropy")]
colnames(metadata_S_sig_test)[2] <- "group" 
metadata_S_sig_test$group <- as.factor(metadata_S_sig_test$group)

res = KwWlx(data = metadata_S_sig_test, i= 3)
res[[1]]
res[[2]]


#sig label
S_sig <- data.frame(
  niche = "Soil",
  group1 = "Village",
  group2 = "Factory",
  y.position = y_max,
  p.adj = res[[2]]$p.adj,
  p.adj.signif = paste0("P = ", res[[2]]$p.adj),
  xmin = c(1),
  xmax = c(2)
)


#visualization
metadata_S_SE <- summarySE(metadata_S, measurevar= "shannon_entropy", groupvars= "Site")
metadata_S_SE$Site <- factor(metadata_S_SE$Site,levels = c("Village","Factory"))
metadata_S_SE$niche <- "Soil"


#
Site_SE <- rbind(metadata_R_SE, metadata_S_SE)
Site_sig <- rbind(R_sig, S_sig)

Site_sig$p.adj.signif <- ifelse(Site_sig$p.adj < 0.05, Site_sig$p.adj.signif, 'ns')

Site_shannon <- rbind(metadata_R, metadata_S)
colnames(Site_shannon)[4] <- "niche"

Site_SE$Site <- factor(Site_SE$Site,levels = c("Village","Factory"))
Site_SE$niche <- factor(Site_SE$niche,levels = c("Root","Soil"))
Site_shannon$Site <- factor(Site_shannon$Site,levels = c("Village","Factory"))
Site_shannon$niche <- factor(Site_shannon$niche, levels = c("Root","Soil"))


#
p_Site <- ggplot(Site_SE, aes(Site, shannon_entropy, fill=Site)) + 
  geom_bar(stat="identity", width = 0.8,color="black") +# 添加color="black"来设置柱子的边框颜色
  facet_grid( ~ niche, scales = "free", space = "free", switch = "x") +
  scale_fill_manual(values=col_Site) +
  xlab("") + ylab("Shannon index ") + 
  #ggtitle("") +
  theme_bw() +
  theme(legend.position="right") + 
  theme(strip.background = element_rect( fill = "white",color = "black",  linewidth = 1,linetype = "solid" ),  
        strip.text = element_text(color = "black", size = 12),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(), 
        axis.line = element_line(color = "black", linewidth = 0.7),
        legend.title = element_text(size=12,color="black"),
        legend.text= element_text(size=12,color="black"),
        axis.text.x = element_text(size=12,color="black", hjust = 0.5, vjust = 1),
        axis.text.y = element_text(size=12,color="black"),
        plot.margin = margin(t = 40, r = 10, b = 20, l = 20),
        axis.title= element_text(size=12,color="black"),
        plot.title =element_text(size=14,color="black", hjust = 0.5, vjust = 1))+
  geom_jitter(data = Site_shannon, aes(Site, shannon_entropy), size=3, width = 0.2, alpha = 0.6) +
  geom_errorbar(data = Site_SE, aes(ymin=shannon_entropy-se, ymax=shannon_entropy+se), width=0.2, colour = "black") +
  scale_y_continuous(limits = c(NA, 11), breaks = seq(0, 10.0, by = 2.5)) # 设置 y 轴的最大值

p1 <-p_Site +  
  stat_pvalue_manual(
    Site_sig,                                          #数据框
    label = "{p.adj.signif}",              #选择的标签，"{p.adj}{p.adj.signif}" 能选择两列作为标签
    y.position = "y.position",        #标签的位置，需要依据成图的数据计算出来。字符刻度  中  每个字符代表一个单位，如 A 就是1， B就是2等等
    xmin = "xmin",                            #标签下的比较线的起始位置，与xmax搭配使用，控制比较线的长度
    xmax = "xmax",                           #
    tip.length = 0.01,                      #比较线的粗细
    size = 4,
    label.fontface = "italic" ,
    coord.flip = F                     #T：翻转标签， F：不翻转，默认不翻转
  ) 

p1
ggsave(filename = "../../figures/Fig.2/Site_shannon.pdf", plot = p1 , width = 5.8, height = 5)





#########
####        Root and Soil |  local Environment(indoor vs outdoor)
##           Root  Village
metadata_R_V <- subset(metadata_R, Site == "Village")
metadata_R_V_SE <- summarySE(metadata_R_V, measurevar= "shannon_entropy", groupvars= "Environment")
metadata_R_V_SE$niche <- "Root"

### Village Environment Significance test
metadata_R_V_sig_test <- metadata_R_V[,c("SampleID","Environment","shannon_entropy")]
colnames(metadata_R_V_sig_test)[2] <- "group" 

res = KwWlx(data = metadata_R_V_sig_test, i= 3)
res[[1]]
res[[2]]


#sig label
RV_sig <- data.frame(
  niche = "Root",
  group1 = "indoor",
  group2 = "outdoor",
  y.position = y_max,
  p.adj = res[[2]]$p.adj,
  p.adj.signif = paste0("P = ", res[[2]]$p.adj),
  xmin = c(1),
  xmax = c(2)
)


##           Soil  Village
metadata_S_V <- subset(metadata_S, Site == "Village")
metadata_S_V_SE <- summarySE(metadata_S_V, measurevar= "shannon_entropy", groupvars= "Environment")
metadata_S_V_SE$niche <- "Soil"

### Village Environment Significance test
metadata_S_V_sig_test <- metadata_S_V[,c("SampleID","Environment","shannon_entropy")]
colnames(metadata_S_V_sig_test)[2] <- "group" 

res = KwWlx(data = metadata_S_V_sig_test, i= 3)
res[[1]]
res[[2]]


#sig label
SV_sig <- data.frame(
  niche = "Soil",
  group1 = "indoor",
  group2 = "outdoor",
  y.position = y_max,
  p.adj = res[[2]]$p.adj,
  p.adj.signif = paste0("P = ", res[[2]]$p.adj),
  xmin = c(1),
  xmax = c(2)
)



# Root and Soil | Village
Village_SE <- rbind(metadata_R_V_SE, metadata_S_V_SE)
Village_sig <- rbind(RV_sig, SV_sig)

Village_sig$p.adj.signif <- ifelse(Village_sig$p.adj < 0.05, Village_sig$p.adj.signif, 'ns')

Village_shannon <- rbind(metadata_R_V, metadata_S_V)
colnames(Village_shannon)[4] <- "niche"

Village_SE$Environment <- factor(Village_SE$Environment,levels = c("indoor","outdoor"))
Village_SE$niche <- factor(Village_SE$niche,levels = c("Root","Soil"))
Village_shannon$Environment <- factor(Village_shannon$Environment,levels = c("indoor","outdoor"))
Village_shannon$niche <- factor(Village_shannon$niche, levels = c("Root","Soil"))


# Village(indoor vs outdoor) | Root and Soil

col_environemt <- c("indoor" = "#99CC00","outdoor" ="#FFCC00")

pv_environment <- ggplot(Village_SE, aes(Environment, shannon_entropy, fill=Environment)) + 
  geom_bar(stat="identity", width = 0.8,color="black") +
  facet_grid( ~ niche, scales = "free", space = "free", switch = "x") +
  scale_fill_manual(values=col_environemt) +
  xlab("") + ylab("Shannon index ") + 
  #ggtitle("") +
  theme_bw() +
  theme(legend.position="right") + 
  theme(strip.background = element_rect( fill = "white",color = "black",  linewidth = 1,linetype = "solid" ),  
        strip.text = element_text(color = "black", size = 12),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(), 
        axis.line = element_line(color = "black", linewidth = 0.7),
        legend.title = element_text(size=12,color="black"),
        legend.text= element_text(size=12,color="black"),
        axis.text.x = element_text(size=12,color="black", hjust = 0.5, vjust = 1),
        axis.text.y = element_text(size=12,color="black"),
        plot.margin = margin(t = 40, r = 10, b = 20, l = 20),
        axis.title= element_text(size=12,color="black"),
        plot.title =element_text(size=14,color="black", hjust = 0.5, vjust = 1))+
  geom_jitter(data = Village_shannon, aes(Environment, shannon_entropy), size=3, width = 0.2, alpha = 0.6) +
  geom_errorbar(data = Village_SE, aes(ymin=shannon_entropy-se, ymax=shannon_entropy+se), width=0.2, colour = "black") +
  scale_y_continuous(limits = c(NA, 11), breaks = seq(0, 10.0, by = 2.5)) 

p2 <-pv_environment +  
  stat_pvalue_manual(
    Village_sig,                                        
    label = "{p.adj.signif}",              
    y.position = "y.position",        
    xmin = "xmin",                           
    xmax = "xmax",                           
    tip.length = 0.01,                      
    size = 4,
    label.fontface = "italic" ,
    coord.flip = F                     
  ) 

p2
ggsave(filename = "../../figures/Fig.2/Village_environment_shannon.pdf", plot = p2 , width = 5.8, height = 5)




#########
####        Root and Soil |  local Environment(indoor vs outdoor)
##          Root  Factory
metadata_R_H <- subset(metadata_R, Site == "Factory")
metadata_R_H_SE <- summarySE(metadata_R_H, measurevar= "shannon_entropy", groupvars= "Environment")
metadata_R_H_SE$niche <- "Root"

# Significance test    Factory Environment
metadata_R_H_sig_test <- metadata_R_H[,c("SampleID","Environment","shannon_entropy")]
colnames(metadata_R_H_sig_test)[2] <- "group" 

res = KwWlx(data = metadata_R_H_sig_test, i= 3)
res[[1]]
res[[2]]


#sig label
RF_sig <- data.frame(
  niche = "Root",
  group1 = "indoor",
  group2 = "outdoor",
  y.position = y_max,
  p.adj = res[[2]]$p.adj,
  p.adj.signif = paste0("P = ", res[[2]]$p.adj),
  xmin = c(1),
  xmax = c(2)
)


##           Soil  Factory
metadata_S_H <- subset(metadata_S, Site == "Factory")
metadata_S_H_SE <- summarySE(metadata_S_H, measurevar= "shannon_entropy", groupvars= "Environment")
metadata_S_H_SE$niche <- "Soil"


# Significance test    Factory Environment
metadata_S_H_sig_test <- metadata_S_H[,c("SampleID","Environment","shannon_entropy")]
colnames(metadata_S_H_sig_test)[2] <- "group" 

res = KwWlx(data = metadata_S_H_sig_test, i= 3)
res[[1]]
res[[2]]


#sig label
SF_sig <- data.frame(
  niche = "Soil",
  group1 = "indoor",
  group2 = "outdoor",
  y.position = y_max,
  p.adj = res[[2]]$p.adj,
  p.adj.signif = paste0("P = ", res[[2]]$p.adj),
  xmin = c(1),
  xmax = c(2)
)


# Root and Soil | Factory
Factory_SE <- rbind(metadata_R_H_SE, metadata_S_H_SE)
Factory_sig <- rbind(RF_sig, SF_sig)

Factory_sig$p.adj.signif <- ifelse(Factory_sig$p.adj < 0.05, Factory_sig$p.adj.signif, 'ns')

Factory_shannon <- rbind(metadata_R_H, metadata_S_H)
colnames(Factory_shannon)[4] <- "niche"

Factory_SE$Environment <- factor(Factory_SE$Environment,levels = c("indoor","outdoor"))
Factory_SE$niche <- factor(Factory_SE$niche,levels = c("Root","Soil"))
Factory_shannon$Environment <- factor(Factory_shannon$Environment,levels = c("indoor","outdoor"))
Factory_shannon$niche <- factor(Factory_shannon$niche, levels = c("Root","Soil"))


# Village(indoor vs outdoor) | Root and Soil

pf_environment <- ggplot(Factory_SE, aes(Environment, shannon_entropy, fill=Environment)) + 
  geom_bar(stat="identity", width = 0.8,color="black") +
  facet_grid( ~ niche, scales = "free", space = "free", switch = "x") +
  scale_fill_manual(values=col_environemt) +
  xlab("") + ylab("Shannon index ") + 
  #ggtitle("") +
  theme_bw() +
  theme(legend.position="right") + 
  theme(strip.background = element_rect( fill = "white",color = "black",  linewidth = 1,linetype = "solid" ),  
        strip.text = element_text(color = "black", size = 12),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(), 
        axis.line = element_line(color = "black", linewidth = 0.7),
        legend.title = element_text(size=12,color="black"),
        legend.text= element_text(size=12,color="black"),
        axis.text.x = element_text(size=12,color="black", hjust = 0.5, vjust = 1),
        axis.text.y = element_text(size=12,color="black"),
        plot.margin = margin(t = 40, r = 10, b = 20, l = 20),
        axis.title= element_text(size=12,color="black"),
        plot.title =element_text(size=14,color="black", hjust = 0.5, vjust = 1))+
  geom_jitter(data = Factory_shannon, aes(Environment, shannon_entropy), size=3, width = 0.2, alpha = 0.6) +
  geom_errorbar(data = Factory_SE, aes(ymin=shannon_entropy-se, ymax=shannon_entropy+se), width=0.2, colour = "black") +
  scale_y_continuous(limits = c(NA, 11), breaks = seq(0, 10.0, by = 2.5)) 

p3 <-pf_environment +  
  stat_pvalue_manual(
    Factory_sig,                                       
    label = "{p.adj.signif}",             
    y.position = "y.position",      
    xmin = "xmin",                           
    xmax = "xmax",                          
    tip.length = 0.01,                      
    size = 4,
    label.fontface = "italic" ,
    coord.flip = F                     
  ) 

p3
ggsave(filename = "../../figures/Fig.2/Factory_environment_shannon.pdf", plot = p3 , width = 5.8, height = 5)









