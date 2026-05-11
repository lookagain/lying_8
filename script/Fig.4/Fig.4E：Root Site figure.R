library(ggplot2)
library(ggpubr)
library(dplyr)
library(patchwork)


# Set Work Path
pwd <- dirname(rstudioapi::getSourceEditorContext()$path)
setwd(pwd)


#读入数据
R_abund <- read.csv("../../data/rf_family_R_top30_res-abund.csv")
R_sig <- read.csv("../../data/rf_family_R_top30_res-abund-wilcox-sig.csv")


##随机森林图
rf_imp <- read.csv("randomForest_importance_family_Root_Site.csv",row.names = 1)
rf_imp <- rf_imp[order(rf_imp$MeanDecreaseAccuracy, decreasing = T),]
rf_imp_top30 <- rf_imp[1:30, ]

rf_imp_top30$phylum <- gsub(".*?p__(.*?);.*", "\\1", rownames(rf_imp_top30))
rf_imp_top30$family <- gsub(".*?f__(.*)", "\\1", rownames(rf_imp_top30))


# sort
rf_imp_top30 <- rf_imp_top30[order(rf_imp_top30$MeanDecreaseAccuracy), ]
family_order <- rf_imp_top30$family
rf_imp_top30$family <- factor(rf_imp_top30$family, levels =family_order )


#set colors
phylum_colors <- c("Proteobacteria" = "#E41A1C",
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

rf_imp_top30$color_R <- phylum_colors[rf_imp_top30$phylum]

p_rf <- ggplot(rf_imp_top30, aes(family, MeanDecreaseAccuracy, fill=phylum, colour = color_R)) + 
  geom_bar(stat="identity",color="black") + 
  geom_errorbar(aes(ymin=MeanDecreaseAccuracy-SD, ymax=MeanDecreaseAccuracy+SD), width = 0.4,colour = "black") + 
  scale_fill_manual(values = phylum_colors) +
  coord_flip() +
  xlab("")+ylab("Mean of decreased accuracy")+
  theme_bw() +
  theme(legend.position="right") + 
  theme(panel.grid = element_blank(),
        panel.border = element_blank(), 
        axis.line = element_line(color = "grey60", linetype = "solid", lineend = "square", linewidth = 0.7),
        legend.title = element_text(size=10,color="black", face = "bold"),
        legend.text= element_text(size=7,color="black"),
        axis.text.x = element_text(size=8,color="black", hjust = 0.5),
        axis.text.y = element_text(size=8,color="black", face = "italic"),
        axis.title= element_text(size=12,color="black"),
        axis.title.x = element_text(hjust = 0.5,vjust = 0),
        plot.margin = margin(10, -5, 10, -15)) +
  guides(colour = "none",
         fill = guide_legend(title = "Phylum",
           keyheight = unit(4, "mm"),  # 设置图例填充项的高度
           keywidth = unit(4, "mm")    # 设置图例填充项的宽度
         ))  

p_rf  

  
###交叉验证图
rf_cv.mean <- read.table("rf_family_Root_Site_cv.mean.txt",sep = "\t",header = 1)

rf_cv <- ggplot(rf_cv.mean, aes(Feature.num, cv.mean.error)) +
  geom_smooth(color = 'black', linewidth = 0.8, method = 'loess', span = 0.1, se = F) +  
  geom_vline(xintercept = 30, linetype = "dashed", color = "black", size = 0.8) +  
  theme(panel.grid = element_blank(),
        panel.background = element_rect(color = 'black', fill = 'transparent'),
        axis.title= element_text(size=8,color="black"),
        axis.text = element_text(size=6,color="black", hjust = 0.5),) +
  scale_x_continuous(breaks = c(0, 30, 100, 200, 300, 400, 500))  +  
  scale_y_continuous(breaks = c(0, 0.05, 0.10, 0.15, 0.20, 0.25, 0.30, 0.35, 0.40), 
                     limits = c(0, 0.40)) +  
  labs(title = '', x = 'Number of feature', y = 'Cross-validation error')

rf_cv  


#将成图转换为子图格式，作为插入图
sub_grob <- ggplotGrob(rf_cv)

xmax <- length(rf_imp_top30$family)/2
ymax <- rf_imp_top30[30,"MeanDecreaseAccuracy"] + rf_imp_top30[30,"SD"]
ymin <- ymax/2 - ymax/18

p <- p_rf + annotation_custom(
  grob = sub_grob,
  xmin = 1, xmax = xmax,  # 控制副图的横坐标范围
  ymin = ymin, ymax = ymax     # 控制副图的纵坐标范围
)

p


#构建手动添加显著性需要的表格
R_abund <- R_abund %>% mutate(y.position = Mean + SE)
R_position <- R_abund[,c("Taxa","y.position")]
R_position <- R_position  %>%
  group_by(Taxa) %>%
  summarise(y.position = max(y.position)) 

R_sig <- R_sig[,c("Taxa","Significance")]
R_sig$group1 <- "Village"
R_sig$group2 <- "Factory"
R_sig <- merge(R_sig,R_position,"Taxa",all = F)
R_sig <- R_sig[order(match(R_sig$Taxa, family_order)), ]

R_abund <- R_abund[order(match(R_abund$Taxa, family_order)), ]

R_sig$num <- c(1:30)
R_sig <- R_sig %>% 
  mutate(min = num - 0.25, max = num + 0.25)

R_sig$Site <- ""

R_sig <- R_sig %>%
  mutate(y.position = y.position + 0.02)


#随机森林的生物标记相对丰度图
color_values = c("#1B9E77","#D95F02","#D95F02")
R_abund$Site <- factor(R_abund$Site, levels = c("Village","Factory"))
R_abund$Taxa <- factor(R_abund$Taxa, levels = family_order)


p2 <- ggplot(R_abund, aes(x = Taxa, y = Mean, color = Site, fill = Site,group = Site)) + 
  scale_color_manual(values = color_values) + 
  scale_fill_manual(values = color_values) + 
  geom_bar(stat = "identity", position = position_dodge()) + 
  geom_errorbar(aes(ymin=Mean-SE, ymax=Mean+SE), width = 0.4,linewidth = 0.5,colour = "black", position = position_dodge(width = 0.8)) + 
  theme_bw() + 
  ylab("Relative abundance") + 
  xlab("") +
  theme(
    legend.title = element_text(size=10,color="black", face = "bold"),
    legend.text= element_text(size=7,color="black"),
    axis.title = element_text(size = 12),
    axis.title.x = element_text(hjust = 0.5,vjust = 0),
    axis.line.x = element_line(color = "grey60", linetype = "solid", lineend = "square", linewidth = 0.7),
    axis.text.x = element_text(size=8,color="black", hjust = 0.5),
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    panel.border = element_blank(), 
    panel.grid.minor.y = element_blank(),
    panel.grid.major.y = element_blank(),
    panel.grid.minor.x = element_blank(), 
    panel.grid.major.x = element_blank(),
    plot.margin = margin(10, -5, 10, -15)
  )  + 
  scale_y_continuous(limits = c(0, 0.28),expand = c(0, 0)) +
  coord_flip() +
  stat_pvalue_manual(
    R_sig,
    label = "Significance",
    y.position = "y.position",
    xmin = "min",
    xmax = "max",
    tip.length = 0.01,
    coord.flip = TRUE,
    size = 2.8,
  ) +
  guides(color = "none",
         fill = guide_legend(
           keyheight = unit(4, "mm"),  # 设置图例填充项的高度
           keywidth = unit(4, "mm")    # 设置图例填充项的宽度
         )) 


p2


p3 <- p + p2 + 
  plot_layout(guides = 'collect')  #合并图例 

p3

ggsave("../../figures/Fig.4/Root_site_combine_figure.pdf", p3, width = 12, height =5.5)



