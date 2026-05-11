#
library(dplyr)
library(reshape2)
library(ggplot2)

# Set Work Path
pwd <- dirname(rstudioapi::getSourceEditorContext()$path)
setwd(pwd)


#
metadata <- read.table("../../data/Metadata_add_Place.txt", sep = "\t",row.names =1 ,header = T, na.strings = "")
sub1 <- metadata[metadata$Area_Microbiome =="Soil",]

phylum <- read.csv("../../data/clean_level2_relative_abundance.csv",check.names = F,row.names = 1)
rownames(phylum) <- gsub(".*p__(.*)", "\\1", rownames(phylum))


##
data1 <- phylum[,colnames(phylum) %in% rownames(sub1)]
data1$rowMeans <- rowMeans(data1, na.rm = TRUE)
data1 <- data1[order(data1$rowMeans,decreasing = T),]
data1$phylum <- rownames(data1)

top10 <- c(
  "Proteobacteria", 
  "Firmicutes", 
  "Bacteroidota", 
  "Actinobacteriota", 
  "Acidobacteriota", 
  "Chloroflexi", 
  "Verrucomicrobiota", 
  "Myxococcota", 
  "Desulfobacterota", 
  "Nitrospirota"
)

data1$phylum <- ifelse(rownames(data1) %in% top10, rownames(data1),"Others")

data2 <- data1 %>%
  group_by(phylum) %>%
  summarise(across(.cols = everything(), ~sum(., na.rm = TRUE)))

data2$phylum <- factor(data2$phylum,levels = c(top10,"Others"))
data2 <- as.data.frame(data2) %>% arrange(phylum)

rownames(data2) <- data2$phylum
data2 <- data2[,colnames(data2) %in% rownames(sub1)]
data2 <- as.data.frame(t(data2))
data2$SampleID <- rownames(data2)

sub2 <- sub1 
sub2$SampleID <- rownames(sub2)
sub2 <- sub2[,c("SampleID", "Place","Plant")]

data3 <- merge(data2,sub2,"SampleID")
rownames(data3) <- data3$SampleID
data3 <- data3[,-1]
abund <- data3 %>%
  group_by(Place, Plant) %>%
  summarise(across(.cols = everything(), ~mean(., na.rm = TRUE)))

rowSums(abund[3:12])
abund$Place <- factor(abund$Place,levels = c("Village indoor", "Village outdoor", "Factory indoor", "Factory outdoor" ))
abund <- abund %>% arrange(Place)
abund$num <- 1:nrow(abund)
sub3 <- abund[,c("num","Place","Plant")]

abund_melt = as.data.frame(melt(abund, id.vars = "num", measure.vars = c(top10,"Others")))
abund_melt <- merge(abund_melt,sub3,by = "num")
colnames(abund_melt) <- c("num","phylum","abundance","Place","Plant")
abund_melt$phylum <- factor(abund_melt$phylum,levels = c(top10,"Others")) 
abund_melt$Place <- factor(abund_melt$Place,levels = c("Village indoor", "Village outdoor", "Factory indoor", "Factory outdoor" ))

abund_melt$Plant <- factor(abund_melt$Plant,levels = unique(abund_melt$Plant))


#colors
phylum_color <- c("Proteobacteria" = "#E41A1C",
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
                  "Patescibacteria" = "#FFCCFF")

color1 <- phylum_color[c(top10,"Others")]

abund_melt$num <- factor(abund_melt$num)
num_plant_map <- abund[,c("num","Plant")]
plant_labels <- num_plant_map$Plant[num_plant_map$num]

p1 <- ggplot(abund_melt, aes(x = Plant, y = 100 * abundance, fill = phylum)) +
  geom_col(position = 'stack', width = 0.8) +
  scale_fill_manual(values = color1, name = "Phylum") +
  facet_grid( ~ Place, scales = "free", space = "free", switch = "x") +
  theme(
    strip.background = element_rect( fill = "white",color = "black",  linewidth = 1,linetype = "solid" ),  
    strip.text = element_text(color = "black", size = 10),  
    panel.spacing = unit(0, "lines"),  
    panel.grid = element_blank(),  
    panel.background = element_rect(fill = "white", color = NA),  
    axis.text.x = element_text(hjust = 1, vjust = 1, color = "black", size = 8, angle = 45),  
    axis.text.y = element_text(hjust = 0.5, color = "black", size = 8), 
    axis.title.y = element_text(size = 12,vjust =3),
    axis.ticks.x = element_blank(), 
    axis.line = element_line(color = "black", linewidth = 0.8),
    legend.title = element_text(size = 12),
    legend.text = element_text( size = 8),
    plot.margin = unit(c(1, 1, 1, 1), "cm")
  ) +
  scale_y_continuous(breaks = c(0, 25, 50, 75, 100), expand = c(0, 0),
                     limits = c(0, 105)) + 
  xlab("") +
  ylab("Relative abundance (%)")

p1
ggsave("../../figures/Fig.4/Soil_Plant_Phylum_composition.pdf", p1, width = 10, height =7)



