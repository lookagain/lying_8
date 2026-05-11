#--导入所需R包#-------
library(psych)
library(reshape2)
library(Hmisc)
library(ggplot2)
library(dplyr)
library(tidyr)
library(igraph)


# Set Work Path
pwd <- dirname(rstudioapi::getSourceEditorContext()$path)
setwd(pwd)

#读入文件
metadata = read.delim("../../data/Metadata.txt",row.names = 1)
metadata$SampleID <- rownames(metadata)
ASV  = read.csv("../../data/ASV Relative Abundance.csv", row.names=1,check.names = F)
taxonomy = read.csv("../../data/ASV_Taxonomy.csv", row.names=1,header = T,check.names = F)
taxonomy[is.na(taxonomy)] <- ""


##
metadata_sub <- metadata[metadata$Area_Microbiome =="Soil",]
Site <- metadata_sub[metadata_sub$Site =="Village",]

#
ASV_sub <- ASV[,colnames(ASV) %in% Site$SampleID]
otutab <- ASV_sub[apply(ASV_sub, 1, function(x) any(x > 0.005)), ]   
taxon <- taxonomy[rownames(taxonomy) %in% rownames(otutab),]


#
otutab <- as.matrix(otutab)
occor <- rcorr(t(otutab), type = 'spearman')
saveRDS(occor, "Soil Village rcorr list.rds")

occor_r = occor$r
occor_p = occor$P
occor_p.adj = p.adjust(occor_p, method = 'fdr')
occor_r[occor_p.adj >0.05|abs(occor_r<0.6)] = 0
diag(occor_r) <- 0
dim(occor_r)

#修剪
r.matrix <- occor_r
remove_rows_cols <- which(rowSums(r.matrix) == 0 & colSums(r.matrix) == 0)
r.matrix <- r.matrix[-remove_rows_cols, -remove_rows_cols]
dim(r.matrix )

igraph <- igraph::graph_from_adjacency_matrix(
  r.matrix,
  mode = "undirected",    
  weighted = TRUE,        
  diag = FALSE           
)


#
A = igraph

#去掉冗余的边
A<-simplify(igraph)

#提取权重
df_weight = E(A)$weight
which(df_weight <=0)

#加入ASV丰度信息
data <- rowSums(otutab)
data1=as.data.frame(data)
df_igraph_size = data1[V(A)$name,] 
df_igraph_size2 = log10(df_igraph_size)
V(A)$Abundance = df_igraph_size2

#加入物种信息,使用不同颜色表示,选择总丰度前6的门，其他归为others
data2=taxon
data2$Phylum <- gsub("^p__", "", data2$Phylum)
data2$ASV <- rownames(data2)
data2$rowSums <- data
result <- data2 %>% group_by(Phylum) %>% summarize(P_sum = sum(rowSums))
result <- result[order(result$P_sum,decreasing = T),]

phylum_top <- result$Phylum[1:6]
phylum_top
data2$Taxon <- data2$Phylum
data2$Taxon[!(data2$Phylum %in% phylum_top)] <- 'Others'
data2 = data2[,"Taxon",drop=F]
df_igraph_col = data2[V(A)$name,]
V(A)$Taxon = as.character(df_igraph_col)

#加一列color 指定节点颜色
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
data2$color <- phylum_color[data2$Taxon]
table(data2$color)
igraph_color = data2[V(A)$name,2]
V(A)$color = as.character(igraph_color)

#生成网络图的结点标签和degree属性；
V(A)$label <- V(A)$name
V(A)$degree <- degree(A)

#计算群体结构（cluster_fast_greedy）
c<- cluster_fast_greedy(A)

#使用默认颜色列表；
V(A)$Modularity <- c$membership

write_graph(A, "../../figures/Fig.3/Soil Village rcorr.graphml", format="graphml")

