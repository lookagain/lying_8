library(ggplot2)
library(ComplexHeatmap)
library(circlize)

# library(stringr)

pwd <- dirname(rstudioapi::getSourceEditorContext()$path)
setwd(pwd)


library(clusterProfiler)
term2gene <- read.csv(file = "../../data/RNA-seq_Locus_ID_GO_term2gene.csv")
term2name <- read.csv(file = "../../data/RNA-seq_Locus_ID_GO_term2name.csv")
term2gene_P <- read.csv(file = "../../data/RNA-seq_Locus_ID_GO_term2gene_P.csv")
term2name_P <- read.csv(file = "../../data/RNA-seq_Locus_ID_GO_term2name_P.csv")


# 
metal_fileList <- list.files(path = "./", pattern = "^metal_log2fc_.*.csv", full.names = TRUE)
metal_list <- lapply(metal_fileList, read.csv)

names(metal_list) <- sapply(metal_fileList, basename)
names(metal_list) <- gsub("metal_log2fc_gt2_q0.05_tpm_scale_", "", names(metal_list))
names(metal_list) <- gsub(".csv", "", names(metal_list))


###BP
BP_GO_list <- list(          cluster1 = subset(metal_list$cluster1, X %in% term2gene_P$gene )$X,
                             cluster2 = subset(metal_list$cluster2, X %in% term2gene_P$gene )$X,
                             cluster3 = subset(metal_list$cluster3, X %in% term2gene_P$gene )$X,
                             cluster4 = subset(metal_list$cluster4, X %in% term2gene_P$gene )$X)

BP_GO_enrich <- compareCluster(BP_GO_list, fun='enricher',TERM2GENE=term2gene ,TERM2NAME=term2name ,pvalueCutoff = 0.05, pAdjustMethod = "BH", qvalueCutoff = 0.05)
BP <- as.data.frame(BP_GO_enrich)
BP <- BP[BP$Description %in% term2name_P$name,]


####
BP_heatmap <- BP[,c("Cluster" ,"ID" ,"Description" ,"qvalue")]
BP_heatmap$name <- paste(BP_heatmap$ID, BP_heatmap$Description, sep = " | ")
BP_heatmap <- BP_heatmap[, c("Cluster" ,"qvalue" ,"name")]
BP_heatmap$log10 <- -log10(BP_heatmap$qvalue)
BP_heatmap <- BP_heatmap[,c(1,4,3)]

cluster1 <- BP_heatmap[BP_heatmap$Cluster =="cluster1",-1]
colnames(cluster1)[1] <- "cluster1"
cluster2 <- BP_heatmap[BP_heatmap$Cluster =="cluster2",-1]
colnames(cluster2)[1] <- "cluster2"
cluster3 <- BP_heatmap[BP_heatmap$Cluster =="cluster3",-1]
colnames(cluster3)[1] <- "cluster3"
cluster4 <- BP_heatmap[BP_heatmap$Cluster =="cluster4",-1]
colnames(cluster4)[1] <- "cluster4"

combine1 <- merge(cluster1 ,cluster2 ,"name",all = T)
combine2 <- merge(cluster3 , cluster4 ,"name",all = T)
combine <- merge(combine1 , combine2 ,"name",all = T)
combine[is.na(combine)] <- 0

rownames(combine) <- combine$name
combine <- combine[,-1]


#排序
sorted_df <- combine[order(
  #第一优先级
  combine$cluster4 == 0,  
  -combine$cluster4,      
  #第二优先级
  combine$cluster3 == 0,  
  -combine$cluster3,      
  #第三优先级
  combine$cluster2 == 0,
  -combine$cluster2,
  #第四优先级
  combine$cluster1 == 0,
  -combine$cluster1
), ]

print(sorted_df)

which(rownames(sorted_df) == "GO:1990961 | xenobiotic detoxification by transmembrane export across the plasma membrane")
nrow(sorted_df)
sorted_df <- sorted_df[c(1:43,45:48,44),]


#
go_order <- Heatmap(
  sorted_df, 
  name = "-log10(qvalue)",
  col = colorRamp2(c(0, 1.3, 9), c("white", "#00CCFF", "#FF0000")),
  cluster_columns = FALSE, 
  cluster_rows = FALSE,
  rect_gp = gpar(col = "#999999", lwd = 1.5),
  show_row_names = TRUE, 
  row_names_gp = gpar(fontsize = 8),
  column_names_gp = gpar(fontsize = 8),
  column_names_rot = 30,
  column_names_centered = TRUE,
  column_names_side = "bottom",
  row_dend_width = unit(4, "cm")
)
go_order

dev.copy(pdf, "metal_4Clusters_logfc_2_qval_0.05_enriched_go_heatmap.pdf", width = 5.5, height = 7.5)
dev.off()

