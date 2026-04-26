#
library(UpSetR)
library(vegan)
library(ComplexHeatmap)
library(circlize)
library(dendextend)


#
pwd <- dirname(rstudioapi::getSourceEditorContext()$path)
setwd(pwd)

##
metal_fileList <- list.files(path = "./", pattern = "Metal_.*sleuth_significant_gene_level.csv", full.names = TRUE)
metal_list <- lapply(metal_fileList, read.csv)

names(metal_list) <- sapply(metal_fileList, basename)
names(metal_list) <- gsub(" sleuth.*", "", names(metal_list))
names(metal_list) <- gsub("Metal_ ", "", names(metal_list))
names(metal_list) <- gsub(" ", "_", names(metal_list))


#
rice_tpm <- read.csv(file = "metal_NIL_Syncom_DGE_sleuth_norm_tpm_genome.csv")
name <- colnames(rice_tpm)[grepl("CK_NIL_|CK_SynCom_|Cd_NIL_|Cd_SynCom_|Mn_NIL_|Mn_SynCom_",colnames(rice_tpm))]


## all DEGs upset
metal_0.05degs_listinput <- list(CK_NIL_vs_Mn_NIL    = subset(metal_list$CK_NIL_vs_Mn_NIL,    qval < 0.05 & abs(log2fc) > 2)$target_id,
                                 CK_NIL_vs_Mn_SynCom = subset(metal_list$CK_NIL_vs_Mn_SynCom, qval < 0.05 & abs(log2fc) > 2)$target_id,
                                 Mn_NIL_vs_Mn_SynCom = subset(metal_list$Mn_NIL_vs_Mn_SynCom, qval < 0.05 & abs(log2fc) > 2)$target_id)

all_id <- unlist(metal_0.05degs_listinput, use.names = F)
all_id <- all_id[!duplicated(all_id)]

p2 <- upset(fromList(metal_0.05degs_listinput), nsets = 3, 
            queries = list(
              list(
                query = elements,  #选择"属于特定集合的元素"
                params = list("CK_NIL_vs_Mn_NIL"),  # 指定目标集合
                color = "#3399CC",  #标记颜色
                active = TRUE
              ),
              list(query = elements, params = list("Mn_NIL_vs_Mn_SynCom"), color = "#FF6B6B", active = T),
              list(query = elements, params = list("CK_NIL_vs_Mn_SynCom"), color = "#CC9933", active = T),
              list(query = intersects, params = list("CK_NIL_vs_Mn_NIL","Mn_NIL_vs_Mn_SynCom"), color = "#33CC00", active = T),
              list(query = intersects, params = list("CK_NIL_vs_Mn_SynCom","CK_NIL_vs_Mn_NIL"), color = "#999999", active = T)
            ),
            order.by = "degree", decreasing = F,
            keep.order = T)  

p2

pdf("Metal_logfc_2_qval_0.05_degs_upset.pdf", width = 7, height = 5)
print(p2)
dev.off()


##                                        DEGs: abs(log2fc) >2 & p <0.05  
DEGs <- list(CK_NIL_vs_CK_SynCom = subset(metal_list$CK_NIL_vs_CK_SynCom, qval < 0.05 & abs(log2fc) > 2)$target_id,
             Cd_NIL_vs_Cd_SynCom = subset(metal_list$Cd_NIL_vs_Cd_SynCom, qval < 0.05 & abs(log2fc) > 2)$target_id,
             CK_NIL_vs_Cd_NIL    = subset(metal_list$CK_NIL_vs_Cd_NIL,    qval < 0.05 & abs(log2fc) > 2)$target_id,
             CK_NIL_vs_Mn_NIL    = subset(metal_list$CK_NIL_vs_Mn_NIL,    qval < 0.05 & abs(log2fc) > 2)$target_id,
             Mn_NIL_vs_Mn_SynCom = subset(metal_list$Mn_NIL_vs_Mn_SynCom, qval < 0.05 & abs(log2fc) > 2)$target_id)
                                 
All_id_log2fc <- unlist(DEGs, use.names = F)
All_id_log2fc <- All_id_log2fc[!duplicated(All_id_log2fc)]



## all DEGs
metal_tpm <- rice_tpm[rice_tpm$Transcript_ID %in% All_id_log2fc, ]
row.names(metal_tpm) <- paste(metal_tpm$Transcript_ID)
metal_tpm_heatmap <- metal_tpm[, grepl("CK_|Cd_|Mn_", colnames(metal_tpm)), drop = FALSE]

metal_tpm_heatmap_scale <- decostand(log(metal_tpm_heatmap + 1, base = 2), "standardize", MARGIN = 1)


#计算聚类簇
metal_tpm_heatmap_scale.cascade <- cascadeKM(metal_tpm_heatmap_scale, 
                                                      inf.gr = 2, sup.gr = 15,
                                                      iter = 100, criterion = "ssi")
plot(metal_tpm_heatmap_scale.cascade, sortg = T)

#####4组
dend = hclust(dist(metal_tpm_heatmap_scale))
dend = color_branches(dend, k = 4)

colscale <- colorRamp2(c(-2, 0, 2), c("#00CCFF", "white", "#FF0000")) 
metal_tpm_heatmap_hclust <- Heatmap(metal_tpm_heatmap_scale, name = "Z-score (log2(TPM+1))",  col = colscale,
                                cluster_columns = F, show_row_names = F,
                                row_names_gp = gpar(fontsize = 8), column_names_gp = gpar(fontsize = 8),
                                clustering_distance_rows = "euclidean",
                                clustering_method_rows = "complete",
                                row_dend_width = unit(2, "cm"),
                                cluster_rows = dend, split = 4)

metal_tpm_heatmap_hclust
dev.copy(pdf, "metal_logfc_2_qval_0.05_degs_norm_tpm_heatmap_scale.pdf", width = 7, height = 8)
dev.off()



##Clusters
metal_tpm_heatmap_hclust <- draw(metal_tpm_heatmap_hclust)
row_dend(metal_tpm_heatmap_hclust)

cluster1 <- row.names(metal_tpm_heatmap_scale)[unlist(row_order(metal_tpm_heatmap_hclust)[1])]
cluster2 <- row.names(metal_tpm_heatmap_scale)[unlist(row_order(metal_tpm_heatmap_hclust)[2])]
cluster3 <- row.names(metal_tpm_heatmap_scale)[unlist(row_order(metal_tpm_heatmap_hclust)[3])]
cluster4 <- row.names(metal_tpm_heatmap_scale)[unlist(row_order(metal_tpm_heatmap_hclust)[4])]

tpm_cluster1_scale <- metal_tpm_heatmap_scale[rownames(metal_tpm_heatmap_scale) %in% cluster1,]
tpm_cluster2_scale <- metal_tpm_heatmap_scale[rownames(metal_tpm_heatmap_scale) %in% cluster2,]
tpm_cluster3_scale <- metal_tpm_heatmap_scale[rownames(metal_tpm_heatmap_scale) %in% cluster3,]
tpm_cluster4_scale <- metal_tpm_heatmap_scale[rownames(metal_tpm_heatmap_scale) %in% cluster4,]

write.csv(tpm_cluster1_scale,"metal_log2fc_gt2_q0.05_tpm_scale_cluster1.csv")
write.csv(tpm_cluster2_scale,"metal_log2fc_gt2_q0.05_tpm_scale_cluster2.csv")
write.csv(tpm_cluster3_scale,"metal_log2fc_gt2_q0.05_tpm_scale_cluster3.csv")
write.csv(tpm_cluster4_scale,"metal_log2fc_gt2_q0.05_tpm_scale_cluster4.csv")





