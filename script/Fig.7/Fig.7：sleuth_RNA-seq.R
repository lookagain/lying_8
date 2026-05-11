# install.packages("BiocManager")
# install.packages("devtools")
# chooseCRANmirror()  
# BiocManager::install("pachterlab/sleuth")
library('sleuth')

#
pwd <- dirname(rstudioapi::getSourceEditorContext()$path)
setwd(pwd)

base_dir <- "../../data/kallisto_quant"
sample_id <- dir(file.path(base_dir))
sample_id

kal_dirs <- sapply(sample_id, function(id) file.path(base_dir, id))
kal_dirs <- gsub("/", "\\\\", kal_dirs)
kal_dirs

s2c <- read.table("../../data/RNA-seq_metadata.txt", header = TRUE, sep='\t',stringsAsFactors=FALSE)
s2c <- dplyr::mutate(s2c, path = kal_dirs)

s2c$group <- paste(gsub("_.*", "_", s2c$sample), s2c$condition, sep="")

print(s2c)


##metal&bacteria
s2c$condition <- factor(s2c$condition ,levels = c("NIL","SynCom"))
levels(s2c$condition)

s2c$metal <- factor(s2c$metal, levels = c("CK","Cd","Mn"))
levels(s2c$metal)

s2c$group <- factor(s2c$group,levels = unique(s2c$group))
levels(s2c$group)


###Construct sleuth object，双因素
so_test = sleuth_prep(s2c, ~ metal + condition)
tpm_norm <- as.data.frame(sleuth_to_matrix(so_test,'obs_norm', 'tpm'))
tpm_raw <- as.data.frame(sleuth_to_matrix(so_test,'obs_raw', 'tpm'))


# Load gene-level annotation
annoData <- read.delim(file="../../data/IRGSP-1.0_representative_annotation_2025-03-19.tsv", header = T, sep = "\t", fill = TRUE, na.strings = "", check.names = F)
names(annoData)
annoData_unique <- annoData[!duplicated(annoData$Transcript_ID), 1:11]
annoData_unique <- annoData_unique[annoData_unique$Transcript_ID != "", ]

tpm_norm$Transcript_ID <- rownames(tpm_norm)
normData = merge(annoData_unique, tpm_norm, by.x = "Transcript_ID", by.y = "Transcript_ID", all.y = T)
normData[is.na(normData)] <- ""

write.csv(normData, file="metal_NIL_Syncom_DGE_sleuth_norm_tpm_genome.csv", row.names = F)
write.csv(tpm_raw, file="metal_NIL_Syncom_DGE_sleuth_raw_tpm.csv")


###两两比较
s2c$group
Metal_comp <- list(c("CK_NIL","Cd_NIL"), c("CK_NIL","Mn_NIL"), c("CK_NIL","CK_SynCom"), c("Cd_NIL","Cd_SynCom"), c("Mn_NIL","Mn_SynCom"),  c("CK_NIL","Mn_SynCom"))

for (xxx in Metal_comp) {
  sleuth_data <- dplyr::filter(s2c, group == xxx[1] | group == xxx[2])
  sleuth_data <- droplevels(sleuth_data) 
  
  so <- sleuth_prep(sleuth_data, ~ group)
  so <- sleuth_fit(so, ~group, fit_name = "full")
  so <- sleuth_fit(so, ~ 1, fit_name = "reduced")
  so_lrt <- sleuth_lrt(so, "reduced", "full")
  
  sleuth_table <- sleuth_results(so_lrt, 'reduced:full', 'lrt', show_all = F)
  sleuth_matrix <- as.data.frame(sleuth_to_matrix(so_lrt, 'obs_norm', 'tpm'))
  
  ordered_columns <- match(sleuth_data$sample, colnames(sleuth_matrix)) 
  sleuth_matrix <- sleuth_matrix[, ordered_columns]
  
  sleuth_matrix_mean <- as.data.frame(lapply(levels(sleuth_data$group), function(lvl) {
    #先根据分组筛选样本名
    samples <- sleuth_data$sample[sleuth_data$group == lvl]
    #再根据样本名提取数据
    rowMeans(sleuth_matrix[, samples, drop = FALSE])
  })
  )
  
  # 显式设置列名
  colnames(sleuth_matrix_mean) <- levels(sleuth_data$group)
  
  sleuth_matrix_mean$log2fc = log2((sleuth_matrix_mean[,2]+0.5)/(sleuth_matrix_mean[,1]+0.5))
  
  sleuth_table = merge.data.frame(sleuth_matrix_mean, sleuth_table, by.x = "row.names", by.y = "target_id")
  names(sleuth_table)[1] <- "target_id"
  
  sleuth_matrix = merge.data.frame(sleuth_matrix_mean, sleuth_matrix, by = "row.names")
  names(sleuth_matrix)[1] <- "target_id"
  
  sleuth_significant <- dplyr::filter(sleuth_table, qval <= 0.05) 
  
  write.csv(sleuth_table, paste("Metal_", xxx[1], "vs", xxx[2],"sleuth_gene_level.csv"),row.names=F)
  write.csv(sleuth_significant, paste("Metal_", xxx[1], "vs", xxx[2],"sleuth_significant_gene_level.csv"),row.names=F)
  write.csv(sleuth_matrix, paste("Metal_", xxx[1],"vs", xxx[2],"sleuth_tpm_norm_gene_level.csv"),row.names=F)    
}





