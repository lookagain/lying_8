#--导入所需R包#-------
library(igraph)
library(EasyStat)
library(dplyr)
library(ggpubr)
library(Rmisc)

# Set Work Path
pwd <- dirname(rstudioapi::getSourceEditorContext()$path)
setwd(pwd)



#读入文件
taxonomy = read.csv("../../data/ASV_Taxonomy.csv", row.names=1,header = T,check.names = F)
taxonomy[is.na(taxonomy)] <- ""


#读入对应的相关性文件
occor <- readRDS("Root Village rcorr list.rds")
occor_r = occor$r
occor_p = occor$P
occor_p.adj = p.adjust(occor_p, method = 'fdr')
occor_r[occor_p.adj >0.05|abs(occor_r<0.6)] = 0
diag(occor_r) <- 0
dim(occor_r)

#去除无关联节点
r.matrix <- occor_r
remove_rows_cols <- which(rowSums(r.matrix) == 0 & colSums(r.matrix) == 0)
r.matrix <- r.matrix[-remove_rows_cols, -remove_rows_cols]
dim(r.matrix )

igraph <- igraph::graph_from_adjacency_matrix(
  r.matrix,
  mode = "undirected",    # 根据你的网络类型选择
  weighted = TRUE,        # 保留相关性值作为权重
  diag = FALSE            # 忽略对角线
)

node1 <- as_data_frame(igraph, what = "vertices")

taxon1 <- taxonomy[rownames(taxonomy) %in% rownames(node1),]

A = igraph
A <- simplify(igraph)

node1$Phylum <- taxon1$Phylum[match(node1$name, rownames(taxon1))]
node1$Phylum <- gsub("^p__", "", node1$Phylum)

node1$Family <- taxon1$Family[match(node1$name, rownames(taxon1))]
node1$Family <- gsub("^f__", "", node1$Family)

sort(table(node1$Phylum), decreasing = TRUE)

node1$degree <- igraph::degree(A, v = igraph::V(A))  # V(A) 表示网络中所有节点

write.csv(node1, "Root_Village_node.csv", row.names = F)




#####
#####Root  Factory
#读入对应的相关性文件
occor <- readRDS("Root Factory rcorr list.rds")
occor_r = occor$r
occor_p = occor$P
occor_p.adj = p.adjust(occor_p, method = 'fdr')
occor_r[occor_p.adj >0.05|abs(occor_r<0.6)] = 0
diag(occor_r) <- 0
dim(occor_r)


#去除无关联节点
r.matrix <- occor_r
remove_rows_cols <- which(rowSums(r.matrix) == 0 & colSums(r.matrix) == 0)
r.matrix <- r.matrix[-remove_rows_cols, -remove_rows_cols]
dim(r.matrix )

igraph <- igraph::graph_from_adjacency_matrix(
  r.matrix,
  mode = "undirected",    # 根据你的网络类型选择
  weighted = TRUE,        # 保留相关性值作为权重
  diag = FALSE            # 忽略对角线
)

node2 <- as_data_frame(igraph, what = "vertices")

taxon2 <- taxonomy[rownames(taxonomy) %in% rownames(node2),]

A = igraph
A <- simplify(igraph)

node2$Phylum <- taxon2$Phylum[match(node2$name, rownames(taxon2))]
node2$Phylum <- gsub("^p__", "", node2$Phylum)

node2$Family <- taxon2$Family[match(node2$name, rownames(taxon2))]
node2$Family <- gsub("^f__", "", node2$Family)

sort(table(node2$Phylum), decreasing = TRUE)

node2$degree <- igraph::degree(A, v = igraph::V(A))  # V(A) 表示网络中所有节点

write.csv(node2, "Root_Factory_node.csv", row.names = F)




#####
#####            Soil  Village
#####
#读入对应的相关性文件
occor <- readRDS("Soil Village rcorr list.rds")
occor_r = occor$r
occor_p = occor$P
occor_p.adj = p.adjust(occor_p, method = 'fdr')
occor_r[occor_p.adj >0.05|abs(occor_r<0.6)] = 0
diag(occor_r) <- 0
dim(occor_r)


#去除无关联节点
r.matrix <- occor_r
remove_rows_cols <- which(rowSums(r.matrix) == 0 & colSums(r.matrix) == 0)
r.matrix <- r.matrix[-remove_rows_cols, -remove_rows_cols]
dim(r.matrix )

igraph <- igraph::graph_from_adjacency_matrix(
  r.matrix,
  mode = "undirected",    # 根据你的网络类型选择
  weighted = TRUE,        # 保留相关性值作为权重
  diag = FALSE            # 忽略对角线
)

node3 <- as_data_frame(igraph, what = "vertices")

taxon3 <- taxonomy[rownames(taxonomy) %in% rownames(node3),]

A = igraph
A <- simplify(igraph)

node3$Phylum <- taxon3$Phylum[match(node3$name, rownames(taxon3))]
node3$Phylum <- gsub("^p__", "", node3$Phylum)

node3$Family <- taxon3$Family[match(node3$name, rownames(taxon3))]
node3$Family <- gsub("^f__", "", node3$Family)

sort(table(node3$Phylum), decreasing = TRUE)

node3$degree <- igraph::degree(A, v = igraph::V(A))  # V(A) 表示网络中所有节点

write.csv(node3, "Soil_Village_node.csv", row.names = F)




#####
#####            Soil  Factory
#####
#读入对应的相关性文件
occor <- readRDS("Soil Factory rcorr list.rds")
occor_r = occor$r
occor_p = occor$P
occor_p.adj = p.adjust(occor_p, method = 'fdr')
occor_r[occor_p.adj >0.05|abs(occor_r<0.6)] = 0
diag(occor_r) <- 0
dim(occor_r)


#去除无关联节点
r.matrix <- occor_r
remove_rows_cols <- which(rowSums(r.matrix) == 0 & colSums(r.matrix) == 0)
r.matrix <- r.matrix[-remove_rows_cols, -remove_rows_cols]
dim(r.matrix )

igraph <- igraph::graph_from_adjacency_matrix(
  r.matrix,
  mode = "undirected",    # 根据你的网络类型选择
  weighted = TRUE,        # 保留相关性值作为权重
  diag = FALSE            # 忽略对角线
)

node4 <- as_data_frame(igraph, what = "vertices")

taxon4 <- taxonomy[rownames(taxonomy) %in% rownames(node4),]

A = igraph
A <- simplify(igraph)

node4$Phylum <- taxon4$Phylum[match(node4$name, rownames(taxon4))]
node4$Phylum <- gsub("^p__", "", node4$Phylum)

node4$Family <- taxon4$Family[match(node4$name, rownames(taxon4))]
node4$Family <- gsub("^f__", "", node4$Family)

sort(table(node4$Phylum), decreasing = TRUE)

node4$degree <- igraph::degree(A, v = igraph::V(A))  # V(A) 表示网络中所有节点

write.csv(node4, "Soil_Factory_node.csv", row.names = F)







#读入文件
R_V <- read.csv("Root_Village_node.csv", check.names = F)
R_F <- read.csv("Root_Factory_node.csv", check.names = F)
S_V <- read.csv("Soil_Village_node.csv", check.names = F)
S_F <- read.csv("Soil_Factory_node.csv", check.names = F)



#-----------------------------------------------------------------------Family degree
##                                           Root    ###
R_V$Site <- "Village"
R_F$Site <- "Factory"

sub1 <- rbind(R_V,R_F)
sort(table(sub1$Family), decreasing = TRUE)

#提取Family-level 名称
nameRF <- names(sort(table(sub1$Family), decreasing = TRUE)[2:12])

subRF <- sub1
subRF$Family <- ifelse(subRF$Family %in% nameRF, subRF$Family, "Others")

nameRF <- c(names(sort(table(sub1$Family), decreasing = TRUE)[2:12]), "Others")

# 初始化空列表，制作显著性表格
sig_results_list <- list()

if (!dir.exists("Tukey检验结果")) {
  dir.create("Tukey检验结果")
}

# 3. 循环遍历每个家族：用seq_along获取索引n（1到11），同时获取家族名称
for (n in seq_along(nameRF)) {
  # 当前家族名称（根据索引n获取）
  family_name <- nameRF[n]
  
  # 3.1 筛选当前家族的数据
  sig1 <- subRF[subRF$Family == family_name, c("name", "Site", "degree")]
  
  # 3.2 重命名列 + 转换因子水平
  colnames(sig1)[2] <- "group" 
  sig1$group <- factor(sig1$group, levels = c("Village", "Factory"))
  
  # 3.3 执行显著性检验
  res <- KwWlx(data = sig1, i = 3)
  res[[2]]
  
  
  # 动态生成文件名：科名称 + 固定后缀，保存到指定文件夹
  file_name <- paste0(
    "Tukey检验结果/",  # 保存到创建的文件夹中（避免杂乱）
    family_name,       # 动态替换为当前家族名称
    "_Root_Site_node_degree_Tukey_result.txt"  # 固定后缀
  )
  
  # 保存文件（制表符分隔，无行名，无引号）
  write.table(
    x = res[[2]],          # 要保存的数据（检验结果）
    file = file_name,      # 动态生成的文件名
    sep = "\t",            # 制表符分隔（方便用Excel打开）
    row.names = FALSE,     # 不保留行名
    col.names = TRUE,      # 保留列名
    quote = FALSE          # 不添加引号（避免数据格式错乱）
  )
  
  
  
  # 3.4 计算y.position（动态获取当前家族数据的y轴最大值，避免遮挡）
  # 假设y轴数据列是"degree"（根据你的代码逻辑），若列名不同可替换
  y_max <- max(sig1$degree, na.rm = TRUE) + 2
  
  # 3.5 动态设置xmin和xmax：第n个家族 → xmin = n-0.5，xmax = n+0.5
  xmin_val <- n - 0.2
  xmax_val <- n + 0.2
  
  # 3.6 生成当前家族的R_sig数据框（含动态xmin/xmax）
  current_R_sig <- data.frame(
    Family = family_name,        # 家族名称（标识用）
    group1 = "Village",          # 比较组1
    group2 = "Factory",          # 比较组2
    y.position = y_max,          # 标注高度（动态计算）
    p.adj = res[[2]] %>% pull(p.adj),  # 校正后p值
    xmin = xmin_val,             # 动态xmin（n-0.5）
    xmax = xmax_val,             # 动态xmax（n+0.5）
    stringsAsFactors = FALSE
  )
  
  # 3.7 将当前结果存入列表
  sig_results_list[[family_name]] <- current_R_sig
}


final_R_sig <- bind_rows(sig_results_list)
final_R_sig$p.adj <- ifelse(final_R_sig$p.adj <0.05, paste0("P = ", final_R_sig$p.adj), 'ns')


##
color_F <- c("Village" = "#8DD3C7","Factory" = "#FFFFB3",
             "Proteobacteria" = "#E41A1C",
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
             "#FFFFB3")



subRF$Family <- factor(subRF$Family, levels = nameRF)
subRF$Site <- factor(subRF$Site, levels = c("Village", "Factory"))
subRF$new <- ifelse(subRF$Family != "Others", subRF$Phylum, "Others")

table(subRF$new)
RP_name <- names(table(subRF$new))[c(4,1,2,3)]
subRF$new <- factor(subRF$new, levels = RP_name)

range(sub1$degree)
# y_max <- max(sub1$degree) + 3
# test <- subRF[subRF$Family == 'Micromonosporaceae',]

set.seed(123)

RP_Family <- ggplot(
  subRF, 
  aes(
    x = Family, 
    y = degree, 
    #fill = Site,        # 用于图例映射（与颜色对应）
    #color = Site,       # 抖动点颜色映射（与 Site 关联）
    group = interaction(Family, Site)  # 明确分组：每个 Family+Site 为一组
  )
) +
  # 箱线图：黑色边框 + 透明填充，每个 Family 下显示 2 个（Village/Factory）
  geom_boxplot(
    aes(color = new),
    fill = NA,                # 透明填充，突出黑色边框
    outlier.shape = NA,
    width = 0.8,              # 箱线图宽度（保留你的设置）
    linewidth = 0.7,          # 边框线宽（保留你的设置）
    #color = "black",           # 关键：所有箱线图边框设为黑色
    position = position_dodge(width = 0.8)  # 分组间距（与点保持一致）
  ) +
  # 抖动点：颜色与 Site 对应，避免重叠
  geom_jitter(
    aes( fill = Site),
    size = 2, 
    alpha = 1,
    color = "black",  # 关键：点的边框固定为黑色
    shape = 21,       # 关键：选择有填充+边框的形状（21是圆形，推荐）
    # 关键：用 position_jitterdodge() 替代 position_dodge()，支持微小抖动避免重叠
    position = position_jitterdodge(
      dodge.width = 0.8,  # 与箱线图的 dodge 间距一致，确保对齐
      jitter.width = 0.3,   # 水平方向无抖动（确保排成垂直线）
      jitter.height = 0.5 # 垂直方向微小抖动（避免点完全重叠，可设为 0）
    )
  ) +
  theme_bw() +
  ggtitle("Root") +
  ylab("Degree") +
  xlab("") +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_blank(),
    axis.line = element_line(color = "black", linewidth = 0.7),
    legend.title = element_text(size=12, color="black", face = "bold"),
    legend.text = element_text(size=12, color="black"),
    axis.text.x = element_text(size=8, color="black", hjust = 0.7, vjust = 0.9, angle = 20),
    axis.text.y = element_text(size=12, color="black"),
    plot.margin = margin(t = 20, r = 10, b = 20, l = 20),
    axis.title = element_text(size=12, color="black", face= "bold"),
    plot.title = element_text(size=14, color="black", face= "bold", hjust = 0.5, vjust = 1)
  )+
  scale_y_continuous(limits = c(NA, 40), breaks = seq(0, 40, by = 10))+
  stat_pvalue_manual(
    final_R_sig,                                          #数据框
    label = "{p.adj}",              #选择的标签，"{p.adj}{p.adj.signif}" 能选择两列作为标签
    y.position = "y.position",        #标签的位置，需要依据成图的数据计算出来。字符刻度  中  每个字符代表一个单位，如 A 就是1， B就是2等等
    xmin = "xmin",                            #标签下的比较线的起始位置，与xmax搭配使用，控制比较线的长度
    xmax = "xmax",                           #
    tip.length = 0.01,                      #比较线的粗细
    size = 3,
    coord.flip = F                     #T：翻转标签， F：不翻转，默认不翻转
  ) 

p3 <- RP_Family +
  scale_color_manual(
    values = color_F,         # 为不同 Site 分配颜色
    name = "Phylum",           
  ) +
  scale_fill_manual(values = color_F)

p3
ggsave("../../figures/Fig.3/Root_Site_Family_degree_plot.pdf",p3, width = 10, height = 5.5)

##



##################
##                                           Soil    ###
S_V$Site <- "Village"
S_F$Site <- "Factory"

sub2 <- rbind(S_V,S_F)
sort(table(sub2$Family), decreasing = TRUE)

#提取Family-level 名称
nameSF <- names(sort(table(sub2$Family), decreasing = TRUE)[c(1:6,8:12)])

subSF <- sub2
subSF$Family <- ifelse(subSF$Family %in% nameSF, subSF$Family, "Others")

nameSF <- c(names(sort(table(sub2$Family), decreasing = TRUE)[c(1:6,8:12)]), "Others")

two_sites_families <- subSF %>%
  group_by(Family) %>%          # 按Family分组
  dplyr::summarise(
    site_count = n_distinct(Site),  # 统计每个Family下唯一的Site数量
    .groups = "drop"               # 取消分组，避免后续警告
  ) %>%
  filter(site_count == 2) %>%    # 筛选有2种Site的Family
  pull(Family)

nameSF <- nameSF[nameSF %in% two_sites_families]


# 初始化空列表，制作显著性表格
sig_results_list <- list()

# 3. 循环遍历每个家族：用seq_along获取索引n（1到11），同时获取家族名称
for (n in seq_along(nameSF)) {
  # 当前家族名称（根据索引n获取）
  family_name <- nameSF[n]
  
  # 3.1 筛选当前家族的数据
  sig1 <- subSF[subSF$Family == family_name, c("name", "Site", "degree")]
  
  # 3.2 重命名列 + 转换因子水平
  colnames(sig1)[2] <- "group" 
  sig1$group <- factor(sig1$group, levels = c("Village", "Factory"))
  
  # 3.3 执行显著性检验
  res <- KwWlx(data = sig1, i = 3)
  res[[2]]
  
  
  # 动态生成文件名：科名称 + 固定后缀，保存到指定文件夹
  file_name <- paste0(
    "Tukey检验结果/",  # 保存到创建的文件夹中（避免杂乱）
    family_name,       # 动态替换为当前家族名称
    "_Soil_Site_node_degree_Tukey_result.txt"  # 固定后缀
  )
  
  # 保存文件（制表符分隔，无行名，无引号）
  write.table(
    x = res[[2]],          # 要保存的数据（检验结果）
    file = file_name,      # 动态生成的文件名
    sep = "\t",            # 制表符分隔（方便用Excel打开）
    row.names = FALSE,     # 不保留行名
    col.names = TRUE,      # 保留列名
    quote = FALSE          # 不添加引号（避免数据格式错乱）
  )
  
  
  
  # 3.4 计算y.position（动态获取当前家族数据的y轴最大值，避免遮挡）
  # 假设y轴数据列是"degree"（根据你的代码逻辑），若列名不同可替换
  y_max <- max(sig1$degree, na.rm = TRUE) + 2
  
  # 3.5 动态设置xmin和xmax：第n个家族 → xmin = n-0.5，xmax = n+0.5
  xmin_val <- n - 0.2
  xmax_val <- n + 0.2
  
  # 3.6 生成当前家族的R_sig数据框（含动态xmin/xmax）
  current_R_sig <- data.frame(
    Family = family_name,        # 家族名称（标识用）
    group1 = "Village",          # 比较组1
    group2 = "Factory",          # 比较组2
    y.position = y_max,          # 标注高度（动态计算）
    p.adj = res[[2]] %>% pull(p.adj),  # 校正后p值
    xmin = xmin_val,             # 动态xmin（n-0.5）
    xmax = xmax_val,             # 动态xmax（n+0.5）
    stringsAsFactors = FALSE
  )
  
  # 3.7 将当前结果存入列表
  sig_results_list[[family_name]] <- current_R_sig
}


final_S_sig <- bind_rows(sig_results_list)
final_S_sig$p.adj <- ifelse(final_S_sig$p.adj <0.05, paste0("P = ", final_S_sig$p.adj), 'ns')


##因为排在第9位的"Planococcaceae"只在一个site中存在，所以这里做一些调整
nameSF <- c(names(sort(table(sub2$Family), decreasing = TRUE)[c(1:6,8:12)]), "Others")

max(subSF$degree[subSF$Family == "Planococcaceae"], na.rm = TRUE) + 2

new_row <- data.frame(
  Family = "Planococcaceae",  # 示例Family名称
  group1 = "Village",           # 固定为Village
  group2 = "Factory",           # 固定为Factory
  y.position = 19,              # 示例y轴位置
  p.adj = "",                # 示例p值
  xmin = 7.8,                   # 示例xmin
  xmax = 8.2,                   # 示例xmax
  stringsAsFactors = FALSE      # 避免自动转因子（关键！）
)

final_S_sig <- rbind(final_S_sig, new_row)
final_S_sig$Family <- factor(final_S_sig$Family, levels = nameSF)

final_S_sig <- final_S_sig %>%  arrange(Family)

#调整xmin和xmax的值，符合排序
final_S_sig$xmin <- c(0.8, 1.8, 2.8, 3.8, 4.8, 5.8, 6.8, 7.8, 8.8, 9.8, 10.8, 11.8)
final_S_sig$xmax <- c(1.2, 2.2, 3.2, 4.2, 5.2, 6.2, 7.2, 8.2, 9.2, 10.2, 11.2, 12.2)


##颜色沿用上面的
coloS_F <- color_F
#

subSF$Family <- factor(subSF$Family, levels = nameSF)
subSF$Site <- factor(subSF$Site, levels = c("Village", "Factory"))
subSF$new <- ifelse(subSF$Family != "Others", subSF$Phylum, "Others")

table(subSF$new)
SP_name <- names(table(subSF$new))[c(4,2,1,3)]
subSF$new <- factor(subSF$new, levels = SP_name)

range(sub2$degree)

# test <- subSF[subSF$Family == "Rhizobiaceae",]
set.seed(123)

SP_Family <- ggplot(
  subSF, 
  aes(
    x = Family, 
    y = degree, 
    group = interaction(Family, Site)  # 明确分组：每个 Family+Site 为一组
  )
) +
  # 箱线图：黑色边框 + 透明填充，每个 Family 下显示 2 个（Village/Factory）
  geom_boxplot(
    aes(color = new),
    fill = NA,                # 透明填充，突出黑色边框
    outlier.shape = NA,
    width = 0.8,              # 箱线图宽度（保留你的设置）
    linewidth = 0.7,          # 边框线宽（保留你的设置）
    position = position_dodge(width = 0.8)  # 分组间距（与点保持一致）
  ) +
  # 抖动点：颜色与 Site 对应，避免重叠
  geom_jitter(
    aes( fill = Site),
    size = 2, 
    alpha = 1,
    color = "black",  # 关键：点的边框固定为黑色
    shape = 21,       # 关键：选择有填充+边框的形状（21是圆形，推荐）
    # 关键：用 position_jitterdodge() 替代 position_dodge()，支持微小抖动避免重叠
    position = position_jitterdodge(
      dodge.width = 0.8,  # 与箱线图的 dodge 间距一致，确保对齐
      jitter.width = 0.3,   # 水平方向无抖动（确保排成垂直线）
      jitter.height = 0.5 # 垂直方向微小抖动（避免点完全重叠，可设为 0）
    )
  ) +
  theme_bw() +
  ggtitle("Soil") +
  ylab("Degree") +
  xlab("") +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_blank(),
    axis.line = element_line(color = "black", linewidth = 0.7),
    legend.title = element_text(size=12, color="black", face = "bold"),
    legend.text = element_text(size=12, color="black"),
    axis.text.x = element_text(size=8, color="black", hjust = 0.7, vjust = 0.9, angle = 20),
    axis.text.y = element_text(size=12, color="black"),
    plot.margin = margin(t = 20, r = 10, b = 20, l = 20),
    axis.title = element_text(size=12, color="black", face= "bold"),
    plot.title = element_text(size=14, color="black", face= "bold", hjust = 0.5, vjust = 1)
  )+
  scale_y_continuous(limits = c(NA, 40), breaks = seq(0, 40, by = 10))+
  stat_pvalue_manual(
    final_S_sig,                                          #数据框
    label = "{p.adj}",              #选择的标签，"{p.adj}{p.adj.signif}" 能选择两列作为标签
    y.position = "y.position",        #标签的位置，需要依据成图的数据计算出来。字符刻度  中  每个字符代表一个单位，如 A 就是1， B就是2等等
    xmin = "xmin",                            #标签下的比较线的起始位置，与xmax搭配使用，控制比较线的长度
    xmax = "xmax",                           #
    tip.length = 0.01,                      #比较线的粗细
    size = 3,
    coord.flip = F                     #T：翻转标签， F：不翻转，默认不翻转
  ) 

p4 <- SP_Family +
  # 颜色配置：抖动点颜色与 Site 对应（使用 coloS_F）
  scale_color_manual(
    values = coloS_F,         # 为不同 Site 分配颜色
    name = "Phylum"
  ) +
  scale_fill_manual(values = coloS_F) 

p4
ggsave("../../figures/Fig.3/Soil_Site_Family_degree_plot.pdf", p4, width = 10, height = 5.5)






#-------------------------------------------------------------------------------------------------

#-------------------------------------------------------total degree
#                        Root
# #Significance test
R_sig_test <- sub1[,c("name","Site","degree")]
colnames(R_sig_test)[2] <- "group"
R_sig_test$group <- factor(R_sig_test$group, levels = c("Village","Factory"))

res = KwWlx(data = R_sig_test, i= 3)
res[[1]]
res[[2]]
#write.table(res[[2]], "Root_Site_node_degree_Tukey_result.txt", sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)

y_max <- max(R_sig_test$degree) + 2

#sig label
R_sig <- data.frame(
  #Site = "",
  group1 = "Village",
  group2 = "Factory",
  y.position = y_max,
  p.adj = res[[2]] %>% pull(p.adj),
  # p.adj.signif = "****",
  xmin = c(1),
  xmax = c(2)
)


#set colors
col_Site <- c("#8DD3C7","#FFFFB3","#FFFFB3")


#
sub1_SE <- summarySE(sub1, measurevar= "degree", groupvars= "Site")
sub1_SE$Site <- factor(sub1_SE$Site, levels = c("Village","Factory"))


P_R_Site <- ggplot(sub1_SE, aes(Site, degree, fill=Site)) + 
  geom_bar(stat="identity", width = 0.8,color="black") +# 添加color="black"来设置柱子的边框颜色
  scale_fill_manual(values=col_Site) +
  xlab("") + ylab("Degree") + 
  ggtitle("Root") +
  theme_bw() +
  theme(legend.position="right",
  ) + 
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        axis.line = element_line(color = "black", linewidth = 0.7),
        legend.title = element_text(size=12,color="black", face = "bold"),
        legend.text= element_text(size=12,color="black"),
        axis.text.x = element_text(size=12,color="black", hjust = 0.5, vjust = 1),
        axis.text.y = element_text(size=12,color="black"),
        plot.margin = margin(t = 20, r = 10, b = 20, l = 20),
        axis.title= element_text(size=12,color="black", face= "bold"),
        plot.title =element_text(size=14,color="black",face= "bold", hjust = 0.5, vjust = 1))+
  geom_jitter(data = sub1, aes(Site, degree), fill = "#99CCFF", color = "#99CCFF", shape = 21, size=3, width = 0.3, alpha = 0.6) +
  geom_errorbar(data = sub1_SE, aes(ymin=degree-se, ymax=degree+se), width=0.2, colour = "black") +
  scale_y_continuous(limits = c(NA, 40), breaks = seq(0, 40, by = 10))+  # 设置 y 轴的最大值
  
  stat_pvalue_manual(
    R_sig,                                          #数据框
    label = "P = {p.adj}",              #选择的标签，"{p.adj}{p.adj.signif}" 能选择两列作为标签
    y.position = "y.position",        #标签的位置，需要依据成图的数据计算出来。字符刻度  中  每个字符代表一个单位，如 A 就是1， B就是2等等
    xmin = "xmin",                            #标签下的比较线的起始位置，与xmax搭配使用，控制比较线的长度
    xmax = "xmax",                           #
    tip.length = 0.01,                      #比较线的粗细
    size = 4,
    label.fontface = "italic" ,
    coord.flip = F                     #T：翻转标签， F：不翻转，默认不翻转
  ) 

P_R_Site

ggsave("../../figures/Fig.3/Root_Site_degree_plot.pdf",P_R_Site, width = 4.5, height = 5)



#-------------------------------------------------------total degree
#                        Soil
S_V$Site <- "Village"
S_F$Site <- "Factory"

sub2 <- rbind(S_V,S_F)
range(sub2$degree)
sub_max <- max(sub2$degree) + 3


#Significance test
S_sig_test <- sub2[,c("name","Site","degree")]
colnames(S_sig_test)[2] <- "group" 
S_sig_test$group <- factor(S_sig_test$group, levels = c("Village","Factory"))

res = KwWlx(data = S_sig_test, i= 3)
res[[1]]
res[[2]]


#sig label
S_sig <- data.frame(
  Site = "",
  group1 = "Village",
  group2 = "Factory",
  y.position = sub_max,
  p.adj = "5e-16",
  p.adj.signif = "****",
  xmin = c(1),
  xmax = c(2)
)


#set colors
col_Site <- c("#8DD3C7","#FFFFB3","#FFFFB3")


#
sub2_SE <- summarySE(sub2, measurevar= "degree", groupvars= "Site")
sub2_SE$Site <- factor(sub2_SE$Site, levels = c("Village","Factory"))


P_S_Site <- ggplot(sub2_SE, aes(Site, degree, fill=Site)) + 
  geom_bar(stat="identity", width = 0.8,color="black") +# 添加color="black"来设置柱子的边框颜色
  scale_fill_manual(values=col_Site) +
  xlab("") + ylab("Degree") + 
  ggtitle("Soil") +
  theme_bw() +
  theme(legend.position="right",
  ) + 
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        axis.line = element_line(color = "black", linewidth = 0.7),
        legend.title = element_text(size=12,color="black", face = "bold"),
        legend.text= element_text(size=12,color="black"),
        axis.text.x = element_text(size=12,color="black", hjust = 0.5, vjust = 1),
        axis.text.y = element_text(size=12,color="black"),
        plot.margin = margin(t = 20, r = 10, b = 20, l = 20),
        axis.title= element_text(size=12,color="black", face= "bold"),
        plot.title =element_text(size=14,color="black",face= "bold", hjust = 0.5, vjust = 1))+
  geom_jitter(data = sub2, aes(Site, degree), fill = "#99CCFF", color = "#99CCFF", shape = 21, size=3, width = 0.3, alpha = 0.6) +
  geom_errorbar(data = sub2_SE, aes(ymin=degree-se, ymax=degree+se), width=0.2, colour = "black") +
  scale_y_continuous(limits = c(NA, 40), breaks = seq(0, 40, by = 10))+  # 设置 y 轴的最大值
  
  stat_pvalue_manual(
    S_sig,                                          #数据框
    label = "P = {p.adj}",              #选择的标签，"{p.adj}{p.adj.signif}" 能选择两列作为标签
    y.position = "y.position",        #标签的位置，需要依据成图的数据计算出来。字符刻度  中  每个字符代表一个单位，如 A 就是1， B就是2等等
    xmin = "xmin",                            #标签下的比较线的起始位置，与xmax搭配使用，控制比较线的长度
    xmax = "xmax",                           #
    tip.length = 0.01,                      #比较线的粗细
    size = 4,
    label.fontface = "italic" ,
    coord.flip = F                     #T：翻转标签， F：不翻转，默认不翻转
  ) 

P_S_Site

ggsave("../../figures/Fig.3/Soil_Site_degree_plot.pdf",P_S_Site, width = 4.5, height = 5)





###########
###                 Phylum 水平的节点分布
#########
sort(table(sub1$Phylum), decreasing = T)
name_R <- c('Proteobacteria', "Actinobacteriota", 'Firmicutes', 'Bacteroidota')
sub1$group <- ifelse(sub1$Phylum %in% name_R, sub1$Phylum, "Others")

sub1_V <- sub1[sub1$Site =="Village",]
sub1_F <- sub1[sub1$Site =="Factory",]

data1_V <- as.data.frame(table(sub1_V$group))
data1_V$group <- "Village"
data1_F <- as.data.frame(table(sub1_F$group))
data1_F$group <- "Factory"

data1 <- rbind(data1_V, data1_F)
colnames(data1) <- c("Phylum", "number", "Site")


#
color_values = c("#8DD3C7","#FFFFB3","#FFFFB3")
data1$Site <- factor(data1$Site, levels = c("Village","Factory"))
data1$Phylum <- factor(data1$Phylum, levels = c(name_R,"Others"))
data1$hjust <- ifelse(data1$Site == "Village", 1.2, -0.5)

p2 <- ggplot(data1, aes(x = Phylum, y = number, fill = Site, group = Site)) + 
  scale_color_manual(values = color_values) + 
  scale_fill_manual(values = color_values) + 
  geom_bar(stat = "identity", ,color="black", position = position_dodge()) + 
  geom_text(label = data1$number, size = 4, colour="black", vjust = -0.8, hjust = data1$hjust)+
  theme_bw() + 
  ggtitle("Root") +
  ylab("Node number") + 
  xlab("") +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_blank(),
    axis.line = element_line(color = "black", linewidth = 0.7),
    legend.title = element_text(size=12,color="black", face = "bold"),
    legend.text= element_text(size=12,color="black"),
    axis.text.x = element_text(size=8,color="black", hjust = 0.5, vjust = 0.7),
    axis.text.y = element_text(size=12,color="black"),
    plot.margin = margin(t = 20, r = 10, b = 20, l = 20),
    axis.title= element_text(size=12,color="black", face= "bold"),
    plot.title =element_text(size=14,color="black",face= "bold", hjust = 0.5, vjust = 1))

p2

ggsave("../../figures/Fig.3/Root_Site_node_number.pdf", p2 , width = 6, height = 5)





#####   Soil   
sort(table(sub2$Phylum), decreasing = T)
name_S <- c('Proteobacteria', "Firmicutes", 'Bacteroidota', 'Actinobacteriota')
sub2$group <- ifelse(sub2$Phylum %in% name_S, sub2$Phylum, "Others")

sub2_V <- sub2[sub2$Site =="Village",]
sub2_F <- sub2[sub2$Site =="Factory",]

data2_V <- as.data.frame(table(sub2_V$group))
data2_V$group <- "Village"
data2_F <- as.data.frame(table(sub2_F$group))
data2_F$group <- "Factory"

data2 <- rbind(data2_V, data2_F)
colnames(data2) <- c("Phylum", "number", "Site")


#
data2$Site <- factor(data2$Site, levels = c("Village","Factory"))
data2$Phylum <- factor(data2$Phylum, levels = c(name_S,"Others"))
data2$hjust <- ifelse(data2$Site == "Village", 1.2, -0.5)

p3 <- ggplot(data2, aes(x = Phylum, y = number, fill = Site, group = Site)) + 
  scale_color_manual(values = color_values) + 
  scale_fill_manual(values = color_values) + 
  geom_bar(stat = "identity", ,color="black", position = position_dodge()) + 
  geom_text(label = data2$number, size = 4, colour="black", vjust = -0.8, hjust = data2$hjust)+
  theme_bw() + 
  ggtitle("Soil") +
  ylab("Node number") + 
  xlab("") +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_blank(),
    axis.line = element_line(color = "black", linewidth = 0.7),
    legend.title = element_text(size=12,color="black", face = "bold"),
    legend.text= element_text(size=12,color="black"),
    axis.text.x = element_text(size=8,color="black", hjust = 0.5, vjust = 0.7),
    axis.text.y = element_text(size=12,color="black"),
    plot.margin = margin(t = 20, r = 10, b = 20, l = 20),
    axis.title= element_text(size=12,color="black", face= "bold"),
    plot.title =element_text(size=14,color="black",face= "bold", hjust = 0.5, vjust = 1))+
  scale_y_continuous(limits = c(NA, 130), breaks = seq(0, 120, by = 30))

p3

ggsave("../../figures/Fig.3/Soil_Site_node_number.pdf", p3 , width = 6, height = 5)






########        Root  Family
###    

sort(table(R_V$Family), decreasing = TRUE)
sort(table(R_F$Family), decreasing = TRUE)

nameR_V <- c(
  "Enterobacteriaceae",
  "Bacillaceae",
  "Pseudomonadaceae",
  "Xanthobacteraceae",
  "Micromonosporaceae",
  "Pseudonocardiaceae",
  "Rhizobiaceae",
  "Paenibacillaceae",
  "Burkholderiaceae",
  "Erwiniaceae"
)

nameR_F <- names(sort(table(R_F$Family), decreasing = TRUE)[c(1, 3:11)])


### Root Village Family-level
sub3 <- R_V[R_V$Family %in% nameR_V,]

data3_V <- as.data.frame(table(sub3$Family))
data3_V$Site <- "Village"
colnames(data3_V)[1:2] <- c("Family", "number")

data3_V$Family <- factor(data3_V$Family, levels = nameR_V)

color_F <- c("Village" = "#8DD3C7","Factory" = "#FFFFB3")

data3_V$Site <- as.factor(data3_V$Site)

p4 <- ggplot(data3_V, aes(x = Family, y = number, fill = Site)) + 
  scale_color_manual(values = color_F) + 
  scale_fill_manual(values = color_F) + 
  geom_bar(stat = "identity", ,color="black", position = position_dodge()) + 
  geom_text(label = data3_V$number, size = 4, colour="black", vjust = -0.8, hjust = 0.5)+
  theme_bw() + 
  ggtitle("Root Village") +
  ylab("Node number") + 
  xlab("") +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_blank(),
    axis.line = element_line(color = "black", linewidth = 0.7),
    legend.title = element_text(size=12,color="black", face = "bold"),
    legend.text= element_text(size=12,color="black"),
    axis.text.x = element_text(size=6,color="black", hjust = 0.7, vjust = 0.9, angle = 20),
    axis.text.y = element_text(size=12,color="black"),
    plot.margin = margin(t = 20, r = 10, b = 20, l = 20),
    axis.title= element_text(size=12,color="black", face= "bold"),
    plot.title =element_text(size=14,color="black",face= "bold", hjust = 0.5, vjust = 1))#+
#scale_y_continuous(limits = c(NA, 130), breaks = seq(0, 120, by = 30))

p4

ggsave("./Root_Village_node_Family_number.pdf", p4 , width = 6, height = 5)











