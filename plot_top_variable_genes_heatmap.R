# Heatmap of Top Variable Genes
# 展示最具变异性的基因表达聚类效果
library(data.table)
library(pheatmap)
library(RColorBrewer)

# 参数：取前N个高变异基因
N <- 100

# 读取表达量数据
expr <- fread('总.CSV', header=TRUE)
gene_names <- expr[[1]]
expr_mat <- as.matrix(expr[,-1, with=FALSE])
rownames(expr_mat) <- gene_names

# 计算log2(FPKM+1)
log_mat <- log2(expr_mat + 1)

# 计算每个基因的方差，选前N个
vars <- apply(log_mat, 1, var, na.rm=TRUE)
top_idx <- order(vars, decreasing=TRUE)[1:N]
top_log_mat <- log_mat[top_idx, ]

# 行Z-score标准化
zscore_mat <- t(scale(t(top_log_mat)))

# 读取分组信息
group_info <- fread('新建 文本文档.txt', header=TRUE)
ann <- data.frame(row.names=group_info$Sample, Group=group_info$Group)

# 保证样本顺序一致
zscore_mat <- zscore_mat[, rownames(ann), drop=FALSE]

# 分组注释颜色，Set1调色板
group_levels <- unique(group_info$Group)
group_colors <- setNames(brewer.pal(max(3, length(group_levels)), 'Set1')[1:length(group_levels)], group_levels)
ann_colors <- list(Group=group_colors)

# 绘制热图
pheatmap(zscore_mat,
         annotation_col=ann,
         annotation_colors=ann_colors,
         color=colorRampPalette(c('#2166AC','white','#B2182B'))(100),
         show_rownames=FALSE,
         show_colnames=TRUE,
         cluster_cols=FALSE,
         cluster_rows=TRUE,
         main=paste('Top', N, 'Most Variable Genes'),
         fontsize=12) 