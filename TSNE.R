library(ggplot2)
library(ggbiplot)
df = read.csv(file = "总.CSV",# 这里读取了网络上的demo数据，将此处换成你自己电脑里的文件
                header = T,    # 指定第一行是列名
                row.names = 1  # 指定第一列是行名
)
#去除方差为0的行
df <- df[apply(df, 1, var)!=0,]
df = t(df) # 对数据进行转置，如果想对基因分组则不用转置

# 读取样本分组数据文件
dfGroup = read.delim("新建 文本文档.txt",
                     header = T,
                     row.names = 1
)

# 自动检查分组信息与样本名是否一致
cat('表达矩阵样本名：', rownames(df), '\n')
cat('分组信息样本名：', rownames(dfGroup), '\n')

# 对齐分组信息顺序
missing_in_group <- setdiff(rownames(df), rownames(dfGroup))
missing_in_data <- setdiff(rownames(dfGroup), rownames(df))
if(length(missing_in_group) > 0) {
  warning(paste('以下样本在分组信息中缺失:', paste(missing_in_group, collapse=',')))
}
if(length(missing_in_data) > 0) {
  warning(paste('以下分组信息在表达矩阵中缺失:', paste(missing_in_data, collapse=',')))
}
dfGroup <- dfGroup[rownames(df), , drop=FALSE]

cat('对齐后分组信息：', dfGroup[,1], '\n')

# PCA计算
pca_result <- prcomp(df,
                     scale=T  # 一个逻辑值，指示在进行分析之前是否应该将变量缩放到具有单位方差
)

# 绘图
ggbiplot(pca_result, 
         var.axes=F,            # 是否为变量画箭头
         obs.scale = 1,         # 横纵比例 
         groups = dfGroup[,1],  # 添加分组信息，为分组文件的第一列
         ellipse = F,           # 是否围绕分组画椭圆
         circle = F)+ 
  geom_text(                      # geom_text一个在图中添加标注的函数
    aes(label=rownames(df)),   # 指定标注的内容为数据框df的行名
    vjust=1.5,            # 指定标记的位置，vjust=1.5 垂直向下1.5个距离。   负数为位置向上标记，正数为位置向下标记
    size=5                # 标记大小
  )

# =====================
# MDS分析与可视化
# =====================
# 计算欧氏距离
mds_dist <- dist(df)
# 进行MDS分析，k=3表示三维
mds_result <- cmdscale(mds_dist, k = 3)
mds_df <- as.data.frame(mds_result)
colnames(mds_df) <- c('MDS1', 'MDS2', 'MDS3')
mds_df$Group <- dfGroup[rownames(mds_df), 1]
mds_df$Sample <- rownames(mds_df)

# 检查MDS分组是否有NA
if(any(is.na(mds_df$Group))) {
  warning('有样本分组信息缺失，请检查分组文件与表达矩阵样本名是否一致！')
  print(mds_df[is.na(mds_df$Group), ])
}

# 绘制三维MDS散点图（plotly）
# 如未安装plotly，请先 install.packages('plotly')
library(plotly)
plot_ly(
  mds_df, x = ~MDS1, y = ~MDS2, z = ~MDS3,
  color = ~Group,
  type = 'scatter3d', mode = 'markers',
  marker = list(size = 6)
) %>%
  layout(title = '3D MDS Plot')

# =====================
# t-SNE分析与可视化
# =====================
# 如未安装Rtsne，请先 install.packages('Rtsne')
library(Rtsne)
set.seed(123) # 保证可重复
# t-SNE输入需为数值型矩阵，且不能有NA
if(any(is.na(df))) stop('表达矩阵含有NA，无法进行t-SNE分析')
tsne_result <- Rtsne(df, dims = 2, perplexity = 3, verbose = TRUE, max_iter = 1000)
tsne_df <- as.data.frame(tsne_result$Y)
colnames(tsne_df) <- c('tSNE1', 'tSNE2')
rownames(tsne_df) <- rownames(df)   # 关键：加上样本名
tsne_df$Group <- dfGroup[rownames(tsne_df), 1]
tsne_df$Sample <- rownames(tsne_df)

# 检查t-SNE分组是否有NA
if(any(is.na(tsne_df$Group))) {
  warning('有样本分组信息缺失，请检查分组文件与表达矩阵样本名是否一致！')
  print(tsne_df[is.na(tsne_df$Group), ])
}

# 绘制t-SNE二维分组散点图（无label，仅图例）
library(ggplot2)
ggplot(tsne_df, aes(x = tSNE1, y = tSNE2, color = Group)) +
  geom_point(size = 4) +
  theme_minimal() +
  labs(title = 't-SNE Plot', x = 't-SNE1', y = 't-SNE2')
