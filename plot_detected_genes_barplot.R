# 统计每个样本FPKM>1的基因数并绘制柱状图
# =====================
library(data.table)
library(ggplot2)

# 读取分组信息
group_info <- fread('新建 文本文档.txt', header=TRUE)

# 读取表达量数据
expr <- fread('总.CSV', header=TRUE)
sample_names <- colnames(expr)[-1]

# 统计每个样本FPKM>1的基因数
count_genes <- sapply(sample_names, function(smp) sum(expr[[smp]] > 1, na.rm=TRUE))
bar_df <- data.frame(Sample=sample_names, Detected=count_genes)

# 合并分组信息（如需分组着色）
bar_df <- merge(bar_df, group_info, by.x='Sample', by.y='Sample', all.x=TRUE)

# 绘制柱状图，配色和主题与plot_density.R一致
p <- ggplot(bar_df, aes(x=Sample, y=Detected, fill=Group)) +
  geom_bar(stat='identity', width=0.7) +
  scale_fill_brewer(palette='Set1') +
  labs(title='Detected Genes per Sample (FPKM > 1)',
       x='Sample', y='Number of Detected Genes') +
  theme_classic(base_size=16) +
  theme(
    axis.text.x=element_text(angle=45, hjust=1, vjust=1, size=12),
    plot.title=element_text(hjust=0.5, face='bold'),
    legend.title=element_text(face='bold'),
    legend.position='right'
  )
print(p) 