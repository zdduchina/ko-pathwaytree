library(data.table)
library(reshape2)
library(ggplot2)

# 1. 读取分组信息
group_info <- read.table("新建 文本文档.txt", header=TRUE, sep="\t", stringsAsFactors=FALSE)

# 2. 读取表达量数据
expr <- fread("总.CSV", header=FALSE)
colnames(expr) <- c("Gene", group_info$Sample)

# 3. 转为长表格
expr_long <- melt(expr, id.vars="Gene", variable.name="Sample", value.name="FPKM")
expr_long <- merge(expr_long, group_info, by="Sample")
expr_long$logFPKM <- log2(as.numeric(expr_long$FPKM) + 1)

# 4. 绘制美观箱线图
p_box <- ggplot(expr_long, aes(x=Sample, y=logFPKM, fill=Group)) +
  geom_boxplot(outlier.size=0.5, outlier.alpha=0.5, width=0.7) +
  scale_fill_brewer(palette="Set1") +
  labs(
    x = "Sample",
    y = expression(log[2]*"(FPKM+1)"),
    title = "Expression Distribution per Sample (Boxplot)"
  ) +
  theme_classic(base_size = 16) +
  theme(
    axis.text.x = element_text(angle=45, hjust=1, vjust=1, size=12),
    plot.title = element_text(hjust=0.5, face="bold"),
    legend.title = element_text(face="bold"),
    legend.position = "right"
  )

print(p_box)