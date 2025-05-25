library(openxlsx)#读取.xlsx文件
library(ggplot2)#柱状图和点状图
library(stringr)#基因ID转换
library(enrichplot)#GO,KEGG,GSEA
library(clusterProfiler)#GO,KEGG,GSEA
library(GOplot)#弦图，弦表图，系统聚类图
library(DOSE)
library(ggnewscale)
library(topGO)#绘制通路网络图
library(circlize)#绘制富集分析圈图
library(ComplexHeatmap)#绘制图例
########GSEA#######
install.packages("C:/Users/zdduc/OneDrive - zju.edu.cn/max_plank_insititu/1.Experiment/0.database/strongyloides_ratti.PRJEB125.WBPS19.CDS_transcripts.fa/data/org.Sratti.eg.db", repos=NULL, type="sources")
library(org.Hnematodes.eg.db)
GO_database <- org.Hnematodes.eg.db
KEGG_database <- "C:/Users/zdduc/OneDrive - zju.edu.cn/zju/投稿文章/dzd/manuscript2/A 2025 new/自建库/KEGGdatabase_hc.Rdata"

info <- read.csv(file = 'HEME-DEG.csv',header = T)

GSEA_input <- info$log2FC
names(GSEA_input) = info$NAME
GSEA_input = sort(GSEA_input, decreasing = TRUE)
GSEA_KEGG <- GSEA(
       geneList = GSEA_input,
       TERM2GENE = pathway2gene,
       TERM2NAME = pathway2name,
       pvalueCutoff = 1 )

GO_TERM2GENE <- read.delim("C:/Users/zdduc/OneDrive - zju.edu.cn/zju/投稿文章/dzd/manuscript2/A 2025 new/自建库/go2gene.txt", header = TRUE, sep = "\t")
GO_TERM2NAME <- read.delim("C:/Users/zdduc/OneDrive - zju.edu.cn/zju/投稿文章/dzd/manuscript2/A 2025 new/自建库/go2name.txt", header = TRUE, sep = "\t")

GSEA_GO <- GSEA(
  geneList = GSEA_input,
  TERM2GENE = GO_TERM2GENE,
  TERM2NAME = GO_TERM2NAME,
  pvalueCutoff = 1
)

ridgeplot(GSEA_KEGG) 
gseaplot2(GSEA_KEGG,1)
gseaplot2(GSEA_KEGG,1:26) ## 30是根据ridgeplot中有30个富集通路得到的
ridgeplot(GSEA_GO) 
gseaplot2(GSEA_GO,1)
gseaplot2(GSEA_GO,1:29)