rm(list = ls())  
options(stringsAsFactors = F)
Sys.setenv(LANGUAGE = "en")
library(WGCNA)
library(FactoMineR)
library(factoextra)  
library(tidyverse) # ggplot2 stringer dplyr tidyr readr purrr  tibble forcats
library(data.table)
fpkm00 <- fread("working-table.CSV",data.table = F)
### 合并具有相同基因名的行
table(duplicated(fpkm00$gene)) #统计重复基因名的基因
gene <- fpkm00$gene ; fpkm00 <- fpkm00[,-1]
fpkm0 <- aggregate(fpkm00, by=list(gene), FUN=sum)
fpkm <- column_to_rownames(fpkm0,"Group.1")
enableWGCNAThreads(nThreads = 0.75*parallel::detectCores()) 
data <- log2(fpkm+1)
### 筛选MAD前5000的基因
keep_data <- data[order(apply(data,1,mad), decreasing = T)[1:5000],]

### 创建datTraits，包含分组、表型等信息
datTraits <- data.frame(row.names = colnames(data),group=colnames(data))
fix(datTraits)
grouptype <- data.frame(group=sort(unique(datTraits$group)),
                        groupNo=1:length(unique(datTraits$group)))
# fix(grouptype)
datTraits$groupNo = "NA"
for(i in 1:nrow(grouptype)){
  datTraits[which(datTraits$group == grouptype$group[i]),'groupNo'] <- grouptype$groupNo[i]}
datTraits
### 转置
datExpr0 <- as.data.frame(t(fpkm))
### 判断数据质量--缺失值
gsg <- goodSamplesGenes(datExpr0,verbose = 3)
gsg$allOK
if (!gsg$allOK){
  # Optionally, print the gene and sample names that were removed:
  if (sum(!gsg$goodGenes)>0)
    printFlush(paste("Removing genes:", paste(names(datExpr0)[!gsg$goodGenes],
                                              collapse = ", ")));
  if (sum(!gsg$goodSamples)>0)
    printFlush(paste("Removing samples:",
                     paste(rownames(datExpr0)[!gsg$goodSamples], collapse = ", ")));
  # Remove the offending genes and samples from the data:
  datExpr0 = datExpr0[gsg$goodSamples, gsg$goodGenes]
}
gsg <- goodSamplesGenes(datExpr0,verbose = 3)
gsg$allOK

if(T){
  #针对样本做聚类树
  sampleTree <- hclust(dist(datExpr0), method = "average")
  par(mar = c(0,5,2,0))
  plot(sampleTree, main = "Sample clustering", sub="", xlab="", cex.lab = 2, 
       cex.axis = 1, cex.main = 1,cex.lab=1)
  ## 若样本有性状、表型，可以添加对应颜色，查看是否聚类合理
  sample_colors <- numbers2colors(as.numeric(factor(datTraits$group)), 
                                  colors = rainbow(length(table(datTraits$group))), 
                                  signed = FALSE)
  ## 绘制样品的系统聚类树及对应性状
  par(mar = c(1,4,3,1),cex=0.8)
  pdf("step1_Sample dendrogram and trait.pdf",width = 8,height = 6)
  plotDendroAndColors(sampleTree, sample_colors,
                      groupLabels = "trait",
                      cex.dendroLabels = 0.8,
                      marAll = c(1, 4, 3, 1),
                      cex.rowText = 0.01,
                      main = "Sample dendrogram and trait" )
  ## Plot a line to show the cut
  # abline(h = 23500, col = "red") #根据实际情况而定
  dev.off()
}

##若存在显著离群点；剔除掉
if(F){
  clust <- cutreeStatic(sampleTree, cutHeight = 23500, minSize = 10) # cutHeight根据实际情况而定
  table(clust)
  keepSamples <- (clust==1)
  datExpr0 <- datExpr0[keepSamples, ]
  datTraits <- datTraits[keepSamples,]
  dim(datExpr0) 
}

##保存数据
datExpr <- datExpr0
nGenes <- ncol(datExpr)
nSamples <- nrow(datExpr)
save(nGenes,nSamples,datExpr,datTraits,file="step1_input.Rdata")
############################### 2.挑选最佳阈值power ###################################
R.sq_cutoff = 0.8  #设置R^2 cut-off，默认为0.85
if(T){
  # Call the network topology analysis function
  #设置power参数选择范围
  powers <- c(seq(1,20,by = 1), seq(22,30,by = 2)) 
  sft <- pickSoftThreshold(datExpr, 
                           networkType = "unsigned",
                           powerVector = powers, 
                           RsquaredCut = R.sq_cutoff,  
                           verbose = 5)
  #SFT.R.sq > 0.8 , slope ≈ -1
  pdf("step2_power-value.pdf",width = 16,height = 12)
  # Plot the results: 寻找拐点，确认最终power取值
  par(mfrow = c(1,2));
  cex1 = 0.9;
  # Scale-free topology fit index as a function of the soft-thresholding power
  plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
       xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n")
  text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
       labels=powers,cex=cex1,col="red")
  # this line corresponds to using an R^2 cut-off of h
  abline(h=R.sq_cutoff ,col="red")
  # Mean connectivity as a function of the soft-thresholding power
  plot(sft$fitIndices[,1], sft$fitIndices[,5],
       xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n")
  text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")
  abline(h=100,col="red")
  dev.off()
}

sft$powerEstimate  #查看估计的最佳power
# power = sft$powerEstimate
power = 15

if(is.na(power)){
  # 官方推荐 "signed" 或 "signed hybrid"
  # 为与原文档一致，故未修改
  type = "unsigned"
  nSamples=nrow(datExpr)
  power = ifelse(nSamples<20, ifelse(type == "unsigned", 9, 18),
                 ifelse(nSamples<30, ifelse(type == "unsigned", 8, 16),
                        ifelse(nSamples<40, ifelse(type == "unsigned", 7, 14),
                               ifelse(type == "unsigned", 6, 12))      
                 )
  )
}
save(sft, power, file = "step2_power_value.Rdata")
##################### 3.一步法构建加权共表达网络，识别基因模块 ####################
load(file = "step1_input.Rdata")
load(file = "step2_power_value.Rdata")
if(T){
  net <- blockwiseModules(
    datExpr,
    power = power,
    maxBlockSize = ncol(datExpr),
    corType = "pearson", #默认为"pearson","bicor"则更能考虑离群点的影响
    networkType = "unsigned",
    TOMType = "unsigned", 
    minModuleSize = 30,    ##越大模块越少
    mergeCutHeight = 0.25, ##越大模块越少
    numericLabels = TRUE, 
    saveTOMs = F,
    verbose = 3
  )
  table(net$colors) 
  # power: 上一步计算的软阈值
  # maxBlockSize:计算机能处理的最大模块的基因数量(默认5000),16G内存可以处理2万个，
  # 计算资源允许的情况下最好放在一个block里面。
  # corType：计算相关性的方法；可选pearson(默认)，bicor。后者更能考虑离群点的影响。
  # networkType:计算邻接矩阵时，是否考虑正负相关性；默认为"unsigned",可选"signed", "signed hybrid"
  # TOMType：计算TOM矩阵时，是否考虑正负相关性；默认为"signed",可选"unsigned"。但是根据幂律转换的邻接矩阵(权重)的非负性，所以认为这里选择"signed"也没有太多的意义。
  # numericLabels: 返回数字而不是颜色作为模块的名字，后面可以再转换为颜色
  # saveTOMs：最耗费时间的计算，可存储起来供后续使用，
  # mergeCutHeight: 合并模块的阈值，越大模块越少,一般为0.25
  # minModuleSize: 每个模块里最少放多少个基因，设定越大模块越少
  # 输出结果根据模块中基因数目的多少，降序排列，依次编号为 `1-最大模块数`。
  # **0 (grey)**表示**未**分入任何模块的基因。
}

## 模块可视化，层级聚类树展示各个模块
if(T){
  # Convert labels to colors for plotting
  moduleColors <- labels2colors(net$colors)
  table(moduleColors)
  # Plot the dendrogram and the module colors underneath
  pdf("step3_genes-modules_ClusterDendrogram.pdf",width = 16,height = 12)
  plotDendroAndColors(net$dendrograms[[1]], moduleColors[net$blockGenes[[1]]],
                      "Module colors",
                      dendroLabels = FALSE, hang = 0.03,
                      addGuide = TRUE, guideHang = 0.05)
  dev.off()
}
save(net, moduleColors, file = "step3_genes_modules.Rdata")
####################### 4.关联基因模块与表型 #####################################
rm(list = ls())  
load(file = "step1_input.Rdata")
load(file = "step2_power_value.Rdata")
load(file = "step3_genes_modules.Rdata")

## 模块与表型的相关性热图
if(T){
  datTraits$group <- as.factor(datTraits$group)
  design <- model.matrix(~0+datTraits$group)
  colnames(design) <- levels(datTraits$group) #get the group
  MES0 <- moduleEigengenes(datExpr,moduleColors)$eigengenes  #Calculate module eigengenes.
  MEs <- orderMEs(MES0)  #Put close eigenvectors next to each other
  moduleTraitCor <- cor(MEs,design,use = "p")
  moduleTraitPvalue <- corPvalueStudent(moduleTraitCor,nSamples)
  textMatrix <- paste0(signif(moduleTraitCor,2),"\n(",
                       signif(moduleTraitPvalue,1),")")
  dim(textMatrix) <- dim(moduleTraitCor)
  
  pdf("step4_Module-trait-relationship_heatmap.pdf",
      width = 2*length(colnames(design)), 
      height = 0.6*length(names(MEs)) )
  par(mar=c(5, 9, 3, 3)) #留白：下、左、上、右
  labeledHeatmap(Matrix = moduleTraitCor,
                 xLabels = colnames(design),
                 yLabels = names(MEs),
                 ySymbols = names(MEs),
                 colorLabels = F,
                 colors = blueWhiteRed(50),
                 textMatrix = textMatrix,
                 setStdMargins = F,
                 cex.text = 0.5,
                 zlim = c(-1,1), 
                 main = "Module-trait relationships")
  dev.off()
  save(design, file = "step4_design.Rdata")
}


### 模块与表型的相关性boxplot图 
if(T){
  mes_group <- merge(MEs,datTraits,by="row.names") 
  library(gplots)
  library(ggpubr)
  library(grid)
  library(gridExtra) 
  draw_ggboxplot <- function(data,Module="Module",group="group"){
    ggboxplot(data,x=group, y=Module,
              ylab = paste0(Module),
              xlab = group,
              fill = group,
              palette = "jco",
              #add="jitter",
              legend = "") +stat_compare_means()
  }
  # 批量画boxplot
  colorNames <- names(MEs)
  pdf("step4_Module-trait-relationship_boxplot.pdf", width = 7.5,height = 1.6*ncol(MEs))
  p <- lapply(colorNames,function(x) {
    draw_ggboxplot(mes_group, Module = x, group = "group")
  })
  do.call(grid.arrange,c(p,ncol=2)) #排布为每行2个
  dev.off()
}


### 基因与模块、表型的相关性散点图
#所有的模块都可以跟基因算出相关系数，所有的连续型性状也可以跟基因算出相关系数， 
#如果跟性状显著相关的基因也跟某个模块显著相关，那么这些基因可能就非常重要。

# 选择离散性状的表型
levels(datTraits$group)
choose_group <- "primed"  

if(T){
  modNames <- substring(names(MEs), 3)
  
  ### 计算模块与基因的相关性矩阵 
  ## Module Membership: 模块内基因表达与模块特征值的相关性
  geneModuleMembership <- as.data.frame(cor(datExpr, MEs, use = "p"))
  MMPvalue <- as.data.frame(corPvalueStudent(as.matrix(geneModuleMembership), nSamples))
  names(geneModuleMembership) <- paste0("MM", modNames)
  names(MMPvalue) <- paste0("p.MM", modNames)
  
  ###  计算性状与基因的相关性矩阵 
  ## Gene significance，GS：比较样本某个基因与对应表型的相关性
  ## 连续型性状
  # trait <- datTraits$groupNo  
  ## 非连续型性状，需转为0-1矩阵, 已存于design中
  trait <- as.data.frame(design[,choose_group])
  
  geneTraitSignificance <- as.data.frame(cor(datExpr,trait,use = "p"))
  GSPvalue <- as.data.frame(corPvalueStudent(as.matrix(geneTraitSignificance),nSamples))
  names(geneTraitSignificance) <- paste0("GS")
  names(GSPvalue) <- paste0("GS")
  
  ### 可视化基因与模块、表型的相关性.
  #selectModule<-c("blue","green","purple","grey")  ##可以选择自己想要的模块
  selectModule <- modNames  ## 全部模块批量作图
  pdf("step4_gene-Module-trait-significance.pdf",width=7, height=1.5*ncol(MEs))
  par(mfrow=c(ceiling(length(selectModule)/2),2)) #批量作图开始
  for(module in selectModule){
    column <- match(module,selectModule)
    print(module)
    moduleGenes <- moduleColors==module
    verboseScatterplot(abs(geneModuleMembership[moduleGenes, column]),
                       abs(geneTraitSignificance[moduleGenes, 1]),
                       xlab = paste("Module Membership in", module, "module"),
                       ylab = "Gene significance for trait",
                       main = paste("Module membership vs. gene significance\n"),
                       cex.main = 1.2, cex.lab = 1.2, cex.axis = 1.2, col = module)
  }
  dev.off()
}

#########################  5. WGCNA可视化：TOMplot  Eigengene-adjacency-heatmap ##################################
rm(list = ls())  
load(file = 'step1_input.Rdata')
load(file = "step2_power_value.Rdata")
load(file = "step3_genes_modules.Rdata")
load(file = "step4_design.Rdata")

if(T){
  TOM=TOMsimilarityFromExpr(datExpr,power=power)
  dissTOM=1-TOM
  ## draw all genes 
  if(T){
    geneTree = net$dendrograms[[1]]
    plotTOM = dissTOM^7
    diag(plotTOM)=NA
    png("step5_TOMplot_Network-heatmap.png",width = 800, height=600) 
    TOMplot(plotTOM,geneTree,moduleColors,
            col=gplots::colorpanel(250,'red',"orange",'lemonchiffon'),
            main="Network heapmap plot")
    dev.off()
  }
  ### draw selected genes to save time...just for test...
  if(F){
    nSelect =0.1*nGenes
    set.seed(123)
    select=sample(nGenes,size = nSelect)
    selectTOM = dissTOM[select,select]
    selectTree = hclust(as.dist(selectTOM),method = "average")
    selectColors = moduleColors[select]
    plotDiss=selectTOM^7
    diag(plotDiss)=NA
    pdf("step5_select_TOMplot_Network-heatmap.pdf",width=8, height=6)
    TOMplot(plotDiss,selectTree,selectColors,
            col=gplots::colorpanel(250,'red',"orange",'lemonchiffon'),
            main="Network heapmap plot of selected gene")
    dev.off()
  }
}


### 模块相关性展示 Eigengene-adjacency-heatmap
if(T){
  MEs = moduleEigengenes(datExpr,moduleColors)$eigengenes
  MET = orderMEs(MEs)
  # 若添加表型数据
  if(T){
    ## 连续型性状
    # MET = orderMEs(cbind(MEs,datTraits$groupNo))
    ## 非连续型性状，需将是否属于这个表型进行0,1数值化，已存于design中
    design
    primed = as.data.frame(design[,3])
    names(primed) = "primed"
    # Add the weight to existing module eigengenes
    MET = orderMEs(cbind(MEs, primed))
  }
  pdf("step5_module_cor_Eigengene-dendrogram.pdf",width = 8,height = 10)
  plotEigengeneNetworks(MET, setLabels="", 
                        marDendro = c(0,4,1,4),  # 留白：下右上左
                        marHeatmap = c(5,5,1,2), # 留白：下右上左
                        cex.lab = 0.8,
                        xLabelsAngle = 90)
  dev.off()
}

#################### 6. 选择感兴趣基因模块进行GO分析 ####################
rm(list = ls())  
load(file = 'step1_input.Rdata')
load(file = "step2_power_value.Rdata")
load(file = "step3_genes_modules.Rdata")
load(file = "step4_design.Rdata")

### 条件设置
OrgDb = "org.Mm.eg.db"  # "org.Mm.eg.db"  "org.Hs.eg.db"
genetype = "SYMBOL"    # "SYMBOL"   "ENSEMBL"
table(moduleColors)
choose_module <- c("black","pink","royalblue")

if(T){
  library(clusterProfiler)
  library(org.Mm.eg.db)
  library(org.Hs.eg.db)
  
  gene_module <- data.frame(gene=colnames(datExpr),
                            module=moduleColors)
  write.csv(gene_module,file = "step6_gene_moduleColors.csv",row.names = F, quote = F) 
  tmp <- bitr(gene_module$gene,fromType = genetype,  # "SYMBOL"   "ENSEMBL"
              toType = "ENTREZID",
              OrgDb = OrgDb )
  gene_module_entrz <- merge(tmp,gene_module, by.x=genetype, by.y="gene")
  
  choose_gene_module_entrz <- gene_module_entrz[gene_module_entrz$module %in% choose_module,]
  
  ###run go analysis
  formula_res <- compareCluster(
    ENTREZID~module,
    data = choose_gene_module_entrz,
    fun = "enrichGO",
    OrgDb = OrgDb,
    ont = "BP",  #One of "BP", "MF", and "CC"  or "ALL"
    pAdjustMethod = "BH",
    pvalueCutoff = 0.25,
    qvalueCutoff = 0.25
  )
  
  ###精简GO富集的结果,去冗余
  lineage1_ego <- simplify( 
    formula_res,
    cutoff=0.5,
    by="p.adjust",
    select_fun=min
  )
  save(gene_module, formula_res, lineage1_ego, file="step6_module_GO_term.Rdata")
  write.csv(lineage1_ego@compareClusterResult,
            file="step6_module_GO_term.csv")
  ### 绘制dotplot图
  dotp <- dotplot(lineage1_ego,
                  showCategory=10,
                  includeAll = TRUE, #将有overlap的结果也展示出来
                  label_format=90)
  ggsave(dotp,filename= "step6_module_GO_term.pdf", #device = cairo_pdf,
         width = 12, 
         height = 15)
}

############################### 7.感兴趣基因模块绘制热图 ######################################
rm(list = ls())  
load(file = 'step1_input.Rdata')
load(file = "step3_genes_modules.Rdata")
table(moduleColors)

module = "royalblue"
### 感兴趣模块画热图 
if(T){
  dat=datExpr[,moduleColors==module]
  library(pheatmap)
  n=t(scale(dat)) #对基因做scale，并转置表达矩阵为行为基因、列为样本形式
  # n[n>2]=2 
  # n[n< -2]= -2
  # n[1:4,1:4]
  
  group_list=datTraits$group
  ac=data.frame(g=group_list)
  rownames(ac)=colnames(n)
  
  pheatmap::pheatmap(n,
                     fontsize = 8,
                     show_colnames =T,
                     show_rownames = F,
                     cluster_cols = T,
                     annotation_col =ac,
                     width = 8,
                     height = 6,
                     angle_col=45,
                     main = paste0("module_",module,"-gene heatmap"),
                     filename = paste0("step7_module_",module,"_Gene-heatmap.pdf"))
  
}

################### 8.感兴趣模块基因导出 VisANT or cytoscape ######################
rm(list = ls())  
load(file = 'step1_input.Rdata')
load(file = "step2_power_value.Rdata")
load(file = "step3_genes_modules.Rdata")
module = "royalblue"
if(T){
  ### 提取感兴趣模块基因名
  gene <- colnames(datExpr) 
  inModule <- moduleColors==module
  modgene <- gene[inModule]
  
  ### 模块对应的基因关系矩阵
  TOM <- TOMsimilarityFromExpr(datExpr,power=power)
  modTOM <- TOM[inModule,inModule]
  dimnames(modTOM) <- list(modgene,modgene)
  
  ### 筛选连接度最大的top100基因
  nTop = 100
  IMConn = softConnectivity(datExpr[, modgene]) #计算连接度
  top = (rank(-IMConn) <= nTop) #选取连接度最大的top100
  filter_modTOM <- modTOM[top, top]
  
  # for visANT
  vis <- exportNetworkToVisANT(filter_modTOM,
                               file = paste("step8_visANTinput-",module,".txt",sep = ""),
                               weighted = T,threshold = 0)
  # for cytoscape
  cyt <- exportNetworkToCytoscape(filter_modTOM,
                                  edgeFile = paste("step8_CytoscapeInput-edges-", paste(module, collapse="-"), ".txt", sep=""),
                                  nodeFile = paste("step8_CytoscapeInput-nodes-", paste(module, collapse="-"), ".txt", sep=""),
                                  weighted = TRUE,
                                  threshold = 0.15,  #weighted权重筛选阈值，可调整
                                  nodeNames = modgene[top], 
                                  nodeAttr = moduleColors[inModule][top])
}