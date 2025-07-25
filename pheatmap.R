
 getwd()

 setwd("C:/Users/zdduc/OneDrive - zju.edu.cn/teacher ma/生信/R数据处理/SRP026668数据作图/")
 getwd()
 library(pheatmap)
############# pheatmap_1.0.12  R version 4.1.2 (2021-11-01)
  data <- read.csv(file = "HEME-DEG.csv",header = T,row.names = 1)
 View(data)
###########将data第一列作为名称
 rownames(data) = data[,1]
###########删除第一列
 data = data [,-7]
 data = data [,-7] 
 data = data [,-7]
 data = data [-61,]
 data = data [-60,]
 data = data [-59,]
#################将（数据框）data改为矩阵格式
  data <- data.matrix(data)
#################这俩没明白那不一样？ 在研究把 一个不行试试另一个
 data <- as.matrix(data)
#################绘制热图
 pheatmap(data,color = blues9,cluster_rows = FALSE,cluster_cols = FALSE,scale = "row",border_color = NA)
#################加载绘制颜色包
 library(RColorBrewer)

 pheatmap(data,color = brewer.pal(3,"Reds"),cluster_rows = TRUE,cluster_cols = FALSE,scale = "row",border_color = NA)
 pheatmap(data,color = brewer.pal(9,"Reds"),cluster_rows = TRUE,cluster_cols = FALSE,scale = "row",border_color = NA)
 pheatmap(data,color = brewer.pal(3,"Reds"),cluster_rows = TRUE,cluster_cols = FALSE,scale = "none",border_color = NA)
 pheatmap(data,color = brewer.pal(3,"Reds"),cluster_rows = TRUE,cluster_cols = FALSE,scale = "none",border_color = NA,fontsize = 9)
 pheatmap(data,color = brewer.pal(3,"Reds"),cluster_rows = TRUE,cluster_cols = FALSE,scale = "none",border_color = NA)
 pheatmap(data,color = brewer.pal(3,"Reds"),cluster_rows = TRUE,cluster_cols = FALSE,scale = "none",border_color = NA,fontsize = 9)
 pheatmap(data,color = brewer.pal(3,"Reds"),cluster_rows = TRUE,cluster_cols = FALSE,scale = "none",border_color = NA,fontsize = 6)
 pheatmap(data,color = brewer.pal(3,"Reds"),cluster_rows = TRUE,cluster_cols = FALSE,scale = "none",border_color = NA,fontsize = 6)
 pheatmap(data,color = brewer.pal(3,"Reds"),cluster_rows = TRUE,cluster_cols = FALSE,scale = "row",border_color = NA,fontsize = 6)
 pheatmap(data,color = brewer.pal(3,"Reds"),cluster_rows = TRUE,cluster_cols = FALSE,scale = "row",border_color = NA,fontsize = 6)
 pheatmap(data,color = brewer.pal(3,"Reds"),cluster_rows = FALSE,cluster_cols = FALSE,scale = "row",border_color = NA,fontsize = 6)
##################编制颜色函数
 b2p1 <- colorRampPalette(c("blue", "red"))
 pheatmap(data,color = b2p1(9),cluster_rows = FALSE,cluster_cols = FALSE,scale = "row",border_color = NA,fontsize = 6)
 b2p1 <- colorRampPalette(c("#00EEEE", "grey","#FF6A6A"))
###################以这个9为模板  这个好看 
 pheatmap(data,color = b2p1(9),cluster_rows = FALSE,cluster_cols = FALSE,scale = "row",border_color = NA,fontsize = 6)
 pheatmap(data,color = b2p1(6),cluster_rows = FALSE,cluster_cols = FALSE,scale = "colmun",border_color = NA,fontsize = 6)
 pheatmap(data,color = b2p1(6),cluster_rows = FALSE,cluster_cols = FALSE,scale = "column",border_color = NA,fontsize = 6)
 pheatmap(data,color = b2p1(6),cluster_rows = FALSE,cluster_cols = TRUE,scale = "column",border_color = NA,fontsize = 6)
 pheatmap(data,color = b2p1(6),cluster_rows = TRUE,cluster_cols = FALSE,scale = "column",border_color = NA,fontsize = 6)
 pheatmap(data,color = b2p1(6),cluster_rows = FALSE,cluster_cols = FALSE,scale = "none",border_color = NA,fontsize = 6)
####################### row normalisation +row  cluster 模板及颜色
 pheatmap(data,color = b2p1(9),cluster_rows = TRUE,cluster_cols = FALSE,scale = "row",border_color = NA,fontsize = 25,cellwidth = 70,cellhight = 4,show_rownames = F)
######################相同的颜色方程跑出来的才是同样的长宽高brewer.pal(3,"Reds") 的颜色代码为"#FEE0D2" "#FC9272" "#DE2D26",放到b2p2 <- colorRampPalette(c("#FEE0D2","#FC9272","#DE2D26"))中可以使用同样的方程
#####颜色使用不同 会导致图形形状不同，要使用同一种颜色公式
######################  no normalisation +row  cluster 模板及颜色
 pheatmap(a,color = brewer.pal(3,"Reds"),cluster_rows = TRUE,cluster_cols = FALSE,scale = "none",border_color = NA,fontsize = 6,cellwidth = 35,cellhight = 7)
######################  保存为pdf格式，背景为透明。
#######创建pdf
 pdf("DAheatmap3.pdf")
#######绘制图片
 pheatmap(d,color = b2p1(16),cluster_rows = TRUE,cluster_cols = FALSE,scale = "row",border_color = NA,fontsize = 6,cellwidth = 50,cellheight = 7,treeheight_row = 50)
#######结束进程
 dev.off()
#######将热图的名字写在左边
library(pheatmap)
library(grid)
library(gtable)

# Modified pheatmap:::heatmap_motor
heatmap_motor <- function (matrix, border_color, cellwidth, cellheight, tree_col, 
    tree_row, treeheight_col, treeheight_row, filename, width, 
    height, breaks, color, legend, annotation_row, annotation_col, 
    annotation_colors, annotation_legend, annotation_names_row, 
    annotation_names_col, main, fontsize, fontsize_row, fontsize_col, 
    hjust_col, vjust_col, angle_col, fmat, fontsize_number, number_color, 
    gaps_col, gaps_row, labels_row, labels_col, ...) 
{
    lo = pheatmap:::lo(coln = labels_col, rown = labels_row, nrow = nrow(matrix), 
        ncol = ncol(matrix), cellwidth = cellwidth, cellheight = cellheight, 
        treeheight_col = treeheight_col, treeheight_row = treeheight_row, 
        legend = legend, annotation_col = annotation_col, annotation_row = annotation_row, 
        annotation_colors = annotation_colors, annotation_legend = annotation_legend, 
        annotation_names_row = annotation_names_row, annotation_names_col = annotation_names_col, 
        main = main, fontsize = fontsize, fontsize_row = fontsize_row, 
        fontsize_col = fontsize_col, angle_col = angle_col, gaps_row = gaps_row, 
        gaps_col = gaps_col, ...)
    res = lo$gt
    mindim = lo$mindim
    if (!is.na(filename)) {
        if (is.na(height)) {
            height = convertHeight(gtable_height(res), "inches", valueOnly = T)
        }
        if (is.na(width)) {
            width = convertWidth(gtable_width(res), "inches", valueOnly = T)
        }
        r = regexpr("\\.[a-zA-Z]*$", filename)
        if (r == -1) 
            stop("Improper filename")
        ending = substr(filename, r + 1, r + attr(r, "match.length"))
        f = switch(ending, pdf = function(x, ...) pdf(x, ...), 
            png = function(x, ...) png(x, units = "in", res = 300, 
                ...), jpeg = function(x, ...) jpeg(x, units = "in", 
                res = 300, ...), jpg = function(x, ...) jpeg(x, 
                units = "in", res = 300, ...), tiff = function(x, 
                ...) tiff(x, units = "in", res = 300, compression = "lzw", 
                ...), bmp = function(x, ...) bmp(x, units = "in", 
                res = 300, ...), stop("File type should be: pdf, png, bmp, jpg, tiff"))
        f(filename, height = height, width = width)
        gt = heatmap_motor(matrix, cellwidth = cellwidth, cellheight = cellheight, 
            border_color = border_color, tree_col = tree_col, 
            tree_row = tree_row, treeheight_col = treeheight_col, 
            treeheight_row = treeheight_row, breaks = breaks, 
            color = color, legend = legend, annotation_col = annotation_col, 
            annotation_row = annotation_row, annotation_colors = annotation_colors, 
            annotation_legend = annotation_legend, annotation_names_row = annotation_names_row, 
            annotation_names_col = annotation_names_col, filename = NA, 
            main = main, fontsize = fontsize, fontsize_row = fontsize_row, 
            fontsize_col = fontsize_col, hjust_col = hjust_col, 
            vjust_col = vjust_col, angle_col = angle_col, fmat = fmat, 
            fontsize_number = fontsize_number, number_color = number_color, 
            labels_row = labels_row, labels_col = labels_col, 
            gaps_col = gaps_col, gaps_row = gaps_row, ...)
        grid.draw(gt)
        dev.off()
        return(gt)
    }
    if (mindim < 3) 
        border_color = NA
    if (!is.na(main)) {
        elem = pheatmap:::draw_main(main, fontsize = 1.3 * fontsize, ...)
        res = gtable_add_grob(res, elem, t = 1, l = 3, name = "main", 
            clip = "off")
    }
    if (!pheatmap:::is.na2(tree_col) & treeheight_col != 0) {
        elem = pheatmap:::draw_dendrogram(tree_col, gaps_col, horizontal = T)
        res = gtable_add_grob(res, elem, t = 2, l = 3, name = "col_tree")
    }
    if (!pheatmap:::is.na2(tree_row) & treeheight_row != 0) {
        elem = pheatmap:::draw_dendrogram(tree_row, gaps_row, horizontal = F)
        res = gtable_add_grob(res, elem, t = 4, l = 1, name = "row_tree")
    }
    elem = pheatmap:::draw_matrix(matrix, border_color, gaps_row, gaps_col, 
        fmat, fontsize_number, number_color)
    res = gtable_add_grob(res, elem, t = 4, l = 3, clip = "off", 
        name = "matrix")
    if (length(labels_col) != 0) {
        pars = list(labels_col, gaps = gaps_col, fontsize = fontsize_col, 
            hjust_col = hjust_col, vjust_col = vjust_col, angle_col = angle_col, 
            ...)
        elem = do.call(pheatmap:::draw_colnames, pars)
        res = gtable_add_grob(res, elem, t = 5, l = 3, clip = "off", 
            name = "col_names")
    }
    if (length(labels_row) != 0) {
        pars = list(labels_row, gaps = gaps_row, fontsize = fontsize_row, 
            ...)
        elem = do.call(pheatmap:::draw_rownames, pars)
        res = gtable_add_grob(res, elem, t = 4, l = 3, clip = "off", 
            name = "row_names")
    }
    if (!pheatmap:::is.na2(annotation_col)) {
        converted_annotation = convert_annotations(annotation_col, 
            annotation_colors)
        elem = pheatmap:::draw_annotations(converted_annotation, border_color, 
            gaps_col, fontsize, horizontal = T)
        res = gtable_add_grob(res, elem, t = 3, l = 3, clip = "off", 
            name = "col_annotation")
        if (annotation_names_col) {
            elem = pheatmap:::draw_annotation_names(annotation_col, fontsize, 
                horizontal = T)
            res = gtable_add_grob(res, elem, t = 3, l = 4, clip = "off", 
                name = "col_annotation_names")
        }
    }
    if (!pheatmap:::is.na2(annotation_row)) {
        converted_annotation = convert_annotations(annotation_row, 
            annotation_colors)
        elem = pheatmap:::draw_annotations(converted_annotation, border_color, 
            gaps_row, fontsize, horizontal = F)
        res = gtable_add_grob(res, elem, t = 4, l = 2, clip = "off", 
            name = "row_annotation")
        if (annotation_names_row) {
            elem = pheatmap:::draw_annotation_names(annotation_row, fontsize, 
                horizontal = F, hjust_col = hjust_col, vjust_col = vjust_col, 
                angle_col = angle_col)
            res = gtable_add_grob(res, elem, t = 5, l = 2, clip = "off", 
                name = "row_annotation_names")
        }
    }
    annotation = c(annotation_col[length(annotation_col):1], 
        annotation_row[length(annotation_row):1])
    annotation = annotation[unlist(lapply(annotation, function(x) !pheatmap:::is.na2(x)))]
    if (length(annotation) > 0 & annotation_legend) {
        elem = pheatmap:::draw_annotation_legend(annotation, annotation_colors, 
            border_color, fontsize = fontsize, ...)
        t = ifelse(is.null(labels_row), 4, 3)
        res = gtable_add_grob(res, elem, t = t, l = 6, b = 5, 
            clip = "off", name = "annotation_legend")
    }
    if (!pheatmap:::is.na2(legend)) {
        elem = pheatmap:::draw_legend(color, breaks, legend, fontsize = fontsize, 
            ...)
        t = ifelse(is.null(labels_row), 4, 3)
        res = gtable_add_grob(res, elem, t = t, l = 5, b = 5, 
            clip = "off", name = "legend")
    }
    return(res)
}

# Modified pheatmap:::lo    
lo <- function (rown, coln, nrow, ncol, cellheight = NA, cellwidth = NA, 
    treeheight_col, treeheight_row, legend, annotation_row, annotation_col, 
    annotation_colors, annotation_legend, annotation_names_row, 
    annotation_names_col, main, fontsize, fontsize_row, fontsize_col, 
    angle_col, gaps_row, gaps_col, ...) 
{
    if (!is.null(coln[1]) | (!pheatmap:::is.na2(annotation_row) & annotation_names_row)) {
        if (!is.null(coln[1])) {
            t = coln
        }
        else {
            t = ""
        }
        tw = strwidth(t, units = "in", cex = fontsize_col/fontsize)
        if (annotation_names_row) {
            t = c(t, colnames(annotation_row))
            tw = c(tw, strwidth(colnames(annotation_row), units = "in"))
        }
        longest_coln = which.max(tw)
        gp = list(fontsize = ifelse(longest_coln <= length(coln), 
            fontsize_col, fontsize), ...)
        coln_height = unit(1, "grobheight", textGrob(t[longest_coln], 
            rot = angle_col, gp = do.call(gpar, gp))) + unit(10, 
            "bigpts")
    }
    else {
        coln_height = unit(5, "bigpts")
    }
    if (!is.null(rown[1])) {
        t = rown
        tw = strwidth(t, units = "in", cex = fontsize_row/fontsize)
        if (annotation_names_col) {
            t = c(t, colnames(annotation_col))
            tw = c(tw, strwidth(colnames(annotation_col), units = "in"))
        }
        longest_rown = which.max(tw)
        gp = list(fontsize = ifelse(longest_rown <= length(rown), 
            fontsize_row, fontsize), ...)
        rown_width = unit(1, "grobwidth", textGrob(t[longest_rown], 
            rot = 0, gp = do.call(gpar, gp))) + unit(10, "bigpts")
    }
    else {
        rown_width = unit(5, "bigpts")
    }
    gp = list(fontsize = fontsize, ...)
    if (!pheatmap:::is.na2(legend)) {
        longest_break = which.max(nchar(names(legend)))
        longest_break = unit(1.1, "grobwidth", 
            textGrob(as.character(names(legend))[longest_break], 
            gp = do.call(gpar, gp)))
        title_length = unit(1.1, "grobwidth", textGrob("Scale", 
            gp = gpar(fontface = "bold", ...)))
        legend_width = unit(12, "bigpts") + longest_break * 1.2
        legend_width = max(title_length, legend_width)
    }
    else {
        legend_width = unit(0, "bigpts")
    }
    if (is.na(main)) {
        main_height = unit(0, "npc")
    }
    else {
        main_height = unit(1.5, "grobheight", textGrob(main, 
            gp = gpar(fontsize = 1.3 * fontsize, ...)))
    }
    textheight = unit(fontsize, "bigpts")
    if (!pheatmap:::is.na2(annotation_col)) {
        annot_col_height = ncol(annotation_col) * (textheight + 
            unit(2, "bigpts")) + unit(2, "bigpts")
        t = c(as.vector(as.matrix(annotation_col)), colnames(annotation_col))
        annot_col_legend_width = unit(1.2, "grobwidth", textGrob(t[which.max(nchar(t))], 
            gp = gpar(...))) + unit(12, "bigpts")
        if (!annotation_legend) {
            annot_col_legend_width = unit(0, "npc")
        }
    }
    else {
        annot_col_height = unit(0, "bigpts")
        annot_col_legend_width = unit(0, "bigpts")
    }
    if (!pheatmap:::is.na2(annotation_row)) {
        annot_row_width = ncol(annotation_row) * (textheight + 
            unit(2, "bigpts")) + unit(2, "bigpts")
        t = c(as.vector(as.matrix(annotation_row)), colnames(annotation_row))
        annot_row_legend_width = unit(1.2, "grobwidth", textGrob(t[which.max(nchar(t))], 
            gp = gpar(...))) + unit(12, "bigpts")
        if (!annotation_legend) {
            annot_row_legend_width = unit(0, "npc")
        }
    }
    else {
        annot_row_width = unit(0, "bigpts")
        annot_row_legend_width = unit(0, "bigpts")
    }
    annot_legend_width = max(annot_row_legend_width, annot_col_legend_width)
    treeheight_col = unit(treeheight_col, "bigpts") + unit(5, 
        "bigpts")
    treeheight_row = unit(treeheight_row, "bigpts") + unit(5, 
        "bigpts")
    if (is.na(cellwidth)) {
        mat_width = unit(1, "npc") - rown_width - legend_width - 
            treeheight_row - annot_row_width - annot_legend_width
    }
    else {
        mat_width = unit(cellwidth * ncol, "bigpts") + length(gaps_col) * 
            unit(4, "bigpts")
    }
    if (is.na(cellheight)) {
        mat_height = unit(1, "npc") - main_height - coln_height - 
            treeheight_col - annot_col_height
    }
    else {
        mat_height = unit(cellheight * nrow, "bigpts") + length(gaps_row) * 
            unit(4, "bigpts")
    }
    gt = gtable(widths = unit.c(treeheight_row, rown_width,  
        mat_width, treeheight_row, legend_width, annot_legend_width), 
        heights = unit.c(main_height, treeheight_col, annot_col_height, 
            mat_height, coln_height), vp = viewport(gp = do.call(gpar, 
            gp)))
    cw = convertWidth(mat_width - (length(gaps_col) * unit(4, 
        "bigpts")), "bigpts", valueOnly = T)/ncol
    ch = convertHeight(mat_height - (length(gaps_row) * unit(4, 
        "bigpts")), "bigpts", valueOnly = T)/nrow
    mindim = min(cw, ch)
    res = list(gt = gt, mindim = mindim)
    return(res)
}

# Modified pheatmap:::draw_rownames      
draw_rownames <- function (rown, gaps, ...) 
{
    coord = pheatmap:::find_coordinates(length(rown), gaps)
    y = unit(1, "npc") - (coord$coord - 0.5 * coord$size)
    res = textGrob(rown, x = unit(-3, "bigpts"), y = y, vjust = 0.5, 
        hjust = 1, gp = gpar(...))
    return(res)
}

assignInNamespace(x="draw_rownames", value=draw_rownames, ns="pheatmap")
assignInNamespace(x="lo", value=lo, ns="pheatmap")
assignInNamespace(x="heatmap_motor", value=heatmap_motor, ns="pheatmap") 

h <- pheatmap(data.random,
        cluster_rows=FALSE,
        cluster_cols=FALSE,
        legend=TRUE,
        color=my_palette,
        breaks=colors,
        show_rownames=TRUE,
        show_colnames=TRUE
)