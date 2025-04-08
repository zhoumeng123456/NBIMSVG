install.packages("ComplexUpset")
install.packages("/Users/zhoumeng/Downloads/UpSetR_1.4.0.tar.gz",repos = NULL, type = "source")

library(ComplexUpset)
library(UpSetR)
library(ggplot2)
library(Seurat)
library(scatterpie)
library(patchwork)
library(RColorBrewer)
library(ggplot2)
library(gridExtra)
library(grid)
library(reshape2)
library(ggplot2)
library(pheatmap)
library(cowplot)
library(magick)
library(ggVennDiagram)





meta_process =function(position1,count1){
  meta=cbind(position1,t(count1))
  colnames(meta)[c(1:2)] = c('x','y')
  meta=as.data.frame(meta)
  meta$x <- (meta$x - min(meta$x))/(max(meta$x) - min(meta$x))
  meta$y <- (meta$y - min(meta$y))/(max(meta$y) - min(meta$y))
  return(meta)
}

pattern_plot2 <- function(pltdat, igene, xy = T, main = F, titlesize = 2, 
                          pointsize = 3, xpand = 0, ypand = 1, title = NULL) {
  if (!xy) {
    xy <- matrix(as.numeric(do.call(rbind, strsplit(as.character(pltdat[, 
                                                                        1]), split = "x"))), ncol = 2)
    rownames(xy) <- as.character(pltdat[, 1])
    colnames(xy) <- c("x", "y")
    pd <- cbind.data.frame(xy, pltdat[, 2:ncol(pltdat)])
  } else {
    pd <- pltdat
  }
  
  # pal <- colorRampPalette(c('seagreen1','orangered')) pal <-
  # colorRampPalette(c('#00274c','#ffcb05')) pal <-
  # colorRampPalette(c('deepskyblue','goldenrod1')) pal <-
  # colorRampPalette(c('deepskyblue','deeppink'))
  pal <- colorRampPalette(c("mediumseagreen", "lightyellow2", "deeppink"))
  gpt <- ggplot(pd, aes(x = x, y = y, color = pd[, igene + 2])) + geom_point(size = pointsize) + 
    # scale_color_gradientn(colours=pal(5))+
    scale_color_gradientn(colours = pal(5)) + scale_x_discrete(expand = c(xpand, 
                                                                          ypand)) + scale_y_discrete(expand = c(xpand, ypand)) + coord_equal() + 
    # labs(title = colnames(pd)[igene+2], x = NULL, y = NULL)+
    #theme_bw()
  theme_void()
  
  if (main) {
    if (is.null(title)) {
      title = colnames(pd)[igene + 2]
    }
    out = gpt + labs(title = title, x = NULL, y = NULL) + theme(legend.position = "none", 
                                                                plot.title = element_text(hjust = 0.5, size = rel(titlesize)))
  } else {
    out = gpt + labs(title = NULL, x = NULL, y = NULL) + theme(legend.position = "none")
  }
  return(out)
}

######upset图############upset图############upset图############upset图############upset图############upset图############
######upset图############upset图############upset图############upset图############upset图############upset图############

##虽然是16，但是导入之后还是15
load("/Users/zhoumeng/Desktop/数据/data1_LIBD_svgene_list_1216_filter.RData")

##并集
uni_set=list("nnSVG-Union"=data1_LIBD_svgene_list_1215_filter[[1]][[5]],"SPARK-Union"=data1_LIBD_svgene_list_1215_filter[[3]][[5]],
     "SPARKX-Union"=data1_LIBD_svgene_list_1215_filter[[2]][[5]],"spVC-Union"=data1_LIBD_svgene_list_1215_filter[[4]][[5]],
     "HEARTSVG-Union"=data1_LIBD_svgene_list_1215_filter[[5]][[5]],Proposed=data1_LIBD_svgene_list_1215_filter[[6]]
     )
upset_uni_set=fromList(uni_set)
##交集
inter_set=list("nnSVG-Inter"=data1_LIBD_svgene_list_1215_filter[[1]][[6]],"SPARK-Inter"=data1_LIBD_svgene_list_1215_filter[[3]][[6]],
               "SPARKX-Inter"=data1_LIBD_svgene_list_1215_filter[[2]][[6]],"spVC-Inter"=data1_LIBD_svgene_list_1215_filter[[4]][[6]],
               "HEARTSVG-Inter"=data1_LIBD_svgene_list_1215_filter[[5]][[6]],Proposed=data1_LIBD_svgene_list_1215_filter[[6]]
)
upset_inter_set=fromList(inter_set)
#整合
integ_set=list("nnSVG-PASTE"=data1_LIBD_svgene_list_1215_filter[[1]][[7]],"SPARK-PASTE"=data1_LIBD_svgene_list_1215_filter[[3]][[7]],
               "SPARKX-PASTE"=data1_LIBD_svgene_list_1215_filter[[2]][[7]],"spVC-PASTE"=data1_LIBD_svgene_list_1215_filter[[4]][[7]],
               "HEARTSVG-PASTE"=data1_LIBD_svgene_list_1215_filter[[5]][[7]],Proposed=data1_LIBD_svgene_list_1215_filter[[6]],DESpace=data1_LIBD_svgene_list_1215_filter[[7]]
)
upset_integ_set=fromList(integ_set)


png("/Users/zhoumeng/Desktop/upset_union.png", width = 8, height = 3, units = "in", res = 300)
upset(upset_uni_set,
      nsets = 20, # 可视化数据集数量
      mainbar.y.max=550,##y轴的上限
      nintersects= 32, # 显示前多少个交集
      main.bar.color = "#91CAE8", # 柱状图颜色
      matrix.color="black", # 集合点的颜色
      sets.bar.color= "#F48892", # 条形图条形的颜色，横着的这个
      set_size.show = T, # 是否在条形图上显示集合大小
      mb.ratio = c(0.6, 0.4), # 矩阵图与主柱图之比
      order.by="freq",
      set_size.scale_max = 3000
)
dev.off()

png("/Users/zhoumeng/Desktop/upset_inter.png", width = 8, height = 3, units = "in", res = 300)
upset(upset_inter_set,
      nsets = 20, # 可视化数据集数量
      mainbar.y.max=1300,##y轴的上限
      nintersects= 32, # 显示前多少个交集
      main.bar.color = "#91CAE8", # 柱状图颜色
      matrix.color="black", # 集合点的颜色
      sets.bar.color= "#F48892", # 条形图条形的颜色，横着的这个
      set_size.show = T, # 是否在条形图上显示集合大小
      mb.ratio = c(0.6, 0.4), # 矩阵图与主柱图之比
      order.by="freq",
      set_size.scale_max = 2000
)
dev.off()

png("/Users/zhoumeng/Desktop/upset_paste.png", width = 8, height = 3, units = "in", res = 300)
upset(upset_integ_set,
      nsets = 20, # 可视化数据集数量
      mainbar.y.max=2200,##y轴的上限
      nintersects= 32, # 显示前多少个交集
      main.bar.color = "#91CAE8", # 柱状图颜色
      matrix.color="black", # 集合点的颜色
      sets.bar.color= "#F48892", # 条形图条形的颜色，横着的这个
      set_size.show = T, # 是否在条形图上显示集合大小
      mb.ratio = c(0.6, 0.4), # 矩阵图与主柱图之比
      order.by="freq",
      set_size.scale_max = 6000
)
dev.off()


img_paths <- c(
  "/Users/zhoumeng/Desktop/upset_union.png",
  "/Users/zhoumeng/Desktop/upset_inter.png",
  "/Users/zhoumeng/Desktop/upset_paste.png"
)

# 读取图片并添加标签
images <- lapply(1:3, function(i) {
  img <- image_read(img_paths[i])
  # 添加标签（位置可调）
  image_annotate(img, 
                 text = LETTERS[i], 
                 size = 50, 
                 color = "black",
                 gravity = "northwest",  # 定位左上角
                 location = "+20+20")    # 偏移量调整
})

# 垂直合并图片
combined <- image_append(do.call(c, images), stack = TRUE)
print(combined)
image_write(combined, "/Users/zhoumeng/Desktop/upset_combined.png")



##########细胞类别########################细胞类别########################细胞类别########################细胞类别##############
##########细胞类别########################细胞类别########################细胞类别########################细胞类别##############

##获得一个矩阵meta，前两列是坐标，后面是细胞类别的分布
#先导入坐标
position4=as.matrix(read.csv("/Users/zhoumeng/Downloads/matrix4_postion_126.csv",row.names = 1))
calpha4=as.matrix(read.csv("/Users/zhoumeng/Downloads/matrix4_celltype_126.csv",row.names = 1))

position3=as.matrix(read.csv("/Users/zhoumeng/Downloads/matrix3_postion_126.csv",row.names = 1))
calpha3=as.matrix(read.csv("/Users/zhoumeng/Downloads/matrix3_celltype_126.csv",row.names = 1))

position2=as.matrix(read.csv("/Users/zhoumeng/Downloads/matrix2_postion_126.csv",row.names = 1))
calpha2=as.matrix(read.csv("/Users/zhoumeng/Downloads/matrix2_celltype_126.csv",row.names = 1))

position1=as.matrix(read.csv("/Users/zhoumeng/Downloads/matrix1_postion_126.csv",row.names = 1))
calpha1=as.matrix(read.csv("/Users/zhoumeng/Downloads/matrix1_celltype_126.csv",row.names = 1))

plot_scatterpie <- function(position, calpha, pie_scale = 0.4) {
  # col_manual <- c(
  #   "Inhib" = "#FFDDDD",
  #   "Oligo" = "#FFFACD",
  #   "OPC" = "#32CD32",
  #   "Excit" = "yellow",
  #   "MicroOligo" ="#B0E0E6" ,
  #   "Astro" = "#8A2BE2",
  #   "EndoMural" = "#FF1493"
  # )
  
  col_manual <- c(
    "Inhib" = "#FFDDDD",
    "Oligo" = "#B0E0E6",
    "OPC" ="#FFFACD" ,
    "Excit" ="#32CD32" ,
    "MicroOligo" = "#FFEB99",
    "Astro" = "#8A2BE2",
    "EndoMural" = "#FF1493"
  )
  focus_types <- c("Astro", "EndoMural", "OPC")
  calpha[, focus_types] <- calpha[, focus_types] * 2  # 关键类型放大
  calpha <- calpha / rowSums(calpha)
  
  oligo_mask <- calpha[, "Oligo"] < 0.75
  oligo_mask1 <- calpha[, "Oligo"] > 0.75
  calpha[oligo_mask, "Oligo"] <- calpha[oligo_mask, "Oligo"] * 0.6  # 缩小到 0.8 倍
  calpha[oligo_mask1, "MicroOligo"] <- calpha[oligo_mask1, "MicroOligo"] * 0.6  # 缩小到 0.8 倍
  calpha[oligo_mask1, "Oligo"] <- calpha[oligo_mask1, "Oligo"] * 1.2  # 缩小到 0.8 倍
  
  calpha <- calpha / rowSums(calpha)
  
  meta <- cbind(position, calpha)
  colnames(meta)[1:2] <- c('x', 'y')
  meta <- as.data.frame(meta)
  
  # 归一化坐标
  meta$x <- (meta$x - min(meta$x)) / (max(meta$x) - min(meta$x))
  meta$y <- (meta$y - min(meta$y)) / (max(meta$y) - min(meta$y))
  
  ggplot() +
    geom_scatterpie(
      data = meta, 
      aes(x = x, y = y),
      cols = colnames(meta)[3:9], 
      color = NA, 
      pie_scale = pie_scale
    ) +
    coord_fixed(ratio = 1) +
    scale_fill_manual(values = col_manual) +
    theme_bw() +
    theme(
      panel.grid = element_blank(),
      # 移除坐标轴相关元素
      axis.line = element_blank(),
      axis.text = element_blank(),
      axis.ticks = element_blank(),
      axis.title = element_blank(),
      panel.border = element_blank(),  # 移除绘图区灰色边框
      legend.text = element_text(size = 14),  # 增加图例文字大小
      plot.tag = element_text(size = 18,face="bold")
    ) +
    labs(x = NULL, y = NULL)  # 移除坐标轴标签
}

enhance_contrast <- function(calpha, power=3) {
  # 使用幂变换增强对比度
  calpha_trans <- calpha^power
  # 重新归一化保持比例总和为1
  calpha_trans / rowSums(calpha_trans)
}

# 使用示例
calpha1 <- enhance_contrast(calpha1, power=4)  # 调节power参数控制增强力度
calpha2 <- enhance_contrast(calpha2, power=4)  # 调节power参数控制增强力度
calpha3 <- enhance_contrast(calpha3, power=4)  # 调节power参数控制增强力度
calpha4 <- enhance_contrast(calpha4, power=4)  # 调节power参数控制增强力度


# 使用不同的 position 和 calpha 绘制四个图
p1 <- plot_scatterpie(position1, calpha1)
p2 <- plot_scatterpie(position2, calpha2)
p3 <- plot_scatterpie(position3, calpha3)
p4 <- plot_scatterpie(position4, calpha4)

print(p1)

combined_plot <- (p1 | p2 | p3 | p4) +
  plot_layout(guides = "collect") +  # 收集所有子图的图例
  plot_annotation(
    title = NULL,
    subtitle = NULL,
    caption = NULL,
    tag_levels = 'A',
    tag_prefix = "",
    tag_sep = "",
    theme = theme(
      plot.tag = element_text(size = 18),  # 放大标签字体
      legend.position = "bottom",  
      legend.box = "horizontal",
      legend.direction = "horizontal",  # 图例水平排列
      legend.text = element_text(size = 14),  # 增大图例文字大小
      legend.key.size = unit(1.5, "cm")  # 增大图例柱子的大小
    )
  ) &
  guides(fill = guide_legend(nrow = 1, title = NULL, keywidth = unit(1.5, "cm")))

print(combined_plot)

ggsave("/Users/zhoumeng/Desktop/cell_type0406.png", plot = combined_plot, width = 16, height =5, dpi = 300)







##绘制每一张切片的
meta1=meta_process(position1,calpha1)
p1=pattern_plot2(meta1,1,xpand=-2,ypand = 2,main = TRUE)






######绘制基因表达##############绘制基因表达##############绘制基因表达##############绘制基因表达##############绘制基因表达########
######绘制基因表达##############绘制基因表达##############绘制基因表达##############绘制基因表达##############绘制基因表达########


###目标1:绘制与细胞类别差异相似的空间表达sv基因
###目标1:绘制与细胞类别差异相似的空间表达sv基因
###目标1:绘制与细胞类别差异相似的空间表达sv基因
#1.找到所有方法共同的基因，在并集里找,因为spv考虑到了细胞类别，所以不考虑进方法
intersv_uni_set=Reduce(intersect,list(nnsvg_uni=data1_LIBD_svgene_list_1215_filter[[1]][[5]],sparkX_uni=data1_LIBD_svgene_list_1215_filter[[2]][[5]],
             spark_uni=data1_LIBD_svgene_list_1215_filter[[3]][[5]],spv_inter=data1_LIBD_svgene_list_1215_filter[[4]][[5]],
             heartsvg_uni=data1_LIBD_svgene_list_1215_filter[[5]][[5]]
             ))
extraour_intersv_uni_set=setdiff(intersv_uni_set, data1_LIBD_svgene_list_1215_filter[[6]])###53个

#2.画出这53个基因中符合细胞类别的
position1=as.matrix(read.csv("/Users/zhoumeng/Downloads/matrix1_postion_105.csv",row.names = 1))
count1=as.matrix(read.csv("/Users/zhoumeng/Downloads/matrix1_count_105.csv",row.names = 1))
position2=as.matrix(read.csv("/Users/zhoumeng/Downloads/matrix2_postion_105.csv",row.names = 1))
count2=as.matrix(read.csv("/Users/zhoumeng/Downloads/matrix2_count_105.csv",row.names = 1))
position3=as.matrix(read.csv("/Users/zhoumeng/Downloads/matrix3_postion_105.csv",row.names = 1))
count3=as.matrix(read.csv("/Users/zhoumeng/Downloads/matrix3_count_105.csv",row.names = 1))
position4=as.matrix(read.csv("/Users/zhoumeng/Downloads/matrix4_postion_105.csv",row.names = 1))
count4=as.matrix(read.csv("/Users/zhoumeng/Downloads/matrix4_count_105.csv",row.names = 1))

#3.筛选出这53个基因
count1=count1[extraour_intersv_uni_set,]
count2=count2[extraour_intersv_uni_set,]
count3=count3[extraour_intersv_uni_set,]
count4=count4[extraour_intersv_uni_set,]

#绘制
meta1=meta_process(position1,count1)
meta2=meta_process(position2,count2)
meta3=meta_process(position3,count3)
meta4=meta_process(position4,count4)

#绘图保存
for(genenum in seq_along(rownames(count1))){
    p1=pattern_plot2(meta1,genenum,xpand=-2,ypand = 2,main = TRUE)##用于绘制
    p2=pattern_plot2(meta2,genenum,xpand=-2,ypand = 2,main = TRUE)
    p3=pattern_plot2(meta3,genenum,xpand=-2,ypand = 2,main = TRUE)
    p4=pattern_plot2(meta4,genenum,xpand=-2,ypand = 2,main = TRUE)
    combined_plot <- arrangeGrob(p1, p2, p3, p4, ncol = 2, nrow = 2)
    ggsave(paste0("/Users/zhoumeng/Desktop/extraour_intersv_uni_set/",genenum,".png"), plot = combined_plot, width = 8, height = 8, dpi = 300)
}


#获取基因表达
#获取基因表达
#获取基因表达
#获取基因表达
#获取基因表达
position1=as.matrix(read.csv("/Users/zhoumeng/Downloads/matrix1_postion_105.csv",row.names = 1))
count1=as.matrix(read.csv("/Users/zhoumeng/Downloads/matrix1_count_105.csv",row.names = 1))
position2=as.matrix(read.csv("/Users/zhoumeng/Downloads/matrix2_postion_105.csv",row.names = 1))
count2=as.matrix(read.csv("/Users/zhoumeng/Downloads/matrix2_count_105.csv",row.names = 1))
position3=as.matrix(read.csv("/Users/zhoumeng/Downloads/matrix3_postion_105.csv",row.names = 1))
count3=as.matrix(read.csv("/Users/zhoumeng/Downloads/matrix3_count_105.csv",row.names = 1))
position4=as.matrix(read.csv("/Users/zhoumeng/Downloads/matrix4_postion_105.csv",row.names = 1))
count4=as.matrix(read.csv("/Users/zhoumeng/Downloads/matrix4_count_105.csv",row.names = 1))

#3.筛选
count1=count1[unique_intersv_our,]
count2=count2[unique_intersv_our,]
count3=count3[unique_intersv_our,]
count4=count4[unique_intersv_our,]

#绘制
meta1=meta_process(position1,count1)
meta2=meta_process(position2,count2)
meta3=meta_process(position3,count3)
meta4=meta_process(position4,count4)

#seq_along(rownames(count1)
for(genenum in seq_along(rownames(count1))){
  p1=pattern_plot2(meta1,genenum,xpand=-2,ypand = 2,main = TRUE)##用于绘制
  p2=pattern_plot2(meta2,genenum,xpand=-2,ypand = 2,main = TRUE)
  p3=pattern_plot2(meta3,genenum,xpand=-2,ypand = 2,main = TRUE)
  p4=pattern_plot2(meta4,genenum,xpand=-2,ypand = 2,main = TRUE)
  combined_plot <- arrangeGrob(p1, p2, p3, p4, ncol = 2, nrow = 2)
  ggsave(paste0("/Users/zhoumeng/Desktop/unique_intersv_our_nodespace/",genenum,".png"), plot = combined_plot, width = 8, height = 8, dpi = 300)
}




#开始绘制
#开始绘制
#开始绘制
###绘制所有方法的交集(除了我们自己的方法)
###绘制所有方法的交集
###绘制所有方法的交集
inter_set=list(nnsvg_inter=data1_LIBD_svgene_list_1215_filter[[1]][[6]],spark_inter=data1_LIBD_svgene_list_1215_filter[[3]][[6]],
               sparkx_inter=data1_LIBD_svgene_list_1215_filter[[2]][[6]],spv_inter=data1_LIBD_svgene_list_1215_filter[[4]][[6]],
               heartsvg_inter=data1_LIBD_svgene_list_1215_filter[[5]][[6]]
               )##即四张切片都选出来了
all_other_inter_sv=Reduce(intersect,inter_set)##共有30个

all_inter_sv=Reduce(intersect,list(data1_LIBD_svgene_list_1215_filter[[6]],all_other_inter_sv))#25个
remain_all_inter_sv=setdiff(all_other_inter_sv,all_inter_sv)#5个

position1=as.matrix(read.csv("/Users/zhoumeng/Downloads/matrix1_postion_126.csv",row.names = 1))
count1=as.matrix(read.csv("/Users/zhoumeng/Downloads/matrix1_count_126.csv",row.names = 1))
position2=as.matrix(read.csv("/Users/zhoumeng/Downloads/matrix2_postion_126.csv",row.names = 1))
count2=as.matrix(read.csv("/Users/zhoumeng/Downloads/matrix2_count_126.csv",row.names = 1))
position3=as.matrix(read.csv("/Users/zhoumeng/Downloads/matrix3_postion_126.csv",row.names = 1))
count3=as.matrix(read.csv("/Users/zhoumeng/Downloads/matrix3_count_126.csv",row.names = 1))
position4=as.matrix(read.csv("/Users/zhoumeng/Downloads/matrix4_postion_126.csv",row.names = 1))
count4=as.matrix(read.csv("/Users/zhoumeng/Downloads/matrix4_count_126.csv",row.names = 1))

#3.筛选
count1=count1[remain_all_inter_sv,]
count2=count2[remain_all_inter_sv,]
count3=count3[remain_all_inter_sv,]
count4=count4[remain_all_inter_sv,]

#绘制
meta1=meta_process(position1,count1)
meta2=meta_process(position2,count2)
meta3=meta_process(position3,count3)
meta4=meta_process(position4,count4)

for(genenum in seq_along(rownames(count1))){
  p1=pattern_plot2(meta1,genenum,xpand=-2,ypand = 2,main = FALSE)##用于绘制
  p2=pattern_plot2(meta2,genenum,xpand=-2,ypand = 2,main = FALSE)
  p3=pattern_plot2(meta3,genenum,xpand=-2,ypand = 2,main = FALSE)
  p4=pattern_plot2(meta4,genenum,xpand=-2,ypand = 2,main = FALSE)
  combined_plot <- arrangeGrob(p1, p2, p3, p4, ncol = 2, nrow = 2)
  ggsave(paste0("/Users/zhoumeng/Desktop/图片结果/1214以后/all_sv_identify/remain5/",rownames(count1)[genenum],".png"), plot = combined_plot, width = 8, height = 8, dpi = 300)
}


#目标2:绘制我们方法找到的独特的空间表达sv基因，用于和刚才的做对比,最好是一个信号强，某些信号弱
#目标2:绘制我们方法找到的独特的空间表达sv基因，出现在了并集，但是没有出现在交集，
#目标2:绘制我们方法找到的独特的空间表达sv基因，也没有出现在despace上
uni_inter_sv_set=Reduce(union,list(nnsvg_inter=data1_LIBD_svgene_list_1215_filter[[1]][[6]],sparkX_inter=data1_LIBD_svgene_list_1215_filter[[2]][[6]],
                                   spark_inter=data1_LIBD_svgene_list_1215_filter[[3]][[6]],spv_inter=data1_LIBD_svgene_list_1215_filter[[4]][[6]],
                                   heartsvg_inter=data1_LIBD_svgene_list_1215_filter[[5]][[6]]))#不包括despace版本，共795个基因

uni_inter_in_our_sv=Reduce(intersect,list(data1_LIBD_svgene_list_1215_filter[[6]],uni_inter_sv_set))#交集的并集被我们识别的560个
extra_otheruni_set=setdiff(data1_LIBD_svgene_list_1215_filter[[6]],uni_inter_sv_set)#1222个


{
  ahc=intersect(data1_LIBD_svgene_list_1215_filter[[6]],uni_inter_sv_set)
  print(paste0("nnsvg_contri: ",length(intersect(data1_LIBD_svgene_list_1215_filter[[1]][[6]],ahc))/2795))
  print(paste0("spark_contri: ",length(intersect(data1_LIBD_svgene_list_1215_filter[[2]][[6]],ahc))/2795))
  print(paste0("spakx_contri: ",length(intersect(data1_LIBD_svgene_list_1215_filter[[3]][[6]],ahc))/2795))
  print(paste0("spv_contri: ",length(intersect(data1_LIBD_svgene_list_1215_filter[[4]][[6]],ahc))/2795))
  print(paste0("heart_contri: ",length(intersect(data1_LIBD_svgene_list_1215_filter[[5]][[6]],ahc))/2795))
  print(paste0("despace_contri: ",length(intersect(data1_LIBD_svgene_list_1215_filter[[7]],ahc))/2795))
}#计算下每个方法的小贡献



#要在别人并集的交集/并集中
other_uni_uni=Reduce(union,list(nnsvg_uni=data1_LIBD_svgene_list_1215_filter[[1]][[5]],sparkx_uni=data1_LIBD_svgene_list_1215_filter[[2]][[5]],
                                spark_uni=data1_LIBD_svgene_list_1215_filter[[3]][[5]],spv_uni=data1_LIBD_svgene_list_1215_filter[[4]][[5]],
                                heartsvg_uni=data1_LIBD_svgene_list_1215_filter[[5]][[5]]))##别人并集的并集
other_uni_inter=Reduce(intersect,list(nnsvg_uni=data1_LIBD_svgene_list_1215_filter[[1]][[5]],sparkx_uni=data1_LIBD_svgene_list_1215_filter[[2]][[5]],
                                      spark_uni=data1_LIBD_svgene_list_1215_filter[[3]][[5]],spv_uni=data1_LIBD_svgene_list_1215_filter[[4]][[5]],
                                      heartsvg_uni=data1_LIBD_svgene_list_1215_filter[[5]][[5]]))#别人并集的交集

#至少出现在了两个方法的并集中
gene_list <- list(nnsvg_uni=data1_LIBD_svgene_list_1215_filter[[1]][[5]],sparkx_uni=data1_LIBD_svgene_list_1215_filter[[2]][[5]],
                  spark_uni=data1_LIBD_svgene_list_1215_filter[[3]][[5]],spv_uni=data1_LIBD_svgene_list_1215_filter[[4]][[5]],
                  heartsvg_uni=data1_LIBD_svgene_list_1215_filter[[5]][[5]])
gene_counts <- table(unlist(gene_list))
genes_in_two_or_more <- names(gene_counts[gene_counts >= 3])


unique_intersv_our=intersect(extra_otheruni_set,other_uni_inter)##16
unique_intersv_our=intersect(extra_otheruni_set,other_uni_uni)#1393个存在于他们并集的并集的，即这些基因都被其他方法识别到过
unique_intersv_our=intersect(extra_otheruni_set,genes_in_two_or_more)


##并集有问题
##并集有问题
##并集有问题
##并集有问题
##并集有问题
other_uni_uni=Reduce(union,list(nnsvg_uni=data1_LIBD_svgene_list_1215_filter[[1]][[5]],spark_uni=data1_LIBD_svgene_list_1215_filter[[3]][[5]],
                                sparkx_uni=data1_LIBD_svgene_list_1215_filter[[2]][[5]],spv_uni=data1_LIBD_svgene_list_1215_filter[[4]][[5]],
                                heartsvg_uni=data1_LIBD_svgene_list_1215_filter[[5]][[5]]))##别人并集的并集
wecant=intersect(other_uni_uni,data1_LIBD_svgene_list_1215_filter[[6]])#1640

unique_intersv_our=intersect(extra_otheruni_set,other_uni_uni)#1080个
unique_sv_our=setdiff(extra_otheruni_set,other_uni_uni)#142个


gene_list <- list(nnsvg_uni=data1_LIBD_svgene_list_1215_filter[[1]][[5]],spark_uni=data1_LIBD_svgene_list_1215_filter[[3]][[5]],
                  sparkx_uni=data1_LIBD_svgene_list_1215_filter[[2]][[5]],spv_uni=data1_LIBD_svgene_list_1215_filter[[4]][[5]],
                  heartsvg_uni=data1_LIBD_svgene_list_1215_filter[[5]][[5]])
gene_counts <- table(unlist(gene_list))
genes_in_two_or_more <- names(gene_counts[gene_counts >= 3])#1140
genes_in_two_or_more <- names(gene_counts[gene_counts >= 2])#1898
genes_in_two_or_more <- names(gene_counts[gene_counts >= 1])#3401

denoise_svgene=intersect(data1_LIBD_svgene_list_1215_filter[[6]],genes_in_two_or_more)#868，1233,1640

noise_svgene=setdiff(genes_in_two_or_more,denoise_svgene)

position1=as.matrix(read.csv("/Users/zhoumeng/Downloads/matrix1_postion_126.csv",row.names = 1))
count1=as.matrix(read.csv("/Users/zhoumeng/Downloads/matrix1_count_126.csv",row.names = 1))
position2=as.matrix(read.csv("/Users/zhoumeng/Downloads/matrix2_postion_126.csv",row.names = 1))
count2=as.matrix(read.csv("/Users/zhoumeng/Downloads/matrix2_count_126.csv",row.names = 1))
position3=as.matrix(read.csv("/Users/zhoumeng/Downloads/matrix3_postion_126.csv",row.names = 1))
count3=as.matrix(read.csv("/Users/zhoumeng/Downloads/matrix3_count_126.csv",row.names = 1))
position4=as.matrix(read.csv("/Users/zhoumeng/Downloads/matrix4_postion_126.csv",row.names = 1))
count4=as.matrix(read.csv("/Users/zhoumeng/Downloads/matrix4_count_126.csv",row.names = 1))

#3.筛选
count1=count1[noise_svgene,]
count2=count2[noise_svgene,]
count3=count3[noise_svgene,]
count4=count4[noise_svgene,]

#绘制
meta1=meta_process(position1,count1)
meta2=meta_process(position2,count2)
meta3=meta_process(position3,count3)
meta4=meta_process(position4,count4)

for(genenum in seq_along(rownames(count1))){
  p1=pattern_plot2(meta1,genenum,xpand=-2,ypand = 2,main = FALSE)##用于绘制
  p2=pattern_plot2(meta2,genenum,xpand=-2,ypand = 2,main = FALSE)
  p3=pattern_plot2(meta3,genenum,xpand=-2,ypand = 2,main = FALSE)
  p4=pattern_plot2(meta4,genenum,xpand=-2,ypand = 2,main = FALSE)
  combined_plot <- arrangeGrob(p1, p2, p3, p4, ncol = 2, nrow = 2)
  ggsave(paste0("/Users/zhoumeng/Desktop/图片结果/1214以后/uniproblem/",rownames(count1)[genenum],".png"), plot = combined_plot, width = 8, height = 8, dpi = 300)
}
#EFEMP1
#ACO1
#ACTA2

##四张空间效应都很弱，但被我们找到了
##四张空间效应都很弱，但被我们找到了##四张空间效应都很弱，但被我们找到了##四张空间效应都很弱，但被我们找到了
##四张空间效应都很弱，但被我们找到了##四张空间效应都很弱，但被我们找到了##四张空间效应都很弱，但被我们找到了

weak_svgene_our=setdiff(data1_LIBD_svgene_list_1215_filter[[6]],other_uni_uni)

position1=as.matrix(read.csv("/Users/zhoumeng/Downloads/matrix1_postion_126.csv",row.names = 1))
count1=as.matrix(read.csv("/Users/zhoumeng/Downloads/matrix1_count_126.csv",row.names = 1))
position2=as.matrix(read.csv("/Users/zhoumeng/Downloads/matrix2_postion_126.csv",row.names = 1))
count2=as.matrix(read.csv("/Users/zhoumeng/Downloads/matrix2_count_126.csv",row.names = 1))
position3=as.matrix(read.csv("/Users/zhoumeng/Downloads/matrix3_postion_126.csv",row.names = 1))
count3=as.matrix(read.csv("/Users/zhoumeng/Downloads/matrix3_count_126.csv",row.names = 1))
position4=as.matrix(read.csv("/Users/zhoumeng/Downloads/matrix4_postion_126.csv",row.names = 1))
count4=as.matrix(read.csv("/Users/zhoumeng/Downloads/matrix4_count_126.csv",row.names = 1))

#3.筛选
count1=count1[weak_svgene_our,]
count2=count2[weak_svgene_our,]
count3=count3[weak_svgene_our,]
count4=count4[weak_svgene_our,]

#绘制
meta1=meta_process(position1,count1)
meta2=meta_process(position2,count2)
meta3=meta_process(position3,count3)
meta4=meta_process(position4,count4)

for(genenum in seq_along(rownames(count1))){
  p1=pattern_plot2(meta1,genenum,xpand=-2,ypand = 2,main = FALSE)##用于绘制
  p2=pattern_plot2(meta2,genenum,xpand=-2,ypand = 2,main = FALSE)
  p3=pattern_plot2(meta3,genenum,xpand=-2,ypand = 2,main = FALSE)
  p4=pattern_plot2(meta4,genenum,xpand=-2,ypand = 2,main = FALSE)
  combined_plot <- arrangeGrob(p1, p2, p3, p4, ncol = 2, nrow = 2)
  ggsave(paste0("/Users/zhoumeng/Desktop/图片结果/1214以后/weakproblem/",rownames(count1)[genenum],".png"), plot = combined_plot, width = 8, height = 8, dpi = 300)
}


##对比多样本方法
#与DEspace
despace_our_svgene=setdiff(data1_LIBD_svgene_list_1215_filter[[6]],data1_LIBD_svgene_list_1215_filter[[7]])
length(intersect(data1_LIBD_svgene_list_1215_filter[[6]],data1_LIBD_svgene_list_1215_filter[[7]]))

despace_remain=intersect(setdiff(data1_LIBD_svgene_list_1215_filter[[7]],data1_LIBD_svgene_list_1215_filter[[6]]),genes_in_two_or_more)
#446


#对比paste方法,以及despace
despace_our_svgene=setdiff(data1_LIBD_svgene_list_1215_filter[[6]],data1_LIBD_svgene_list_1215_filter[[7]])
despace_our_svgene=setdiff(data1_LIBD_svgene_list_1215_filter[[6]],data1_LIBD_svgene_list_1215_filter[[2]][[7]])
despace_our_svgene=setdiff(data1_LIBD_svgene_list_1215_filter[[6]],data1_LIBD_svgene_list_1215_filter[[5]][[7]])
despace_our_svgene=setdiff(data1_LIBD_svgene_list_1215_filter[[6]],data1_LIBD_svgene_list_1215_filter[[1]][[7]])


nnsvg_paste=intersect(data1_LIBD_svgene_list_1215_filter[[1]][[7]],data1_LIBD_svgene_list_1215_filter[[1]][[5]])#69
nnsvg_paste=intersect(data1_LIBD_svgene_list_1215_filter[[1]][[7]],data1_LIBD_svgene_list_1215_filter[[1]][[6]])#38
nnsvg_paste=intersect(data1_LIBD_svgene_list_1215_filter[[1]][[7]],uni_inter_sv_set)#86


spark_paste=intersect(data1_LIBD_svgene_list_1215_filter[[3]][[7]],uni_inter_sv_set)#194
spark_paste=intersect(data1_LIBD_svgene_list_1215_filter[[3]][[7]],data1_LIBD_svgene_list_1215_filter[[6]])#138
spark_paste=intersect(data1_LIBD_svgene_list_1215_filter[[3]][[7]],all_other_inter_sv)#194
remain_spark=setdiff(data1_LIBD_svgene_list_1215_filter[[3]][[7]],spark_paste)

hhh=intersect(remain_spark,noise_svgene)

count1=count1[hhh,]
count2=count2[hhh,]
count3=count3[hhh,]
count4=count4[hhh,]

meta1=meta_process(position1,count1)
meta2=meta_process(position2,count2)
meta3=meta_process(position3,count3)
meta4=meta_process(position4,count4)

for(genenum in seq_along(rownames(count1))){
  p1=pattern_plot2(meta1,genenum,xpand=-2,ypand = 2,main = FALSE)##用于绘制
  p2=pattern_plot2(meta2,genenum,xpand=-2,ypand = 2,main = FALSE)
  p3=pattern_plot2(meta3,genenum,xpand=-2,ypand = 2,main = FALSE)
  p4=pattern_plot2(meta4,genenum,xpand=-2,ypand = 2,main = FALSE)
  combined_plot <- arrangeGrob(p1, p2, p3, p4, ncol = 2, nrow = 2)
  ggsave(paste0("/Users/zhoumeng/Desktop/图片结果/1214以后/spark_paste/spark_paste/",rownames(count1)[genenum],".png"), plot = combined_plot, width = 8, height = 8, dpi = 300)
}


###绘制总图
###绘制总图
###绘制总图
easy_set=c("SCGB1D2","GFAP","KRT19")
interproblemset=c("GABBR2","PPP2CA","CAP2")
uniproblemset=c("ACER3","B3GAT2","NR2F2")
weakproblemset=c("ALKBH5","C16orf72","HNRNPA0")
total_gene_choose=c("SCGB1D2","GFAP","KRT19","GABBR2","CAP2","PPP2CA","ACER3","B3GAT2","NR2F2","ALKBH5","C16orf72","HNRNPA0")

count1=count1[total_gene_choose,]
count2=count2[total_gene_choose,]
count3=count3[total_gene_choose,]
count4=count4[total_gene_choose,]

#绘制
meta1=meta_process(position1,count1)
meta2=meta_process(position2,count2)
meta3=meta_process(position3,count3)
meta4=meta_process(position4,count4)

pattern_plot2 <- function(pltdat, igene, xy = TRUE, main = FALSE, titlesize = 2, 
                          pointsize = 3, xpand = 0, ypand = 1, title = NULL, legend = TRUE) {
  if (!xy) {
    xy <- matrix(as.numeric(do.call(rbind, strsplit(as.character(pltdat[, 1]), split = "x")), ncol = 2))
                 rownames(xy) <- as.character(pltdat[, 1])
                 colnames(xy) <- c("x", "y")
                 pd <- cbind.data.frame(xy, pltdat[, 2:ncol(pltdat)])
  } else {
    pd <- pltdat
  }
  
  # 提取目标列数据并计算相对百分比（0%~100%）
  target_col <- pd[, igene + 2]
  min_val <- min(target_col, na.rm = TRUE)
  max_val <- max(target_col, na.rm = TRUE)
  pd$rel_value <- (target_col - min_val) / (max_val - min_val) * 100  # 转换为百分比
  
  # 定义颜色渐变（固定范围为0~100）
  pal <- colorRampPalette(c("mediumseagreen", "lightyellow2", "deeppink"))
  gpt <- ggplot(pd, aes(x = x, y = y, color = rel_value)) + 
    geom_point(size = pointsize) + 
    scale_color_gradientn(
      colours = pal(5),
      limits = c(0, 100),  # 固定颜色范围为0%~100%
      labels = scales::percent_format(scale = 1)  # 标签显示为百分比（如50%）
    ) +
    scale_x_discrete(expand = c(xpand, ypand)) + 
    scale_y_discrete(expand = c(xpand, ypand)) + 
    coord_equal() + 
    theme_bw() + 
    labs(color = "Relative Expression Level ")  # 修改图例标题
    
  # 调整图例位置和样式
  if (legend) {
    gpt <- gpt + 
      theme(legend.position = "bottom",
            legend.justification = "center",
            legend.text = element_text(
              size = 8,
              margin = margin(t = 0, b = 0, unit = "cm"),  # 文字左移3mm
              hjust = 0,                                             # 文字左对齐
              vjust = 0.5                                            # 保持垂直居中
            )) +  # 确保图例居中
      guides(color = guide_colorbar(
        direction = "horizontal",
        barwidth = unit(5, "cm"),
        barheight = unit(0.3, "cm"),
        title.position = "left",
        label.position = "bottom"
      ))
  } else {
    gpt <- gpt + theme(legend.position = "none")
  }
  
  if (main) {
    out <- gpt + labs(title = title %||% colnames(pd)[igene + 2], x = NULL, y = NULL) + 
      theme(plot.title = element_text(hjust = 0.5, size = rel(titlesize)))
  } else {
    out <- gpt + labs(title = NULL, x = NULL, y = NULL)
  }
  
  return(out)
}

combined_list=list()

for(genenum in c(1:12)){
p1 <- pattern_plot2(meta1, genenum, xpand = 0.01, ypand = 0.01, main = FALSE, pointsize = 1,legend = FALSE)
p2 <- pattern_plot2(meta2, genenum, xpand = 0.01, ypand = 0.01, main = FALSE, pointsize = 1,legend = FALSE)
p3 <- pattern_plot2(meta3, genenum, xpand = 0.01, ypand = 0.01, main = FALSE, pointsize = 1,legend = FALSE)
p4 <- pattern_plot2(meta4, genenum, xpand = 0.01, ypand = 0.01, main = FALSE, pointsize = 1,legend = FALSE)

combined_list[[genenum]] <- wrap_elements((p1 | p2) / (p3 | p4)) + 
  labs(title = total_gene_choose[genenum]) + 
  theme(
    plot.title = element_text(size = 16, face = "bold", hjust = 0.5, vjust = 0),
  )
}

print(combined_list[[genenum]])


combined_plot <- 
  (wrap_elements(combined_list[[1]]|combined_list[[2]]|combined_list[[3]]) + 
      labs(tag = "A") + 
      theme(plot.tag = element_text(size = 20, face = "bold", margin = margin(0)),
            plot.tag.position = c(0.01, 1))) /
  (wrap_elements(combined_list[[4]]|combined_list[[5]]|combined_list[[6]]) + 
     labs(tag = "B") + 
     theme(plot.tag = element_text(size = 20, face = "bold", margin = margin(0)),
           plot.tag.position = c(0.01, 1))) /
  (wrap_elements(combined_list[[7]]|combined_list[[8]]|combined_list[[9]]) + 
     labs(tag = "C") + 
     theme(plot.tag = element_text(size = 20, face = "bold", margin = margin(0)),
           plot.tag.position = c(0.01, 1))) /
  (wrap_elements(combined_list[[10]]|combined_list[[11]]|combined_list[[12]]) + 
     labs(tag = "D") + 
     theme(plot.tag = element_text(size = 20, face = "bold", margin = margin(0)),
           plot.tag.position = c(0.01, 1))) +
  plot_layout(heights = c(1.2, 1.2,1.2,1.2,0.2)) 

print(combined_plot)

ggsave(paste0("/Users/zhoumeng/Desktop/图片结果/1214以后/data_express/combined1.png"), plot = combined_plot, width = 15, height = 24, dpi = 300)














# 合并图形（图例统一显示0%~100%）
combined_plot <- ((p1 | p2) / (p3 | p4)) + 
  plot_layout(guides = "collect", heights = c(1, 1)) +  
  plot_annotation(
    title = total_gene_choose[genenum],
    theme = theme(
      plot.title = element_text(hjust = 0.5, size = 12,face = "bold"),
      plot.tag = element_text(size = 14,face = "bold")) 
  ) &
  theme(legend.position = "bottom",
        plot.margin = margin(1, 1, 1, 1, "mm"))

ggsave(paste0("/Users/zhoumeng/Desktop/图片结果/1214以后/data_express/",total_gene_choose[genenum],".png"), plot = combined_plot, width = 4, height = 4, dpi = 300)



img_paths <-"/Users/zhoumeng/Desktop/图片结果/1214以后/data_express/"
images <- lapply(1:12, function(genenum) {
  img <- image_read(paste0(img_paths,total_gene_choose[genenum],".png"))
  # 添加标签（位置可调）
})

image_plots <- lapply(images, function(img) ggdraw() + draw_image(img))

combined_plot <- (
  wrap_elements(full = (image_plots[[1]] | image_plots[[2]] | image_plots[[3]])) /
    wrap_elements(full = (image_plots[[4]] | image_plots[[5]] | image_plots[[6]]))/
    wrap_elements(full = (image_plots[[7]] | image_plots[[8]] | image_plots[[9]]))/
  wrap_elements(full = (image_plots[[10]] | image_plots[[11]] | image_plots[[12]]))
) +
  plot_layout(guides = "collect",heights = c(1.2,1.2,1.2,1.2,0.3))+
  plot_annotation(
    title = NULL,
    subtitle = NULL,
    caption = NULL,
    tag_levels = list(c("A", "B","C","D")),  # 手动指定大块标签
    tag_prefix = "",
    tag_sep = "",
    theme = theme(plot.tag = element_text(size = 12, face = "bold"))  # 调整标签样式
  )

print(combined_plot)
# 查看效果
print(final_plot)
ggsave(paste0("/Users/zhoumeng/Desktop/图片结果/1214以后/data_express/combined.png"), plot = combined_plot, width = 12, height = 16, dpi = 300)





get_legend <- function(p) {
  tmp <- ggplotGrob(p)
  legend <- tmp$grobs[[which(sapply(tmp$grobs, function(x) x$name) == "guide-box")]]
  return(legend)
}
# 定义pattern_plot2函数
pattern_plot2 <- function(pltdat, igene, xy = TRUE, main = FALSE, titlesize = 2, 
                          pointsize = 3, xpand = 0, ypand = 1, title = NULL, legend = TRUE) {
  if (!xy) {
    xy <- matrix(as.numeric(do.call(rbind, strsplit(as.character(pltdat[, 1]), split = "x")), ncol = 2))
    rownames(xy) <- as.character(pltdat[, 1])
    colnames(xy) <- c("x", "y")
    pd <- cbind.data.frame(xy, pltdat[, 2:ncol(pltdat)])
  } else {
    pd <- pltdat
  }
  
  # 提取目标列数据并计算相对百分比（0%~100%）
  target_col <- pd[, igene + 2]
  min_val <- min(target_col, na.rm = TRUE)
  max_val <- max(target_col, na.rm = TRUE)
  pd$rel_value <- (target_col - min_val) / (max_val - min_val) * 100  # 转换为百分比
  
  # 定义颜色渐变（固定范围为0~100）
  pal <- colorRampPalette(c("mediumseagreen", "lightyellow2", "deeppink"))
  gpt <- ggplot(pd, aes(x = x, y = y, color = rel_value)) + 
    geom_point(size = pointsize) + 
    scale_color_gradientn(
      colours = pal(5),
      limits = c(0, 100),  # 固定颜色范围为0%~100%
      labels = scales::percent_format(scale = 1)  # 标签显示为百分比（如50%）
    ) +
    scale_x_discrete(expand = c(xpand, ypand)) + 
    scale_y_discrete(expand = c(xpand, ypand)) + 
    coord_equal() + 
    theme_bw() + 
    labs(color = "Relative Expression Level ")+
    theme(
      legend.position = "bottom",
      legend.key.width = unit(1, "cm"),
      legend.key.height = unit(0.8, "cm"),
      legend.spacing.x = unit(0.5, "cm"),# 增加色条高度
      legend.title = element_text(size = 10, margin = margin(b = 5)),  # 标题下边距
      legend.text = element_text(
        size = 8,
        margin = margin(l=-0.5,t = 1, b = 1),  # 上下对称边距
        vjust = 0.5  # 关键参数：垂直居中
      ),
      legend.box.spacing = unit(0.2, "cm")  # 控制图例容器间距
    )
  
  # 是否显示图例
  if (!legend) {
    gpt <- gpt + theme(legend.position = "none")
  }
  
  if (main) {
    out <- gpt + labs(title = title %||% colnames(pd)[igene + 2], x = NULL, y = NULL) + 
      theme(plot.title = element_text(hjust = 0.5, size = rel(titlesize)))
  } else {
    out <- gpt + labs(title = NULL, x = NULL, y = NULL)
  }
  
  return(out)
}

legend_plot <- pattern_plot2(meta1, 1, legend = TRUE)
legend <- get_legend(legend_plot)  # 提取图例
print(legend)

ggsave(paste0("/Users/zhoumeng/Desktop/图片结果/1214以后/data_express/legend.png"), plot = legend, width = 6, height = 4, dpi = 300)









###绘制vn图###绘制vn图###绘制vn图
###绘制vn图###绘制vn图###绘制vn图
###绘制vn图###绘制vn图###绘制vn图###绘制vn图###绘制vn图###绘制vn图
###绘制vn图###绘制vn图###绘制vn图
max_length <- max(sapply(data1_LIBD_svgene_list_1215_filter$SPARKX, length))  # 找到最长向量的长度
gene_matrix <- do.call(cbind, lapply(data1_LIBD_svgene_list_1215_filter$SPARKX, function(x) c(x, rep(NA, max_length - length(x)))))
colnames(gene_matrix) <- c("section1","section2","section3","section4","section5","section6","section7")
# 导出为 txt 文件
write.table(gene_matrix[,1:4], file = "/Users/zhoumeng/Desktop/图片结果/1214以后/gene_list_for_sparkx_1217.txt", sep = "\t", row.names = FALSE, quote = FALSE, na = "")


max_length <- max(sapply(data1_LIBD_svgene_list_1215_filter$spVC, length))  # 找到最长向量的长度
gene_matrix <- do.call(cbind, lapply(data1_LIBD_svgene_list_1215_filter$spVC, function(x) c(x, rep(NA, max_length - length(x)))))
colnames(gene_matrix) <- c("section1","section2","section3","section4","section5","section6","section7")
# 导出为 txt 文件
write.table(gene_matrix[,1:4], file = "/Users/zhoumeng/Desktop/图片结果/1214以后/gene_list_for_spVC_1217.txt", sep = "\t", row.names = FALSE, quote = FALSE, na = "")


max_length <- max(sapply(data1_LIBD_svgene_list_1215_filter$nnSVG, length))  # 找到最长向量的长度
gene_matrix <- do.call(cbind, lapply(data1_LIBD_svgene_list_1215_filter$nnSVG, function(x) c(x, rep(NA, max_length - length(x)))))
colnames(gene_matrix) <- c("section1","section2","section3","section4","section5","section6","section7")
# 导出为 txt 文件
write.table(gene_matrix[,1:4], file = "/Users/zhoumeng/Desktop/图片结果/1214以后/gene_list_for_nnSVG_1217.txt", sep = "\t", row.names = FALSE, quote = FALSE, na = "")

max_length <- max(sapply(data1_LIBD_svgene_list_1215_filter$SPARK, length))  # 找到最长向量的长度
gene_matrix <- do.call(cbind, lapply(data1_LIBD_svgene_list_1215_filter$SPARK, function(x) c(x, rep(NA, max_length - length(x)))))
colnames(gene_matrix) <- c("section1","section2","section3","section4","section5","section6","section7")
# 导出为 txt 文件
write.table(gene_matrix[,1:4], file = "/Users/zhoumeng/Desktop/图片结果/1214以后/gene_list_for_SPARK_1217.txt", sep = "\t", row.names = FALSE, quote = FALSE, na = "")

max_length <- max(sapply(data1_LIBD_svgene_list_1215_filter$HEARTSVG, length))  # 找到最长向量的长度
gene_matrix <- do.call(cbind, lapply(data1_LIBD_svgene_list_1215_filter$HEARTSVG, function(x) c(x, rep(NA, max_length - length(x)))))
colnames(gene_matrix) <- c("section1","section2","section3","section4","section5","section6","section7")
# 导出为 txt 文件
write.table(gene_matrix[,1:4], file = "/Users/zhoumeng/Desktop/图片结果/1214以后/gene_list_for_HEARTSVG_1217.txt", sep = "\t", row.names = FALSE, quote = FALSE, na = "")





###画细胞类别分类rigion图：
###画细胞类别分类rigion图：
###画细胞类别分类rigion图：

spot.region <- c(rep(c(rep(1,20),rep(2,12)),16),rep(c(rep(3,12),rep(1,8),rep(4,12)),4),rep(rep(3:4,each=16),12))

pal2 <- c(rgb(200, 50, 54, maxColorValue = 200), 
          rgb(224, 190, 56, maxColorValue = 230),
          rgb(50, 183, 195, maxColorValue = 255), 
          rgb(131, 31, 138, maxColorValue = 180))
pal2 <- setNames(pal2, c("1", "2", "3", "4"))

spot.coor <- as.data.frame(expand.grid(y=0:31,x=0:31)[,2:1])

gg <- ggplot(spot.coor, aes(x = x, y = y))
pl <- gg + geom_point(size = 3, 
                      aes(color = as.factor(spot.region))) +
  scale_color_manual(values=pal2) +
  theme(legend.text=element_text(size=15),        
        legend.position = "top", 
        panel.border = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), 
        axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.title.y = element_blank(),
        axis.text.y = element_blank(), 
        axis.ticks = element_blank(),
        panel.background = element_blank()) +
  guides(color = guide_legend(title = "Region", 
                              title.theme = element_text(size = 10),
                              override.aes = list(size = 7)))

ggsave(paste0("figures/region.png"), pl, width = 5, height = 5)


print(pl)

zero_ratio <- colSums(gene_matrix == 0) / nrow(gene_matrix)
barplot(zero_ratio, main = "Proportion of Zero Counts per Sample", 
        xlab = "Samples", ylab = "Zero Count Ratio", col = "skyblue", las = 2)



##绘制真实的聚类结果
##绘制真实的聚类结果
##绘制真实的聚类结果
##绘制真实的聚类结果
load("/Users/zhoumeng/Desktop/数据/layer_barcode_1218.RData")

custom_colors <- c("#54beaa", "#fccccb", "#bdb5e1", "#b0d992", "#f9d580", "#99b9e9","#e3716e","#eca680","#7ac7e2","#f3deb7")
for(datanum in c(1:4)){
position1=as.data.frame(read.csv(paste0("/Users/zhoumeng/Downloads/matrix",datanum,"_postion_126.csv"),row.names = 1))
position1$cluster=layer_barcode[[datanum]]
plot=ggplot(position1, aes(x = x, y = y, color = as.factor(cluster))) +
  geom_point(size = 3) +
  scale_color_manual(values = custom_colors, name = "Cluster") 
ggsave(paste0("/Users/zhoumeng/Desktop/细胞整合结果/1215/true_cluster_sec",datanum,".png"), plot = plot, width = 6, height = 6, dpi = 300)
}






###画introduction三种固定基因表达模式的图
###画introduction三种固定基因表达模式的图
###画introduction三种固定基因表达模式的图
###画introduction三种固定基因表达模式的图
###画introduction三种固定基因表达模式的图
#使用真实的坐标，但是以给出方法的
#使用真实的空间域划分
sim_create2 <- function(gene_size =100,svgene_size=0.1,sv_mark=c(1,1),no_sv_mark = c(0,0),inf_size=0,phi=100,etamean=2,
                       xspace="linear",yspace="linear",seed=1,cell_dist=rep(1,6),domainnum=1){
  ###
  ###
  set.seed(seed)
  position=as.matrix(read.csv("/Users/zhoumeng/Downloads/matrix3_postion_126.csv",row.names = 1))
  x_coords <- position[,1]
  y_coords <- position[,2]
  z_type=layer_barcode[[3]]

  location <- as.matrix(data.frame(x = x_coords, y = y_coords,z=z_type))
  
  x <- location[,1]
  y <- location[,2]
  
  x <- x-mean(x)
  x <- x/sd(x)
  y <- y-mean(y)
  y <- y/sd(y)
  
  x <- switch(xspace, "focal" = exp(-x^2/2),
              "period" = cos(2*pi*x),
              "linear" = x,
              stop("Invalid xspace!"))
  
  y <- switch(yspace,
              "focal" = exp(-y^2/2),
              "period" = cos(2*pi*y),
              "linear" = y,
              stop("Invalid yspace!"))
  
  kern_coord<-cbind(x,y)
  
  npoints <- nrow(location)
  rownames(location) = paste('spot', 1:npoints, sep = '')
  expres_marx = as.data.frame(matrix(NA, nrow = npoints, ncol = gene_size))
  rownames(expres_marx) = paste('spot', 1:npoints, sep = '')
  colnames(expres_marx) = paste('gene', 1:gene_size, sep = '')
  
  sv_points=svgene_size*gene_size
  sv_gene <- c(1:sv_points)##设定前多少个为SV基因
  no_sv_gene <- setdiff(1:gene_size, sv_gene)##设定其余的基因为非SV基因
  
  eta <- rnorm(gene_size,mean = 2,sd = 0.5)
  
  cell_matrix <- matrix(NA,nrow = npoints,ncol = 6)
  
  for(i in 1:npoints){
    if(z_type[i]=="L1"){
      cell_matrix[i,] <- rdirichlet(1,alpha = c(1,1,1,1,1,1))
    }else if(z_type[i]=="L2"){
      cell_matrix[i,] <- rdirichlet(1,alpha = c(1,3,5,7,9,11))
    }else if(z_type[i]=="L3"){
      cell_matrix[i,] <- rdirichlet(1,alpha = c(14,12,10,8,6,4))
    }else if(z_type[i]=="L4"){
      cell_matrix[i,] <- rdirichlet(1,alpha = c(1,4,4,4,4,1))
    }else{
      cell_matrix[i,] <- rdirichlet(1,alpha = c(18,16,14,12,10,8))
    }
  }
  
  cell_mark <- rnorm(6,mean=0,sd=1)
  
  for(i in sv_gene){
    for(t in 1:npoints){
      
      mu = exp(eta[i]+sum(kern_coord[t,]*sv_mark)+sum(cell_matrix[t,]*cell_mark))
      
      expres_marx[t,i] <- rnbinom(1,mu=mu,size=phi)
    }
  }
  for(i in no_sv_gene){
    for(t in 1:npoints){
      
      mu= exp(eta[i]+sum(kern_coord[t,]*no_sv_mark)+sum(cell_matrix[t,]*cell_mark))
      expres_marx[t,i] <- rnbinom(1,mu=mu,size=phi)
    }
  }
  
  ##添加零膨胀率。
  expres_marx <- as.matrix(expres_marx)
  total.size <- npoints*gene_size
  zero_size<-floor(npoints*gene_size*inf_size)
  zeroin<-sample(c(1:total.size),zero_size)
  expres_marx[zeroin]<-0
  location <- as.matrix(location)
  spe <- list(expres_marx,location[,1:2],cell_matrix)
  
  return(spe)
}

position=as.matrix(read.csv("/Users/zhoumeng/Downloads/matrix3_postion_126.csv",row.names = 1))

#绘制三种情形
shinySRT1 <- SRTsim_shiny()
simSRT1 <- Shiny2SRT(shinySRT1)
count=simSRT1@simCounts
count=as.matrix(count)
meta1=meta_process(position,count)
g=pattern_plot2(meta1,1,xpand=-2,ypand = 2,main = FALSE)
ggsave(paste0("/Users/zhoumeng/Desktop/图片结果/1214以后/spatialpattern/period.png"), plot = g, width = 6, height = 6, dpi = 150)




#绘制对应基因
#绘制对应基因
#绘制对应基因
#绘制对应基因
count=as.matrix(read.csv("/Users/zhoumeng/Downloads/matrix3_count_126.csv",row.names = 1))
position=as.matrix(read.csv("/Users/zhoumeng/Downloads/matrix3_postion_126.csv",row.names = 1))
count=count[c("AGRN","KRT19"),]
meta1=meta_process(position,count)
g=pattern_plot2(meta1,2,xpand=-2,ypand = 2,main = FALSE)
ggsave(paste0("/Users/zhoumeng/Desktop/图片结果/1214以后/spatialpattern/KRT19.png"), plot = g, width = 6, height = 6, dpi = 150)


other_uni_uni=Reduce(union,list(nnsvg_uni=data1_LIBD_svgene_list_1215_filter[[1]][[5]],spark_uni=data1_LIBD_svgene_list_1215_filter[[3]][[5]],
                                sparkx_uni=data1_LIBD_svgene_list_1215_filter[[2]][[5]]))


complexgene=setdiff(Reduce(intersect,list(data1_LIBD_svgene_list_1215_filter[[6]],spv_uni=data1_LIBD_svgene_list_1215_filter[[4]][[5]],
                                          heartsvg_uni=data1_LIBD_svgene_list_1215_filter[[5]][[5]])),
                    Reduce(intersect,list(nnsvg_uni=data1_LIBD_svgene_list_1215_filter[[1]][[5]],spark_uni=data1_LIBD_svgene_list_1215_filter[[3]][[5]]
                                      ,sparkx_uni=data1_LIBD_svgene_list_1215_filter[[2]][[5]])))

complexgene=setdiff(data1_LIBD_svgene_list_1215_filter[[6]],
                    Reduce(union,list(nnsvg_uni=data1_LIBD_svgene_list_1215_filter[[1]][[5]],spark_uni=data1_LIBD_svgene_list_1215_filter[[3]][[5]]
                                      ,sparkx_uni=data1_LIBD_svgene_list_1215_filter[[2]][[5]])))

position1=as.matrix(read.csv("/Users/zhoumeng/Downloads/matrix1_postion_126.csv",row.names = 1))
count1=as.matrix(read.csv("/Users/zhoumeng/Downloads/matrix1_count_126.csv",row.names = 1))
position2=as.matrix(read.csv("/Users/zhoumeng/Downloads/matrix2_postion_126.csv",row.names = 1))
count2=as.matrix(read.csv("/Users/zhoumeng/Downloads/matrix2_count_126.csv",row.names = 1))
position3=as.matrix(read.csv("/Users/zhoumeng/Downloads/matrix3_postion_126.csv",row.names = 1))
count3=as.matrix(read.csv("/Users/zhoumeng/Downloads/matrix3_count_126.csv",row.names = 1))
position4=as.matrix(read.csv("/Users/zhoumeng/Downloads/matrix4_postion_126.csv",row.names = 1))
count4=as.matrix(read.csv("/Users/zhoumeng/Downloads/matrix4_count_126.csv",row.names = 1))

dim(position4)
#3.筛选
count1=count1[complexgene,]
count2=count2[complexgene,]
count3=count3[complexgene,]
count4=count4[complexgene,]

count1=count1[data1_LIBD_svgene_list_1215_filter[[6]],]
count2=count2[data1_LIBD_svgene_list_1215_filter[[6]],]
count3=count3[data1_LIBD_svgene_list_1215_filter[[6]],]
count4=count4[data1_LIBD_svgene_list_1215_filter[[6]],]


#绘制
meta1=meta_process(position1,count1)
meta2=meta_process(position2,count2)
meta3=meta_process(position3,count3)
meta4=meta_process(position4,count4)

for(genenum in seq_along(rownames(count1))){
  p1=pattern_plot2(meta1,genenum,xpand=-2,ypand = 2,main = FALSE)##用于绘制
  p2=pattern_plot2(meta2,genenum,xpand=-2,ypand = 2,main = FALSE)
  p3=pattern_plot2(meta3,genenum,xpand=-2,ypand = 2,main = FALSE)
  p4=pattern_plot2(meta4,genenum,xpand=-2,ypand = 2,main = FALSE)
  combined_plot <- arrangeGrob(p1, p2, p3, p4, ncol = 2, nrow = 2)
  ggsave(paste0("/Users/zhoumeng/Desktop/图片结果/1214以后/proposedidentify/",rownames(count1)[genenum],".png"), plot = combined_plot, width = 8, height = 8, dpi = 300)
}




##绘制miscell图
##绘制miscell图
##绘制miscell图
##绘制miscell图
plot_miscell<-function(name1="miscell",name2="",seeddomain=c(1:50),title = "Linear pattern"){
  len=length(seeddomain)
  ##绘制devolution箱形图
  result_F1 = matrix(NA,nrow = 8*len,ncol = 2)
  result_F1[,1]=rep(c(100, 80, 60, 40, 20, 10, 5,1), each = len)
  for(seed1 in seeddomain){
    file_path <- paste0("/Volumes/Elements/", name1, "/data_miscell_",name2,"seed", seed1, ".RData")
    load(file_path) 
    for(i in c(1:8)){
      total_result=result_total[[i]][[1]]
      total_result <- apply(total_result, 1, max)
      thrs <- BayFDR(total_result, 0.05/2/length(total_result))
      total_result[total_result > thrs] <- 1
      total_result[total_result <= thrs] <- 0
      a=mean(total_result[1:500],na.rm=TRUE)
      b=mean(total_result[501:5000],na.rm=TRUE)
      F1=2*a /(1+9 * b + a)
      result_F1[seed1+(i-1)*len,2]=F1
    }
  }
  
  custom_colors <- c("#54beaa", "#fccccb", "#bdb5e1", "#b0d992", "#f9d580", "#99b9e9","#e3716e","#eca680","#7ac7e2","#f3deb7")
  colnames(result_F1) <- c("concentration_parameter", "F1")
  result_F1 <- as.data.frame(result_F1)
  result_F1$concentration_parameter <- factor(result_F1$concentration_parameter, 
                                              levels = c(100, 80, 60, 40, 20, 10, 5, 1), 
                                              ordered = TRUE)
  gg=ggplot(result_F1, aes(x = concentration_parameter, y = F1, fill = factor(concentration_parameter))) +
    geom_boxplot(color = "black", outlier.shape = 16, outlier.colour = "red") +
    facet_grid(scales = "free_y") +
    scale_fill_manual(values = custom_colors) +  
    scale_x_discrete(labels = c(100, 80, 60, 40, 20, 10, 5, 1)) +  # 确保 x 轴标签正确
    labs(title = title, x = expression("Concentration parameter" ~ alpha[0]), y = "F1 value",fill =NULL ) +
    theme_minimal()+
    theme(
      plot.title = element_text(hjust = 0.5),  # 让标题居中
      )
  return(gg)
}

p1=plot_miscell()
print(p1)
ggsave(paste0("/Users/zhoumeng/Desktop/图片结果/1214以后/miscell/linear_miscell.png"), plot = p1, width = 5, height = 4, dpi = 300)

p2=plot_miscell(name1="miscell",name2="focal_",seeddomain=c(1:32),title = "Focal pattern")
print(p1)
ggsave(paste0("/Users/zhoumeng/Desktop/图片结果/1214以后/miscell/focal_miscell.png"), plot = p1, width = 10, height = 4, dpi = 300)

p3=plot_miscell(name1="miscell",name2="period_",seeddomain=c(1:17),title = "Period pattern")

combined_plot <- ((p1 | p2 | p3)) +
  plot_layout(guides = "collect") + 
  theme(legend.position = "right") & 
  guides(fill = guide_legend(ncol = 1))

print(combined_plot)

ggsave(paste0("/Users/zhoumeng/Desktop/图片结果/1214以后/miscell/miscell.png"), plot = combined_plot, width = 10, height = 5, dpi = 300)



result1=normal_create(inf_size = 0.1,domainnum = 5,kernel="exp1",svgene_size=0.5,seed=5,gene_size=20,tau1=1,tau2=1,num_cores =22,simu_type=1)
meta1=meta_process(result1[[2]],t(result1[[1]]))
meta2=meta_process(result2[[2]],t(result2[[1]]))
meta3=meta_process(position3,count3)
meta4=meta_process(position4,count4)

for(genenum in seq_along(rownames(count1))){
  p1=pattern_plot2(meta1,genenum,xpand=-2,ypand = 2,main = FALSE)##用于绘制
  p2=pattern_plot2(meta2,genenum,xpand=-2,ypand = 2,main = FALSE)
  p3=pattern_plot2(meta3,genenum,xpand=-2,ypand = 2,main = FALSE)
  p4=pattern_plot2(meta4,genenum,xpand=-2,ypand = 2,main = FALSE)
  combined_plot <- arrangeGrob(p1, p2, p3, p4, ncol = 2, nrow = 2)
  ggsave(paste0("/Users/zhoumeng/Desktop/图片结果/1214以后/proposedidentify/",rownames(count1)[genenum],".png"), plot = combined_plot, width = 8, height = 8, dpi = 300)
}


for(genenum in c(1:20)){
  p1=pattern_plot2(meta1,genenum,xpand=-3,ypand = 3,main = FALSE)
  print(p1)
}

result1=normal_create(inf_size = 0.5,domainnum = 5,kernel="exp",svgene_size=0.5,seed=1,gene_size=20,tau1=9,tau2=1,num_cores =22)
result2=normal_create(inf_size = 0.5,domainnum = 5,kernel="exp",svgene_size=0.5,seed=2,gene_size=20,tau1=4,tau2=1,num_cores =22)
meta1=meta_process(result1[[2]],t(result1[[1]]))
meta2=meta_process(result2[[2]],t(result2[[1]]))

genenum=3
p2=pattern_plot2(meta2,genenum,xpand=-3,ypand = 3,main = FALSE)
print(p2)

p1=pattern_plot2(meta1,genenum,xpand=-3,ypand = 3,main = FALSE)
print(p1)



##绘制对应空间模式的图
##绘制对应空间模式的图
##绘制对应空间模式的图
##绘制对应空间模式的图
pattern_create <- function(gene_size =100,svgene_size=0.1,sv_mark=c(1,1),no_sv_mark = c(0,0),inf_size=0,phi=100,etamean=2,
                          xspace="linear",yspace="linear",seed=1,cell_dist=rep(1,6),domainnum=1,eta_sd=0.5){
  ###
  ###
  set.seed(seed)
  x_coords <- rep(0:31, each = 32)
  y_coords <- rep(0:31, times = 32)
  if(domainnum==1){
    z_type <- c(rep(rep(1:2,each=16),9),rep(1,96),rep(c(rep(3,12),rep(4,20)),10),rep(3,128),rep(rep(3:4,each=16),6))#横平竖直的分成了四个区域
  }else if(domainnum==2){
    z_type <- rep(5,1024)
  }else if(domainnum==3){
    z_type <- c(rep(1,512),rep(2,512))
  }else if(domainnum==4){
    z_type <- c(rep(rep(1:2,each=16),16),rep(3,512))
  }else if(domainnum==5){
    z_type <- c(rep(c(rep(1,20),rep(2,12)),16),rep(c(rep(3,12),rep(1,8),rep(4,12)),4),rep(rep(3:4,each=16),12))
  }else if(domainnum==6){
    z_type <- c(rep(c(rep(2,12),rep(1,20)),16),rep(c(rep(2,12),rep(3,20)),6),rep(3,320))
  }else if(domainnum==7){
    #判断是否是3号的偏差影响太大
    z_type <- c(rep(c(rep(2,12),rep(1,20)),16),rep(c(rep(2,12),rep(3,20)),6),rep(5,320))
  }else if(domainnum==8){
    z_type <- c(rep(c(rep(2,12),rep(1,20)),16),rep(c(rep(2,12),rep(3,20)),6),rep(4,320))
  }else if(domainnum==9){
    z_type <- c(rep(c(rep(1,20),rep(2,12)),16),rep(c(rep(3,12),rep(1,8),rep(4,12)),8),rep(rep(3:4,each=16),8))
  }else if(domainnum==10){
    z_type <- c(rep(c(rep(1,20),rep(2,12)),16),rep(c(rep(3,4),rep(1,16),rep(4,12)),4),rep(rep(3:4,each=16),12))
  }else if(domainnum==11){
    #
    z_type <- c(rep(rep(2:3,each=16),16),rep(1,512))
  }
  
  
  location <- as.matrix(data.frame(x = x_coords, y = y_coords,z=z_type))
  
  x <- location[,1]
  y <- location[,2]
  
  x <- x-mean(x)
  x <- x/sd(x)
  y <- y-mean(y)
  y <- y/sd(y)
  
  
  x <- switch(xspace, "focal" = exp(-x^2/2),
              "period" = cos(2*pi*x),
              "linear" = x,
              "focalext"=exp(-(x-1)^2/2),
              "periodext" =cos(2*pi*x+1),
              "exp"=exp(-x/2),
              "sigmoid"=1/(1+exp(-x)),
              "polynomial"=0.5 * (x + 1) * (x - 0.8) * (x - 1.6),
              "polynomial1"=-0.5*(x^3)+0.3*x,
              "polynomial2"=0.2*(x^4)-0.5*(x^2)+1.2,
              "polynomial3"=0.4*(x^3)-0.2*x^2+0.1*x,
              "polynomial4"=0.15*(x^4)-0.1*x^2+0.7,
              "polynomial5"=-0.05*(x^4)+0.4*(x^3)-0.3*(x^2)+0.2*x,
              "polynomial6"=0.25 * x^3 + 0.1 * x^2 - 0.15 * x + 0.3,
              "polynomial7"=-0.3 * x^3 + 0.6 * x^2 - 0.2 * x,
              "polynomial8"=-0.5 * x^3 - 0.2 * x^2 - 0.1 * x,
              "focalext1"=exp(-(x-0.5)^2/2),
              "sigmoid1"=1/(1+exp(-x/2)),
              stop("Invalid xspace!"))
  
  y <- switch(yspace,
              "focal" = exp(-y^2/2),
              "period" = cos(2*pi*y),
              "linear" = y,
              "focalext"=exp(-(y-1)^2/2),
              "periodext" =cos(2*pi*y+1),
              "exp"=exp(-y/2),
              "sigmoid"=1/(1+exp(-y)),
              "polynomial"=0.5 * (y + 1) * (y - 0.8) * (y - 1.6),
              "polynomial1"=-0.5*(y^3)+0.3*y,
              "polynomial2"=0.2*(y^4)-0.5*(y^2)+1.2,
              "polynomial3"=0.4*(y^3)-0.2*y^2+0.1*y,
              "polynomial4"=0.15*(y^4)-0.1*y^2+0.7,
              "polynomial5"=-0.05*(y^4)+0.4*(y^3)-0.3*(y^2)+0.2*y,
              "polynomial6"=0.25 * y^3 + 0.1 * y^2 - 0.15 * y + 0.3,
              "polynomial7"=-0.3 * y^3 + 0.6 * y^2 - 0.2 * y,
              "polynomial8"=-0.5 * y^3 - 0.2 * y^2 - 0.1 * y,
              "focalext1"=exp(-(y-0.5)^2/2),
              "sigmoid1"=1/(1+exp(-y/2)),
              stop("Invalid yspace!"))
  
  kern_coord<-cbind(x,y)
  
  npoints <- nrow(location)
  rownames(location) = paste('spot', 1:npoints, sep = '')
  expres_marx = as.data.frame(matrix(NA, nrow = npoints, ncol = gene_size))
  rownames(expres_marx) = paste('spot', 1:npoints, sep = '')
  colnames(expres_marx) = paste('gene', 1:gene_size, sep = '')
  
  sv_points=svgene_size*gene_size
  sv_gene <- c(1:sv_points)##设定前多少个为SV基因
  no_sv_gene <- setdiff(1:gene_size, sv_gene)##设定其余的基因为非SV基因
  
  eta <- rnorm(gene_size,mean = 2,sd = eta_sd)
  
  cell_matrix <- matrix(NA,nrow = npoints,ncol = 6)
  
  for(i in 1:npoints){
    if(z_type[i]==1){
      cell_matrix[i,] <- rdirichlet(1,alpha = c(1,1,1,1,1,1))
    }else if(z_type[i]==2){
      cell_matrix[i,] <- rdirichlet(1,alpha = c(1,3,5,7,9,11))
    }else if(z_type[i]==3){
      cell_matrix[i,] <- rdirichlet(1,alpha = c(14,12,10,8,6,4))
    }else if(z_type[i]==4){
      cell_matrix[i,] <- rdirichlet(1,alpha = c(1,4,4,4,4,1))
    }else{
      cell_matrix[i,] <- rep(1/6,6)
    }
  }
  
  #cell_mark <- rnorm(6,mean=0,sd=1)
  cell_mark <- rep(0,6)
  
  for(i in sv_gene){
    for(t in 1:npoints){
      mu = exp(eta[i]+sum(kern_coord[t,]*sv_mark)+sum(cell_matrix[t,]*cell_mark))
      expres_marx[t,i] <- rnbinom(1,mu=mu,size=phi)
    }
  }
  for(i in no_sv_gene){
    for(t in 1:npoints){
      mu= exp(eta[i]+sum(kern_coord[t,]*no_sv_mark)+sum(cell_matrix[t,]*cell_mark))
      expres_marx[t,i] <- rnbinom(1,mu=mu,size=phi)
    }
  }
  
  ##添加零膨胀率。
  expres_marx <- as.matrix(expres_marx)
  total.size <- npoints*gene_size
  zero_size<-floor(npoints*gene_size*inf_size)
  zeroin<-sample(c(1:total.size),zero_size)
  expres_marx[zeroin]<-0
  location <- as.matrix(location)
  spe <- list(expres_marx,location[,1:2],cell_matrix)
  
  return(spe)
}

"polynomial4"
"polynomial6"
result1=pattern_create(gene_size = 20,sv_mark = c(1,1),svgene_size=0.5,inf_size = 0.01,phi =15,xspace = "polynomial",yspace ="polynomial",seed =2,domainnum = 2,eta_sd=0.5)
meta1=meta_process(result1[[2]],t(result1[[1]]))
for(genenum in c(1:10)){
  p1=pattern_plot2(meta1,genenum,xpand=0.01,ypand =0.01,main = FALSE)
  print(p1)
}

#10\10、9、3（0.3）、7(poly6)、10(poly)
p1=pattern_plot2(meta1,10,xpand=0.01,ypand =0.01,main = FALSE)
ggsave(paste0("/Users/zhoumeng/Desktop/sigmoid.png"), plot = p1, width = 4, height = 4, dpi = 300)

normal_create <- function(gene_size =100,svgene_size=0.1,inf_size=0,kernel="exp"
                          ,seed=1,cell_dist=rep(1,6),domainnum=1,tau1=1,tau2=1,num_cores=8){
  ###
  ###
  set.seed(seed)
  x_coords <- rep(0:31, each = 32)
  y_coords <- rep(0:31, times = 32)
  if(domainnum==1){
    z_type <- c(rep(rep(1:2,each=16),9),rep(1,96),rep(c(rep(3,12),rep(4,20)),10),rep(3,128),rep(rep(3:4,each=16),6))#横平竖直的分成了四个区域
  }else if(domainnum==2){
    z_type <- rep(1,1024)
  }else if(domainnum==3){
    z_type <- c(rep(1,512),rep(2,512))
  }else if(domainnum==4){
    z_type <- c(rep(rep(1:2,each=16),16),rep(3,512))
  }else if(domainnum==5){
    z_type <- c(rep(c(rep(1,20),rep(2,12)),16),rep(c(rep(3,12),rep(1,8),rep(4,12)),4),rep(rep(3:4,each=16),12))
  }else if(domainnum==6){
    z_type <- c(rep(c(rep(2,12),rep(1,20)),16),rep(c(rep(2,12),rep(3,20)),6),rep(3,320))
  }else if(domainnum==7){
    #判断是否是3号的偏差影响太大
    z_type <- c(rep(c(rep(2,12),rep(1,20)),16),rep(c(rep(2,12),rep(3,20)),6),rep(5,320))
  }else if(domainnum==8){
    z_type <- c(rep(c(rep(2,12),rep(1,20)),16),rep(c(rep(2,12),rep(3,20)),6),rep(4,320))
  }else if(domainnum==9){
    z_type <- c(rep(c(rep(1,20),rep(2,12)),16),rep(c(rep(3,12),rep(1,8),rep(4,12)),8),rep(rep(3:4,each=16),8))
  }else if(domainnum==10){
    z_type <- c(rep(c(rep(1,20),rep(2,12)),16),rep(c(rep(3,4),rep(1,16),rep(4,12)),4),rep(rep(3:4,each=16),12))
  }else if(domainnum==11){
    #
    z_type <- c(rep(rep(2:3,each=16),16),rep(1,512))
  }
  
  
  location <- as.matrix(data.frame(x = x_coords, y = y_coords,z=z_type))
  
  x <- location[,1]
  y <- location[,2]
  
  x <- x-mean(x)
  x <- x/sd(x)
  y <- y-mean(y)
  y <- y/sd(y)
  
  #产生一个协方差，为n✖️n
  kernel_matrix=matrix(NA,nrow = 1024,ncol = 1024)
  for (i in 1:1024) {
    for (j in 1:1024) {
      kernel_matrix[i, j] <- sqrt((x[i] - x[j])^2 + (y[i] - y[j])^2)  # 欧式距离公式
    }
  }
  base_corr=tau2*diag(rep(1,1024))
  
  #生成对应的协方差矩阵
  if(kernel=="exp"){
    kernel_matrix=exp(-kernel_matrix^2/2)*tau1+base_corr
  }else if(kernel == "cos"){
    kernel_matrix=cos(2*pi*kernel_matrix)*tau1+base_corr
  }
  
  kernel_matrix=round(kernel_matrix,3)
  kern_coord<-cbind(x,y)
  
  npoints <- nrow(location)
  rownames(location) = paste('spot', 1:npoints, sep = '')
  
  expres_marx = as.data.frame(matrix(NA, nrow = npoints, ncol = gene_size))
  rownames(expres_marx) = paste('spot', 1:npoints, sep = '')
  colnames(expres_marx) = paste('gene', 1:gene_size, sep = '')
  
  sv_points=svgene_size*gene_size
  sv_gene <- c(1:sv_points)##设定前多少个为SV基因
  no_sv_gene <- setdiff(1:gene_size, sv_gene)##设定其余的基因为非SV基因
  
  cell_matrix <- matrix(NA,nrow = npoints,ncol = 6)
  rownames(cell_matrix)=paste('spot', 1:npoints, sep = '')
  for(i in 1:npoints){
    if(z_type[i]==1){
      cell_matrix[i,] <- rdirichlet(1,alpha = c(1,1,1,1,1,1))
    }else if(z_type[i]==2){
      cell_matrix[i,] <- rdirichlet(1,alpha = c(1,3,5,7,9,11))
    }else if(z_type[i]==3){
      cell_matrix[i,] <- rdirichlet(1,alpha = c(14,12,10,8,6,4))
    }else if(z_type[i]==4){
      cell_matrix[i,] <- rdirichlet(1,alpha = c(1,4,4,4,4,1))
    }else{
      cell_matrix[i,] <- rdirichlet(1,alpha = c(18,16,14,12,10,8))
    }
  }
  
  #cell_mark <- rnorm(6,mean=0,sd=1)
  cell_mark <- rep(0,6)
  nor_mean=rep(1,1024)
  for(t in c(1:npoints)){
    nor_mean[t]=sum(cell_matrix[t,]*cell_mark)
  }
  
  kernel_matrix=make_positive_definite(kernel_matrix)
  
  # for(i in sv_gene){
  #   expres_marx[,i]=mvrnorm(n = 1, mu = nor_mean, Sigma = kernel_matrix)
  # }
  # for(i in no_sv_gene){
  #   expres_marx[,i]=mvrnorm(n = 1, mu = nor_mean, Sigma = base_corr)
  # }  
  num_cores <- min(num_cores, parallel::detectCores() - 2)
  cl <- makeCluster(num_cores)
  clusterEvalQ(cl,{
    library(MASS)
  })
  
  sv_results <- parLapply(cl, sv_gene, function(i) {
    mvrnorm(n = 1, mu = nor_mean, Sigma = kernel_matrix)
  })
  
  no_sv_results <- parLapply(cl, no_sv_gene, function(i) {
    mvrnorm(n = 1, mu = nor_mean, Sigma = base_corr)
  })
  
  stopCluster(cl)
  
  # Assign the results back to expres_marx
  for (i in seq_along(sv_gene)) {
    expres_marx[, sv_gene[i]] <- sv_results[[i]]
  }
  
  for (i in seq_along(no_sv_gene)) {
    expres_marx[, no_sv_gene[i]] <- no_sv_results[[i]]
  }
  
  expres_marx[expres_marx<0]=0
  expres_marx=round(expres_marx)
  ##添加零膨胀率。
  expres_marx <- as.matrix(expres_marx)
  total.size <- npoints*gene_size
  zero_size<-floor(npoints*gene_size*inf_size)
  zeroin<-sample(c(1:total.size),zero_size)
  expres_marx[zeroin]<-0
  location <- as.matrix(location)
  
  spe <- list(expres_marx,location[,1:2],cell_matrix)
  
  return(spe)
}


result1=normal_create(inf_size = 0.01,domainnum = 5,kernel="exp",svgene_size=0.5,seed=1,gene_size=20,tau1=9,tau2=1,num_cores =22)



x=seq(0,1,0.01)
y=0.25 * x^3 + 0.1 * x^2 - 0.15 * x + 0.3
plot(x,y)





name1=c("linear","focal","period","complex","linear_gene","focal_gene","period_gene","complex_gene")
img_paths <-"/Volumes/Elements/supple_0323/"
images <- lapply(1:8, function(genenum) {
  img <- image_read(paste0(img_paths,name1[genenum],".png"))
})


##使用图片合并
name2=c("Linear","Focal","Period","Complex","AQP4","COX6C","CAMK2N1","AGR2")
image_plots <- lapply(seq_along(images), function(i) {
  ggdraw() + 
    draw_image(images[[i]]) +  # 用索引提取图像
    draw_label(
      name2[i],               # 用索引提取对应的标题
      x = 0.5, y = 0.95, 
      hjust = 0.5, vjust = 1,
      size = 10, 
      color = "black",
      fontface = "bold"
    )
})

combined_plot <- (
  wrap_elements(full = (image_plots[[1]] | image_plots[[2]] | image_plots[[3]]|image_plots[[4]])) /
    wrap_elements(full = (image_plots[[5]] |image_plots[[6]] | image_plots[[7]]|image_plots[[8]] ))
) +
  plot_layout(guides = "collect",heights = c(0.8,0.8))+
  plot_annotation(
    title = NULL,
    subtitle = NULL,
    caption = NULL,
    tag_levels = list(c("A", "B")),  # 手动指定大块标签
    tag_prefix = "",
    tag_sep = "",
    theme = theme(plot.tag = element_text(size = 14, face = "bold"))  # 调整标签样式
  )

print(combined_plot)
# 查看效果

ggsave(paste0("/Users/zhoumeng/Desktop/图片结果/1214以后/data_express/combined_fig2.png"), plot = combined_plot, width = 8, height = 5.5, dpi = 300)



##绘制vn图##绘制vn图
##绘制vn图##绘制vn图
##绘制vn图##绘制vn图##绘制vn图
##绘制vn图##绘制vn图##绘制vn图
library(VennDiagram)

name1=c("nnsvg_venn_0404","spark_venn_0404","spvc_venn_0404","sparkx_venn_0404","heartsvg_venn_0404")
img_paths <-"/Volumes/Elements/supple_0323/"
images <- lapply(1:5, function(genenum) {
  img <- image_read(paste0(img_paths,name1[genenum],".png"))
})
name2=c("nnSVG","SPARK","spVC","SPARKX","HEARTSVG")
image_plots <- lapply(seq_along(images), function(i) {
  ggdraw() + 
    draw_image(images[[i]]) +  # 用索引提取图像
    draw_label(
      name2[i],               # 用索引提取对应的标题
      x = 0.5, y = 0.95, 
      hjust = 0.5, vjust = 1,
      size = 10, 
      color = "black",
      fontface = "bold"
    )
})
combined_plot <- 
  (image_plots[[1]] | image_plots[[2]] | image_plots[[5]]|image_plots[[4]] |image_plots[[3]] ) +
  plot_layout(guides = "collect")+
  plot_annotation(
    title = NULL,
    subtitle = NULL,
    caption = NULL
  )

print(combined_plot)
ggsave(paste0("/Users/zhoumeng/Desktop/图片结果/1214以后/data_express/combined_venn.png"), plot = combined_plot, width = 7, height = 5, dpi = 300)






name1=c("linearfocal","focalperiod","linearperiod","ZINNGP","sigmoid","polynomial1","polynomial2","polynomial3","polynomial4")
img_paths <-"/Users/zhoumeng/Desktop/图片结果/1214以后/special pattern/"
images <- lapply(1:9, function(genenum) {
  img <- image_read(paste0(img_paths,name1[genenum],".png"))
})

##使用图片合并
name2=c("Linear-Focal","Focal-Period","Linear-Period","ZINNGP","Sigmoid","Polynomial1","Polynomial2","Polynomial3","Polynomial4")
image_plots <- lapply(seq_along(images), function(i) {
  ggdraw() + 
    draw_image(images[[i]]) +  # 用索引提取图像
    draw_label(
      name2[i],               # 用索引提取对应的标题
      x = 0.5, y = 1.06, 
      hjust = 0.5, vjust = 1,
      size = 9, 
      fontface = "bold"
    )+
    theme(
      plot.margin = margin(10, 0, 10, 0)  # 上、右、下、左边距分别为20、20、40、20（增加下方边距，避免标题遮挡）
    )
})

combined_plot <- 
  (wrap_elements(full = (image_plots[[4]] | image_plots[[1]] | image_plots[[2]])) /
    wrap_elements(full = (image_plots[[3]] |image_plots[[5]] | image_plots[[6]]))/
    wrap_elements(full = (image_plots[[7]] |image_plots[[8]] | image_plots[[9]])))+
  plot_layout(guides = "collect")+
  plot_annotation(
    title = NULL,
    subtitle = NULL,
    caption = NULL,
    tag_levels = list(c("A", "B","C")),  # 手动指定大块标签
    tag_prefix = "",
    tag_sep = "",
    theme = theme(plot.tag = element_text(size = 12, face = "bold")
                  )  # 调整标签样式
  )

print(combined_plot)
ggsave(paste0("/Users/zhoumeng/Desktop/图片结果/1214以后/special pattern/combined_extend_pattern.png"), plot = combined_plot, width = 6, height =6, dpi = 300)

print(image_plots[[4]])










install.packages("ggvenn")
library(ggvenn)

venn_list=list()
namelist=c("nnSVG","SPARKX","SPARK","spVC","HEARTSVG")
for(i in c(1:5)){
  data <- list(
    "Section 1"=data1_LIBD_svgene_list_1215_filter[[i]][[1]],
    "Section 2"=data1_LIBD_svgene_list_1215_filter[[i]][[2]],
    "Section 3"=data1_LIBD_svgene_list_1215_filter[[i]][[3]],
    "Section 4"=data1_LIBD_svgene_list_1215_filter[[i]][[4]]
  )
  
  colors <- c(A = "#E41A1C", B = "#377EB8", C = "#4DAF4A", D = "#FFD700") # 红、蓝、绿、黄
  venn_plot <- ggvenn(
    data,
    fill_color = colors,
    stroke_size = 0,
    set_name_size = 0,
    show_percentage = FALSE
  ) + 
    coord_fixed(ratio = 5/4) +   # 拉伸纵向（上下拉伸）
    ggtitle(namelist[i]) +  # 添加标题
    theme(plot.title = element_text(size = 18, face = "bold", hjust = 0.5,vjust = -0.5),
          text = element_text(size = 6) ) 
  
  venn_list[[i]]=venn_plot
}


combined_plot <- (wrap_elements(full = (venn_list[[5]] | venn_list[[2]] | venn_list[[3]]|venn_list[[4]]|venn_list[[1]]))) 
print(combined_plot)

ggsave(paste0("/Users/zhoumeng/Desktop/图片结果/1214以后/data_express/combined_venn.png"), plot = combined_plot, width = 16, height =4, dpi = 300)













