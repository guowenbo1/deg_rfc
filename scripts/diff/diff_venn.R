#!/usr/bin/env Rscript

suppressPackageStartupMessages(library(tidyverse, quietly=TRUE))
suppressPackageStartupMessages(library(argparser, quietly=TRUE))



argv <- arg_parser("维恩图分析") %>%
  add_argument("--feature", help="feature table ",default = "all_DEGs_list.xls") %>%
  add_argument("--pdf", help = "韦恩图路径", default = "venn_or_flower.pdf") %>%  
  add_argument("--output", help = "韦恩图路径", default = "venn_or_flower.png") %>%
  add_argument("--height", help = "height, 图片高度",default = 6, type = "numeric") %>%
  add_argument("--width", help = "width, 图片宽度", default = 6, type = "numeric") %>%
  add_argument("--titile", help = "图片标题", default = "DEGs") %>%
  parse_args()

library(tidyverse)
library(VennDiagram)

myPalette <- c("#66C2A5","#DA9870","#AB98C8","#DF92B6","#ADCF60","#3787BC","BE979C")

data <- read_tsv(argv$feature)
fig1 <- argv$pdf
fig <- argv$output
title <- argv$title
height <- argv$height
width <- argv$width

outdir <- dirname(fig)

lst<-as.list(data)

for(i in 1:length(data)) {
lst[[i]]<-na.omit(lst[[i]])
}

if(length(lst)<=5 & length(lst) >1){
  venn.diagram(lst, resolution = 300, imagetype = "png", alpha = rep(0.8, length(lst)),
              fill = myPalette[1:length(lst)], cat.fontface=4, fontfamily=3,
              col = myPalette[1:length(lst)],
              main = title,
              main.cex = 2, main.fontface = 2, main.fontfamily = 3,
              margin = 0.2,
              filename = fig)
 p<-venn.diagram(lst, resolution = 300, alpha = rep(0.8, length(lst)),
              fill = myPalette[1:length(lst)],
              col = myPalette[1:length(lst)],
              main = title,
              main.cex = 3,
              margin = 0.2,
              filename = NULL)
pdf(file= fig1 )
grid.draw(p)
dev.off()  
} else if (length(lst) < 2){
    data<-as.data.frame(lst)
    data$value <- rep(1, time =dim(data)[1])
    png(file= fig)
    pie(data[1,]$value, labels = dim(data)[1], col = "#66C2A5", border = NA, main = names(data)[1])
    dev.off()
    pdf(file= fig1)
    pie(data[1,]$value, labels = dim(data)[1], col = "#66C2A5", border = NA, main = names(data)[1])
    dev.off()

}else {


library(stringi)
data_name <- read.delim(argv$feature, header = T,sep = '\t', stringsAsFactors = F, check.names = F)

sample_id <- colnames(data_name)
otu_id <- unique(data_name[,1])
otu_id <- otu_id[otu_id != '']
core_otu_id <- otu_id
otu_num <- length(otu_id)

for (i in 2:ncol(data_name)) {
  otu_id <- unique(data_name[,i])
  otu_id <- otu_id[otu_id != '']
  core_otu_id <- intersect(core_otu_id, otu_id)
  otu_num <- c(otu_num, length(otu_id))
}
core_num <- length(core_otu_id)

#下文中绘制椭圆的draw.ellipse()、draw.circle()等命令，需要plotrix包实现
library(plotrix)

#定义备选颜色(20个颜色)
ellipse_col <- c('#6181BD4E','#F348004E','#64A10E4E','#9300264E','#464E044E','#049a0b4E','#4E0C664E','#D000004E','#FF6C004E','#FF00FF4E','#c7475b4E','#00F5FF4E','#BDA5004E','#A5CFED4E','#f0301c4E','#2B8BC34E','#FDA1004E','#54adf54E','#CDD7E24E','#9295C14E')


flower_plot <- function(sample, otu_num, core_otu, start, a, b, r, ellipse_col, circle_col) {
  par( bty = 'n', ann = F, xaxt = 'n', yaxt = 'n', mar = c(1,1,1,1))
  plot(c(0,10),c(0,10),type='n')
  n   <- length(sample)
  deg <- 360 / n
  res <- lapply(1:n, function(t){
    draw.ellipse(x = 5 + cos((start + deg * (t - 1)) * pi / 180),
                 y = 5 + sin((start + deg * (t - 1)) * pi / 180),
                 col = ellipse_col[t],
                 border = ellipse_col[t],
                 a = a, b = b, angle = deg * (t - 1))
    text(x = 5 + 2.5 * cos((start + deg * (t - 1)) * pi / 180),
         y = 5 + 2.5 * sin((start + deg * (t - 1)) * pi / 180),
         otu_num[t])
    
    if (deg * (t - 1) < 180 && deg * (t - 1) > 0 ) {
      text(x = 5 + 3.3 * cos((start + deg * (t - 1)) * pi / 180),
           y = 5 + 3.3 * sin((start + deg * (t - 1)) * pi / 180),
           sample[t],
           srt = deg * (t - 1) - start,
           adj = 1,
           cex = 1
      )
    } else {
      text(x = 5 + 3.3 * cos((start + deg * (t - 1)) * pi / 180),
           y = 5 + 3.3 * sin((start + deg * (t - 1)) * pi / 180),
           sample[t],
           srt = deg * (t - 1) + start,
           adj = 0,
           cex = 1
      )
    }     
  })
  draw.circle(x = 5, y = 5, r = r, col = circle_col, border = NA)
  text(x = 5, y = 5, paste('Core:', core_otu))
}

#调用上述函数作图，视情况修改参数
png(file = fig, width = 1000, height = 1000, res = 100, units = 'px')
flower_plot(sample = sample_id, otu_num = otu_num, core_otu = core_num,
            start = 90, a = 0.5, b = 2, r = 1, ellipse_col = ellipse_col, circle_col = 'white')
dev.off()

pdf(file = fig1)
flower_plot(sample = sample_id, otu_num = otu_num, core_otu = core_num,
            start = 90, a = 0.5, b = 2, r = 1, ellipse_col = ellipse_col, circle_col = 'white')
dev.off()
}



