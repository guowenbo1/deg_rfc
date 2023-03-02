#nsets需要自己改一下，表示有多少个样本
#text.scale中的参数表示（交集标题，交集标签，设置标题，设置标签，设置样本名称，条上的数字）的字体大小，样本多时修改样本名的大小

suppressPackageStartupMessages(library(tidyverse, quietly=TRUE))
suppressPackageStartupMessages(library(argparser, quietly=TRUE))
argv <- arg_parser("维恩图分析") %>%
  add_argument("--feature", help="feature table ",default = "otutab.txt") %>%
  add_argument("--pdf", help = "韦恩图路径", default = "venn.pdf") %>%  
  add_argument("--output", help = "韦恩图路径", default = "venn.png") %>%
  add_argument("--height", help = "height, 图片高度",default = 14, type = "numeric") %>%
  add_argument("--width", help = "width, 图片宽度", default = 16, type = "numeric") %>%
  parse_args()

feature <- colnames(data)[[1]]
fig1 <- argv$pdf
fig <- argv$output
title <- argv$title
height <- argv$height
width <- argv$width

outdir <- dirname(fig)

library(UpSetR)
otu <- read.delim(argv$feature, header = T, row.names = 1, sep = '\t', stringsAsFactors = F, check.names = F)
otu[otu > 0] <- 1
png(file= fig, width = width, height = height, res = 300, units = 'in')
upset(otu, nsets = as.numeric(length(names(otu))),point.size = 2, line.size = 0.3,
      mainbar.y.label = "Intersections size", sets.x.label = "Set size", 
      text.scale = c(1, 0.8, 1, 0.75, 1,0.75),order.by = "freq",
      show.numbers = "no",
      matrix.color = "#333333",main.bar.color = "#999999",
      sets.bar.color = "#666666",
      shade.color = "#CCCCCC",att.color = "gay")
dev.off()

pdf(file= fig1,onefile=F,width = width, height = height)
upset(otu, nsets = as.numeric(length(names(otu))),point.size = 2, line.size = 0.3,
      mainbar.y.label = "Intersections size", sets.x.label = "Set size", 
      text.scale = c(1, 0.8, 1, 0.75, 1,0.75),order.by = "freq",
      show.numbers = "no",
      matrix.color = "#333333",main.bar.color = "#999999",
      sets.bar.color = "#666666",
      shade.color = "#CCCCCC",att.color = "gay")
dev.off()