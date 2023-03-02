library(argparse)
library(ggplot2)
library(tidyverse)
library(ggrepel)
library(data.table)
library(RColorBrewer)

parser <- ArgumentParser(description='statistics and visualization of DEG')
parser$add_argument('-i','--input',type='character',nargs="*",
        help="A_vs_B_diff.xls *< one or more files>")
parser$add_argument('-p','--pdf',type='character',
        default='all_diff_plot.pdf',help="output pdf")
parser$add_argument('-o','--out',type='character',
         default='all_diff_plot.png',help="output png")

args <- parser$parse_args()
if(is.null(args$input)){
        parser$print_help()
        q()
}

file <- unique(args$input)
df = NULL
for (i in 1:length(file)) {
  filename = basename(file[i])
  message(paste0("parsing: ",i," ",file[i]))
  data <- fread(file[i],header=T,fill=TRUE,sep="\t") %>% select("gene_id","gene_name","logFC","level","pvalue","FDR","sampleA","sampleB")
  data$cluster <- paste0(data$sampleA,"_vs_",data$sampleB)
  data$cluster_num <-  rep (i, times = dim(data)[1])
  df <- rbind(df, data)
}

#添加显著性标签：
# df$label <- ifelse(df$FDR<0.01,"FDR<0.01","0.01<=FDR<0.05")

#获取每个cluster中表达差异最显著的10个基因；
top10sig = NULL
for (i in unique(df$cluster)) {
  message(paste0("parsing: ",i," ",file[i]))
  top10 <- filter(df,cluster== i) %>% distinct(gene_id,.keep_all = T) %>% top_n(10,abs(logFC))
  top10sig <- rbind(top10sig, top10)
}

#新增一列，将Top10的差异基因标记为2，其他的标记为1；
# df$size <- case_when(!(df$gene_id %in% top10sig$gene_id)~ 1,
#                      df$gene_id %in% top10sig$gene_id ~ 2)
# #提取非Top10的基因表格；
# dt <- filter(df,size==1)

#根据log2FC区间确定背景柱长度：
#绘制背景柱：
## y轴上半部分阴影
bar_up <- df %>% filter(level =="Increased") %>% group_by(cluster)  %>% top_n(1,abs(logFC)) %>% 
  transmute(logFC=c(logFC+0.5),cluster_num,cluster) %>% as.data.frame()
# y轴下半部分阴影
bar_down <- df  %>% filter(level=="Decreased") %>% group_by(cluster)  %>% top_n(1,abs(logFC)) %>% 
  transmute(logFC=c(logFC-0.5),cluster_num,cluster)

# 生成用于绘制方块的数据
tile <-data.frame(cluster_num=bar_down$cluster_num,y=0,cluster=bar_down$cluster)

# 定义差异基因的颜色，RColorBrewer包中的红蓝
red <- "#A8152B"
blue <- "#009966"
# 定义五颜六色方块的颜色
tile_color <- brewer.pal(length(file),"Set2")

# 定义柱形阴影宽度
width_bar <- 0.8
# 绘图
p <- ggplot() +
  # 首先绘制柱形阴影
  geom_col(bar_up,mapping=aes(cluster_num,logFC),
           fill = "#dcdcdc",alpha = 0.6,
           size=0.85,width = width_bar ) +
  geom_col(bar_down,mapping=aes(cluster_num,logFC),
           fill = "#dcdcdc",alpha = 0.6,
           size=0.85,width = width_bar)+
  # 绘制差异基因红蓝圆点
  geom_jitter(data = df,aes(cluster_num,logFC,color=level),
              size=0.85,width = 0.4,alpha= .5)+
  scale_color_manual(values = c(blue,red)) +
  # 再次绘制上下调变化倍数最大的基因
  geom_jitter(data = top10sig,aes(cluster_num,logFC,color=level),
              size=1.2,width = 0.4,alpha= .7)+
  scale_color_manual(values = c(blue,red)) +
  # 为上下调变化倍数最大的基因添加标签
  geom_text_repel(
    data=top10sig,
    aes(x=cluster_num,y=logFC,label=gene_name),
    force = 2,size=2, max.overlaps = 100,
    arrow = arrow(length = unit(0.008:0.0008, "npc"),
                  type = "open", ends = "last")
  )+
  # 绘制五颜六色的方块
  geom_tile(tile,mapping=aes(cluster_num,y,fill=cluster),
            height=1*1.8,
            color = "black",
            #fill = brewer.pal(length(files), "Set2"),
            fill = tile_color,
            alpha = 0.6,
            #show.legend = F,
            width=width_bar) +
  # x,y轴标题
  labs(x="",y="logFC")+
  geom_text(data=tile,
            aes(x=cluster_num,y=0,label=cluster),
            size =4.5,
            color ="black") +
  # 设置绘图主题
  theme_minimal()+
  theme(
    axis.title = element_text(size = 13,
                              color = "black",
                              face = "bold"),
    axis.line.y = element_line(color = "black",
                               size = 1.2),
    axis.line.x = element_blank(),
    axis.text.x = element_blank(),
    panel.grid = element_blank(),
    legend.position = "right",
    legend.direction = "vertical",
    #legend.justification = c(1,0),
    legend.text = element_text(size = 15),
  )
  # 添加差异基因筛选的阈值
 # + annotate("text", x=1, y=15, 
 #         label=paste0("|log2FC|>=",1,"\np-adj<0.05"))



width <- 2.5*length(file)


ggsave(args$pdf, height = 5, width = width)
ggsave(args$out, height = 5, width = width)