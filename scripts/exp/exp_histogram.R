library(data.table)
library(ggplot2)

args <- commandArgs(T)

if(length(args) == 2) {
  width = 10
  height = 10
} else if(length(args) == 100) {
  width = as.numeric(args[3])
  height = as.numeric(args[4])
} else {
  writeLines("
------------------------------------------------------
Usage1:Rscript *.R FPKM_table out_dir
Usage2:Rscript *.R FPKM_table out_dir width height
------------------------------------------------------
")
  q()
}

FPKM_table <- args[1]
outdir <- args[2]

# outdir <- strsplit("Expression/Histogram/",'/')[[1]]

data <- fread(FPKM_table)
samples <- names(data)[-1]
alldata = NULL
for(i in samples) {
  # print(i)
  data_sample <- data[[i]]
  # cut默认left = T, 左开右闭
  cut_data <- table(cut(data_sample, breaks = c(0, 0.01,  0.1, 1, 10, 100, 1000, Inf)))
  dd <- as.data.frame(cut_data)
  names(dd) <- c('Expression', 'Count')
  arrange_level = c('0~0.01', '0.01~0.1', '0.1~1', '1~10', '10~100', '100~1000', '>1000')
  dd$`Expression` = arrange_level
  dd$`Expression` = factor(dd$`Expression`, levels =  arrange_level)
  
  dd$lab <- as.character(sprintf("%0.2f", 100 * dd$`Count` / sum(dd$`Count`)))
  # dd$lab <- paste(dd$`Count`, "\n", paste0(dd$lab, "%"))
  dd$lab <- paste0(dd$`Count`,"(",dd$lab,"%",")")
  dd$plan <- rep(i,time=dim(dd)[1])

  alldata <- rbind(alldata, dd)
  
  p <- ggplot(data = dd, mapping=aes(x = `Expression`, y = `Count`, label = `lab`)) + 
    geom_bar(stat = "identity",
             width = 0.5,
             fill = "#add8e6") +
    geom_text(mapping=aes(x = `Expression`, y = `Count`, label = `lab`, vjust = -0.1)) + 
    theme_bw() + 
    theme(# element_blank(),
          panel.background = element_rect(fill = 'white'),
          panel.border = element_rect(colour = 'black', fill=NA),
          axis.text.x = element_text(size = 8, angle = 60, hjust = 1),
          axis.text.y = element_text(size = 8)) +
    ggtitle(i)

  ggsave(paste0(outdir,'/',i,"_histogram" ,".png"), p, width = width, height = height)
  ggsave(paste0(outdir,'/',i,"_histogram" ,".pdf"), p, width = width, height = height)
}


p2<- ggplot(alldata,aes(x=plan,y=Count,fill=Expression))+
  geom_col(position="dodge",width=0.6,key_glyph=draw_key_rect)+
  # scale_fill_manual(values=c('#A6CEE3','#B2DF8A','#FB9A99','#FDC67E','#FDBF6F','#FF7F00','#FF942B'))+
  scale_fill_brewer("FPKM.Interval",label=c('0~0.01', '0.01~0.1', '0.1~1', '1~10', '10~100', '100~1000', '>1000'),palette="Set2")+
  labs(y="Count")+labs(x="Sample ID") + 
  ggtitle("FPKM Interval Distribution")+
  theme(plot.title = element_text(hjust = 0.5))+
  theme_bw() + 
  theme(plot.title = element_text(hjust = 0.5))+
  theme(# element_blank(),
        panel.background = element_rect(fill = 'white'),
        panel.border = element_rect(colour = 'black', fill=NA),
        axis.text.x = element_text(size = 8, angle = 60, hjust = 1),
        axis.text.y = element_text(size = 8))


  ggsave(paste0(outdir,'/',"_histogram_all" ,".png"), p2, width = width+2, height = height-1)
  ggsave(paste0(outdir,'/',"_histogram_all" ,".pdf"), p2, width = width+2, height = height-1)
