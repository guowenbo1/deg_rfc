library(data.table)
library(ggplot2)

args <- commandArgs(T)

if(length(args) == 2) {
  width = 10
  height = 10
} else if(length(args) == 4) {
  width = as.numeric(args[3])
  height = as.numeric(args[4])
} else {
  writeLines("
------------------------------------------------------
Usage1:Rscript *.R Count_table out_dir
Usage2:Rscript *.R Count_table out_dir width height
------------------------------------------------------
")
  q()
}

FPKM_table <- args[1]
out_dir <- args[2]

data <- fread(FPKM_table)
samples <- names(data)[-1]
for(i in samples) {
  # print(i)
  data_sample <- data[[i]]
  # cut默认left = T, 左开右闭
  cut_data <- table(cut(data_sample, breaks = c(0, 0.01,  0.1, 1, 10, 100, 1000, Inf)))
  dd <- as.data.frame(cut_data)
  names(dd) <- c('Expression Level', 'Count')
  arrange_level = c('0~0.01', '0.01~0.1', '0.1~1', '1~10', '10~100', '100~1000', '>1000')
  dd$`Expression Level` = arrange_level
  dd$`Expression Level` = factor(dd$`Expression Level`, levels =  arrange_level)
  
  dd$lab <- as.character(sprintf("%0.2f", 100 * dd$`Count` / sum(dd$`Count`)))
  dd$lab <- paste(dd$`Count`, "\n", paste0(dd$lab, "%"))
  
  p <- ggplot(data = dd, mapping=aes(x = `Expression Level`, y = `Count`, label = `lab`)) + 
    geom_bar(stat = "identity",
             width = 0.5,
             fill = "#add8e6") +
    geom_text(mapping=aes(x = `Expression Level`, y = `Count`, label = `lab`, vjust = -0.1)) + 
    theme_bw() + 
    theme(# element_blank(),
          panel.background = element_rect(fill = 'white'),
          panel.border = element_rect(colour = 'black', fill=NA),
          axis.text.x = element_text(size = 8, angle = 60, hjust = 1),
          axis.text.y = element_text(size = 8)) +
    ggtitle(i)
  if(!dir.exists(out_dir)) {
    dir.create(out_dir, recursive = T)
  }
  ggsave(paste0(out_dir, "/", i,"_histogram" ,".png"), p, width = width, height = height)
  ggsave(paste0(out_dir, "/", i,"_histogram" ,".pdf"), p, width = width, height = height)
}
