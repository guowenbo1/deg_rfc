library(tidyverse)

myPalette <- c("#FF4040","#228B22","#984EA3","#FF7F00","#6959CD","#F781BF","#458B74","#A52A2A","#8470FF","#8B4513","#556B2F","#CD5B45","#483D8B","#EEC591","#8B0A50","#696969", "#8B6914","#20B2AA", "#8B636C","#36648B","#9ACD32","#8B3A62","#68838B","#CDBE70","#D3D3D3","#EEA2AD","#53868B","#EE9572","#FFA500","#8B4789","#548B54","#F4A460","#3A5FCD","#B0E0E6","#6495ED","#EECBAD","#CD69C9","#436EEE","#8B8B00","#8B7E66","#CD853F","#8B7B8B","#00C5CD","#008B45","#D2B48C","#CDCD00","#836FFF","#7D26CD","#006400","slateblue2","slategray3","violetred2","#00868B","palegreen1","#7A8B8B","palevioletred2","powderblue","orangered2","lightgoldenrod2","deepskyblue", "cyan1","midnightblue","seashell2","yellowgreen","tomato3","violetred4","magenta","chocolate2","darkorange3","olivedrab2","#8B3626","#0000FF")


index <- snakemake@wildcards[["index"]]
group <- snakemake@wildcards[["group"]]

data <- read_tsv(snakemake@input[["matrix"]])
metadata <- read_tsv(snakemake@input[["metadata"]])

gather_cols <- colnames(data)[-1]
data <- data %>% rename(GeneID = `Gene ID`) %>% 
  gather_(key_col = "SampleID", value_col = index, gather_cols = gather_cols) %>%
  left_join(metadata)


# plot

p <- ggplot(data, aes(x = log10(.data[[index]]), color = .data[[group]],fill = .data[[group]])) +
  geom_density(alpha=.4) + ylab("density") + xlab(paste0("log10(", index, ")")) +
  theme_bw() + 
  scale_color_manual(values=myPalette, name = NULL) + 
  scale_fill_manual(values=myPalette, name = NULL) + 
  theme(
    panel.grid = element_blank(),
    axis.text.x = element_text(face="bold", color="black"),
    axis.text.y = element_text(face="bold",  color="black"),
    axis.title.x = element_text(face="bold", color="black"),
    axis.title.y = element_text(face="bold",color="black"),    
    legend.title = element_blank()
  )
  
  
ggsave(snakemake@output[["fig"]][1], p, width = 8, height = 6)
ggsave(snakemake@output[["fig"]][2], p, width = 8, height = 6)