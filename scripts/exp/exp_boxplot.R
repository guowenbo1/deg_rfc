library(tidyverse)
library(tidyr)
library(data.table)
library(ggplot2)

myPalette <- c("#FF4040","#228B22","#984EA3","#FF7F00","#6959CD","#F781BF","#458B74","#A52A2A","#8470FF","#8B4513","#556B2F","#CD5B45","#483D8B","#EEC591","#8B0A50","#696969", "#8B6914","#20B2AA", "#8B636C","#36648B","#9ACD32","#8B3A62","#68838B","#CDBE70","#D3D3D3","#EEA2AD","#53868B","#EE9572","#FFA500","#8B4789","#548B54","#F4A460","#3A5FCD","#B0E0E6","#6495ED","#EECBAD","#CD69C9","#436EEE","#8B8B00","#8B7E66","#CD853F","#8B7B8B","#00C5CD","#008B45","#D2B48C","#CDCD00","#836FFF","#7D26CD","#006400","slateblue2","slategray3","violetred2","#00868B","palegreen1","#7A8B8B","palevioletred2","powderblue","orangered2","lightgoldenrod2","deepskyblue", "cyan1","midnightblue","seashell2","yellowgreen","tomato3","violetred4","magenta","chocolate2","darkorange3","olivedrab2","#8B3626","#0000FF")

index <- snakemake@wildcards[["index"]]
group <- snakemake@wildcards[["group"]]

data <- read_tsv(snakemake@input[["matrix"]])
metadata <- read_tsv(snakemake@input[["metadata"]])

setnames(data, "Gene ID", "GeneID")
setnames(metadata, "group", "Group")

data <- data %>% melt(id.vars = "GeneID", 
                      variable.name = "SampleID",
                      value.name = "FPKM") %>% 
  merge(metadata[, 1:2], all.x = T)


# boxplot (remove fpkm = 0)
# myPalette <- c("#FF4040", "#228B22", "#FF7F00", "#6959CD")

p1 <- ggplot(data, aes(SampleID, log10(FPKM) ,color = Group)) + 
  geom_boxplot(width = 0.5) +
  theme_bw() + scale_color_manual(values=myPalette) +
  theme(
    panel.grid = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1)
    ) + 
  xlab("SampleID") +
  ylab(paste0("log10(", toupper(index), ")"))

p2 <- ggplot(data, aes(Group, log10(FPKM), color = Group)) + 
  geom_boxplot() +
  theme_bw() + scale_color_manual(values=myPalette) +
  theme(
    panel.grid = element_blank(), 
    axis.text.x = element_text(angle = 45, hjust = 1)
    ) + 
  xlab("GroupID") +
  ylab(paste0("log10(", toupper(index), ")"))

ggsave(snakemake@output[["fig"]][1], p1, width = 8, height = 6)
ggsave(snakemake@output[["fig"]][2], p1, width = 8, height = 6)
ggsave(snakemake@output[["fig1"]][1], p2, width = 8, height = 6)
ggsave(snakemake@output[["fig1"]][2], p2, width = 8, height = 6)