library(tidyr)
library(data.table)
library(pheatmap)
library(RColorBrewer)
library(argparser)
library(ggplot2)

height <- snakemake@params[["height"]]
width <- snakemake@params[["width"]]
data <- fread(snakemake@input[["matrix"]])

metadata <- fread(snakemake@input[["metadata"]]) %>% as.data.frame

# rownames(data) <- data[["Gene ID"]]
# data <- data[,-1]
# data <- data[rowSums(data)>0, ]

# sim <- cor(data,method="pearson")
# sim2 <- as.data.frame(sim) %>% rownames_to_column(" ")
# write_tsv(sim2, snakemake@output[[1]])
#color <- colorRampPalette(c("blue", "white", "red"))(256)
# color <- colorRampPalette(brewer.pal(11,'RdBu')[10:2])(256)
# p <- pheatmap(sim, color = color, border_color = NaN)

dt <- data[, -1]
dt <- dt[rowSums(dt) > 0, ]
sim <- cor(dt, method = "pearson") %>% as.data.frame
write.table(sim, snakemake@output[[1]], sep = "\t", quote = F)
# plot heatmap
row.anno <- data.frame(metadata$group, row.names = metadata$SampleID)
colnames(row.anno) <- "group"
#myPalette <- c("#FF4040", "#228B22", "#FF7F00", "#6959CD")
#ann_colors = list(group = myPalette)
#names(ann_colors$group) <- levels(row.anno$group)

p <- pheatmap(sim, border_color = NaN, annotation_row = row.anno, 
              color = colorRampPalette(brewer.pal(11,'Reds')[2:9])(256),
              cluster_col=FALSE,
              cluster_row=FALSE,
              display_numbers = TRUE
              # treeheight_col = 0
              #annotation_colors = ann_colors
              )

p2 <- pheatmap(sim, border_color = NaN, annotation_row = row.anno, 
              color = colorRampPalette(brewer.pal(11,'Reds')[2:9])(256),
              cluster_col=TRUE,
              cluster_row=TRUE,
              display_numbers = TRUE
              # treeheight_col = 0
              #annotation_colors = ann_colors
              )

p3 <- pheatmap(sim, border_color = NaN, annotation_row = row.anno, 
              color = colorRampPalette(brewer.pal(11,'Reds')[2:9])(256),
              cluster_col=FALSE,
              cluster_row=FALSE,
              display_numbers = FALSE
              # treeheight_col = 0
              #annotation_colors = ann_colors
              )

p4 <- pheatmap(sim, border_color = NaN, annotation_row = row.anno, 
              color = colorRampPalette(brewer.pal(11,'Reds')[2:9])(256),
              cluster_col=TRUE,
              cluster_row=TRUE,
              display_numbers = FALSE
              # treeheight_col = 0
              #annotation_colors = ann_colors
              )

ggsave(snakemake@output[["fig"]][1], p, height = height, width = width)
ggsave(snakemake@output[["fig"]][2], p, height = height, width = width)
ggsave(snakemake@output[["fig1"]][1], p2, height = height, width = width)
ggsave(snakemake@output[["fig1"]][2], p2, height = height, width = width)

ggsave(snakemake@output[["fig3"]][1], p3, height = height, width = width)
ggsave(snakemake@output[["fig3"]][2], p3, height = height, width = width)
ggsave(snakemake@output[["fig4"]][1], p4, height = height, width = width)
ggsave(snakemake@output[["fig4"]][2], p4, height = height, width = width)