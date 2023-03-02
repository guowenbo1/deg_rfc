library(tidyverse)
library(pheatmap)
library(RColorBrewer)

contrast <- snakemake@wildcards[["contrast"]]

diff_lst <- read_tsv(snakemake@input[["diff_lst"]], col_names = FALSE) %>% filter(X1==contrast)
group <- diff_lst[[2]]

fpkm <- read_tsv(snakemake@input[["fpkm"]])
colnames(fpkm)[1] <- "gene_id"
fpkm <- fpkm %>% column_to_rownames("gene_id")

metadata <- read_tsv(snakemake@input[["metadata"]]) %>%   select_("SampleID", group) %>% .[.[["group"]] %in% diff_lst[3:4],]


deg <- read_tsv(snakemake@input[["deg"]])
deg <- filter(deg, level %in% c("Increased", "Decreased")) %>% .[["gene_id"]]
fpkm <- fpkm[rownames(fpkm) %in% deg, metadata$SampleID]

annotation_col <- metadata %>% column_to_rownames("SampleID")
# fpkm <- fpkm[which(rowSums(fpkm) > 0),]
p <- pheatmap(fpkm,scale = "row",fontsize_row = 10,fontsize_col=10,border_color = NA,
         color <- colorRampPalette(brewer.pal(11,'RdBu')[10:2])(256),
         annotation_col = annotation_col, angle_col = "45",show_rownames = FALSE, main = contrast)

num=length(fpkm[,1])


if (num>=100){
  height1 <- 10  
  }else if (num>=70) {
    height1 <- 9
  }else if (num>=40) {
    height1 <- 9   
  }else if (num>=30) {
    height1 <- 8
  }else if (num>=20) {
    height1 <- 6
  }else {
    height1 <- 5
  }


if (num>=100){
  height2 <- 0.015*num+9
  }else if (num>=80) {
    height2 <- 10  
  }else if (num>=60) {
    height2 <- 9.5
  }else if (num>=40) {
    height2 <- 9   
  }else if (num>=30) {
    height2 <- 8
  }else if (num>=20) {
    height2 <- 6
  }else {
    height2 <- 5
  }

if (num>=100){
  size <- 1
  }else if (num>=80) {
    size <- 3
  }else if (num>=70) {
    size <- 3.5      
  }else if (num>=60) {
    size <- 4    
  }else if (num>=40) {
    size <- 5 
  }else if (num>=30) {
    size <- 6     
  }else if (num>=20) {
    size <- 7.5
  }else {
    size <- 8
  }

p2 <- pheatmap(fpkm,scale = "row",fontsize_row = size,fontsize_col=10,border_color = NA,
         color <- colorRampPalette(brewer.pal(11,'RdBu')[10:2])(256),
         annotation_col = annotation_col, angle_col = "45",show_rownames = T, main = contrast)

# ggsave(snakemake@output[["fig1"]], p, height = height1, width = 8)
ggsave(snakemake@output[["fig2"]], p2, height = height2, width = 8, limitsize = FALSE)
