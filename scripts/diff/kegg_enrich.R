library(clusterProfiler)
library(data.table)
library(ggplot2)
library(viridis)
library(tidyr)
library(dplyr)
library(readr)
library(RColorBrewer)
label <- snakemake@wildcards[["label"]]
contrast <- snakemake@wildcards[["contrast"]]

kegg <- fread(snakemake@input[["kegg_anno"]])
deg <- fread(snakemake@input[["deg"]])

if (label == "all") {
  # gene_list <- filter(deg, level %in% c("Increased","Decreased")) %>% .[["gene_id"]]
  gene_list <- deg[level %in% c("Increased","Decreased"),]$gene_id
}else if (label == "up") {
  # gene_list <- filter(deg, level== "Increased") %>% .[["gene_id"]]
  gene_list <- deg[level %in% c("Increased"),]$gene_id
}else if (label == "down") {
  # gene_list <- filter(deg, level== "Decreased") %>% .[["gene_id"]]
  gene_list <- deg[level %in% c("Decreased"),]$gene_id
}
if(label == "all"){
  color <- colorRampPalette(brewer.pal(9,'Blues')[3:9])(256)
}else if (label == "down"){
  color <- colorRampPalette(brewer.pal(9,'Greens')[3:9])(256)
}else if (label == "up"){
  color <- colorRampPalette(brewer.pal(9,'Reds')[3:9])(256)
}
kegg_enrich <- enricher(gene_list, pvalueCutoff = 1, qvalueCutoff = 1, 
                        TERM2GENE = kegg[,c("level3_pathway_id","gene_id")],
                        TERM2NAME = kegg[,c("level3_pathway_id","level3_pathway_name")])

res <- as.data.table(kegg_enrich)
fwrite(res, snakemake@output[["res"]], sep = "\t")

res_sig <- res %>% arrange(pvalue) %>% dplyr::slice(1:20)
res_sig <- arrange(res_sig,Count)

# 宽度
max_length <- 0
for(i in res_sig$Description){
  max_length <- ifelse(max_length > nchar(as.character(i)), max_length, nchar(as.character(i)))
}
my_width <- 0.1 * max_length + 6

res_sig <- tidyr::unite(res_sig, "des", ID, Description,sep = ":", remove = FALSE)
res_sig$des <- factor(res_sig$des, levels = res_sig$des)
res_sig$shape <- rep(" ",times = length(rownames(res_sig)))
res_sig$shape[which(res_sig$pvalue < 0.05)] <- "*"
#res_sig$Description <- factor(res_sig$Description, levels = rev(res_sig$Description))

ggplot(res_sig, aes(x= des, y = Count, fill = -log10(pvalue))) + geom_col(colour="#666666") +
  theme_bw() + coord_flip()  + 
  xlab("") + ylab("Num of Genes") + ggtitle(title) +
  labs(title= paste("KEGG pathway", contrast, label)) +
  theme(panel.grid = element_blank()) + 
  theme(axis.title.x =element_text(size=14), axis.title.y=element_text(size=14),legend.title=element_text(size=13),axis.text=element_text(size=12,color  = "black")) + 
  scale_fill_gradientn(colours=color) +
  geom_text(aes(label = shape, vjust = 0.5, hjust = -0.2))
ggsave(snakemake@output[["barplot"]][1], width = my_width, height = 6.5)
ggsave(snakemake@output[["barplot"]][2], width = my_width, height = 6.5)

num_level3 <- kegg[, .N, by = level3_pathway_id]
a<-as.data.table(num_level3)
dt <- res_sig
names(a) <- c("level3_pathway_id","N")
dt <- merge(dt, a, by.x = "ID", by.y = "level3_pathway_id", all.x = T)
dt$Ratio <- dt$Count/dt$N
dt<- arrange(dt,Count)
dt$des <- factor(dt$des, levels = dt$des)

ggplot(dt, aes(x= des, y = Ratio)) + 
  geom_point(aes(color = -log10(pvalue), size = Count)) +
  theme_bw() + coord_flip() +  xlab("")  +
  ylab("Ratio") + 
  labs(title= paste("KEGG pathway", contrast, label)) +
  theme(axis.title.x =element_text(size=14), axis.title.y=element_text(size=14),legend.title=element_text(size=13),axis.text=element_text(size=12,color  = "black")) + 
  scale_color_gradientn(colours=color) + 
  # 扩大点的大小 
  scale_size_continuous(range=c(3,8)) +
  # 固定图例位置
  guides(color = guide_colorbar(order = 0),size = guide_legend(order = 1))

ggsave(snakemake@output[["dotplot"]][1], width = my_width, height = 6.5)
ggsave(snakemake@output[["dotplot"]][2], width = my_width, height = 6.5)
