# library(clusterProfiler)
# library(tidyverse)
# library(viridis)
# library(tidyr)
library(clusterProfiler)
library(data.table)
library(ggplot2)
library(viridis)
library(tidyr)
library(dplyr)
library(readr)

# label <- snakemake@wildcards[["label"]]
# contrast <- snakemake@wildcards[["contrast"]]
# # title <- paste0(contrast, "/", label)
# go <- fread(snakemake@input[["go_anno"]])
# deg <- fread(snakemake@input[["deg"]])
Args <- commandArgs(TRUE)
label = Args[1]
contrast = Args[2]
go = Args[3]
deg = Args[4]
output = Args[5]
barplot = Args[6]
barplot2 = Args[7]
dotplot = Args[8]
dotplot2 = Args[9]

go<-fread(go)
deg<-fread(deg)

go2ont <- go[,c(3,5)] %>% unique()

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

go_enrich <- enricher(gene_list, pvalueCutoff = 1, qvalueCutoff = 1, 
                      TERM2GENE = go[,c(3,1)],TERM2NAME = go[,c(3,6)])


res <- as.data.table(go_enrich) %>% na.omit() %>% left_join(go2ont, by = c("ID"="GOID"))

fwrite(res, output, sep = "\t")
# write_tsv(res, snakemake@output[["res"]])

# res_sig <- res %>% group_by(ONTOLOGY) %>% arrange(ONTOLOGY, desc(Count)) %>% 
#   dplyr::slice(1:10)

res_sig <- res %>% group_by(ONTOLOGY) %>% arrange(ONTOLOGY, p.adjust) %>% dplyr::slice(1:10)
res_sig2 <- arrange(res_sig,desc(Count))

res_sig2$Description <- factor(res_sig2$Description, levels = rev(res_sig2$Description))
# res_sig$ID <- factor(res_sig$ID, levels = rev(res_sig$ID))
p <- ggplot(res_sig2, aes(x= Description, y = Count, fill = -log10(p.adjust))) + geom_col(colour="#666666") +
  theme_bw() + coord_flip() + facet_grid(ONTOLOGY~.,scales="free") +   
  xlab("") + ylab("Num of Genes") + ggtitle(title) +
  labs(title= paste("GO term", contrast, label)) +
  theme(panel.grid = element_blank()) + 
  theme(axis.title.x =element_text(size=14), axis.title.y=element_text(size=14),legend.title=element_text(size=13),axis.text=element_text(size=12,color  = "black")) + 
  scale_fill_gradientn(colours=viridis_pal(option = "D")(20)[5:20]) 

ggsave(barplot, p, width = 16, height = 10)
ggsave(barplot2, p, width = 16, height = 10)

# a<- res_sig %>% separate(GeneRatio, c("n1", "n2"), "/")
# a$n1 <- as.numeric(a$n1)
# a$n2 <- as.numeric(a$n2)
# a$GeneRatio <- a$n1/a$n2

gonumber <- go[, .N, by = GOID]
dt <- merge(res_sig2, gonumber, by.x = "ID", by.y = "GOID", all.x = T)
dt$Ratio <- dt$Count/dt$N
dt<- arrange(dt,ONTOLOGY,Ratio)
dt$Description <- factor(dt$Description, levels = dt$Description)
p2<-ggplot(dt, aes(x= Description, y = Ratio)) + 
  geom_point(aes(color = -log10(p.adjust), size = Count)) +
  theme_bw() + coord_flip() + 
  facet_grid(ONTOLOGY~.,scales="free") + 
  xlab("") + ylab("Ratio") + ggtitle(title) +
  labs(title= paste("GO term", contrast, label)) +
  theme(panel.grid = element_blank()) + geom_point(aes(size = Count),shape=21) +
  theme(axis.title.x =element_text(size=15), axis.title.y=element_text(size=15),legend.title=element_text(size=15),axis.text=element_text(size=14,color  = "black")) + 
  scale_color_gradientn(colours=viridis_pal(option = "D")(20)[5:20]) 

ggsave(dotplot, p2, width = 18, height = 10)
ggsave(dotplot2, p2, width = 18, height = 10)

