library(dplyr)
library(readr)
library(GenomicRanges)
library(ggbio)
library(RColorBrewer)
library(patchwork)
args <- commandArgs(T)

if(length(args) == 9) {
  width = 10
  height = 8
} else if(length(args) == 11) {
  width = as.numeric(args[9])
  height = as.numeric(args[10])
} else {
  writeLines("
------------------------------------------------------------------------------------------------------------
Usage:Rscript *.R gene_fpkm_anno sample-metadata chrom_info M2_vs_B1 M2_vs_B1_diff circle_pdf  circle_png
exp: Rscript ggbio.R  gene_fpkm.xls  sample-metadata.txt chrom_info.csv M2_vs_B1 M2_vs_B1_diff_diffanno.xls circle.pdf  circle.png
------------------------------------------------------------------------------------------------------------
")
  q()
}
# ext_file <- "/nfs2/public2/User/lixiaofei/projects/mRNA/20210720-HT20210621236-chenguangxia/genome_file/chr_info.xls"
fpkm_table <- args[1]
sample_table <- args[2]
chrom_info <- args[3]
ext_file <- args[4]
comparesname <- args[5]
diff_table <- args[6]
ggbio_num <- args[7]
circle_pdf <- args[8]
circle_png <- args[9]

ext_table<-read_tsv(ext_file)
# data<-read_tsv("gene_fpkm.xls")
fpkm<-read_tsv(fpkm_table)
fpkm<-merge(fpkm, ext_table, all.x=TRUE, by="gene_id")
# anno<-read_tsv("hsa_geneid2description.xls")
# aa<-merge(a,anno,all.x=T)

data <- fpkm %>% select(chr,Start,End,Strand,gene_name,gene_id,everything())

group <- read_tsv(sample_table)

comparesname <- comparesname
compares <- strsplit(comparesname,'_vs_')[[1]]

for (i in group$group) {
    data[[i]] <- log2(rowMeans(data[,group$SampleID[which(i == group$group)]])+1)
}

chro <- read.csv(chrom_info,header = F)
names(chro) <- c("chr","start","end")
chr<-as.data.frame(chro[,1])
names(chr) <- "chr"
mRNA<-merge(chr,data,all.x=T)
head(mRNA)
# chro<- arrange(chro,chr)
# chro$chr <- factor(chro$chr, levels = chro$chr)

gmrna <- GRanges(mRNA)
head(gmrna)
gchro <- GRanges(chro)
head(gchro)
seqlengths(gchro) <- chro$end
seqlengths(gchro)
seqlengths(gmrna) <- seqlengths(gchro)

num <- as.numeric(ggbio_num)
#绘制外圈
p <- ggplot() + layout_circle(gchro, geom = "ideo", fill = rainbow(num), radius = 30, trackWidth = 1)
#添加刻度
p <- p + layout_circle(gchro, geom = "scale", size = 1, radius = 31, trackWidth = 0.5)
#添加染色体名称
p <- p + layout_circle(gchro, geom = "text", aes(label = seqnames), vjust = 0,size = 3, radius = 32, trackWidth = 2)
#绘制实验组圈图
p1 <- p + layout_circle(gmrna, geom = "bar", aes_string(y = compares[1]),size=0.1,color = "#E41A1C", fill = "#E41A1C",radius = 25,trackWidth = 5)
#绘制对照组圈图
p2 <- p1 + layout_circle(gmrna, geom = "bar",aes_string(y= compares[2]),size=0.1, color = "#377EB8", fill = "#377EB8",radius = 20,trackWidth = 5)

#读入差异基因，使用logFC绘制圈图
diff <- read_tsv(diff_table)
diff <- merge(diff, ext_table, all.x=TRUE, by="gene_id")
# diff <- diff[level %in% c("Increased","Decreased"),]
diff <- diff[which(diff$level!="nonsignificant"),]
# anno<- fread(gene_chr_anno)
# diffanno<-merge(diff,anno,all.x=T)
diff1<-diff[,c("chr","Start",'End','logFC')]

# setdiff(chro$chr,unique(diff1$chr))
d1 <- chro[which(chro$chr %in% setdiff(chro$chr,unique(diff1$chr))),]
d1$logFC <- rep(0,time=nrow(d1))
names(d1) <- names(diff1)
d2<-rbind(diff1,d1)
d3<-merge(chr,d2,all.x=T)

gmrna <- GRanges(d3)
seqlengths(gmrna) <- seqlengths(gchro)
p3 <- p2 + layout_circle(gmrna, geom = "bar",aes(y=logFC),size=0.1, color = "#4DAF4A", radius = 15,trackWidth = 5)

#添加图例
col<-c("#4DAF4A", "#377EB8","#E41A1C" )
names(col) <- c("DE-mRNA",compares[2],compares[1])
# pdf(circle_pdf)
p4<-p3+ggplot() + annotate("point", x=1,y=1:3,shape=15, size=3.5,color=col) +annotate("text", x=1.01, y=1:3, label=names(col), hjust=0) +  xlim(0.99, 1.1) + ylim(-2, 10) + theme_void()
# dev.off()
ggsave(circle_pdf,p4,width = width, height = height)
ggsave(circle_png,p4,width = width, height = height)