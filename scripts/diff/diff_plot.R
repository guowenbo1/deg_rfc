#! /nfs1/public2/User/wzc/miniconda3/envs/rnaseq/bin/Rscript
library(tidyverse)
library(ggsci)

Args <- commandArgs(TRUE)

method = Args[1]
contrast = Args[2]
rdata = Args[3]
type = Args[4]
logfc_thres = as.numeric(Args[5])
value_thres = as.numeric(Args[6])
out_v = Args[7]
out_ma = Args[8]

data <- read_tsv(rdata)

dataup <- data[which(data$level=="Increased"), ]
datadown <- data[which(data$level=="Decreased"), ]
datano <- data[which(data$level=="nonsignificant"), ]

uplable = paste0("Increased:", nrow(dataup))
downlable = paste0("Decreased:", nrow(datadown))
nonsignificant = paste0("nonsignificant", nrow(datano))

dataup$sig <- uplable
datadown$sig <- downlable
datano$sig <- "FALSE"
datanew <- rbind(datano, dataup, datadown)


if(nrow(dataup) > 0 & nrow(datadown) > 0) {
    manual_color <- c("#009966", "grey", "#A8152B")
}
if(nrow(dataup) > 0 & nrow(datadown) == 0) {
    manual_color <- c( "#A8152B", "grey")
}
if(nrow(dataup) == 0 & nrow(datadown) > 0) {
    manual_color <- c("#009966", "grey")
}
if(nrow(dataup) == 0 & nrow(datadown) == 0) {
    manual_color <- c("grey")
}


ggplot(datanew,aes_string(x = 'logFC', y = paste0('-log10(', type, ')'), colour = 'sig')) +
    geom_point(alpha = 0.9, size = 1)+
	ggtitle(contrast) +
    scale_color_manual("", breaks=c(uplable,downlable,NA), values = manual_color)+
	geom_vline(xintercept = c(-logfc_thres, logfc_thres),lty = 4,col = "grey",lwd = 0.6)+
    geom_hline(yintercept = -log10(as.numeric(value_thres)),lty = 4,col = "grey",lwd = 0.6)+
    ggtitle(contrast) +
    theme_bw() + 
	guides(color=guide_legend(override.aes = list(size=1.5)))+
	theme(panel.grid=element_blank(),legend.title = element_blank(),
    axis.text.x = element_text(face="bold", color="black"),
    axis.text.y = element_text(face="bold",  color="black"),
    axis.title.x = element_text(face="bold", color="black"),
    axis.title.y = element_text(face="bold",color="black"))+
    labs(x="log2(fold change)",y=paste0('-log10(', type, ')'))
ggsave(out_v, height = 5, width = 6.25)


if (method == "DESeq2") {
  sum <- sum(data$baseMean)
  ggplot(datanew) +
    geom_point(aes(x = log2(baseMean*10^6/sum), y = logFC, colour = sig), size = 1,alpha=0.8)+
    xlab("logCPM") +
	ggtitle(contrast) +
	guides(color=guide_legend(override.aes = list(size=1.5)))+
    scale_color_manual("", breaks=c(uplable,downlable,NA), values = manual_color)+
    theme_bw() + theme(panel.grid=element_blank(),legend.title = element_blank(),
    axis.text.x = element_text(face="bold", color="black"),
    axis.text.y = element_text(face="bold",  color="black"),
    axis.title.x = element_text(face="bold", color="black"),
    axis.title.y = element_text(face="bold",color="black"))
}else {
  ggplot(datanew) +
    geom_point(aes(x = logCPM, y = logFC, colour = sig), size = 1,alpha=0.8)+
    scale_color_manual("", breaks=c(uplable,downlable,NA), values = manual_color)+
	guides(color=guide_legend(override.aes = list(size=1.5)))+
    ggtitle(contrast) +
    theme_bw() + theme(panel.grid=element_blank(),legend.title = element_blank(),
    axis.text.x = element_text(face="bold", color="black"),
    axis.text.y = element_text(face="bold",  color="black"),
    axis.title.x = element_text(face="bold", color="black"),
    axis.title.y = element_text(face="bold",color="black"))
}
ggsave(out_ma, height = 5, width = 6.25)
