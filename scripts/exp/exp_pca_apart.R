library(tidyverse)
library(ggrepel)
library(factoextra)
args <- commandArgs(T)

exp_file <- args[1]
sample_file <- args[2]
label <- args[3]
case <- args[4]
out_file <- args[5]
out_normal_png <- args[6]
out_normal_pdf <- args[7]
out_label_png <- args[8]
out_label_pdf <- args[9]
out_ellipse_png <- args[10]
out_ellipse_pdf <- args[11]

myPalette <- c("#FF4040","#228B22","#984EA3","#FF7F00","#6959CD","#F781BF","#458B74","#A52A2A","#8470FF","#8B4513","#556B2F","#CD5B45","#483D8B","#EEC591","#8B0A50","#696969", "#8B6914","#20B2AA", "#8B636C","#36648B","#9ACD32","#8B3A62","#68838B","#CDBE70","#D3D3D3","#EEA2AD","#53868B","#EE9572","#FFA500","#8B4789","#548B54","#F4A460","#3A5FCD","#B0E0E6","#6495ED","#EECBAD","#CD69C9","#436EEE","#8B8B00","#8B7E66","#CD853F","#8B7B8B","#00C5CD","#008B45","#D2B48C","#CDCD00","#836FFF","#7D26CD","#006400","slateblue2","slategray3","violetred2","#00868B","palegreen1","#7A8B8B","palevioletred2","powderblue","orangered2","lightgoldenrod2","deepskyblue", "cyan1","midnightblue","seashell2","yellowgreen","tomato3","violetred4","magenta","chocolate2","darkorange3","olivedrab2","#8B3626","#0000FF")

index <- label
group <- "group"

data <- read_tsv(exp_file) %>% column_to_rownames("Gene ID")
metadata <- read_tsv(sample_file)
compares <- strsplit(case, "_vs_")
included_group <- metadata[which(metadata$group %in% compares[[1]]), ]


data <- data[rowSums(data)>0, ] %>% select(included_group$SampleID)
data <- t(data)
data <- data[ , which(apply(data, 2, var) != 0)]

pca <- prcomp(data, center = T, scale. = T)

pca_x <- as.data.frame(pca$x)
plotdata <- pca_x %>% rownames_to_column("SampleID") %>%
  left_join(metadata) %>% select(-fq1, -fq2)

write_tsv(plotdata, out_file)

eig.val <- get_eigenvalue(pca)
e1 <- eig.val[1,2]
e2 <- eig.val[2,2]


p1 <- ggplot(plotdata, aes_string(x = "PC1", y = "PC2", color = group)) +
  geom_point(size = 3.5) +
  xlab(paste("PC1: ", round(e1, digits = 2), "%")) +
  ylab(paste("PC2: ", round(e2, digits = 2), "%")) +
  scale_color_manual(values=myPalette, name = NULL) +
  theme_bw() +
  ggtitle(paste(index, "PCA")) +
    theme(panel.grid=element_blank(), 
      axis.text.x = element_text(face="bold", color="black"),
      axis.text.y = element_text(face="bold",  color="black"),
      axis.title.x = element_text(face="bold", color="black"),
      axis.title.y = element_text(face="bold",color="black"), 
      plot.title = element_text(lineheight=2.5, face="bold",hjust=0.5))

ggsave(out_normal_png, p1, width = 6.5, height = 5)
ggsave(out_normal_pdf, p1, width = 6.5, height = 5)

p2 <- ggplot(plotdata, aes_string(x = "PC1", y = "PC2", color = group)) +
    geom_point(size = 3.5) +
    xlab(paste("PC1: ", round(e1,digits = 2), "%")) +
    ylab(paste("PC2: ", round(e2,digits = 2), "%")) +
    scale_color_manual(values=myPalette, name = NULL) +
    theme_bw() + geom_text_repel(aes(label = SampleID)) +
    ggtitle(paste(index, "PCA"))+
    theme(panel.grid=element_blank(), 
      axis.text.x = element_text(face="bold", color="black"),
      axis.text.y = element_text(face="bold",  color="black"),
      axis.title.x = element_text(face="bold", color="black"),
      axis.title.y = element_text(face="bold",color="black"), 
      plot.title = element_text(lineheight=2.5, face="bold",hjust=0.5))

ggsave(out_label_png, p2, width = 6.5, height = 5)
ggsave(out_label_pdf, p2, width = 6.5, height = 5)

p3 <- ggplot(plotdata, aes_string(x = "PC1", y = "PC2", color = group)) +
    geom_point(size = 3.5) +
    xlab(paste("PC1: ", round(e1, digits = 2), "%")) +
    ylab(paste("PC2: ", round(e2, digits = 2), "%")) +
    scale_color_manual(values=myPalette, name = NULL) +
    theme_bw() + stat_ellipse(level = 0.75) +
    ggtitle(paste(index, "PCA"))+
    theme(panel.grid=element_blank(), 
      axis.text.x = element_text(face="bold", color="black"),
      axis.text.y = element_text(face="bold",  color="black"),
      axis.title.x = element_text(face="bold", color="black"),
      axis.title.y = element_text(face="bold",color="black"), 
      plot.title = element_text(lineheight=2.5, face="bold",hjust=0.5))

ggsave(out_ellipse_png, p3, width = 6.5, height = 5)
ggsave(out_ellipse_pdf, p3, width = 6.5, height = 5)