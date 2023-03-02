
Args <- commandArgs(TRUE)

read_distribution = Args[1]
out_pdf = Args[2]
out_png = Args[3]

library(ggplot2)
library(RColorBrewer)

get_name <- strsplit(read_distribution,'/')[[1]][4]
sample <- strsplit(get_name,'[.]')[[1]][1]
data<-read.table(read_distribution, header = FALSE, sep = "\t")
data<-as.matrix(data)
# as.numeric(data[,2])
exonratio<-(as.numeric(data[2,2])+as.numeric(data[3,2])+as.numeric(data[4,2]))/as.numeric(data[1,2])
intronratio<-as.numeric(data[5,2])/as.numeric(data[1,2])
interratio<-(as.numeric(data[6,2])+as.numeric(data[7,2]))/as.numeric(data[1,2])
exonratio<-round(exonratio, 4)*100
interratio<-round(interratio, 4)*100
intronratio<-round(intronratio, 4)*100
number <- c( exonratio , interratio , intronratio)
Type <- c( paste("exon", " ",  exonratio, "%", sep = ""),
            paste("intergenic", " ", interratio, "%", sep = ""),
            paste("intron", " ", intronratio, "%", sep = ""))
df <- data.frame(Type, number)
bp<- ggplot(df, aes(x="", y=number, fill=Type))+labs(x = "",y ="") +
  geom_bar(width = 1, stat = "identity")

blank_theme <- theme_minimal()+
  theme(
  axis.title.x = element_blank(),
  axis.title.y = element_blank(),
  panel.border = element_blank(),
  panel.grid=element_blank(),
  axis.ticks = element_blank(),
  plot.title=element_text(size=14, face="bold",hjust = 0.5)
  )
pie <- bp +  coord_polar("y", start=0)
p <- pie + scale_fill_brewer(palette="Set2",direction = -1) +  blank_theme +
  theme(axis.text.x=element_blank())+
    ggtitle(sample)
ggsave(out_pdf, p, width = 6, height = 6)
ggsave(out_png, p, width = 6, height = 6)

# exon<-(as.numeric(data[2,2])+as.numeric(data[3,2])+as.numeric(data[4,2]))
# intron<-as.numeric(data[5,2])
# inter<-(as.numeric(data[6,2])+as.numeric(data[7,2]))

# count <- c( exon , inter , intron)
# type <- c('exon','inter','intron')
# df <- data.frame(type = type, count = count)

# p <- ggplot(data = df, mapping = aes(x = '', y = count, fill = type)) + geom_bar(stat = 'identity', position = 'stack', width = 1)

# pie <- bp + coord_polar(theta = 'y') + labs(x = '', y = '', title = '')  + theme(axis.text = element_blank()) + scale_fill_brewer(palette="Set2",direction = -1)

# pdf(file= out_pdf)

# dev.off()

# pie(count, labels=type, col=cols)