library(argparse)
library(data.table)
library(tidyverse)
library(forcats)
library(ggprism)
library(labeling)

parser <- ArgumentParser(description='statistics and visualization of DEG')
parser$add_argument('-i','--input',type='character',nargs="*",
        help="DESeq2.DE_results.filter.DEG.tsv *< one or more files>")
parser$add_argument('-o','--output',type='character',
        default='.',help="output directory,[./]")
parser$add_argument('-n','--name',type='character',
        help="output file prefix,[DEG.stat]")
parser$add_argument('-t','--title',type='character',nargs="*",
        default="Statistics of All Differentially Expressed Genes",
        help="plot title,[\"Statistics of All Differentially Expressed Genes\"]")

args <- parser$parse_args()
if(is.null(args$input)){
        parser$print_help()
        q()
}

file <- unique(args$input)
# print(file)
write(paste0("Total files Numbers: ",n_distinct(file)),stderr())
dir.create(args$output,showWarnings = F)
output <- normalizePath(args$output)
prefix <- ifelse(is.null(args$name),"DEG.stat",args$name) 

alldata = NULL
for (i in 1:length(file)) {
    filename = basename(file[i])
    message(paste0("parsing: ",i," ",file[i]))
    data <- fread(file[i],header=T,fill=TRUE,sep="\t") %>% .[,"level"]
    data$plan <- str_split(filename,'_diff',2,simplify = T)[1]
    alldata <- rbind(alldata, data)
}

theme_set(theme_prism(border=TRUE)+
          theme(
               panel.border = element_blank(),
               axis.line.y = element_line(color='black',size=.6),
               axis.ticks.x = element_blank(),
               axis.ticks.y = element_line(size=.6),
               plot.title = element_text(face = 2,size = 20,hjust = 0.5),
               axis.text.y = element_text(face = 1),
               axis.title.x = element_blank(),
               axis.title.y = element_text(face = 2,size = 15)
          )
)

# write.csv(alldata,"alldata.csv")

df <- alldata %>%
	count(plan,level) %>%
	filter(level %in% c("Increased","Decreased"))
df$level <- factor(df$level,levels=c("Increased","Decreased"))
df[df$level=="Increased",2] <- "Increased"
df[df$level=="Decreased",2] <- "Decreased"

group <- unique(df$plan)

for (i in group){
  a<- df[which(df$plan == i),]
  a
  if (dim(a)[1] !=2){
    if (a$level %in% "Increased"){
      tmp <- data.frame(plan=i,level="Decreased",n=as.numeric(0))
      df <- rbind(df,tmp)
    }    
    else{
      tmp <- data.frame(plan=i,level="Increased",n=as.numeric(0))
      df <- rbind(df,tmp)  
    }
  }
}

max <- max(df$n)*1.2
tmp = with(df, labeling::extended(0, range(n)[2]*1.2, m = 6))
lm = tmp[c(1,length(tmp))]

p <- ggplot(df,aes(x=plan,y=n,fill=level))+
	geom_col(position="dodge",width=0.6,color="black",key_glyph=draw_key_rect)+
	scale_fill_manual(values=c("#A8152B","#009966"))+
	scale_y_continuous(expand=c(0,0),limits = lm,breaks = tmp,guide = "prism_minor",
                       labels=function(x) format(x, big.mark = "'", scientific = FALSE)
               		  )+
	expand_limits(y=max)+
	labs(y="Gene Numbers")+
	ggtitle(paste0(args$title,collapse=" "))+
	guides(colour = guide_legend(override.aes = list(size = 2)))

if(n_distinct(df$plan) >= 6){
    p <- p + theme(axis.text.x = element_text(face = 2,size = 10,angle = 45,hjust = 1, vjust = 1))
    p <- p + geom_text(aes(label=prettyNum(n,big.mark=","),y=n),
             		size = 4, colour = 'black',angle=90,
             		hjust= -.2,vjust = .5,position=position_dodge(0.5)
             )
} else {
    p <- p + theme(axis.text.x = element_text(face = 2,size = 10,angle = 0,hjust = .5, vjust = .5))
    p <- p + geom_text(aes(label=prettyNum(n,big.mark=","),y=n),
             		size = 4, colour = 'black',
             		hjust=0.5,vjust = -.5,position=position_dodge(0.6)
             )
}                   
ggsave(paste0(output,"/",prefix,'.barplot.pdf'),p,width = 30, height = 20,units = "cm")
print(paste0("Saved to ",output,"/",prefix,'.barplot.pdf'))
ggsave(paste0(output,"/",prefix,'.barplot.png'),p,width = 30, height = 20,units = "cm",dpi=1000)
print(paste0("Saved to ",output,"/",prefix,'.barplot.png'))
df <- df %>%
	setnames(c("Plan","State","Number")) %>%
	mutate_if(is.numeric,prettyNum,big.mark=",") %>%
	arrange(Plan,State) 
fwrite(df,paste0(output,"/",prefix,'.tsv'),sep="\t")
print(paste0("Saved to ",output,"/",prefix,'.tsv'))
