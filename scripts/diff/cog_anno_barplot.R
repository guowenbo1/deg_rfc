library(dplyr)
library (tidyr)
library(ggplot2)
library(data.table)
library(argparser)

argv <- arg_parser("COG_anno_barplot") %>%
  add_argument("--contrast", help="level table",default = "A_vs_B") %>%
  add_argument("--anno", help="level table",default = "COG_annotations.txt") %>%
  add_argument("--deg", help="level table",default = "A_vs_B_DEGs_list.xls") %>%
  add_argument("--out", help="outdir",default = "COG") %>%
  parse_args()

contrast<- argv$contrast
cog_anno<-fread(argv$anno)
deg_list<-fread(argv$deg)
names(deg_list)<-"gene_ID"

cog<- cog_anno %>% filter(gene_ID %in% deg_list$gene_ID)

fwrite(cog,paste0(argv$out,"/",contrast,"_COG_annotations.xls"),sep="\t")

if (dim(cog)[1] == 0){
  pdf(paste0(argv$out,"/",contrast,"_COG_anno.pdf"))
  plot(1:5, 1:5, xlim = c(0,6), ylim = c (0,6),axes = F,type = "n",xlab="", ylab="")+text(x = 3, y = 3, labels = "result doesn't exist")
  dev.off()
  png(paste0(argv$out,"/",contrast,"_COG_anno.png"))
  plot(1:5, 1:5, xlim = c(0,6), ylim = c (0,6),axes = F,type = "n",xlab="", ylab="")+text(x = 3, y = 3, labels = "result doesn't exist")
  dev.off()
}
if (dim(cog)[1] != 0){
cog$COGFunctionalCategory<-substr(cog$COGFunctionalCategory,start=1,stop=1)
cog$anno<-paste(cog$COGFunctionalCategory,cog$COGFunctionalCategoryDescription,sep=":")

classnumber <- cog[, .N, by = anno]
names(classnumber) <- c("COGFunctional","genenumber")

data<-separate(classnumber, COGFunctional, c("number", "class"), sep = ":", remove = FALSE)

data<- arrange(data,genenumber)
data$COGFunctional <- paste(data$COGFunctional,"(",data$genenumber,")",sep ="")

p<-ggplot(data=data, aes(x=number, y=genenumber, fill=COGFunctional)) +
   geom_bar(stat="identity", width=0.8) +
   # scale_fill_manual(values = CPCOLS) + 
   theme_test() + 
   # scale_x_discrete(labels=labels) +
   scale_x_discrete() +
   xlab("") + 
   ylab("Num of matched genes") + 
   theme(axis.text=element_text(face = "plain", color="#000000")) + 
   # theme(axis.text.x=element_text(angle=70, hjust = 1)) + 
   theme(legend.title=element_blank(),legend.position = "right")+
   theme(legend.title=element_blank(),legend.text = element_text(size = 6),legend.key.size = unit(10, "pt"),legend.position = "right") +
   guides(fill = guide_legend( ncol = 1, byrow = TRUE)) +
   theme(plot.margin = unit(c(2,5,1,2), "cm"))

ggsave(paste0(argv$out,"/",contrast,"_COG_anno.pdf"),p,width = 15, height = 8)
ggsave(paste0(argv$out,"/",contrast,"_COG_anno.png"),p, width = 15, height = 8)
}