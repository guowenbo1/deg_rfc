library(clusterProfiler) # 富集分析R包
library(stringr) # 标签换行
library(AnnotationDbi)
library(org.Hs.eg.db) # 参考物种的基因注释数据库hg19
library(DOSE)
library(ggplot2) # 绘图
library(ggrepel) # 标签相关
# 将gene symbol转换为Entrez ID,防止分析出错
deg_path <- "/Users/guowenbo/Library/CloudStorage/OneDrive-Personal/生物信息/ngsdata/论文/medsci/table/deg_genename.csv"
deg_gene <- read.csv(deg_path)
id_list <- mapIds(org.Hs.eg.db,deg_gene$name,"ENTREZID","SYMBOL")
# 去除转换不成功的结果，即id=NA的情况
id_list <- na.omit(id_list)
# GO富集分析
go <- enrichGO(gene = id_list, # Entrez ID列表
               OrgDb = org.Hs.eg.db, # 指定物种数据库
               keyType = "ENTREZID", # 指定给定的名称类型
               ont = "ALL", # 可选,BP(生物学过程)/CC(细胞组分)/MF(分子功能)/ALL(同时指定)
               pAdjustMethod = "BH", # P值校正方法,还可以是fdr
               pvalueCutoff = 0.05,qvalueCutoff = 0.2, # p/q值阈值
               readable = T # 将ID转换为symbol
)
go.res <- data.frame(go) # 将GO结果转为数据框，方便后续分析（不转也可以，看个人习惯）
write.csv(go.res,"Table_GO_result.csv",quote = F) # 输出GO富集分析结果
# 绘制GO富集分析条形图，结果默认按qvalue升序，分别选出前十的term进行绘图即可
goBP <- subset(go.res,subset = (ONTOLOGY == "BP"))[1:10,]
goCC <- subset(go.res,subset = (ONTOLOGY == "CC"))[1:10,]
goMF <- subset(go.res,subset = (ONTOLOGY == "MF"))[1:10,]
go.df <- rbind(goBP,goCC,goMF)
# 使画出的GO term的顺序与输入一致
go.df$Description <- factor(go.df$Description,levels = rev(go.df$Description))
# 绘图
go_bar <- ggplot(data = go.df, # 绘图使用的数据
                 aes(x = Description, y = Count,fill = ONTOLOGY))+ # 横轴坐标及颜色分类填充
  geom_bar(stat = "identity",width = 0.9)+ # 绘制条形图及宽度设置
  coord_flip()+theme_bw()+ # 横纵坐标反转及去除背景色
  scale_x_discrete(labels = function(x) str_wrap(x,width = 50))+ # 设置term名称过长时换行
  labs(x = "GO terms",y = "GeneNumber",title = "Barplot of Enriched GO Terms")+ # 设置坐标轴标题及标题
  theme(axis.title = element_text(size = 13), # 坐标轴标题大小
        axis.text = element_text(size = 11), # 坐标轴标签大小
        plot.title = element_text(size = 14,hjust = 0.5,face = "bold"), # 标题设置
        legend.title = element_text(size = 13), # 图例标题大小
        legend.text = element_text(size = 11), # 图例标签大小
        plot.margin = unit(c(0.5,0.5,0.5,0.5),"cm")) # 图边距
ggsave(go_bar,filename = "GO_Barplot.pdf",width = 9,height = 7)


# KEGG富集分析
kegg <- enrichKEGG(gene = id_list, 
                   organism = "hsa",keyType = "kegg", 
                   pAdjustMethod = "BH",pvalueCutoff = 0.05,qvalueCutoff = 0.2,
                   minGSSize = 10,maxGSSize = 500,use_internal_data = F)
kegg_enrich <- kegg
res <- as.data.frame(kegg_enrich)
write.csv(res, 'kegg_enrich.csv')

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
  labs(title= 'KEGG Patthway') +
  theme(panel.grid = element_blank()) + 
  theme(axis.title.x =element_text(size=14), axis.title.y=element_text(size=14),legend.title=element_text(size=13),axis.text=element_text(size=12,color  = "black")) + 
  scale_fill_gradientn(colours=color) +
  geom_text(aes(label = shape, vjust = 0.5, hjust = -0.2))
ggsave('kegg_pathway_barplot.png', width = my_width, height = 6.5)

library("pathview")
pathview(gene.data  = deg_gene,
         pathway.id = "hsa04080",
         species    = "hsa",
         limit      = list(gene=max(abs(deg_gene)), cpd=1))


