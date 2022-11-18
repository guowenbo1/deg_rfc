library(limma)
library(dplyr)
library(edgeR)
data_path <- '/Users/guowenbo/Library/CloudStorage/OneDrive-Personal/生物信息/ngsdata/论文/medsci/gse203206/GSE203206_log_Subramaniam.ADRC_brain.counts.tsv'
gse15222_path <- "/Users/guowenbo/Library/CloudStorage/OneDrive-Personal/生物信息/ngsdata/论文/medsci/code/gse15222.csv"
df <- read.table(data_path,header=T)

list <- c(rep("AD", 39), rep("CON",8)) %>% factor(., levels = c("AD", "CON"), ordered = F)
list <- model.matrix(~factor(list)+0)
colnames(list) <- c("AD", "CON")

df.fit <- lmFit(df, list)
df.matrix <- makeContrasts(CON - AD , levels = list)
fit <- contrasts.fit(df.fit, df.matrix)
fit <- eBayes(fit)
tempOutput <- topTable(fit,number=100000)
write.csv(tempOutput,'GSE203206.csv',row.names = FALSE)
