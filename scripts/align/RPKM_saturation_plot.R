library(tidyverse)


sample <- snakemake@wildcards[["sample"]]

data <- read_tsv(snakemake@input[["rawCount"]]) %>% .[, 7:26]

myfunc <- function(x){
  for(i in seq_along(x)){
    if (x[i]>0) {
      x[i] <- 1
    }
  }
  return(x)
}


data <- map_df(data,myfunc) %>% colSums(.) %>% as.data.frame(.)

colnames(data) <- "genenum"

data <- data %>% rownames_to_column("seq") %>% mutate(seq2 = str_replace(seq, "%", ""))

data$seq2 <- as.numeric(data$seq2)

p <- ggplot(data, aes(seq2, genenum)) +  geom_point() + geom_line(aes(group = 1)) +
  theme_bw() + xlab("Mapped reads(%)") + ylab("Detected genes Number") + ggtitle(paste(sample, "saturation curve")) +
  theme(panel.grid = element_blank(), legend.title = element_blank())

ggsave(snakemake@output[["fig"]][1], p, width = 8, height = 5)
ggsave(snakemake@output[["fig"]][2], p, width = 8, height = 5)