library(tidyverse)
library(GEOquery)
library(ggpubr)


data3 <- read_csv("C:/Users/mc/Documents/data3.csv")


#You've performed differential gene expression analyses on a sample. The data is in data3.csv. 
#How many genes can be classified as significantly different at an false discovery rate of 0.05? (2 pts)

#add a new column of the pvalues adjusted by a false discovery rate
data3 <- data3 %>%
  mutate(adjust.p = p.adjust(pvalue,method = 'fdr'))

#filter for significance based on the adjusted p-values
sig <- data3 %>%
  filter(adjust.p < 0.05)
not.sig <- data3 %>%
  filter(adjust.p >= .05)


#3628 genes are significantly different



#Is there a difference in the log2FoldChange.rna between genes that are significantly different and 
#the ones that are not? 


#perform wilcox as n > 30 and shape isn't normal to compare statistical means
hist(sig$log2FoldChange.rna)
hist(not.sig$log2FoldChange.rna)
wilcox.test(sig$log2FoldChange.rna, not.sig$log2FoldChange.rna)



#Does gene length impact RNA siP?

#filter based on gene lengths
long_genes <- data3 %>%
  filter(Length > 5000)

short_genes <- data3 %>%
  filter(Length < 1000)

#create a histogram to visualize the data
hist(short_genes$RNA.siP.mean.tpm, xlim = c(0,1100), breaks = 1000)
hist(long_genes$RNA.siP.mean.tpm, xlim = c(0,1100), breaks = 1000)

#non normal so conduct wilcox test
wilcox.test(short_genes$RNA.siP.mean.tpm, long_genes$RNA.siP.mean.tpm)

#p-value of 2.2e-16 indicates relationship

