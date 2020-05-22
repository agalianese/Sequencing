## Questions

#libraries
library(GenomicRanges)
library(rtracklayer)


#read in yeast information
yeast.gff <- readGFFAsGRanges("C:/Users/mc/Documents/saccharomyces_cerevisiae.gff")
#turn it into a data frame to work with
yeast.df <- as.data.frame(yeast.gff)
```


#### 1. How many genes are on the positive strand versus negative strand?
```{r}
#use data frame version of the code
yeast.df %>%
  #filter so type is gene
  filter(type == 'gene') %>%
  #group by the type of strand, removing null
  group_by(strand) %>%
  #tally the groups
  tally()

```
#### 2. How many genes have introns?
```{r}

yeast.df %>%
  #filtered for introns
  filter(type == 'intron') %>%
  #took out all its information
  pull(Name) %>%
  #found unique readings
  unique() %>%
  #found how many
  length()

`````
#### 3. What fraction of the genes are on the mitochondria?
```{r}

#number of mitochondrial genes
prob_mito <- yeast.df %>%
  #filtered for genes on the mitochondria
  filter(seqnames == 'chrmt', type == 'gene') %>%
  #counted how many
  tally()

#number of total genes
prob_geno <- yeast.df %>%
  #filtered for just genes
  filter(type == 'gene') %>%
  #tallied
  tally()

#calculates probability
print(prob_mito / prob_geno)

```
#### 4. Find the longest gene on chromosome II
```{r}
yeast.df %>%
  #filter so looking at 2nd chromosome
  filter(seqnames == 'chrII') %>%
  #arrange in descending order of width
  arrange(desc(width)) %>%
  #longest gene
  .[2,]

````
#### 5. Which chromosome has the most number of tRNA genes? 
```{r}
yeast.df %>%
  #filter for only tRNA gene
  filter(type == 'tRNA_gene') %>%
  #group by chromosomes
  group_by(seqnames) %>%
  #tally the tRNA genes for each
  tally() %>%
  #arrange in descending order
  arrange(desc(n)) %>%
  #print the chromosome
  .[1,]

````
