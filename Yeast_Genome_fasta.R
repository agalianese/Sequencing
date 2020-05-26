```{r}
library(Biostrings)
library(rtracklayer)
library(tidyverse)
library(ggplot2)
library(ggpubr)
library(ggseqlogo)


###Given input code ####
yeast.fasta <- readDNAStringSet("yeast.fasta")
yeast.gff <- readGFFAsGRanges("saccharomyces_cerevisiae.gff.txt")

yeast.cds.annot <- yeast.gff[yeast.gff$type == "CDS"]
yeast.cds.all.seq <- yeast.fasta[yeast.cds.annot]
names(yeast.cds.all.seq) <- yeast.cds.annot$Name
cds.on.neg.strand <- strand(yeast.cds.annot) == "-"

yeast.cds.all.seq[cds.on.neg.strand] <- reverseComplement(yeast.cds.all.seq[cds.on.neg.strand])


cds.names <- unique(names(yeast.cds.all.seq))
list.of.combd.seq <- list()

for(dummy in cds.names){
  cds.dna <- yeast.cds.all.seq[names(yeast.cds.all.seq) == dummy]
  collapsed.cds <- unlist(cds.dna)
  
  list.of.combd.seq[[dummy]] <- collapsed.cds
}

all.yeast.cds.seq <- DNAStringSet(list.of.combd.seq)


yeast.cds <- yeast.gff[yeast.gff$type == "CDS"]

# Get a list of genes without any introns
genes.wo.introns <- yeast.cds %>%
  as.data.frame() %>%
  group_by(Name) %>%
  tally() %>%
  filter(n == 1) %>%
  pull(Name)

yeast.cds.no.introns <- yeast.cds[yeast.cds$Name %in% genes.wo.introns]

# Get the sequences of aCDS
yeast.cds.seq <- yeast.fasta[yeast.cds.no.introns]

# Add gene names to the new DNAStringSet
names(yeast.cds.seq) <- yeast.cds.no.introns$Name

# Make sure all sequences are in the right orientation

yeast.cds.seq[strand(yeast.cds.no.introns) == "-"] <- reverseComplement(yeast.cds.seq[strand(yeast.cds.no.introns) == "-"])

```

# Create forward and reverse primers that are 30bp long for each CDS in the yeast genome
```{r}
#get the first 30bp in the all.yeasts
forward.yeast.30 <- subseq(all.yeast.cds.seq, start = 1, end = 30)
#get the last 30 bp in all.yeasts
reverse.yeast.30 <- subseq(all.yeast.cds.seq, start = width(all.yeast.cds.seq) - 29, end = width(all.yeast.cds.seq))
#take the reverse complement
reverse.yeast.30 <- reverseComplement(reverse.yeast.30)

#print output
print(forward.yeast.30)
print(reverse.yeast.30)

```


#Which CDS in the yeast genome have long stretches (15 bp) of same nucleotide? 
#Which bases are repeated more often?
```{r}
#create the patterns to be matched
bp.a = "AAAAAAAAAAAAAAA"
bp.t = 'TTTTTTTTTTTTTTT'
bp.g = 'GGGGGGGGGGGGGGG'
bp.c = 'CCCCCCCCCCCCCCC'

#use vcountPattern to find the amount of times the pattern occurs in the cds region
#of the yeast genome, then sum and print them
a.count = sum(vcountPattern(pattern = bp.a, all.yeast.cds.seq, max.mismatch = 0))
print(a.count)
t.count = sum(vcountPattern(pattern = bp.t, all.yeast.cds.seq, max.mismatch = 0))
print(t.count)
g.count = sum(vcountPattern(pattern = bp.g, all.yeast.cds.seq, max.mismatch = 0))
print(g.count)
c.count = sum(vcountPattern(pattern = bp.c, all.yeast.cds.seq, max.mismatch = 0))
print(c.count)

#add the number of long reads together
long_read.count = a.count + t.count + g.count + c.count
print(long_read.count)
``



### Which amino acids are most abundant in the C-terminal ends (last 5 residues) 
#of yeast proteins?
```{r}
#translate into amino acids
aa.yeast <- translate(all.yeast.cds.seq)

#take the last 5 residues
aa.yeast.5 <- subseq(aa.yeast, start = width(aa.yeast) - 5, end = width(aa.yeast))

#find the frequencies and output into a chart
alphabetFrequency(aa.yeast.5, collapse = T) 

#exluding *, K is the most abundant with 3403 instances. 



#Using the genome annotation file, create a histogram of gene lengths of all `tRNA_gene`
````{r}
gff.df <- as.data.frame(yeast.gff)
table(gff.df$type)
hist_trna <- gff.df %>%
  filter(type == 'tRNA_gene')

hist(hist_trna$width)

````


#Which chromosome has the most arginine amino acids and which one has the
#highest proportion of arginine amino acids?
```{r}

#translate the yeast fasta into its amino acid code
aa__cds_yeast = translate(yeast.fasta)


#count the instances of Arginine (R)
arg_count_n = letterFrequency(aa_yeast, letters = 'R')


#count the instances of Arginine (R) and report it as a probability
arg_count_prob = letterFrequency(aa_yeast, letters = 'R', as.prob = T)

#create a data frame of each chromosome and its corresponding arginine amount
arg_df <- data.frame('Chr' = c(1:17),
                     "Arginine Count" = arg_count_n,
                     "Arginine Prob" = arg_count_prob)

#arrange in descending order and report which chromosome has the highest number of arginine
arg_chr_count <- arg_df %>%
  arrange(desc(R)) %>%
  .[1,]

#arrange in descending order and report which chromosome has the highest proportion of arginine
arg_chr_prop <- arg_df %>%
  arrange(desc(R.1)) %>%
  .[1,]
```




