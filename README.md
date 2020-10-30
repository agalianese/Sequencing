# Sequencing
This is for different types of sequencing in which I have experience. 

RNA_Seq.R has an example of some RNA Seq data where after differential gene expression analyses was performed on a sample, I found how many genes can be classified as significantly different at an false discovery rate of 0.05. I found if there was a difference in the log2FoldChange.rna between genes that are significantly different and the ones that are not, and found whether gene length impacted RNA siP values. 

Rosalind_Seq.py are from problems found at rosalind.info, a platform for learning bioinformatics through problem solving. The problems below deal with analyzing sequencing data using python, while additional problems dealing with bioinformatic problems can be found at https://github.com/agalianese/Python_Code.

S_Cerevisiae_gff.R is the analyzation of a yeast genome. I analyzed how many genes are on the positive strand versus negative strand, how many genes have introns, what fraction of the genes are on the mitochondria, the longest gene on chromosome II, and which chromosome has the most number of tRNA genes. 

Yeast_Genome_fasta.R is a yeast genome along with its annotations. I created forward and reverse primers that are 30bp long for each CDS in the yeast genome, found which CDS in the yeast genome have long stretches (15 bp) of same nucleotide and which bases are repeated more often, identified which amino acids are most abundant in the C-terminal ends (last 5 residues of an amino acid chain) in yeast proteins, used the genome annotation file to create a histogram of gene lengths of all tRNA_genes, and identified which chromosome has the most arginine amino acids and which one has the highest proportion of arginine amino acids. 
