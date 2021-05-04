# Load the required packages
require(ggplot2)
require(ggseqlogo)

# Read input
clusters <- read.table('clusters.txt', header=TRUE, sep="\t")

# For each cluster, make a list of sequences
sequences = list()
for (i in unique(clusters$cluster)){
  sequences[[i]] <- clusters[clusters$cluster==i,]$CDR3
}
sequences <- sequences[lengths(sequences) != 0]

# Plot sequence logos using ggseqlogo package
ggseqlogo(sequences, ncol=2)