# Install and load the necessary packages
library(phangorn)
library(ape)
library(ggtree)
library(gridExtra)
library(treeio)

library(ggplot2)

primates <- read.phyDat(file = "data/total_msa.fasta", format = "fasta")
dist_ml = dist.ml(primates)

tree_upgma_ml  <- upgma(dist_ml)
p = msaplot(p=ggtree(tree_upgma_ml), fasta="data/total_msa.fasta", window=c(1, 28000))
print(p)
