library(phangorn)
library(ape)
library(ggtree)
library(gridExtra)

primates <- read.phyDat(file = "data/reduce_msa.fasta", format = "fasta")

dist_ml = dist.ml(primates)
dist_poly = dist.p(primates)
dist_hamming = dist.hamming(primates)

tree_upgma_ml  <- upgma(dist_ml)
tree_upgma_poly  <- upgma(dist_poly)
tree_upgma_hamming  <- upgma(dist_hamming)

tree_nj_ml  <- NJ(dist_ml)
tree_nj_poly  <- NJ(dist_poly)
tree_nj_hamming  <- NJ(dist_hamming)


print(tree_upgma_ml)
print(tree_upgma_poly)
print(tree_upgma_hamming)
#ml_model <- modelTest(primates, model=c("JC", "F81", "K80", "HKY", "SYM", "GTR"), control = pml.control(trace = 0))

# Optimize tree topology, branch lengths, and substitution model
#optimized_tree <- pml(tree = tree_upgma_ml, data = primates, model = ml_model$model)

