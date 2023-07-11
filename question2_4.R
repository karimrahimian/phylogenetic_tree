library(phangorn)
library(ape)
library(ggtree)
library(gridExtra)

primates <- read.phyDat(file = "data/reduce_msa.fasta", format = "fasta")

dist_ml = dist.ml(primates)

tree_upgma_ml  <- upgma(dist_ml)
tree_nj_ml  <- NJ(dist_ml)

optimizer_tree_upgma <- pml(tree = tree_upgma_ml, data = primates)
optimizer_tree_nj <- pml(tree = tree_nj_ml, data = primates)

optimized_tree_upgma <- optim.pml(optimizer_tree_upgma, data = primates, model = "GTR")
optimized_tree_nj <- optim.pml(optimizer_tree_nj, data = primates, model = "GTR")


#ml_model <- modelTest(primates, model=c("JC", "F81", "K80", "HKY", "SYM", "GTR"), control = pml.control(trace = 0))

# Optimize tree topology, branch lengths, and substitution model
#optimized_tree <- pml(tree = tree_upgma_ml, data = primates, model = ml_model$model)

