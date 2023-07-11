library(phangorn)
library(ape)
library(ggtree)
library(gridExtra)

primates <- read.phyDat(file = "data/reduce_msa.fasta", format = "fasta")

dist_ml = dist.ml(primates)
bn
tree_upgma_ml  <- upgma(dist_ml)
tree_nj_ml  <- NJ(dist_ml)

optimizer_tree_upgma <- pml(tree = tree_upgma_ml, data = primates)
optimizer_tree_nj <- pml(tree = tree_nj_ml, data = primates)

fit_tree_upgma <- optim.pml(optimizer_tree_upgma, data = primates, model = "GTR")
fit_tree_nj <- optim.pml(optimizer_tree_nj, data = primates, model = "GTR")
par(mfrow = c(2, 1))
bs <- bootstrap.pml(fit_tree_upgma, bs=100, optNni=TRUE,
                    control = pml.control(trace = 0))
plotBS(midpoint(fit_tree_upgma$tree), bs, p = 50, type="p", main="Standard bootstrap For UPGMA",cex=0.4)

bs <- bootstrap.pml(fit_tree_upgma, bs=100, optNni=TRUE,
                    control = pml.control(trace = 0))
plotBS(midpoint(fit_tree_nj$tree), bs, p = 50, type="p", main="Standard bootstrap For NJ",cex=0.4)

for (i in 1:length(optimized_tree$edge)) {
  branch <- optimized_tree$edge[i, ]
  support_count <- sum(sapply(bootstrap_trees, function(tree) edge.exists(tree, branch)))
  support_percentage <- (support_count / num_replicates) * 100
  branch_support <- c(branch_support, support_percentage)
}


#ml_model <- modelTest(primates, model=c("JC", "F81", "K80", "HKY", "SYM", "GTR"), control = pml.control(trace = 0))

# Optimize tree topology, branch lengths, and substitution model
#optimized_tree <- pml(tree = tree_upgma_ml, data = primates, model = ml_model$model)

