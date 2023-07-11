library(phangorn)
library(ape)
library(ggtree)
library(gridExtra)
primates <- read.phyDat(file = "data/total_msa.fasta", format = "fasta")

dist_ml = dist.ml(primates)
dist_poly = dist.p(primates)
dist_hamming = dist.hamming(primates)

tree_upgma_ml  <- upgma(dist_ml)
tree_upgma_poly  <- upgma(dist_poly)
tree_upgma_hamming  <- upgma(dist_hamming)

tree_nj_ml  <- NJ(dist_ml)
tree_nj_poly  <- NJ(dist_poly)
tree_nj_hamming  <- NJ(dist_hamming)

par(mfrow = c(3, 2))
rooted = FALSE

if (rooted == TRUE){
  fontSize = 0.3
  plot(tree_upgma_ml, main="UPGMA using Ml ",which.plots=2,cex=fontSize)
  plot(tree_nj_ml, main="NJ using ML",which.plots=2,cex=fontSize)
  
  plot(tree_upgma_poly, main="UPGMA using polymorphism",which.plots=2,cex=fontSize)
  plot(tree_nj_poly, main="NJ using polymorphism",which.plots=2,cex=fontSize)
  
  plot(tree_upgma_hamming, main="UPGMA using Hamming",which.plots=2,cex=fontSize)
  plot(tree_nj_hamming, main="NJ using hamming",which.plots=2,cex=fontSize)
  #ggtree_obj <- ggtree(tree_nj_ml)
  #plot(ggtree_obj)
  
}else{
  
  plot(tree_upgma_ml, "unrooted", main="UPGMA using ML",cex = 0.3)
  plot(tree_nj_ml, "unrooted", main="NJ using ML",cex = 0.3)
  
  plot(tree_upgma_poly, "unrooted", main="UPGMA using Polymorphism",cex = 0.3)
  plot(tree_nj_poly, "unrooted", main="NJ using Polymorphism",cex = 0.3)
  
  plot(tree_upgma_hamming, "unrooted", main="UPGMA using Hamming",cex = 0.3)
  plot(tree_nj_hamming, "unrooted", main="NJ using hamming",cex = 0.3)
  
}

#mt <- modelTest(primates, model=c("JC", "F81", "K80", "HKY", "SYM", "GTR"), control = pml.control(trace = 0))

# Optimize tree topology, branch lengths, and substitution model
#optimized_tree <- pml(tree = tree_nj_ml, data = primates, model = ml_model$model)
