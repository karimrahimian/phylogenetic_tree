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
tree_nj_ml  <- NJ(dist_ml)

plot1 =  ggtree(tree_upgma_ml, size=0.1,layout="circular", branch.length='none')+geom_tiplab(aes(angle=angle), color='blue',size = 3)+geom_point(aes(shape=isTip, color=isTip), size=1)+ggtitle("label")
plot2 =  ggtree(tree_nj_ml, size=0.1,layout="circular", branch.length='none')+geom_tiplab(aes(angle=angle), color='blue',size = 3)+geom_point(aes(shape=isTip, color=isTip), size=1) + geom_text(aes(label=node), hjust=-.3,size=2)

#grid.arrange(plot1, plot2, plot1, plot2, ncol=2,nrow=2)
grid.arrange(plot2, ncol=1, nrow=1)




