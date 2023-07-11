library(phangorn)
library(ape)
library(ggtree)

library(circlize)
library(ComplexHeatmap)

primates <- read.phyDat(file = "data/total_msa.fasta", format = "fasta")
primates_name = names(primates)

dist_ml = as.matrix(dist.ml(primates))
dist_poly = as.matrix(dist.p(primates))
dist_hamming = as.matrix(dist.hamming(primates))

rownames(dist_ml) = primates_name
colnames(dist_ml) = primates_name

rownames(dist_poly) = primates_name
colnames(dist_poly) = primates_name

rownames(dist_hamming) = primates_name
colnames(dist_hamming) = primates_name

figure = 1
if (figure ==1){
  heatmap( dist_ml)
}else if (figure==2){
  heatmap(dist_poly)
}else if (figure==3){
  heatmap(dist_hamming)
}
