library(phangorn)
library(ape)
library(ggtree)
library(gridExtra)

primates <- read.phyDat(file = "data/total_msa.fasta", format = "fasta")

mt <- modelTest(fasta, model=c("JC", "F81", "K80", "HKY", "SYM", "GTR"), 
                control = pml.control(trace = 0))

write.csv(mt, file = "result.csv", row.names = FALSE)
print(mt)