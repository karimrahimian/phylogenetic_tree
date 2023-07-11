#install.packages("Biostrings")
library(Biostrings)

findDuplicateTaxonName<-function (fasta_data){
  print ("Duplicate Taxon name  : ")
  taxon_names <- names(fasta_data)
  duplicate_taxon_names <- taxon_names[duplicated(taxon_names)]
  print(duplicate_taxon_names)
}

findDuplicateSequence<- function (fasta_data){
  print ("Duplicate Sequence are  : ")
  taxon_names <- names(fasta_data)
  duplicate_taxon_names <- taxon_names[duplicated(taxon_names)]
  duplicate_sequences <- fasta_data[duplicated(fasta_data)]
  print(duplicate_sequences)
}

findInvalidCharacterInSequence<- function (fasta_data){
  print ("Invalid taxon sequences : ")
  taxon_names <- names(fasta_data)
  sequences <- as.character(fasta_data)
  valid_pattern <- "^[ACGT-]+$"
  invalid_sequences <- sequences[!grepl(valid_pattern, sequences)]
  print(invalid_sequences)
}

findInvalidCharacterInTaxonName<- function (fasta_data){
  print ("Invalid taxon taxon : ")
  header_names <- names(fasta_data)

  taxon_names <- sub("/.*", "", header_names)
  valid_pattern <- "^[A-Za-z0-9_\\s]+$"
  invalid_taxon_names <- vector()
  
  for (taxon_name in taxon_names) {
    if (!grepl(valid_pattern, taxon_name)) {
      invalid_taxon_names <- c(invalid_taxon_names, taxon_name)
    }
  }
  
  print(invalid_taxon_names)
}

selectSubsetOfData<-function(fasta_data){
  
  fasta_data <- unique(fasta_data)
  taxon_names <- names(fasta_data)
  
  new_taxon_names <- sapply(strsplit(taxon_names, "|",fixed=TRUE), function(x) x[1])
  new_taxon_names <- sapply(strsplit(new_taxon_names, "/",fixed=TRUE), function(x) paste0(x[3]))
  
  subset_indices <- c()
  size = length(fasta_data)
  for (i in 1:165){
    subset_indices = c(subset_indices,i)
  }
  
  subset_sequences <- fasta_data[subset_indices]
  subset_taxon_names <- new_taxon_names[subset_indices]
  
  sequence_names <- names(subset_sequences)
  sequences <- as.character(subset_sequences)
  names(sequences) <- subset_taxon_names
  dna_set <- DNAStringSet(sequences)

  writeXStringSet(dna_set, file = "data/total_msa.fasta")
  #write.fasta(sequences, names = sequence_names ,file.out = "data/reduce_msa.fasta")
  
}

#filename = "data/gisaid_hcov-19_2023_06_30_09.fasta"
filename = "data/total_wuhan_msa.fasta"
fasta_data <- readDNAStringSet(filename)

findInvalidCharacterInSequence(fasta_data)
findInvalidCharacterInTaxonName(fasta_data)
findDuplicateTaxonName(fasta_data)
findDuplicateSequence(fasta_data)
selectSubsetOfData(fasta_data)

