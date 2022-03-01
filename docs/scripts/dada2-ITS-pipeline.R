rm(list=ls())

################################
### Dada2 ITS pipeline (1.8) ###
################################

# Link: https://benjjneb.github.io/dada2/ITS_workflow.html


# load packages
library(dada2); packageVersion('dada2')
library(ShortRead); packageVersion('ShortRead')
library(Biostrings); packageVersion('Biostrings')
library(tidyverse); packageVersion('tidyverse')

# gather data
path <- 'D:/ehil/Analysis_NGS01607_ITS/Trimmed_reads'
list.files(path)

fnFs <- sort(list.files(path, pattern = 'R1_trimmed.fastq.gz', full.names = TRUE))
fnRs <- sort(list.files(path, pattern = 'R2_trimmed.fastq.gz', full.names = TRUE))


# identify primers
FWD <- 'CCTACGGGNGGCWGCAG'
REV <- 'GACTACHVGGGTATCTAATCC'

allOrients <- function(primer) {
  # Create all orientations of the input sequence
  require(Biostrings)
  dna <- DNAString(primer)  # The Biostrings works w/ DNAString objects rather than character vectors
  orients <- c(Forward = dna, Complement = complement(dna), Reverse = reverse(dna), 
               RevComp = reverseComplement(dna))
  return(sapply(orients, toString))  # Convert back to character vector
}

FWD.orients <- allOrients(FWD)
REV.orients <- allOrients(REV)
FWD.orients

fnFs.filtN <- file.path(path, "filtN", basename(fnFs)) # Put N-filterd files in filtN/ subdirectory
fnRs.filtN <- file.path(path, "filtN", basename(fnRs))
filterAndTrim(fnFs, fnFs.filtN, fnRs, fnRs.filtN, maxN = 0, multithread = FALSE)

primerHits <- function(primer, fn) {
  # Counts number of reads in which the primer is found
  nhits <- vcountPattern(primer, sread(readFastq(fn)), fixed = FALSE)
  return(sum(nhits > 0))
}
rbind(FWD.ForwardReads = sapply(FWD.orients, primerHits, fn = fnFs.filtN[[1]]), 
      FWD.ReverseReads = sapply(FWD.orients, primerHits, fn = fnRs.filtN[[1]]), 
      REV.ForwardReads = sapply(REV.orients, primerHits, fn = fnFs.filtN[[1]]), 
      REV.ReverseReads = sapply(REV.orients, primerHits, fn = fnRs.filtN[[1]]))