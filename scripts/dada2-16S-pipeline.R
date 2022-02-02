rm(list=ls())

################################
### Dada2 16S pipeline (1.8) ###
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


