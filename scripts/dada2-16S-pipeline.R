rm(list=ls())

################################
### Dada2 16S pipeline (1.8) ###
################################

# Link: https://benjjneb.github.io/dada2/tutorial_1_8.html

# load packages
library(dada2); packageVersion('dada2')
library(ShortRead); packageVersion('ShortRead')
library(Biostrings); packageVersion('Biostrings')
library(tidyverse); packageVersion('tidyverse')

# gather data
path <- 'D:/ehil/Analysis_NGS01607_16S/Original_reads/'
path <- 'D:/ehil/Analysis_NGS01607_16S/Trimmed_reads/'
list.files(path)

fnFs <- sort(list.files(path, pattern = 'R1.fastq.gz', full.names = TRUE))
fnRs <- sort(list.files(path, pattern = 'R2.fastq.gz', full.names = TRUE))


fnFs <- sort(list.files(path, pattern = 'R1_trimmed.fastq.gz', full.names = TRUE))
fnRs <- sort(list.files(path, pattern = 'R2_trimmed.fastq.gz', full.names = TRUE))

sample.names <- sapply(strsplit(basename(fnFs), "_"), `[`, 1)

plotQualityProfile(fnFs[1:2])

plotQualityProfile(fnRs[1:2])

# Place filtered files in filtered/ subdirectory
filtFs <- file.path(path, "filtered", paste0(sample.names, "_F_filt.fastq.gz"))
filtRs <- file.path(path, "filtered", paste0(sample.names, "_R_filt.fastq.gz"))

out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs, truncLen=c(230,220),
                     maxN=0, maxEE=c(2,2), truncQ=2, rm.phix=TRUE,
                     compress=TRUE, multithread=FALSE) # On Windows set multithread=FALSE
head(out)

errF <- learnErrors(filtFs, multithread=TRUE)
errR <- learnErrors(filtRs, multithread=TRUE)

plotErrors(errF, nominalQ=TRUE)
plotErrors(errR, nominalQ=TRUE)

derepFs <- derepFastq(filtFs, verbose=TRUE)
derepRs <- derepFastq(filtRs, verbose=TRUE)

# Name the derep-class objects by the sample names
names(derepFs) <- sample.names
names(derepRs) <- sample.names
