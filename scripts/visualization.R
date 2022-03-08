rm(list=ls())

library(tidyverse)
library(readxl)
library(ampvis2)

taxmat <- read.csv("D:/archive/ehil/Analysis_NGS01607_16S/Trimmed_reads/mothur/final.tx.1.cons.taxonomy", sep = "") %>% 
  separate(col = Taxonomy, sep = ";", into = c('Kingdom', 'Phylum', 'Class', 'Order', 'Family', 'Genus')) %>% 
  select(-Size)

otumat <- read.csv("D:/archive/ehil/Analysis_NGS01607_16S/Trimmed_reads/mothur/final.tx.shared", sep = "") %>% 
  select(-c(1, 3))

metmat <- read_xlsx("D:/Dokumente/Projekte/2019-EHIL/data/metadata.xlsx", sheet = 2)

amp_tax <- taxmat %>% 
  column_to_rownames('OTU') %>% 
  mutate(OTU = rownames(.))

amp_met <- metmat %>% 
  select(SampleID, Control_Treatment, Treatment, Location, Year, Field_row)

amp_otu <- setNames(data.frame(t(otumat[ , - 1])), otumat[ , 1]) %>% 
  mutate(OTU = rownames(.))
  

amp <- amp_load(
  otutable = amp_otu,
  taxonomy = amp_tax,
  metadata = amp_met
  )

amp_heatmap(
  amp, 
  facet_by = 'Location', 
  group_by = 'SampleID',
  tax_aggregate = 'Genus', 
  plot_values = F,
  tax_show = 30, 
)

amp_ordinate(
  amp,
  type = 'nmds',
  distmeasure = 'bray', 
  sample_color_by = 'Location', 
  transform = 'total',
  sample_colorframe = 'Treatment',
)
