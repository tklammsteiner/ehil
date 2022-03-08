rm(list = ls())

library(tidyverse)
library(reshape2)
library(readxl)

setwd("O:PhD/Projekte/2019-EHIL/")
list.files()

metadata <- read_xlsx("data/metadata.xlsx", sheet = 2)
