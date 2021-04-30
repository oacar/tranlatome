library(tidyverse)

# Read Input files -----------------
# Read costanzo 2016 dataset downloaded on Mon, Nov  5 2018 from http://thecellmap.org/costanzo2016/
# SGA_NxN for nonessential-nonessential network
# SGA_ExN_NxE for essential-nonessential & nonessential-essential network
# SGA_ExE for essential-essential network

filename_SGA_NxN <- "analysis/data/raw_data/SGA_NxN.txt"
filename_SGA_ExN_NxE <- "analysis/data/raw_data/SGA_ExN_NxE.txt"
filename_SGA_ExE <- "analysis/data/raw_data/SGA_ExE.txt"

SGA_NxN <- read_delim(filename_SGA_NxN, delim = "	") %>% mutate(data_source = 'NxN')
SGA_ExN_NxE <- read_delim(filename_SGA_ExN_NxE, delim = "	") %>% mutate(data_source = 'ExN_NxE')
SGA_ExE <- read_delim(filename_SGA_ExE, delim = "	") %>% mutate(data_source = 'ExE')


# Mutate systematic names---------
# Original input files has extensions on SGD systematic names (like 'YAL002_sn273')
# This makes finding genes in the future analysis difficult, so I remove those parts after '_' character
# for both query and array id columns

SGA_NxN$`Array Strain ID` <- sapply(strsplit(SGA_NxN$`Array Strain ID`, "_"), "[", 1)
SGA_NxN$`Query Strain ID` <- sapply(strsplit(SGA_NxN$`Query Strain ID`, "_"), "[", 1)
SGA_ExN_NxE$`Array Strain ID` <- sapply(strsplit(SGA_ExN_NxE$`Array Strain ID`, "_"), "[", 1)
SGA_ExN_NxE$`Query Strain ID` <- sapply(strsplit(SGA_ExN_NxE$`Query Strain ID`, "_"), "[", 1)
SGA_ExE$`Array Strain ID` <- sapply(strsplit(SGA_ExE$`Array Strain ID`, "_"), "[", 1)
SGA_ExE$`Query Strain ID` <- sapply(strsplit(SGA_ExE$`Query Strain ID`, "_"), "[", 1)

# Get the data frame from edited input data
net.df <- bind_rows(SGA_NxN, SGA_ExN_NxE, SGA_ExE)
#following alleles has 2 temperature reported while supplementary methods of Costanzo 2016 says they used particular temperatures and reported only 1 according to quality. Thus I stick with their explanation and remove 30 degrees data
supplist<- c('vac8-supp1','mtg1-supp1','tfb1-6-supp1','mob2-11-supp1')
net.df.filtered=net.df %>% filter(!(`Query allele name`=='med6-ts'&`Arraytype/Temp`=='DMA30')&!(`Query allele name`%in%supplist&`Arraytype/Temp`=='DMA30')&!(`Query allele name`%in%supplist&`Arraytype/Temp`=='TSA30'))
#save data
saveRDS(net.df.filtered,'analysis/data/derived_data/SGA_data_combined.rds.gz',compress = 'gzip')



