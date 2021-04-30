library(tidyverse)
# this script reads SGA_data_combined data frame which contains costanzo 2016 data
# and creates a data frame containing SGD systematic gene names and the Allele names used in the experiments


net.df <- readRDS('analysis/data/derived_data/SGA_data_combined.rds.gz')

# Take query data and array data as separate data frames
# give meaningful column names to both
# combine two data frames, remove suppressor mutations, extract unique rows and save
q.data <- net.df[, c(1, 2)] %>% distinct()
a.data <- net.df[, c(3, 4)] %>% distinct()
colnames(q.data) <- c("Systematic gene name", "Allele Gene name")
colnames(a.data) <- c("Systematic gene name", "Allele Gene name")

strain_ids <- bind_rows(q.data, a.data) %>%
  distinct() %>%
  filter(grepl("supp", `Allele Gene name`) == F)


essetial.query <- net.df %>%
  filter(data_source=='ExE') %>%
  select(`Query allele name`) %>%
  pull()
essetial.array <- net.df %>%
  filter(data_source=='ExE') %>%
  select(`Array allele name`) %>%
  pull()
exn.query <- net.df %>%
  filter(data_source == "ExN_NxE" & (`Arraytype/Temp` == "TSA26" | `Arraytype/Temp` == "TSA30")==F) %>%
  select(`Query allele name`) %>%
  distinct() %>%
  pull()
essential.alleles <- unique(c(essetial.query,essetial.array,exn.query))

exp.number.data <- strain_ids %>%
  mutate(maincat = ifelse(`Allele Gene name`%in%essential.alleles,'essential','nonessential'))

write_csv(exp.number.data, "analysis/data/derived_data/strain_ids_with_experiment_count_all.csv")


