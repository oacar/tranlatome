library(tidyverse)
library(readxl)
library(ggsignif)
library(glue)
# Read Input files -----------------
con <- DBI::dbConnect(RMariaDB::MariaDB(), groups = "mariaDB")
translation <- DBI::dbGetQuery(con, "SELECT * FROM aaron.scer_orfs_translation_info")
DBI::dbDisconnect(con)

pgs <- data.table::fread("/home/acwach/YeastTest/named_transient_orfs")
colnames(pgs) <- c("V1", "orf_name")
overlappingorfs <- read_csv("analysis/data/raw_data/overlappingorfs.csv")
# Read costanzo 2016 dataset downloaded on Mon, Nov  5 2018 from http://thecellmap.org/costanzo2016/
# SGA_NxN for nonessential-nonessential network
# SGA_ExN_NxE for essential-nonessential & nonessential-essential network
# SGA_ExE for essential-essential network

filename_SGA_NxN <- "analysis/data/raw_data/SGA_NxN.txt"
filename_SGA_ExN_NxE <- "analysis/data/raw_data/SGA_ExN_NxE.txt"
filename_SGA_ExE <- "analysis/data/raw_data/SGA_ExE.txt"

SGA_NxN <- read_delim(filename_SGA_NxN, delim = "	") %>% mutate(data_source = "NxN")
SGA_ExN_NxE <- read_delim(filename_SGA_ExN_NxE, delim = "	") %>% mutate(data_source = "ExN_NxE")
SGA_ExE <- read_delim(filename_SGA_ExE, delim = "	") %>% mutate(data_source = "ExE")


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


# Create a dataframe with systematic gene names and the Allele names used in the experiments and whether gene is essential or nonessential
q.data <- net.df[, c(1, 2)] %>% distinct()
a.data <- net.df[, c(3, 4)] %>% distinct()
colnames(q.data) <- c("Systematic gene name", "Allele Gene name")
colnames(a.data) <- c("Systematic gene name", "Allele Gene name")

# filter out any suppressor mutations
strain_ids <- bind_rows(q.data, a.data) %>%
    distinct() %>%
    filter(grepl("supp", `Allele Gene name`) == F)


# find which orfs are essential and nonessential
essetial.query <- net.df %>%
    filter(data_source == "ExE") %>%
    select(`Query allele name`) %>%
    pull()
essetial.array <- net.df %>%
    filter(data_source == "ExE") %>%
    select(`Array allele name`) %>%
    pull()
exn.query <- net.df %>%
    filter(data_source == "ExN_NxE" & (`Arraytype/Temp` == "TSA26" | `Arraytype/Temp` == "TSA30") == F) %>%
    select(`Query allele name`) %>%
    distinct() %>%
    pull()
essential.alleles <- unique(c(essetial.query, essetial.array, exn.query))

# assign essential and nonessential to each gene
exp.number.data <- strain_ids %>%
    mutate(maincat = ifelse(`Allele Gene name` %in% essential.alleles, "essential", "nonessential"))

# Calculate genetic interaction ratio
# pgs<-pgs[pgs$orf_name %in% overlappingorfs$orf_name==FALSE,]
# net.df <- readRDS("analysis/data/derived_data/SGA_data_combined.rds.gz")
net_df_significant <- filter(net.df, `P-value` <= 0.05 & `Query Strain ID` %in% overlappingorfs$orf_name == FALSE & `Array Strain ID` %in% overlappingorfs$orf_name == FALSE)
# net_df_significant_sl <- filter(net_df_significant, `Genetic interaction score (ε)` <= -0.2)#& `Double mutant fitness standard deviation`<=0.1)
exp.data <- exp.number.data %>% mutate(group = ifelse(`Systematic gene name` %in% pgs$orf_name, "proto-gene", maincat))

# get YER175W-A genetic interactions at epsilon < -0.2
net_df_significant %>%
    filter(`P-value` <= 0.05, `Genetic interaction score (ε)` <= -0.2) %>%
    filter(`Query Strain ID` == "YER175W-A" | `Array Strain ID` == "YER175W-A")

## Calculate strong (eps<-0.2) and lethal (eps<-0.35) interactions -----
strong_orf_list <- net_df_significant %>%
    filter(`P-value` <= 0.05, `Genetic interaction score (ε)` <= -0.2) %>%
    select(`Query Strain ID`, `Array Strain ID`) %>%
    as.list() %>%
    unlist() %>%
    unique()
lethal_orf_list <- net_df_significant %>%
    filter(`P-value` <= 0.05, `Genetic interaction score (ε)` <= -0.35) %>%
    select(`Query Strain ID`, `Array Strain ID`) %>%
    as.list() %>%
    unlist() %>%
    unique()

exp.data <- exp.data %>%
    mutate(
        strong_interaction = ifelse(`Systematic gene name` %in% strong_orf_list, TRUE, FALSE),
        lethal_interaction = ifelse(`Systematic gene name` %in% lethal_orf_list, TRUE, FALSE)
    ) %>%
    filter(`Systematic gene name` %in% overlappingorfs$orf_name == F)

exp_data_nones <- exp.data %>% filter(group != "essential")
# check numbers
translation %>%
    filter(gene_systematic_name %in% pgs$orf_name[pgs$orf_name %in% exp_data_nones$`Systematic gene name`]) %>%
    select(orf_class) %>%
    group_by(orf_class) %>%
    count()

# 3901/(3901+772)
table(exp_data_nones$group, exp_data_nones$strong_interaction)
prop.table(table(exp_data_nones$group, exp_data_nones$strong_interaction), 1)
table(exp_data_nones$group, exp_data_nones$lethal_interaction)
prop.table(table(exp_data_nones$group, exp_data_nones$lethal_interaction), 1)
## Compare interaction ratios of transient orfs and nonessential orfs with fisher.test ----
strong_interaction_pval <- table(exp_data_nones$group, exp_data_nones$strong_interaction) %>% fisher.test() # 0.04652
lethal_interaction_pval <- table(exp_data_nones$group, exp_data_nones$lethal_interaction) %>% fisher.test() # 0.02636

percent_plot_data <- exp_data_nones %>%
    select(group, strong_interaction) %>%
    group_by(group, strong_interaction) %>%
    summarise(str_count = n()) %>%
    mutate(freq = str_count / sum(str_count)) %>%
    filter(strong_interaction == TRUE) %>%
    mutate(cat = "strong") %>% # %>% gather()
    bind_rows(
        exp_data_nones %>% select(group, lethal_interaction) %>% group_by(group, lethal_interaction) %>% summarise(leth_count = n()) %>% mutate(freq = leth_count / sum(leth_count)) %>% filter(lethal_interaction == TRUE) %>% mutate(cat = "lethal")
    ) %>%
    mutate(cat = factor(cat, levels = c("strong", "lethal"))) %>%
    mutate(group = factor(group, levels = c("proto-gene", "nonessential")))


percent_plot <- percent_plot_data %>%
    ggplot(aes(x = cat, y = freq, fill = group)) +
    geom_bar(stat = "identity", position = position_dodge(width = 0.95)) +
    theme_bw() +
    theme(
        axis.title.x = element_text(size = 8), axis.text.x = element_text(size = 10), axis.text.y = element_text(size = 8), axis.title.y = element_text(size = 12),
        legend.text = element_text(size = 8), legend.position = "bottom", legend.title = element_blank(), legend.margin = margin(t = -0.25, unit = "in"),
        legend.spacing.x = unit(0.1, "in")
    ) +
    scale_fill_manual(labels = c("Transient ORFs", "Nonessential\nGenes"), values = c("#673ab7", "#9e9e9e")) +
    scale_x_discrete(labels = c(paste0("\u03B5", "< -0.2"), paste0("\u03B5", "<-0.35"))) +
    ylab("Percent with at least one\ninteraction at given threshold") +
    scale_y_continuous(labels = scales::percent) +
    xlab("") +
    ggsignif::geom_signif(
        y_position = c(1.05, 0.85), xmin = c(0.8, 1.8), xmax = c(1.2, 2.2),
        annotation = c(glue::glue("p={format.pval(strong_interaction_pval$p.value,2)}"), glue::glue("p={format.pval(lethal_interaction_pval$p.value,2)}")), tip_length = .1
    )
ggsave(plot = percent_plot, filename = "Figure6D.pdf", width = 4, height = 4)
