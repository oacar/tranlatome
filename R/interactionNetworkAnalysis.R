## Inputs -------------

library(tidyverse)
library(glue)
library(igraph)
library(ggsignif)
library(cowplot)
library(readxl)

pgs <- readr::read_csv('analysis/data/raw_data/transient_annotated_intergenic_041421')
overlappingorfs <- read_csv("analysis/data/raw_data/overlappingorfs.csv")
net.df <- readRDS("analysis/data/derived_data/SGA_data_combined.rds.gz")
net_df_significant <- filter(net.df, `P-value` <= 0.05 & `Query Strain ID` %in% overlappingorfs$orf_name == FALSE & `Array Strain ID` %in% overlappingorfs$orf_name == FALSE)
net_df_significant_sl <- filter(net_df_significant, `Genetic interaction score (ε)` <= -0.2)#& `Double mutant fitness standard deviation`<=0.1)
exp.data <- read_csv("analysis/data/derived_data/strain_ids_with_experiment_count_all.csv")%>% mutate(group=ifelse(`Systematic gene name`%in%pgs$orf_name,'proto-gene',maincat))
strain_ids_smf <- read_excel("analysis/data/raw_data/strain_ids_and_single_mutant_fitness.xlsx")



## Calculate strong (eps<-0.2) and lethal (eps<-0.35) interactions -----
strong_orf_list <- net.df %>% filter(`P-value`<=0.05,`Genetic interaction score (ε)` <= -0.2) %>% select(`Query Strain ID`, `Array Strain ID`) %>% as.list() %>% unlist() %>% unique()
lethal_orf_list <- net.df %>% filter(`P-value`<=0.05,`Genetic interaction score (ε)` <= -0.35) %>% select(`Query Strain ID`, `Array Strain ID`) %>% as.list() %>% unlist() %>% unique()
exp.data <- exp.data %>% mutate(strong_interaction=ifelse(`Systematic gene name`%in%strong_orf_list,TRUE,FALSE),
                          lethal_interaction=ifelse(`Systematic gene name`%in%lethal_orf_list,TRUE,FALSE))

exp_data_nones=exp.data %>% filter(group!='essential')

## Compare interaction ratios of transient orfs and nonessential orfs with fisher.test ----
strong_interaction_pval <- table(exp_data_nones$group,exp_data_nones$strong_interaction) %>% fisher.test() #0.04652
lethal_interaction_pval <- table(exp_data_nones$group,exp_data_nones$lethal_interaction) %>% fisher.test() #0.02636

percent_plot_data <- exp_data_nones %>%
  select(group,strong_interaction) %>%
  group_by(group,strong_interaction) %>%
  summarise(str_count=n()) %>%
  mutate(freq = str_count / sum(str_count)) %>%
  filter(strong_interaction==TRUE) %>% mutate(cat='strong') %>% #%>% gather()
  bind_rows(
    exp_data_nones %>% select(group,lethal_interaction) %>% group_by(group,lethal_interaction) %>% summarise(leth_count=n()) %>% mutate(freq = leth_count / sum(leth_count)) %>% filter(lethal_interaction==TRUE) %>% mutate(cat='lethal')
  ) %>% mutate(cat=factor(cat,levels=c('strong','lethal')))%>%
  mutate(group = factor(group, levels = c("proto-gene", "nonessential")))


percent_plot <-percent_plot_data %>%
  ggplot(aes(x=cat,y=freq,fill=group))+
  geom_bar(stat = "identity", position = position_dodge(width = 0.95)) + theme_classic() + theme(
    axis.title.x = element_text(size = 8), axis.text.x = element_text(size = 10), axis.text.y = element_text(size = 8),axis.title.y=element_text(size=12),
    legend.text = element_text(size = 8), legend.position = "right", legend.title = element_blank(), legend.margin = margin(t=-0.25,unit='in'),
    legend.spacing.x = unit(0.1, 'in')
  ) +
  scale_fill_manual(labels = c("Transient ORFs", "Nonessential\nGenes"), values = c("#1CBDC2", "#EF3D23")) +
  scale_x_discrete(labels = c(expression(epsilon * "< -0.2"), expression(epsilon * "<-0.35"))) +
  ylab("Ratio with at least one\ninteraction at given threshold") + scale_y_continuous(labels = scales::percent) + xlab("") +
  geom_signif(
    y_position = c(1.05, 0.85), xmin = c(0.8, 1.8), xmax = c(1.2, 2.2),
    annotation = c(glue("p={format.pval(strong_interaction_pval$p.value,2)}"),glue("p={format.pval(lethal_interaction_pval$p.value,2)}")), tip_length = .1
  )
ggsave(plot = percent_plot,filename = 'figures/Figure7B.pdf',width = 4,height = 4)


#compare smf to wild type fitness of 1 -----
strain_ids_smf$smf <- strain_ids_smf[, 4] %>%
  pull() %>%
  as.numeric()
strain_ids_joined <- exp.data %>%
  inner_join(strain_ids_smf) %>%
  pivot_wider(names_from = group, values_from = smf)

wt_pvalue <- strain_ids_joined$`proto-gene`[is.na(strain_ids_joined$`proto-gene`) == F] %>% t.test(mu = 1) #p=0.05663

#compare nonessential smf to transient smf distributions ---------
figure7a <- strain_ids_joined %>%
  # filter(subcat != "essential") %>%
  ggplot() +
  geom_histogram(aes(x = `proto-gene`), fill = "#21BDC2", color = "#21BDC2", binwidth=0.01) +
  geom_density(aes(x = nonessential), color = "#EF4024") +
  theme_classic() +
  scale_x_continuous(name = "Single Mutant Fitness", expand = c(0, 0), limits = c(0.5, NA)) +
  scale_y_continuous(name = "Counts", expand = c(0, 0), limits=c(0,20)) +
  theme(axis.title = element_text(size = 9), axis.text = element_text(size = 6))
ggsave("figures/Figure7A.pdf", plot=figure7a, width = 2, height = 2)

##Number of interactions histogram -----
colname <- colnames(exp.data)[2]
pgs_allele <- exp.data %>%
  filter(group=='proto-gene') %>%
  select(`Allele Gene name`) %>%
  pull()
allnet <- graph_from_data_frame(net_df_significant_sl[,c(2,4,6)],directed=FALSE)
allnet <- simplify(allnet)
allnet.deg<-degree(allnet)



#if(plot_histogram){
allnet.deg.pgs <- allnet.deg[names(allnet.deg)%in%pgs_allele]
ggplot()+geom_histogram(binwidth=1,aes(x=allnet.deg.pgs),fill="#21BDC2")+theme_classic()+xlab('Number of interactions at  < -0.2')+
  theme(axis.title = element_text(size = 9), axis.text = element_text(size = 6))+ylab('Count')
ggsave('figures/Figure7D.pdf',width=2,height = 2)

## get yer175w-a ego graph to csv -----
induced_subgraph(allnet,ego(allnet,order=1,nodes='yer175w-a')[[1]])%>%get.edgelist()%>%as.data.frame()%>%
  write_csv('analysis/data/derived_data/yer175w-a_ego.csv', col_names = FALSE)




## Interaction Density -----
thr <- -0.2
allnet <- graph_from_data_frame(net_df_significant[,c(2,4,6)],directed=FALSE)
adj <- as_adjacency_matrix(allnet,type='both',attr="Genetic interaction score (ε)") %>% as.matrix()
adj_nonw <- as_adjacency_matrix(allnet,type='both') %>% as.matrix()
colname <- colnames(exp.data)[2]

exp.data$essint <- NA
exp.data$nonesint <- NA
exp.data$esscount <- NA
exp.data$nonescount <- NA
for( i in seq_along(exp.data$`Allele Gene name`)){
  name <- exp.data$`Allele Gene name`[i]
  if(name%in%colnames(adj)){
    interactions <- adj[name,]
    int_nonw <- adj_nonw[name,]
    lethal_int_name <- names(interactions[interactions<=thr])
    exp.data$essint[i] <- sum(lethal_int_name%in%exp.data$`Allele Gene name`[exp.data$maincat=='essential'])
    exp.data$nonesint[i] <- sum(lethal_int_name%in%exp.data$`Allele Gene name`[exp.data$maincat=='nonessential'])
    #names(int_nonw[int_nonw==1])
    exp.data$esscount[i] <- sum(names(int_nonw[int_nonw==1])%in%exp.data$`Allele Gene name`[exp.data$maincat=='essential'])
    exp.data$nonescount[i] <- sum(names(int_nonw[int_nonw==1])%in%exp.data$`Allele Gene name`[exp.data$maincat=='nonessential'])
  }

}

exp.data <- exp.data %>%
  dplyr::mutate(essint_density = essint / esscount, nonesint_density = nonesint / nonescount)

pe12 = coin::independence_test(group~essint,exp.data %>% mutate(group = as.factor(group)) %>% filter(group!='essential'))
pn12 = coin::independence_test(group~nonesint,exp.data %>% mutate(group = as.factor(group)) %>% filter(group!='essential'))
#1.028e-05 nones~pgs for nones ints
#0.0165 nones~pgs for ess ints
#p-value = 0.00124 nones~pgs for nones ints
#p-value = 0.01571 nones~pgs for ess ints
exp.data$group <- factor(exp.data$group,levels = c('proto-gene','nonessential','essential'))
plt <- exp.data %>%
  select(group, essint_density, nonesint_density) %>%
  gather(type, int_density, -group) %>%
  ggplot(aes(x = type, y = int_density, fill = group)) + # geom_boxplot()
  stat_summary(fun.y = mean, geom = "bar", position = position_dodge(width = 1)) +
  stat_summary(fun.data = mean_se, geom = "errorbar", position = position_dodge(width = 1), width = 0.5) +
  ylab("Interaction density") +
  xlab("") +
  scale_x_discrete(labels = c("Interactions with\nEssential genes", "Interactions with\nNonssential genes")) +
  theme_classic() + theme(
    axis.title.x = element_text(size = 8), axis.text.x = element_text(size = 12), axis.text.y = element_text(size = 10),axis.title.y = element_text(size = 12),
    legend.text = element_text(size = 8), legend.position = "bottom", legend.title = element_blank(), legend.margin = margin(t=-0.25,unit='in'),
    legend.spacing.x = unit(0.1, 'in')
  ) +
  scale_fill_manual(labels = c("Transient ORFs", "Nonessential\nGenes", "Essential\nGenes"), values = c("#1CBDC2", "#EF3D23", "#FAA51A")) +
  geom_signif(
    annotations = c(paste0("p= ", formatC(coin::pvalue(pe12), digits = 2)), paste0("p= ", formatC(coin::pvalue(pn12), digits = 1))),
    y_position = c(0.07, 0.05), xmin = c(0.65, 1.65), xmax = c(1, 2), tip_length = 0.005, textsize = 4
  )
ggsave('figures/Figure7C.pdf',plot= plt , width = 4,height = 3)

