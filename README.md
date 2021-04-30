# Figure 7 creation scripts
R/ folder contains 3 scripts. 

- `R/create_SGA_data_combined.R` uses raw data from Costanzo 2016, downloaded from thecellmap.org and combines them into single dataframe

- `R/create_strain_ids.R` creates a dataframe containing systematic names, allele names and nonessential/essential categorization of all the ORFs in the Costanzo 2016

- `R/interactionNetworkAnalysis.R` uses files created by above scripts as well as overlapping_orfs.csv file in `analysis/data/raw_data/` folder to create figures 7A to 7D and data for 7E

Figure 7E was plotted using Cytoscape and then all figures are then combined using inkscape.
