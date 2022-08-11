# Introduction
This repository contains all the code used to calculate the network scores as seen in the paper (...). 

# R version
All code was written for R version 4.0.3 (2020-10-10) -- "Bunny-Wunnies Freak Out" on a x86_64-pc-linux-gnu (64-bit) environment. The code will probably work on other versions of R from at least R version 3.6.3 and on other platforms, but this is not guaranteed. 

# Required libraries
You will need the following libraries:
```R
library(plyr)
library(dplyr)
library(data.table) # Faster and more memory-friendly version of data frames. 
library(gtools) # For making permutations. 
library(igraph) # Network visualization and calculation. 
library(jubilee) # For the jubilee.mcsapply function.
library(matrixcalc) # Hadamard product. 
library(parallel)
library(regexPipes)
library(stringi)
library(stringr)
library(tibble)
library(tidyr)
library(miceadds) # For the source.all function. 
```
# Required inputs
To calculate network scores, you will need the following input files (all .RDS format):
- A **pData** file.
- One or more **disease-network** files. Each disease network file corresponds to one sample (cell line, patient sample, etc.) and contains the expression data for the sample.
- A **reference protein-protein–interaction** (PPI) **network** file.
- A **drug-target data** file.
Further details are given in the following subsections.

## pData
## Disease network
## Reference PPI network
## Drug-target data

# How to calculate network scores
The basic function for calculating network scores will 
```R
calcNetworkScores(i, 
                  pData = pData, 
                  ref_network_list = ref_network_list, 
                  drug_target_data = drug_target_data, 
                  disease_network_folder = disease_network_folder, 
                  results_dir = results_dir, 
                  data_prefix = data_prefix, 
                  deg_type = deg_type, 
                  weighted = T, 
                  export_file = T, 
                  return_value = F, 
                  drug_combos = T, 
                  drug_drug_ixns = F, 
                  calc_z_score = T, 
                  allow_drugs_no_target = allow_drugs_no_target, 
                  temp_dir = temp_dir, 
                  num_cores = num_cores)
```
