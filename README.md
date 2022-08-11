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
- One or more **disease-network** files. Each disease-network file corresponds to one sample (cell line, patient sample, etc.) and contains the expression data for the sample. ***The disease-network files have a special nomenclature that MUST be followed.***
- A **reference protein-protein–interaction** (PPI) **network** file.
- A **drug-target data** file.
Further details are given in the following subsections.

## pData
The pData is an R **data.frame** that lists all the samples × drug treatments that you want to have scored, and it MUST contain at least the following three columns: **Sample**, **Drugs**, and **Dataset**. 
Additional columns are allowed, and the order of the columns is not important, but at the bare minimum those three columns must be present and named as such. 
Each row consists of one experimental observation: one sample (cell line, patient sample, etc.) treated with one or more drugs. 

- The **Sample** column is a **character vector** containing the identifier of the sample treated with the drugs. The identifier will be part of the disease-network file name; see the next subsection for how to name these files.
- The **Drugs** column is a **character vector** in which all the drugs used to treat the corresponding sample are listed in a single string, the individual drugs separated with an underscore (_).
- The **Dataset** column is a **character vector** containing an arbitrary identifier of the data set the sample is a part of.

Below is an example of a pData data frame:
| Sample | Drugs | Dataset |
| ---         |     ---      |          --- |
| A549   | carboplatin_paclitaxel    | NCI-60    |
| A549     | erlotinib       | NCI-60      |
| A549   | oxaliplatin_doxorubicin     | NCI-60    |
| UO-31     | imatinib       | NCI-60      |
| UO-31   | arsenic trioxide_paclitaxel     | NCI-60    |
| UO-31     | docetaxel       | NCI-60      |
| TCGA-2G-AAF1   | bevacizumab     | TCGA    |
| TCGA-2G-AAF1     | oxaliplatin_irinotecan      | TCGA      |

## Disease network

## Reference PPI network
## Drug-target data
The drug-target data is an R **list** in which item is a data.frame. 

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
