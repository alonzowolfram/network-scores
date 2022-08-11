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
The pData file consists of an R **data.frame** that lists all the samples × drug treatments that you want to have scored, and it MUST contain at least the following three columns: **Sample**, **Drugs**, and **Dataset**. 
Additional columns are allowed, and the order of the columns is not important, but at the bare minimum those three columns must be present and named as such. 
Each row consists of one experimental observation: one sample (cell line, patient sample, etc.) treated with one or more drugs. 

- The **Sample** column is a **character vector** containing the identifier of the sample treated with the drugs. The identifier will be part of the disease-network file name; see the next subsection for how to name these files.
- The **Drugs** column is a **character vector** in which all the drugs used to treat the corresponding sample are listed in a single string, the individual drugs separated with an underscore (_).
- The **Dataset** column is a **character vector** containing an arbitrary identifier of the data set the sample is a part of.

Below is an example of a pData data frame:
| Sample | Drugs | Dataset | CellGrowthPercent |
| ---         |     ---      |          --- |     ---      |
| A549   | carboplatin_paclitaxel    | NCI-60    | 27   |  
| A549     | erlotinib       | NCI-60      | 102   |
| A549   | oxaliplatin_doxorubicin     | NCI-60    | 84    |
| UO-31     | imatinib       | NCI-60      | -30   |
| UO-31   | arsenic trioxide_paclitaxel     | NCI-60    | 50    |
| UO-31     | docetaxel       | NCI-60      | 95   |
| TCGA-2G-AAF1   | bevacizumab     | TCGA    | 43    |
| TCGA-2G-AAF1     | oxaliplatin_irinotecan      | TCGA      | 22    |

## Disease network(s)
***(Please make sure to read the note on file nomenclature at the end of this section!)***

Each disease-network file corresponds to one sample (cell line, patient sample, etc.) and contains differential gene expression (DGE) data and disease-module membership data for each gene (whether each gene is part of the sample's disease module). It consists of an R **list** that contains at least the following items that must be named as such: 1) **Map** and 2) **Intensity**. Additional items are allowed, but will not be used by the program. 

- The $**Map** item is a data frame of dimensions *n* × *m*, where *n* is the number of genes measured for the sample, and *m* >= 4, with at least the following four columns:
  - **GeneSymbol**. This is a character vector with the HUGO gene symbol for each gene. 
  - **RankChange**. This is an integer or numeric vector with the change in rank between normal and tumor tissue for each gene. 
  - **DiseaseModuleGene**. This is a character vector indicating whether each gene is part of the disease module ("Y") or not ("N"). 
  - **STRING_id**. This is a character vector with the corresponding STRING identifier for each gene. 
- The $**Intensity** item is a data frame.

We are aware that the **Intensity** data frame could probably be merged into the **Map** data frame, and the disease-network file could be reduced to simply a data frame instead of a list. This will be addressed in future updates. 

> ***Important note on nomenclature!*** The disease-network files MUST use the following format for the filename: 
> 
> `[sample name]_disease_network_[DEG type].rds`
> 
> where `[sample name]` matches exactly the appropriate value of the **Sample** column in the pData data frame
> 
> and `[DEG type]` is an arbitrary indicator of how the differentially expressed genes (DEGs) were calculated. For example, I might put `2SD` to indicate that the genes whose rank change was >= 2SD away from the mean rank change were chosen as DEGs.

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
