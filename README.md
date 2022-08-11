# Introduction
This repository contains all the code used to calculate the network scores as seen in our forthcoming paper (details to be provided as soon as it is published). 

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

## Required input 1: pData
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

## Required input 2: Disease network(s)
***(Please make sure to read the note on file nomenclature at the end of this section!)***

Each disease-network file corresponds to one sample (cell line, patient sample, etc.) and contains differential gene expression (DGE) data and disease-module membership data for each gene (whether each gene is part of the sample's disease module). It consists of an R **list** that contains at least the following items that must be named as such: 1) **Map** and 2) **Intensity**. Additional items are allowed, but will not be used by the program. 

- The $**Map** item is a data frame of dimensions *n* × *m*, where *n* is the number of genes measured for the sample, and *m* >= 4, with at least the following four columns:
  - **GeneSymbol**. This is a character vector with the HUGO gene symbol for each gene. 
  - **RankChange**. This is an integer or numeric vector with the change in rank between normal and tumor tissue for each gene. 
  - **DiseaseModuleGene**. This is a character vector indicating whether each gene is part of the disease module ("Y") or not ("N"). 
  - **STRING_id**. This is a character vector with the corresponding STRING identifier for each gene. 
- The $**Intensity** item is a data frame containing the _raw_ intensities of each measured gene. These intensities will come from the appropriate raw microarray or RNA-seq files. 

We are aware that the **Intensity** data frame could probably be merged into the **Map** data frame, and the disease-network file could be reduced to simply a data frame instead of a list. This will be addressed in future updates. 

> ***Important note on nomenclature!*** The disease-network files MUST use the following format for the filename: 
> 
> `[sample name]_disease_network_[DEG type].rds`
> 
> where `[sample name]` matches exactly the appropriate value of the **Sample** column in the pData data frame
> 
> and `[DEG type]` is an arbitrary indicator of how the differentially expressed genes (DEGs) were calculated. For example, I might put `2SD` to indicate that the genes whose rank change was >= 2SD away from the mean rank change were chosen as DEGs.

## Required input 3: Reference PPI network
The reference PPI network file should contain an R **list** with the following items:
- The $**Map** item is a data frame of dimensions *n* × *m*, where *n* is the number of genes measured for the sample, and *m* >= 4, with at least the following four columns:
  - **GeneSymbol**. This is a character vector with the HUGO gene symbol for each gene. 
  - **STRING_id**. This is a character vector with the corresponding STRING identifier for each gene. 
- The $**Network** item is an iGraph object containing the PPI network.
- The $**Diameter** item is a numeric scalar consisting of the diameter of the $**Network** item.

## Required input 4: Drug-target data
The drug-target data is an R **list** in which each item corresponds to a single drug and is a data frame. The names of the items should match exactly the names as they appear in the pData $**Drugs** column. For example, if one of the treatments in the pData data frame is **oxaliplatin_fluorouracil**, there should be at least two items in the drug-target–data list, named **oxaliplatin** and **fluorouracil**, respectively. If, for example, you named the fluorouracil entry 5-fluorouracil, this would not work; the names much match up exactly. (Alternatively, you could change the pData entry to **oxaliplatin_5-fluorouracil**.)

Each data frame should contain at least the following three columns, **Drug**, **Target**, and **Effect**, named as such.
- The **Drug** column is a **character vector** containing the name of the drug. Because each item in the list corresponds to only one drug, every entry in this column should be the same. 
- The **Target** column is a **character vector** containing the drug's targets. ***These should be HUGO gene names.***
- The **Effect** column is a **numeric vector** the directional effect of the drug on the corresponding target. The values of these will be –1 (inhibition or downregulation), 0 (neutral or unknown), and +1 (activation or upregulation). 

Below is a sample data frame:
| Drug | Target | Effect |
| ---         |     ---      |          --- |
| sunitinib   | CSF1    | –1    |
| sunitinib   | CSF1R       | –1     |
| sunitinib  | FLT1     | –1    |
| sunitinib     | FLT3       | –1      |
| sunitinib  | FLT4     | –1    |
| sunitinib    | KDR       | –1      |
| sunitinib   | KIT     | –1    |
| sunitinib     | RET      | –1      |

# How to calculate network scores
The basic function for calculating network scores is `calcNetworkScores()`, found in the file `calcNetworkScores.R`. Its usage is as follows, with explanation of the parameters after the code:
```R
calcNetworkScores(
                  pData = pData,
                  i,  
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
                  num_cores = num_cores
                  )
```
- **pData**. This is the path to the pData file.
- **i**. This is an integer telling `calcNetworkScores()` which row of the pData to calculate the network scores for. 
- **drug_target_data**. This is the path to the drug-target–data file.
- **disease_network_folder**. This is the path to the folder containing the disease-network file(s).
- **results_dir**. This is the path to the folder in which the results files will be placed.
- **data_prefix**. This is an arbitrary string denoting the data set the results are a part of. You can name it whatever you want, but it is recommended to pick a descriptive name. 
- **deg_type**. This is an arbitrary string denoting how the DEGs were generated. You can name it whatever you want, but it is recommended to pick a descriptive name that reflects the methodology used to generate the DEGs. For example, I might put `2SD` to indicate that the genes whose rank change was >= 2SD away from the mean rank change were chosen as DEGs.
- **weighted**. This is a Boolean value indicating whether or not to weight the following network scores: Z-score, centrality.
- **export_file**. This is a Boolean value indicating whether or not to export (CSV) files of the results.
- **return_value**. This is a Boolean value indicating whether or not the function should return the scores. 
- **drug_combos**. This is a Boolean value indicating whether or not to calculate the separation score. 
- **drug_drug_ixns**. This is a Boolean value indicating whether or not to calculate drug-drug interactions.
- **calc_z_score**. This is a Boolean value indicating whether or not to to calculate Z-scores.
- **allow_drugs_no_target**. This is a Boolean value indicating whether or not to allow drugs that do not have corresponding drug-target data. In this case, these drugs will contribute 0 to the network scores. If `FALSE`, if a drug combination contains at least one drug without corresponding drug-target data, the entire row will be skipped.
- **temp_dir**. This is the path to the directory in which to store temporary files.
- **num_cores**. This is an integer indicating how many cores to use for parallelization. 
