calcNetworkScores = function(i, pData, ref_network_list, drug_target_data, disease_network_folder, results_dir, data_prefix, deg_type, weighted = F, export_file = T, return_value = F, drug_combos = F, drug_drug_ixns = F, calc_z_score = F, allow_drugs_no_target = F, temp_dir, num_cores) {
  #######################################################################################################
  # i (integer): Number corresponding to a row in the pData data frame. 
  # pData (data.frame): Data frame containing the following columns: Drugs, Sample, Dataset. 
  # data_prefix (character): Character string containing the prefix to be affixed onto output files. 
  # ref_network_list (list): List containing the following elements: Network, Map. Network is an igraph network containing the full PPI interactome, and Map is a data.frame that acts as a dictionary mapping node names to other identifiers. 
  # data_dir (character): Character string containing the path to the input data.
  # results_dir (character): Character string containing the path to the output data. 
  # disease_network_folder (character): Character string containing the folder name (not full path!) with the disease networks. 
  # drug_target_data (list): List in which each element is a data.frame for one (1) drug. 
  # weighted (Boolean): Flag indicating whether to weight the scores.
  # export_file (Boolean): Flag indicating whether to export a file containing the results of row i. 
  # return_value (Boolean): Flag indicating whether to return the results of row i.
  # drug_combos (Boolean): Flag indicating whether to calculate scores applicable only to combinations of drugs.
  # drug_drug_ixns (Boolean): Flag indicating whether to calculate scores for drug-drug interactions.
  # allow_drugs_no_target (Boolean): Flag indicating whether to allow drugs that have no target data. If true, drugs without targets will simply not contribute to the final network scores. 
  #######################################################################################################
  
  print(paste0("Calculating the score for row ", i, "."))
  print(paste0("Weighted: ", weighted))
  print(paste0("Export file: ", export_file))
  print(paste0("Return value: ", return_value))
  print(paste0("Calculate Z-scores: ", calc_z_score))
  print(paste0("Drug combos: ", drug_combos))
  print(paste0("Drug-drug interactions: ", drug_drug_ixns))
  print(paste0("Allow drugs without target data: ", allow_drugs_no_target))
  
  ######################
  ## SET UP VARIABLES ##
  ######################
  if(!(i %in% 1:nrow(pData))) {
    print(paste0("This row is out of bounds for the pData! Skipping."))
    return(NULL)
  }
  drugs = pData[i, "Drugs"] %>% unlist() %>% as.character %>% strsplit("_") %>% unlist()
  drugs_oneline = pData[i, "Drugs"] %>% unlist()
  # Check if there are more than 1 drug passed as arguments.
  num_drugs = length(unique(drugs)) # num_drugs = 1 if there is only one drug in the "combination."
  sample = pData[i, "Sample"] %>% unlist() %>% as.character()
  dataset = pData[i, "Dataset"] %>% unlist() %>% as.character()
  
  if(length(drugs) < 1 | length(sample) < 1) {
    print("Please supply both a valid sample and drug treatment.")
    return(NULL)
  }
  if(length(dataset) < 1)  {
    print("No valid dataset supplied, defaulting to NA.")
    dataset = NA
  }
  
  # Check if the temporary directory has been created. If not, create it. 
  #if(!dir.exists(temp_dir)) dir.create(temp_dir)
    
  # Check if this one has been done already. 
  file_name = paste0(results_dir, "/network_scores_", data_prefix, "_", sample, "_", drugs_oneline, "_", deg_type, ".csv") # Originally written to results_dir. 
  if(file.exists(file_name)) {
    print(paste0("The sample x drug combination consisting of sample ", sample, " and drug(s) ", drugs, " has already been done! Skipping to the next row in the pData."))
    return(NULL)
  } 
  
  ##########################
  ## LOAD DISEASE NETWORK ##
  ##########################
  # Load the appropriate network. 
  disease_network_file_name = paste0(disease_network_folder, "/", sample, "_disease_network_", deg_type, ".rds")
  if(!file.exists(disease_network_file_name)) {
    print(paste0("The disease network ", disease_network_file_name, " does not exist! Skipping to the next row in the pData."))
    return(NULL)
  }
  
  possibleError = tryCatch(
    {
      disease_network_list = readRDS(disease_network_file_name)
    }, error=function(cond) {
      message(paste0("There was a problem in calcNetworkScore() loading the network for sample ", sample, ": "))
      message(cond)
      cond
    }
  )
  
  if(inherits(possibleError, "error") | !exists("disease_network_list")) { # Network could not be loaded correctly. 
    score_table = handleError(sample, i, "Network could not be loaded correctly.", export_file)
    
    # Return the value if indicated.
    if(return_value) return(score_table) 
      
    return(NULL)
  }
  
  disease_map = disease_network_list$Map
  disease_nodes = disease_map %>% 
    .[.$DiseaseModuleGene=="Y","STRING_id"] %>% # dplyr version, removed to allow use on TACC Stampede2: dplyr::filter(DiseaseModuleGene=="Y") %>% 
    #.$STRING_id %>% # dplyr::select(STRING_id) %>% 
    #unlist %>% 
    as.character
  disease_intensity = disease_network_list$Intensity
  
  ###################################
  ## LOAD REFERENCE (FULL) NETWORK ##
  ###################################
  ref_network = ref_network_list$Network
  ref_network_map = ref_network_list$Map
  ref_network_diameter = ref_network_list$Diameter
  
  ###############################
  ## GET TARGETS FOR EACH DRUG ##
  ###############################
  possibleError = tryCatch( # https://stackoverflow.com/a/8094059, https://stackoverflow.com/a/12195574
    {
      drug_targets = lapply(drugs, FUN = getTargets, 
                            #mc.cores = num_cores,
                            disease_map = disease_map, 
                            sample = sample, 
                            drug_target_data = drug_target_data)
      names(drug_targets) = drugs
    }, error=function(cond) {
      message(paste0("There was a problem retrieving the targets for the drugs ", drugs, ": "))
      message(cond)
      cond
    }
  )
  
  if(inherits(possibleError, "error") | sum(is.na(drug_targets)) > 0) {
    score_table = handleError(sample, i, paste0("Problem retrieving the targets for the drugs ", drugs, "."), export_file)
    
    # Return the value if indicated.
    if(return_value) return(score_table) 
    
    return(NULL)
  }
  
  # If the allow_drugs_no_target flag is false, check if any of the elements of the list drug_targets is null or has no data.
  if(!allow_drugs_no_target) {
    drug_targets_null = 0
    for(element in names(drug_targets)) {
      if(is.null(drug_targets[[element]])) {drug_targets_null = drug_targets_null + 1}
      else if(nrow(drug_targets[[element]]) < 1) {drug_targets_null = drug_targets_null + 1}
    }
    if(drug_targets_null > 0) {
      print(paste0("One or more of the drugs in ", paste(drugs, collapse = ", "), " has no target data. Skipping this row."))
      return(NULL)
    }
  }
  
  #######################
  ## SCORE CALCULATION ##
  #######################
  ##########################
  ## 1. CALCULATE Z-SCORE ##
  ##########################
  if(calc_z_score) {
    print(paste0("Calculating z-score for row ", i, "."))
    
    # Get the degree of every node in the network. 
    degree_all_nodes = igraph::degree(ref_network)
    # Get the number of vertices in the network. 
    order_ref_network = gorder(ref_network)
    
    # Convert the proteins of the current disease module to ENSEMBL.
    disease_module_ensembl = disease_map %>% 
      .[.$DiseaseModuleGene=="Y","STRING_id"] %>% 
      #.$STRING_id %>% 
      #unlist %>% 
      as.character() # Dplyr version, removed to allow use on TACC Stampede2 servers: disease_module_ensembl = disease_map %>% dplyr::filter(DiseaseModuleGene=="Y") %>% dplyr::select(STRING_id) %>% unlist %>% as.character()
    disease_module_dist = degree_all_nodes[names(degree_all_nodes) %in% disease_module_ensembl]
    order_disease_module = length(disease_module_dist)
    
    z_score = lapply(drugs, FUN = calcZScore, 
                     #mc.cores = min(length(drugs), num_cores),
                     drug_target_data = drug_target_data, # drug_target_data contains the data on all drugs, not just the ones for this row in the pData (sample x drug [x drug] combination), while drug_targets contains the data on just the drugs for this row in the pData. 
                     ref_network = ref_network,
                     ref_network_map = ref_network_map,
                     ref_network_diameter = ref_network_diameter,
                     disease_map = disease_map, 
                     disease_module_dist = disease_module_dist, 
                     disease_module_ensembl = disease_module_ensembl, 
                     degree_all_nodes = degree_all_nodes, 
                     order_ref_network = order_ref_network, 
                     order_disease_module = order_disease_module,
                     num_cores = num_cores)
    z_score = do.call(rbind, z_score)
    z_score = ifelse(is.finite(z_score), z_score, NA)
    # Add up all the Z-scores. 
    # Previously we normalized by the number of drugs to get the final drug-network overlap score, but come to think of it, I'm not sure that makes sense. 
    z_score_final = apply(z_score, 2, sum, na.rm = T) #/ length(drugs)
    z_score_final_unweighted = z_score_final[1]
    z_score_final_weighted_normalized = z_score_final[2]
    z_score_final_weighted_unnormalized = z_score_final[3]
    # For troubleshooting: 
    # drug = drugs[1]
    # Zscore_test = calcZScore(drug = drugs[1], drug_target_data = drug_target_data, ref_network = ref_network, ref_network_map = ref_network_map, disease_map = disease_map, disease_module_dist = disease_module_dist, disease_module_ensembl = disease_module_ensembl, degree_all_nodes = degree_all_nodes, order_ref_network = order_ref_network, order_disease_module = order_disease_module)
  } else {
    z_score_final_unweighted = NA
    z_score_final_weighted_normalized = NA
    z_score_final_weighted_unnormalized = NA
  }
  
  #############################
  ## 2. CALCULATE CENTRALITY ##
  #############################
  print(paste0("Calculating centrality for row ", i, "."))
  
  centrality_scores = lapply(drugs, FUN = calcCentrality, 
                             #mc.cores = num_cores,
                             ref_network = ref_network, 
                             ref_network_map = ref_network_map,
                             disease_map = disease_map, 
                             drug_targets = drug_targets
                             )
  centrality_scores = do.call(rbind, centrality_scores)
  centrality_scores_sums = apply(centrality_scores, 2, sum, na.rm=T) # calcCentrality() calculates the degree centrality for a single drug; each row in centrality_scores corresponds to ONE drug, hence our summing the rows up. 
  
  centrality_final_unweighted = centrality_scores_sums[1]
  centrality_final_weighted_normalized = centrality_scores_sums[2]
  centrality_final_weighted_unnormalized = centrality_scores_sums[3]
  #bridging_centrality_weighted = centrality_scores_sums[4]
  #bridging_centrality_unweighted = centrality_scores_sums[5]
  
  ################################
  ## 3. CALCULATE WINTHER SCORE ##
  ################################
  print(paste0("Calculating WINTHER score for row ", i, "."))
  
  winther_score = sapply(drugs, FUN = calcWINTHERscore, # Originally jubilee.mcsapply()
                         #mc.cores = num_cores,
                         drug_targets = drug_targets, 
                         ref_network = ref_network, 
                         ref_network_map = ref_network_map, 
                         disease_map = disease_map,
                         intensity = disease_intensity
                         )
  winther_score = ifelse(is.finite(winther_score), winther_score, 0)
  # Add up all the WINTHER scores for individual drugs. 
  # Previously we normalized by the number of drugs to get the final drug-network overlap score, but come to think of it, I'm not sure that makes sense. 
  winther_score_final = sum(winther_score, na.rm=T) #/ length(drugs)
  
  ###################################
  ## 4. CALCULATE SEPARATION SCORE ##
  ###################################
  # Separation score = drug-drug overlap.
  if(num_drugs < 2 | !drug_combos) {
    # There is only 1 unique drug or we're not calculating combination-specific scores. Do not calculate the overlap between drugs.
    sAB_final = NA
  } else {
    print(paste0("Calculating drug-drug overlap for row ", i, "."))
    
    # There are >= 2 unique drugs in the combination.
    # Create a list of drug pairs. 
    drug_pairs = createListOfDrugPairs(drugs, "combination")
    
    # For each drug pair AB, calculate the overlap (distance; sAB) between the targets of the pair.
    # The more positive the sAB, the greater the separation/distance between the targets of the pair.
    sABs = sapply(1:nrow(drug_pairs), FUN = calcSeparationScore, # Originally jubilee.mcsapply()
                  #mc.cores = num_cores,
                  drug_pairs = drug_pairs, 
                  drug_targets = drug_targets, 
                  ref_network = ref_network, 
                  ref_network_map = ref_network_map,
                  disease_map = disease_map
                  )
    # Add up all the sABs to get the final sAB. Normalize by the number of drug pairs. (In other words, calculate the average sAB for all pairwise drug combinations.)
    sAB_final = sum(sABs, na.rm=T) / nrow(drug_pairs)
  }
  
  #########################################
  ## 5. CALCULATE DRUG-DRUG INTERACTIONS ##
  #########################################
  if(num_drugs < 2 | !drug_drug_ixns) {
    # There is only 1 unique drug. Do not calculate the chemical interactions between drugs.
    dd_chem_ixns_final = NA
  } else {
    print(paste0("Calculating drug-drug interactions for row ", i, "."))
    
    # There are >= 2 unique drugs in the combination.
    # So we will need to generate a list of all the _PERMUTATIONS_ of drug A and drug B. 
    drug_perms = createListOfDrugPairs(drugs, "permutation")
    
    dd_chem_ixns_1 = jubilee.mcsapply(1:nrow(drug_perms), FUN = calcDrugDrugChemIxns, drug_perms = drug_perms, map = map, ref_network = ref_network, sample = sample, mc.cores = num_cores)
    dd_chem_ixns_2 = jubilee.mcsapply(1:nrow(drug_perms), FUN = calcCotreatments, drug_perms = drug_perms, map = map, ref_network = ref_network, sample = sample, mc.cores = num_cores)
    dd_chem_ixns_1 = ifelse(is.finite(dd_chem_ixns_1), dd_chem_ixns_1, 0)
    dd_chem_ixns_2 = ifelse(is.finite(dd_chem_ixns_2), dd_chem_ixns_2, 0)
    dd_chem_ixns = dd_chem_ixns_1 + dd_chem_ixns_2
    dd_chem_ixns_final = sum(dd_chem_ixns, na.rm=T)
  }
  
  ############
  ## OUTPUT ##
  ############
  score_table = data.frame(
    ZScoreUnweighted = z_score_final_unweighted,
    ZScoreWeightedNormalized = z_score_final_weighted_normalized, 
    ZScoreWeightedUnnormalized = z_score_final_weighted_unnormalized, 
    WINTHERscore = winther_score_final,
    CentralityUnweighted = centrality_final_unweighted,
    CentralityWeightedNormalized = centrality_final_weighted_normalized,
    CentralityWeightedUnnormalized = centrality_final_weighted_unnormalized,
    sAB = sAB_final, 
    DrugDrugInteractions = dd_chem_ixns_final,
    Sample = sample, 
    Dataset = dataset, 
    Drugs = drugs_oneline)
  
  # Save to CSV file if indicated.
  if(export_file) write.table(score_table, file_name, row.names = F, col.names = F, sep = ',')
  
  # Write to database if indicated. 
  
  # Return the value if indicated.
  if(return_value) return(score_table)
  
}

