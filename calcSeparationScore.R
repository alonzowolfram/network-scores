calcSeparationScore = function(j, drug_pairs, drug_targets, ref_network, ref_network_map, disease_map) {
  drug1 = drug_pairs[j,1]
  drug2 = drug_pairs[j,2]
  drug_targets_1 = drug_targets[[drug1]]
  drug_targets_2 = drug_targets[[drug2]]
    
  # For each target of drug1, calculate the shortest distance between it and another target of drug1. (I.e., calculate the shortest distance WITHIN the targets of drug1.)
  
  # b1 = all the drug targets for drug1, but their ENSEMBL names. 
  # From calcZScore.R, where b1 is named drug_targets_j_ensembl. 
  b1 = ref_network_map %>% 
    .[.$GeneSymbol %in% drug_targets_1$Target,"STRING_id"] %>% # dplyr version: dplyr::filter(GeneSymbol %in% drug_targets_j) %>% 
    # dplyr version: dplyr::select(STRING_id) %>% 
    # unlist %>% 
    as.character 
  
  if(nrow(drug_targets_1)==0 | nrow(drug_targets_2)==0 | is.null(drug_targets_1) | is.null(drug_targets_2)) return(NA)

  # Get all shortest distances from target k to other targets within the targets of drug 1.
  shortest_distances_1 = jubilee.mcsapply(drug_targets_1$Target, FUN = getShortestDistance, ref_network = ref_network, ref_network_map = ref_network_map, disease_map = disease_map, b = b1, within = T, mc.cores = num_cores)
  # Remove non-numeric values (e.g. Inf, NA, NaN) and take the mean.
  shortest_distances_1 = shortest_distances_1[is.finite(shortest_distances_1)]
  mean_shortest_distance_1 = ifelse(length(shortest_distances_1)==0, 0, mean(shortest_distances_1))
    
  # Get the mean shortest distance within the targets of drug 2. 
  # For each target of drug2, calculate the shortest distance between it and another target of drug2. (I.e., calculate the shortest distance WITHIN the targets of drug2.)
    
  # b2 = all the drug targets for drug2, but their ENSEMBL names. 
  # From calcZScore.R, where b2 is named drug_targets_j_ensembl. 
  b2 = ref_network_map %>% 
    .[.$GeneSymbol %in% drug_targets_2$Target,"STRING_id"] %>% # dplyr version: dplyr::filter(GeneSymbol %in% drug_targets_j) %>% 
    # dplyr version: dplyr::select(STRING_id) %>% 
    # unlist %>% 
    as.character 
    
  # Get all shortest distances from target k to other targets within the targets of drug 2. 
  shortest_distances_2 = jubilee.mcsapply(drug_targets_2$Target %>% as.character, FUN = getShortestDistance, ref_network = ref_network, ref_network_map = ref_network_map, b = b2, within = T, mc.cores = num_cores)
  # Remove non-numeric values (e.g. Inf, NA, NaN) and take the mean. 
  shortest_distances_2 = shortest_distances_2[is.finite(shortest_distances_2)]
  mean_shortest_distance_2 = ifelse(length(shortest_distances_2)==0, 0, mean(shortest_distances_2))
    
  # Calculate the shortest distances BETWEEN the targets of drug1 and drug2.
  shortest_distances_1_2 = jubilee.mcsapply(drug_targets_1$Target %>% as.character, FUN = getShortestDistance, ref_network = ref_network, ref_network_map = ref_network_map, disease_map = disease_map, b = b2, within = F, mc.cores = num_cores)
  shortest_distances_2_1 = jubilee.mcsapply(drug_targets_2$Target %>% as.character, FUN = getShortestDistance, ref_network = ref_network, ref_network_map = ref_network_map, disease_map = disease_map, b = b1, within = F, mc.cores = num_cores)
  # Remove non-numeric values (e.g. Inf, NA, NaN) and take the mean. 
  shortest_distances_1_2 = shortest_distances_1_2[is.finite(shortest_distances_1_2)]
  shortest_distances_2_1 = shortest_distances_2_1[is.finite(shortest_distances_2_1)]
  shortest_distances = c(shortest_distances_1_2, shortest_distances_2_1)
  mean_shortest_distance_1_2 = ifelse(length(shortest_distances)==0, 0, mean(shortest_distances))
    
  # Calculate the resulting separation (sAB) of the drug pair. 
  sAB = mean_shortest_distance_1_2 - 0.5*(mean_shortest_distance_1 + mean_shortest_distance_2)

  return(sAB)
}

getShortestDistance = function(target, ref_network, ref_network_map, disease_map, b, within) {
  # Convert the target to its corresponding Ensembl protein (peptide) ID.
  target_ensembl = ref_network_map[ref_network_map$GeneSymbol==target,"STRING_id"] %>% as.character
  target_ensembl = ifelse(length(target_ensembl) < 1, NA, target_ensembl)
  if(!(target_ensembl %in% V(ref_network)$name)) return(NA)
  
  if(within) {
    # Remove the target in question from the list of targets to calculate the distance from (otherwise, the distance would just be 0.) 
    c = b[b != target_ensembl] 
  } else {
    c = b
  }
  
  shortest_distance = min(distances(ref_network, v=V(ref_network)[name==target_ensembl], to=V(ref_network)[name %in% c], weights=NA)) # shortestDistance(target, drug_targets[[drug1]])
  
  return(shortest_distance)
}

## OLD CODE
#calcSeparationScore = function(j, drug_pairs, drug_targets, ref_network, ref_network_map, disease_map) {
#  drug1 = drug_pairs[j,1]
#  drug2 = drug_pairs[j,2]
#  drug_targets_1 = drug_targets[[drug1]]
#  drug_targets_2 = drug_targets[[drug2]]
#  
#  # For each target of drug1, calculate the shortest distance between it and another target of drug1. (I.e., calculate the shortest distance WITHIN the targets of drug1.)
#  
#  # drug_targets_1 = drug_targets_1[drug_targets_1$Target != target,,drop=F] # Remove the target in question from the list of targets to calculate the distance from (otherwise, the distance would just be 0.)  
#  # Convert the drug_targets_1 to their Ensembl protein (peptide) IDs.
#  # a and b _should_ be the same length, because in getTargets() we removed all the targets that do not have corresponding IDs (i.e., a should not be > b.)
#  a1 = as.character(unique(drug_targets_1$Target)) %>% .[order(.)]
#  b1 = ref_network_map[ref_network_map$GeneSymbol %in% a1 & ref_network_map$GeneSymbol %in% disease_map$GeneSymbol,] # Originally ref_network_map[ref_network_map$GeneSymbol %in% a1,]. See notes in getTargets.R. 2022-06-26.
#  b1 = b1[order(b1$GeneSymbol),"STRING_id"]
#  names(b1) = a1
#  
#  if(nrow(drug_targets_1)==0 | nrow(drug_targets_2)==0 | is.null(drug_targets_1) | is.null(drug_targets_2)) return(NA)
#  
#  # Get all shortest distances from target k to other targets within the targets of drug 1.
#  shortest_distances_1 = jubilee.mcsapply(drug_targets_1$Target, FUN = getShortestDistance, ref_network = ref_network, ref_network_map = ref_network_map, disease_map = disease_map, b = b1, within = T, mc.cores = num_cores)
#  # Remove non-numeric values (e.g. Inf, NA, NaN) and take the mean.
#  shortest_distances_1 = shortest_distances_1[is.finite(shortest_distances_1)]
#  mean_shortest_distance_1 = ifelse(length(shortest_distances_1)==0, 0, mean(shortest_distances_1))
#  
#  # Get the mean shortest distance within the targets of drug 2. 
#  # For each target of drug2, calculate the shortest distance between it and another target of drug2. (I.e., calculate the shortest distance WITHIN the targets of drug2.)
#  
#  # Convert the drug_targets_2 to their Ensembl protein (peptide) IDs.
#  # a and b _should_ be the same length, because in getTargets() we removed all the targets that do not have corresponding IDs (i.e., a should not be > b.)
#  a2 = as.character(unique(drug_targets_2$Target)) %>% .[order(.)] # Originally as.character(unique(drug_targets_2$Target)) %>% .[order(.)]. See notes in getTargets.R, 2022-06-26. 
#  b2 = ref_network_map[ref_network_map$GeneSymbol %in% a2 & ref_network_map$GeneSymbol %in% disease_map$GeneSymbol,] # Originally ref_network_map[ref_network_map$GeneSymbol %in% a2,]. See notes in getTargets.R. 2022-06-26.
#  b2 = b2[order(b2$GeneSymbol),"STRING_id"]
#  names(b2) = a2
#  
#  # Get all shortest distances from target k to other targets within the targets of drug 2. 
#  shortest_distances_2 = jubilee.mcsapply(drug_targets_2$Target %>% as.character, FUN = getShortestDistance, ref_network = ref_network, ref_network_map = ref_network_map, b = b2, within = T, mc.cores = num_cores)
#  # Remove non-numeric values (e.g. Inf, NA, NaN) and take the mean. 
#  shortest_distances_2 = shortest_distances_2[is.finite(shortest_distances_2)]
#  mean_shortest_distance_2 = ifelse(length(shortest_distances_2)==0, 0, mean(shortest_distances_2))
#  
#  # Calculate the shortest distances BETWEEN the targets of drug1 and drug2.
#  shortest_distances_1_2 = jubilee.mcsapply(drug_targets_1$Target %>% as.character, FUN = getShortestDistance, ref_network = ref_network, ref_network_map = ref_network_map, disease_map = disease_map, b = b2, within = F, mc.cores = num_cores)
#  shortest_distances_2_1 = jubilee.mcsapply(drug_targets_2$Target %>% as.character, FUN = getShortestDistance, ref_network = ref_network, ref_network_map = ref_network_map, disease_map = disease_map, b = b1, within = F, mc.cores = num_cores)
#  # Remove non-numeric values (e.g. Inf, NA, NaN) and take the mean. 
#  shortest_distances_1_2 = shortest_distances_1_2[is.finite(shortest_distances_1_2)]
#  shortest_distances_2_1 = shortest_distances_2_1[is.finite(shortest_distances_2_1)]
#  shortest_distances = c(shortest_distances_1_2, shortest_distances_2_1)
#  mean_shortest_distance_1_2 = ifelse(length(shortest_distances)==0, 0, mean(shortest_distances))
#  
#  # Calculate the resulting separation (sAB) of the drug pair. 
#  sAB = mean_shortest_distance_1_2 - 0.5*(mean_shortest_distance_1 + mean_shortest_distance_2)
#  
#  return(sAB)
#}