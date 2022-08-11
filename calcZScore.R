calcZScore = function(drug, drug_target_data, ref_network, ref_network_map, ref_network_diameter, disease_map, disease_module_dist, disease_module_ensembl, degree_all_nodes, order_ref_network, order_disease_module, num_cores, d_X_Y_rand_dist_size = 1000) {
  # This function takes a drug (as an argument), a list of the drug's targets, and a protein- or gene-interaction network and returns a score measuring the overlap of the drug's targets with the network's vertices. 
  # calcZScore() -> numeric
  
  # Get the targets of the drug.
  # If there are no drug targets, return NA for the Z-score. 
  if(is.null(drug_target_data[[drug]]) | nrow(drug_target_data[[drug]]) < 1) return(NA)
  # Otherwise...
  
  # Convert the targets of the current drug(s) to ENSEMBL.
  drug_targets_j = drug_target_data[[drug]] %>% 
    .$Target %>% # dplyr::select(Target) %>% 
    unlist %>% 
    as.character #drug_target_data %>% dplyr::filter(Drug==drug) %>% dplyr::select(Target) %>% unlist
  drug_targets_j_ensembl = ref_network_map %>% 
    .[.$GeneSymbol %in% drug_targets_j,"STRING_id"] %>% # dplyr version: dplyr::filter(GeneSymbol %in% drug_targets_j) %>% 
    # dplyr version: dplyr::select(STRING_id) %>% 
    # unlist %>% 
    as.character
  drug_targets_j_dist = degree_all_nodes[names(degree_all_nodes) %in% drug_targets_j_ensembl]
  
  # Create random reference distribution.
  # Each reference distribution corresponds to 1 drug x disease module (and each patient sample/cell line has 1 and only 1 disease module), so we will calculate it here. 
  d_X_Y_rand_dist = mclapply(1:d_X_Y_rand_dist_size, generateRandomNetworkProximity, 
                             disease_module_dist = disease_module_dist, 
                             drug_targets_j_dist = drug_targets_j_dist, 
                             disease_module_ensembl = disease_module_ensembl, 
                             drug_targets_j_ensembl = drug_targets_j_ensembl, 
                             ref_network = ref_network, 
                             degree_all_nodes = degree_all_nodes, 
                             order_ref_network = order_ref_network, 
                             order_disease_module = order_disease_module, 
                             mc.cores = num_cores) %>% unlist # !!!!!!! Change number of cores if memory is an issue.
  # Replace all Inf values with the diameter. https://igraph.org/r/doc/diameter.html
  # ref_network_diameter = diameter(ref_network, unconnected = F) Now pre-calculated and passed as an argument to calcZScore.
  d_X_Y_rand_dist = d_X_Y_rand_dist %>% replaceNonnumeric(replacement = ref_network_diameter)
  
  # Calculate z-scores.
  d = distances(ref_network, v = V(ref_network)[name %in% drug_targets_j_ensembl], to = V(ref_network)[name %in% disease_module_ensembl]) %>% 
    replaceNonnumeric(replacement = ref_network_diameter) %>% 
    apply(2, min, na.rm = T) %>% 
    sum %>% 
    prod((1/order_disease_module))
  mu = mean(d_X_Y_rand_dist, na.rm = T)
  sigma = sd(d_X_Y_rand_dist, na.rm = T)
  z_score = (d - mu) / sigma
  
  # Calculate weights for z-scores. 
  weights = lapply(drug_targets_j, FUN = calcIndividualTargetWeight, 
                   #mc.cores = num_cores,
                   ref_network = ref_network,
                   disease_map = disease_map, 
                   drug_target_data_j = drug_target_data[[drug]])
  weights = do.call(rbind, weights)
  weight_sums = apply(weights, 2, sum, na.rm = T)
  
  z_scores = c(z_score, hadamard.prod(rep(z_score,length(weight_sums)), weight_sums))
  return(z_scores) # z_scores[1] = unweighted z-score; z_scores[2] = weighted, normalized z-score; z_scores[3] = weighted, unnormalized z-score.
  
}

calcIndividualTargetWeight = function(target, ref_network, disease_map, drug_target_data_j) {
  # First, determine whether the target is even represented in the network. 
  # Convert the target to its corresponding Ensembl protein (peptide) ID.
  target_ensembl = disease_map[disease_map$GeneSymbol==target,"STRING_id"] %>% as.character
  target_ensembl = ifelse(length(target_ensembl) < 1, NA, target_ensembl)
  if(!(target_ensembl %in% V(ref_network)$name)) return(c(NA, NA))
  
  # If the target is represented in the network:
  # For each target, determine 
  # 1) how the drug modulates the target,
  # 2) how the target is dysregulated between tumor and normal tissue,
  # 3) the extent of target dysregulation, measured by log fold-change,
  # 4) and the distance from the target to the disease module. 
    
  # 1) Get the direction the drug affects the target. 
  drug_effect = sum(drug_target_data_j[drug_target_data_j$Target==target,"Effect"]) # See getTargets.R. Sum because there are sometimes conflicting directions (i.e., it will have both -1 and 1, maybe because different studies found different things.) This is just the fastest way to resolve the issue ... cutting the Gordian knot, as it were. 
    
  # 2) Get the direction of dysregulation of the target between normal and tumor tissue. 
  # We will need the list of DEGs as calculated by limma (degs.)
  # Get the rank change of the target.
  target_rank_change = ifelse(length(disease_map[disease_map$GeneSymbol==target,]$RankChange)==0, 0, disease_map[disease_map$GeneSymbol==target,]$RankChange) 
    
  # 3) Get the extent of target dysregulation, measured by rank change.
  # We already got it above ... target_rank_change.
    
  # 4) Distance from the target to the CLOSEST disease gene (see https://www.nature.com/articles/ncomms10331.)
  # "The z-score is obtained by comparing the observed distance to a reference distance distribution between a randomly selected group of proteins of matching size and degree distribution as the disease proteins and drug targets in the human interactome." (https://www.nature.com/articles/s41467-019-09186-x#Sec7)
  # Remember, the "from" vertices are the rows, and the "to" vertices are the columns in the results from the igraph::distances() function. 
  # Since we want the shortest distance to each disease node, we'll take the minimum of each column. 
  
    
  # All right, now we can put it all together.
  # The final score for this target will be 0 if the rank change is 0, and -1 * target_rank_change * drug_effect otherwise. 
  # (07/05/2020) Initially, the final score was -1 * target_rank_change * drug_effect. However, an increase in the rank of a drug going from normal to tumor would show up as a negative number, because lower numbers indicate higher ranks. For example, if gene X were ranked 500 in normal tissue and 20 in tumor tissue, the rank would have increased, but 20 - 500 = -480. 
  # Therefore, if the drug upregulates the target, then drug_effect will be positive. If the target is upregulated going from normal to tumor, target_rank_change will be negative. If we want the score to be negative if drug regulation and tumor dysregulation move in the same direction, all we have to do is multiply drug_effect by (weighted) target_rank_change ... positive x negative = negative. Do not multiply by that extra -1.
  # (05/16/2021) But now, we're doing z-scores, and we want the distance to be smaller (i.e. z-score to be more negative) ... so we're back to multiplying by -1 again!! Hurrah! 
  # (06/24/2021) Actually ... the z-scores are already negative, so to keep them negative, the weights should be positive if the drug effect is as desired. ... 
  if(target_rank_change != 0) {
    # The drug regulates the gene in the opposite direction. Good. The score will be positive.
    weight_normalized = (target_rank_change/max(abs(disease_map$RankChange))) * drug_effect # Weight the target rank change.
    weight_unnormalized = (target_rank_change/abs(target_rank_change)) * drug_effect
  } else {
    # The gene is not dysregulated in tumor vs normal. No way to tell what the drug "should" be doing, so we will neither increase nor decrease the score.
    weight_normalized = NA
    weight_unnormalized = NA
  }
    
  return(c(weight_normalized, weight_unnormalized))
}