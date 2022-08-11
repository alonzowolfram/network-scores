calcCentrality = function(drug, drug_targets, ref_network, ref_network_map, disease_map, syngenet = F) {
  # calcCentrality() takes a SINGLE drug and a reference network (the full PPI network, not a disease network!) and returns the degree centrality measure of that drug. 
  
  # Get the targets of the drug.
  # If there are no drug targets, return NA for the centrality scores. 
  if(is.null(drug_targets[[drug]]) | nrow(drug_targets[[drug]]) < 1) return(c(NA, NA, NA))
  
  # Otherwise, set up a vector to hold individual target scores and then calculate centrality for each target.
  drug_j_targets = drug_targets[[drug]]$Target %>% as.character
  
  centrality_scores = lapply(drug_j_targets, FUN = calcIndividualTargetCentrality, 
                             #mc.cores = num_cores,
                             ref_network = ref_network, 
                             ref_network_map = ref_network_map, 
                             disease_map = disease_map, 
                             drug_targets_j = drug_targets[[drug]], 
                             syngenet = syngenet
                             )
  centrality_scores = do.call(rbind, centrality_scores)
  
  centrality_scores_sums = apply(centrality_scores, 2, sum, na.rm = T)
  centrality_score_unweighted = ifelse(is.finite(centrality_scores_sums[1]), centrality_scores_sums[1], 0)
  centrality_score_weighted_normalized = ifelse(is.finite(centrality_scores_sums[2]), centrality_scores_sums[2], 0)
  centrality_score_weighted_unnormalized = ifelse(is.finite(centrality_scores_sums[3]), centrality_scores_sums[3], 0)
  
  return(c(centrality_score_unweighted, centrality_score_weighted_normalized, centrality_score_weighted_unnormalized))
}

calcIndividualTargetCentrality = function(target, drug_targets_j, ref_network, ref_network_map, disease_map, syngenet) {
  # calcIndividualTargetCentrality() takes a target and a disease network and returns the degree centrality measure of that target. 
  # Make sure the target is even represented in the network and has a corresponding Ensembl ID.
  # If not, the score should be neither positive nor negative.
  if(nrow(ref_network_map[ref_network_map$GeneSymbol==target,])==0) { # | is.na(ref_network_map[ref_network_map$GeneSymbol==target,"STRING_id"])
    centrality_score_raw = NA
    centrality_score_weighted = NA
    centrality_score_unweighted = NA
    return(c(centrality_score_raw, centrality_score_weighted, centrality_score_unweighted))
  } 
  
  # Convert the target to its corresponding Ensembl protein (peptide) ID.
  target_ensembl = ref_network_map[ref_network_map$GeneSymbol==target,"STRING_id"]
  target_ensembl = ifelse(length(target_ensembl) < 1, NA, target_ensembl)
  
  # Make sure this target is present in the network (why are we having to do this again???)
  if(!(target_ensembl %in% V(ref_network)$name)) {
    centrality_score_raw = NA
    centrality_score_weighted = NA
    centrality_score_unweighted = NA
    return(c(centrality_score_raw, centrality_score_weighted, centrality_score_unweighted))
  }
      
  if(!syngenet) {
    centrality_score_raw = igraph::degree(ref_network, v = target_ensembl, normalized = T)
  } else {
    # https://kateto.net/networks-r-igraph
    # Calculate the betweenness, closeness, and page-rank centrality of target "target_ensembl."
    # Should we estimate these metrics (i.e. use a path-length cutoff)? If we don't, the algorithm would consider all nodes in the network, and that would take forever. However, there is no estimate_ version for page_rank().
    # https://stackoverflow.com/questions/41753929/how-long-does-it-take-for-igraph-r-to-estimate-network-centrality-measures-for-a
    bc = igraph::betweenness(ref_network, v = target_ensembl, directed = FALSE, normalized = TRUE) # If we are going to compare network scores across different studies (e.g. a BRCA data set and a CRC data set, we should probably normalize.)
    cc = igraph::closeness(ref_network, vids = target_ensembl, normalized = TRUE)
    prc = igraph::page_rank(ref_network, vids = target_ensembl, directed = FALSE)$vector[1]
    # If betweenness centrality doesn't work, we can try degree(). 
    centrality_score_raw = mean(c(bc, cc, prc), na.rm = T) # prc
  }
    
  # Determine whether the direction that the drug affects the target and the direction the target is dysregulated in tumor vs normal tissue are the same.
  # If they are the same, the score should be NEGATIVE and weighted by the absolute hub score as calculated above.
  # If they are opposite, the score should be POSITIVE and weighted by the absolute hub score as calculated above.
  # We will need the list of DEGs (map.)
  # Get the rank change of the target.
  target_rank_change = ifelse(length(disease_map[disease_map$GeneSymbol==target,]$RankChange)==0, 0, disease_map[disease_map$GeneSymbol==target,]$RankChange) 
    
  # Get the direction the drug affects the target. 
  drug_effect = sum(drug_targets_j[drug_targets_j$Target==target,"Effect"], na.rm=T) # See getTargets.R. 

  # Calculate the final centrality score based on the criteria above. 
  # The final score for this target will be 0 if the rank change is 0, and -1 * target_rank_change * drug_effect otherwise. 
  # (07/05/2020) Initially, the final score was -1 * target_rank_change * drug_effect. However, an increase in the rank of a drug going from normal to tumor would show up as a negative number, because lower numbers indicate higher ranks. For example, if gene X were ranked 500 in normal tissue and 20 in tumor tissue, the rank would have increased, but 20 - 500 = -480. 
  # Therefore, if the drug upregulates the target, then drug_effect will be positive. If the target is upregulated going from normal to tumor, target_rank_change will be negative. If we want the score to be negative if drug regulation and tumor dysregulation move in the same direction, all we have to do is multiply drug_effect by (weighted) target_rank_change ... positive x negative = negative. Do not multiply by that extra -1.
  if(target_rank_change != 0) {
    # The drug regulates the gene in the opposite direction. Good. The score will be positive.
    centrality_score_weighted = (target_rank_change/max(abs(disease_map$RankChange))) * drug_effect * centrality_score_raw # Weight the target rank change. Initially -1 * (target_rank_change/max(abs(map$RankChange))) * drug_effect
    centrality_score_unweighted = (target_rank_change/abs(target_rank_change)) * drug_effect * centrality_score_raw # Initially -1 * (target_rank_change/abs(target_rank_change)) * drug_effect. 
  } else {
    # The gene is not dysregulated in tumor vs normal. No way to tell what the drug "should" be doing, so we will neither increase nor decrease the score.
    centrality_score_weighted = NA
    centrality_score_unweighted = NA
  }
  
  return(c(centrality_score_raw, centrality_score_weighted, centrality_score_unweighted))
  
}

calcBridgingCentrality = function(target, drug_targets, ref_network, ref_network_map, disease_map) {
  # calcBridgingCentrality() takes a target and a disease network and returns the bridging centrality measure of that target.
  # Make sure the target is even represented in the network and has a corresponding Ensembl ID. 
  if(nrow(ref_network_map[ref_network_map$GeneSymbol==target,])==0) { # | is.na(ref_network_map[ref_network_map$GeneSymbol==target,"STRING_id"])
    # If not, the score should be neither positive nor negative.
    bridging_centrality_score_weighted = 0
    bridging_centrality_score_unweighted = 0
  } else {
    # Convert the target to its corresponding Ensembl protein (peptide) ID.
    target_ensembl = ref_network_map[ref_network_map$GeneSymbol==target,"STRING_id"]
    
    if(!(target_ensembl %in% V(ref_network)$name)) {
      bridging_centrality_score_weighted = 0
      bridging_centrality_score_unweighted = 0
    } else {
      # Calculate bridging centrality.
      # Bridging centrality = betweenness centrality * bridging coefficient. 
      # Calculate bridging coefficient.
      bc = betweenness(ref_network, v = target_ensembl, directed = FALSE, normalized = TRUE) # If we are going to compare network scores across different studies (e.g. a BRCA data set and a CRC data set, we should probably normalize.)
      degree = degree(ref_network, target_ensembl, normalized = TRUE)
      neighbors = neighbors(ref_network, target_ensembl)
      inv_deg_neighbors = 1/degree(ref_network, neighbors, normalized = TRUE)
      br_coeff = (1/degree) / sum(inv_deg_neighbors) 
      brc = bc * br_coeff
      
      # Determine whether the direction that the drug affects the target and the direction the target is dysregulated in tumor vs normal tissue are the same.
      # If they are the same, the score should be NEGATIVE and weighted by the absolute hub score as calculated above.
      # If they are opposite, the score should be POSITIVE and weighted by the absolute hub score as calculated above.
      # We will need the list of DEGs (map.)
      # Get the rank change of the target.
      target_rank_change = ifelse(length(disease_map[disease_map$GeneSymbol==target,]$RankChange)==0, 0, disease_map[disease_map$GeneSymbol==target,]$RankChange)
      
      # Get the direction the drug affects the target. 
      drug_effect = sum(drug_targets[drug_targets$Target==target,"Direction"]) # See getTargets.R. 
      
      # Calculate the final centrality score based on the criteria above. 
      # The final score for this target will be 0 if the rank change is 0, and -1 * target_rank_change * drug_effect otherwise. 
      # (07/05/2020) Initially, the final score was -1 * target_rank_change * drug_effect. However, an increase in the rank of a drug going from normal to tumor would show up as a negative number, because lower numbers indicate higher ranks. For example, if gene X were ranked 500 in normal tissue and 20 in tumor tissue, the rank would have increased, but 20 - 500 = -480. 
      # Therefore, if the drug upregulates the target, then drug_effect will be positive. If the target is upregulated going from normal to tumor, target_rank_change will be negative. If we want the score to be negative if drug regulation and tumor dysregulation move in the same direction, all we have to do is multiply drug_effect by (weighted) target_rank_change ... positive x negative = negative. Do not multiply by that extra -1.
      if(target_rank_change != 0) {
        # The drug regulates the gene in the opposite direction. Good. The score will be positive.
        bridging_centrality_score_weighted = (target_rank_change/max(abs(disease_map$RankChange))) * drug_effect * brc # Weight the target rank change. Initially -1 * (target_rank_change/max(abs(map$RankChange))) * drug_effect
        bridging_centrality_score_unweighted = (target_rank_change/abs(target_rank_change)) * drug_effect * brc # Initially -1 * (target_rank_change/abs(target_rank_change)) * drug_effect. 
      } else {
        # The gene is not dysregulated in tumor vs normal. No way to tell what the drug "should" be doing, so we will neither increase nor decrease the score.
        bridging_centrality_score_weighted = 0
        bridging_centrality_score_unweighted = 0
      }
    } # End else() target_ensembl is in the network.
    
  } # End else() the target is represented. 
  
  return(c(bridging_centrality_score_weighted, bridging_centrality_score_unweighted))
}

calcCentralityTest = function(target, drug_targets, ref_network, ref_network_map, disease_map, cutoff = 3) {
  # calcCentralityTest() takes a target and a disease network and returns the degree, estimated closeness, and estimated betweenness of that target. 
  
  # Make sure the target is even represented in the network and has a corresponding Ensembl ID. 
  if(nrow(ref_network_map[ref_network_map$GeneSymbol==target,])==0) { # | is.na(ref_network_map[ref_network_map$GeneSymbol==target,"STRING_id"])
    # If not, the score should be neither positive nor negative.
    centrality_score_weighted = 0
    centrality_score_unweighted = 0
  } else {
    # Convert the target to its corresponding Ensembl protein (peptide) ID.
    target_ensembl = ref_network_map[ref_network_map$GeneSymbol==target,"ref_network_map"]
    
    # Make sure this target is present in the network (why are we having to do this again???)
    if(!(target_ensembl %in% V(ref_network)$name)) {
      centrality_score_weighted = 0
      centrality_score_unweighted = 0
    } else {
      # https://kateto.net/networks-r-igraph
      # Benchmark the calculation of the betweenness, closeness, and degree of target "target_ensembl."
      # Betweenness.
      sys_time_bc = system.time({bc = estimate_betweenness(ref_network, vids = target_ensembl, directed = FALSE, normalized = TRUE, cutoff = cutoff)})
      # Closeness.
      sys_time_cc = system.time({cc = estimate_closeness(ref_network, vids = target_ensembl, normalized = TRUE, cutoff = cutoff)})
      # Degree.
      sys_time_degree = system.time({degree = degree(ref_network, v = target_ensembl, normalized = TRUE)})
      
    } # End else() the target is present in the network. 
    
  } # End else() ... the target is present in the network. 
  
  return(c(centrality_score_weighted, centrality_score_unweighted))
}