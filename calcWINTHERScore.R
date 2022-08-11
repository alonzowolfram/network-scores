# WINTHER score, per Rodon et al (Nat Med, 2019.)
# This function calculates the WINTHER score for a given drug. 
calcWINTHERscore = function(drug, drug_targets, ref_network, ref_network_map, disease_map, intensity) {
  # Get the targets of the drug.
  # If there are no drug targets, return NA for the WINTHER score. 
  if(is.null(drug_targets[[drug]]) | nrow(drug_targets[[drug]]) < 1) return(NA)

  # Set up a vector to hold individual target scores.
  drug_j_targets = drug_targets[[drug]]$Target %>% as.character
  Fxi = sapply(drug_j_targets, FUN = calcIndividualFxgi, # Originally jubilee.mcsapply()
                         #mc.cores = num_cores,
                         disease_map = disease_map, 
                         drug_targets_j = drug_targets[[drug]], 
                         ref_network = ref_network, 
                         ref_network_map = ref_network_map
                         )
  names(Fxi) = drug_j_targets
  # WINTHER score = 100 * sum(F^x_gi * I_g/|M_i|),
  # where F^x_g,i is as given above, I_g is the intensity of gene g *in set M_i*, and |M_i| is the size of the set of genes for which F_g,i != 0. 
  M_i = Fxi %>% .[. != 0 & is.finite(.)] %>% unlist() # Edit 08/06/2021: This will have some NAs, because not all of the targets are represented on the network, especially for data sets like CTD and TRRUST with more drug targets. So we'll have to remove those. Also interesting to note: those length issues where the length of one vector doesn't match the length of another stem from the fact that when Fxi gets subsetted, the names of the entries with NA also become NA themselves. So if, e.g., gene BEX5 in Fxi has a value of NA, the name changes from BEX5 to NA in M_i. I don't know why. T_T But let's see if that actually affects the hadamard product. 
  size_M_i = length(M_i)
  I_M_i = intensity %>% .[rownames(.) %in% names(M_i),"Intensity"]  # Vector of all I_g's. Edit 08/06/2021: A check reveals that the intensities in I_M_i are indeed in the same order as M_i. Whew. 
  S_i = sum(hadamard.prod(M_i, I_M_i)/size_M_i, na.rm=T) / 1000 # Originally multiplied by 100, but that was per Rodon et al, where their intensities were much smaller. 
  # Edit 02/09/2021: originally sum(hadamard.prod(Fxi, I)/size_M_i), but it should be sum(hadamard.prod(M_i, I)/size_M_i), right? That's why we were getting those "the length of this vector doesn't match the length of that vector" warnings, right?
  # Edit 08/06/2021: the length of the vectors DOES change the value of the hadamard product. Of course. T_T But now that they should match, everything should be hunky dory. 
    
  return(S_i)
}

# Function for calculating the concordance-adjusted fold change (although here it won't be fold change but change in rank.)
# This calculates for a single gene g for a single drug x in individual i. (There is no argument for which drug, because this function is called only by calcWINTHERscore, which already is called only for a single drug.)
calcIndividualFxgi = function(target, drug_targets_j, ref_network, ref_network_map, disease_map) {
  # First, determine whether the target is even represented in the network. 
  # Convert the target to its corresponding Ensembl protein (peptide) ID.
  target_ensembl = ref_network_map[ref_network_map$GeneSymbol==target,"STRING_id"]
  target_ensembl = ifelse(length(target_ensembl) < 1, NA, target_ensembl)
  if(!(target_ensembl %in% V(ref_network)$name)) return(NA)
    
  # For each target, determine 
  # 1) how the drug modulates the target,
  # 2) how the target is dysregulated between tumor and normal tissue,
  # and 3) the extent of target dysregulation, measured by log fold-change.
    
  # 1) Get the direction the drug affects the target. 
  drug_effect = sum(drug_targets_j[drug_targets_j$Target==target,"Effect"], na.rm=T) # See getTargets.R. Sum because there are sometimes conflicting directions (i.e., it will have both -1 and 1, maybe because different studies found different things.) This is just the fastest way to resolve the issue ... cutting the Gordian knot, as it were. 
    
  # 2) Get the direction of dysregulation of the target between normal and tumor tissue. 
  # We will need the list of DEGs as calculated by limma (degs.)
  # Get the rank change of the target.
  target_rank_change = ifelse(length(disease_map[disease_map$GeneSymbol==target,]$RankChange)==0, 0, disease_map[disease_map$GeneSymbol==target,]$RankChange) 
    
  # 3) Get the extent of target dysregulation, measured by rank change.
  # We already got it above ... target_rank_change.
    
  # All right, now we can put it all together.
  # The final score for this target will be 0 if the rank change is 0, and -1 * target_rank_change * drug_effect otherwise. 
  # (07/05/2020) Initially, the final score was -1 * target_rank_change * drug_effect. However, an increase in the rank of a drug going from normal to tumor would show up as a negative number, because lower numbers indicate higher ranks. For example, if gene X were ranked 500 in normal tissue and 20 in tumor tissue, the rank would have increased, but 20 - 500 = -480. 
  # Therefore, if the drug upregulates the target, then drug_effect will be positive. If the target is upregulated going from normal to tumor, target_rank_change will be negative. If we want the score to be negative if drug regulation and tumor dysregulation move in the same direction, all we have to do is multiply drug_effect by (weighted) target_rank_change ... positive x negative = negative. Do not multiply by that extra -1.
    
  if(target_rank_change != 0) {
    # Direction-adjusted change (Fxgi) = f_g,i * d_g,x,
    # where f_g,i = log fold-change = target_rank_change, and d_gx = direction of drug effect = drug_effect. 
    # Fxgi = ifelse(f_gi * d_gx < 0, 0, f_gi * d_gx)
    #Fxgi = ifelse(target_rank_change * drug_effect < 0, 0, target_rank_change * drug_effect)
    Fxgi = target_rank_change * drug_effect
  } else {
    # The gene is not dysregulated in tumor vs normal. No way to tell what the drug "should" be doing, so we will neither increase nor decrease the score.
    Fxgi = NA
  }
  
  return(Fxgi)
  
}