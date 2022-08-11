calcDrugDrugChemIxns = function(i, drug_perms, map, network, sample) {
  # calcDrugDrugChemIxns() is a function that takes a pair of drugs, their respective targets, and a disease protein-protein network and returns a score based on whether either drug affects the action of the other drug.
  # calcDrugDrugChemIxns(character, character, data.frame, iGraph) -> numeric
  
  # For now, we won't consider interactions of 3+ drugs at a time, such as this one: "10,10-bis(4-pyridinylmethyl)-9(10H)-anthracenone inhibits the reaction [Valproic Acid inhibits the reaction [Kainic Acid results in increased expression of FOS protein]]"
  
  # The steps:
  # 1) From the CTD data, get all the interactions of drug A that are of the form "drug A [and not an analog of it] inhibits/promotes the reaction [drug B results in increased/decreased expression/activity of X protein/mRNA]."
  # 2) Determine whether the reaction "drug B results in increased/decreased expression/activity of X protein/mRNA" is a good or bad one.
  # 3) If the reaction is a good one, then the score will be positive if drug A promotes it and negative if drug A inhibits it. If the reaction is a bad one, then the score will be positive if drug A inhibits it and negative if drug A promotes it. 
  
  drugs = drug_perms[i,]
  drug1 = drugs[1]
  drug2 = drugs[2]
  
  # Load the chemical-gene interactions from CTD (saved as a CSV file in /Data/CTD/ALMANAC.)
  possibleError = tryCatch(
    {
      # ctd_all_drugs should already exist in the global environment. 
      # Keep only the interactions for the drugs "drug1" and "drug2" and found in humans.
      regex1 = paste0("^", drug1, "$")
      regex2 = paste0("^", drug2, "$")
      ctd_subset = ctd_all_drugs[ctd_all_drugs$OrganismID=="9606" & (ctd_all_drugs$ChemicalName %in% c(drug1, drug2) | ctd_all_drugs$Input %in% c(drug1, drug2) | base::grepl(regex1, ctd_all_drugs$Input, ignore.case = T) | base::grepl(regex2, ctd_all_drugs$Input, ignore.case = T)),,drop=F]
    }, error=function(cond) {
      message(paste0("For calcDrugDrugChemIxns(), there was a problem retrieving the chemical-gene interactions for the drugs ", drug1, " and ", drug2, ":"))
      message(cond)
      cond
    }
  )
  
  if(inherits(possibleError, "error") | !exists("ctd_subset") | nrow(ctd_subset) < 1) {
    print(paste0("For calcDrugDrugChemIxns(): skipping the drug combination consisting of ", drug1, " and ", drug2, "."))
  
  # Write the error to the appropriate file. 
   sample_error_file = paste0(results_dir, "NetworkScore/Logs/", sample, "_log.txt")
  fileConn = file(sample_error_file)
  writeLines(paste0("For calcDrugDrugChemIxns(), there was a problem calculating the drug-drug interactions for the drug combination consisting of ", paste(drugs, sep="", collapse=" "), "."), fileConn)
  close(fileConn)
  
  return(0)
  } else {
    # Make sure drug1 and drug2 are the exact names used in CTD (get them from ChemicalName column.)
    drug1 = ctd_subset[ctd_subset$OrganismID=="9606" & (ctd_subset$ChemicalName==drug1| ctd_subset$Input==drug1 | base::grepl(regex1, ctd_subset$Input, ignore.case = T) | base::grepl(regex1, ctd_subset$ChemicalName, ignore.case = T)),"ChemicalName"]
    if(!is.atomic(drug1)) {
      drug1 = drug1[complete.cases(drug1),][[1]][1] #drug1[complete.cases(drug1)][1]
    } else {
      drug1 = drug1[complete.cases(drug1)][1]
    }
    
    drug2 = ctd_subset[ctd_subset$OrganismID=="9606" & (ctd_subset$ChemicalName==drug2| ctd_subset$Input==drug2 | base::grepl(regex2, ctd_subset$Input, ignore.case = T)| base::grepl(regex2, ctd_subset$ChemicalName, ignore.case = T)),"ChemicalName"]
    if(!is.atomic(drug2)) {
      drug2 = drug2[complete.cases(drug2),][[1]][1]
    } else {
      drug2 = drug2[complete.cases(drug2)][1]
    }
    
    # Create the regular expressions to be used to filter the ctd_subset.
    regex1 = paste0("^", drug1, " (?!\\banalog\\b)(promotes|inhibits) the reaction") # The regular expression that identifies all the drug-drug interactions. It's important that this one is first!!! (Think why ...)
    # regex = paste0("(", drug1, "|", drug2, ")", " (?!\\banalog\\b)(promotes|inhibits) the reaction") -- we don't have to do this because we enumerate all the permutations, not combinations.
    
    regex2 = paste0(drug2, " (?!\\banalog\\b)") # This will be set to negative, because we do NOT want the analogs. Just the drug itself. 
    regex3 = "(susceptibility to|co-treated|modified form)" # This will be set to negative, because we do NOT want co-treatments (in this script at least--see calcCotreatments.R) or modified forms of proteins. 
    regex4 = "results in (de|in)creased (activity|expression) of" # This will be set to positive, because we want to make sure that this phrase is present. 
    regex5 = paste0("\\[", drug2, " (promotes|inhibits)") # This will be set to negative, because we do NOT want something like the situation "Dactinomycin inhibits the reaction [Tretinoin promotes the reaction [IFNG protein results in increased expression of IRF1 protein]]" which is one of the results that comes up when drug1 is dactinomycin and drug2 is tretinoin.
    # We will use the regexPipes library to do this. 
    interactions = ctd_subset$Interaction %>% regexPipes::grep(regex1, perl = T, value = T) %>% regexPipes::grep(regex2, perl = T, value = T) %>% regexPipes::grep(regex3, value = T, invert = T) %>% regexPipes::grep(regex4, value = T) %>% regexPipes::grep(regex5, value = T, invert = T)
    
    # interactions_table = data.frame(Target=character(), InhibitionOrPromotion=character()) # To be filled with results of 2) and 3).
    # I think this is already accomplished in the overall_reaction_scores table.
    
    if(length(interactions)==0) {
      # No interactions between drug1 and drug2. 
      return(0)
    } else {
      # 2) Now take all those reactions generated in step 1) and determine whether they are good or bad.
      # Because of the possibility of duplicates, we will set this up as a data frame and then remove duplicate rows.
      reactions = mclapply(interactions, FUN = determineRxnGoodOrBad, map = map, mc.cores = num_cores)
      reactions = as.data.frame(do.call(rbind, reactions))
      colnames(reactions) = c("Target", "DrugDirection", "TargetRankChange")  
      reactions$DrugDirection = as.integer(as.character(reactions$DrugDirection))
      reactions$TargetRankChange = as.integer(as.character(reactions$TargetRankChange))
      
      # 3) This step is kind of in parallel with step 2). Determine from the reactions generated in step 1) whether drug A inhibits or promotes it. 
      drugA_effects = jubilee.mcsapply(interactions, FUN = determineDrugARole, mc.cores = num_cores)
      
      # Add drugA_effects to the reactions data frame _and THEN_ remove the duplicate rows from the reactions data frame. (If we removed the duplicate rows first, the number of rows in the reactions data frame would not match the length of the drugA_effects vector.)
      reactions$DrugAEffects = drugA_effects
      reactions = reactions %>% distinct()
      # Remove targets not on the network. 
      reactions = reactions[reactions$Target %in% map$GeneSymbol,,drop=F]
      
      if(nrow(reactions)==0) {
        return(0)
      } else {
        # Now calculate the overall reaction scores, i.e. the final assessment of the reaction in square brackets.
        reaction_scores = jubilee.mcsapply(1:nrow(reactions), FUN = calcOverallRxnScore, reactions = reactions, map = map, mc.cores = num_cores)
      
        # 4) Cbind the reaction_scores and drugA_effects vectors into an overall_reaction_scores data frame. 
        overall_reaction_scores = data.frame(RxnScore=reaction_scores, DrugAEffects=reactions$DrugAEffects) # Do NOT use DrugAEffects = drugA_effects, because if there are duplicate rows in reactions, then length(drugA_effects) != length(reactions$DrugAEffects).
      
        # 5) The final score will be sum(overall_reaction_scores[i,1] * overall_reaction_scores[i,2]).
        return(sum(hadamard.prod(overall_reaction_scores[,1], overall_reaction_scores[,2])))
      } # End else there are targets left on the network. 
    } # End else there are interactions between the drugs.
  
  } # End else loading the CTD chem-gene ixns file was successful.
  
}

determineRxnGoodOrBad = function(interaction, map) {
  # Activity takes priority over expression, so check if there's the phrase "increased/decreased activity of".
  if(base::grepl("(in|de)creased activity of", interaction)) {
    # There is activity. The direction will be the direction of the change in activity. 
    direction = ifelse(base::grepl("decreased activity of", interaction), -1, 1)
    
    # Extract the target. 
    # Split the string by spaces. The last two words should be the name of the gene/protein and "mRNA/protein".
    split_string = strsplit(interaction, " ")[[1]]
    target = ifelse(base::grepl("(mRNA|protein)",split_string[length(split_string)]), split_string[length(split_string) - 1], split_string[length(split_string)])
    
  } else {
    # No activity. The direction will be the direction of the change in expression.
    direction = ifelse(base::grepl("decreased expression of", interaction), -1, 1)
    
    # Extract the target. 
    # Split the string by spaces. The last two words should be the name of the gene/protein and "mRNA/protein".
    split_string = strsplit(interaction, " ")[[1]]
    target = ifelse(base::grepl("(mRNA|protein)",split_string[length(split_string)]), split_string[length(split_string) - 1], split_string[length(split_string)])
    
  }
  
  # The drug direction is the effect that _drug B_ has on the target, not drug A. 
  drug_direction = ifelse(base::grepl("decreased", interaction), -1, 1)
  
  # Compare that to how the target is dysregulated in tumor vs normal tissue. 
  # We will need the list of DEGs in the "map" variable.
  # Get the rank change of the target.
  target_rank_change = ifelse(length(map[map$GeneSymbol==target,]$RankChange)==0, 0, map[map$GeneSymbol==target,]$RankChange) 
  
  return(c(target, drug_direction, target_rank_change))
}

determineDrugARole = function(interaction) {
  drugA_effect = ifelse(base::grepl("promotes", interaction), 1, -1)
  return(drugA_effect)
}

calcOverallRxnScore = function(i, reactions, map) {
  # If the reaction is good, its assessment score should be positive.
  # If the reaction is bad, its assessment score should be negative.
  # (07/05/2020) Initially, the final score was -1 * target_rank_change * drug_effect. However, an increase in the rank of a drug going from normal to tumor would show up as a negative number, because lower numbers indicate higher ranks. For example, if gene X were ranked 500 in normal tissue and 20 in tumor tissue, the rank would have increased, but 20 - 500 = -480. 
  # Therefore, if the drug upregulates the target, then drug_effect will be positive. If the target is upregulated going from normal to tumor, target_rank_change will be negative. If we want the score to be negative if drug regulation and tumor dysregulation move in the same direction, all we have to do is multiply drug_effect by (weighted) target_rank_change ... positive x negative = negative. Do not multiply by that extra -1.
  drug_direction = reactions[i,]$DrugDirection
  target_rank_change = reactions[i,]$TargetRankChange
  affected_reaction_score = drug_direction * (target_rank_change/max(abs(map$RankChange))) # Originally -1 * drug_direction * (target_rank_change/max(abs(map$RankChange)))
  
  return(affected_reaction_score)
}