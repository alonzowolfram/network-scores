calcCotreatments = function(i, drug_perms, map, network, sample) {
  # calcDrugDrugChemIxns() is a function that takes a pair of drugs, their respective targets, and a disease protein-protein network and returns a score based on whether either drug affects the action of the other drug.
  # calcDrugDrugChemIxns(character, character, data.frame, iGraph) -> numeric
  
  # For now, we won't consider interactions of 3+ drugs at a time, such as this one: "10,10-bis(4-pyridinylmethyl)-9(10H)-anthracenone inhibits the reaction [Valproic Acid inhibits the reaction [Kainic Acid results in increased expression of FOS protein]]"
  
  # The steps:
  # 1) From the CTD data, get all the interactions of drug A that are of the form "[drug A [and not an analog of it] co-treated with drug B [and not an analog of it]] results in increased/decreased expression/activity of X protein/mRNA."
  # 2) Determine whether the reaction "[drug A co-treated with drug B] results in increased/decreased expression/activity of X protein/mRNA" is a good or bad one.
  
  drugs = drug_perms[i,]
  drug1 = drugs[1]
  drug2 = drugs[2]
  
  # Load the chemical-gene interactions from CTD (saved as a CSV file in /Data/CTD/ALMANAC.)
  possibleError = tryCatch(
    {
      # ctd_all_drugs should already be loaded in the global environment. 
      # Keep only the interactions for the drugs "drug1" and "drug2" and found in humans.
      regex1 = paste0("^", drug1, "$")
      regex2 = paste0("^", drug2, "$")
      ctd_subset = ctd_all_drugs[ctd_all_drugs$OrganismID=="9606" & (ctd_all_drugs$ChemicalName %in% c(drug1, drug2) | ctd_all_drugs$Input %in% c(drug1, drug2) | base::grepl(regex1, ctd_all_drugs$Input, ignore.case = T) | base::grepl(regex2, ctd_all_drugs$Input, ignore.case = T)),,drop=F]
    }, error=function(cond) {
      message(paste0("For calcCotreatments(), there was a problem retrieving the chemical-gene interactions for the drugs ", drug1, " and ", drug2, ":"))
      message(cond)
      cond
    }
  )
  
  if(inherits(possibleError, "error") | !exists("ctd_subset") | nrow(ctd_subset) < 1) {
    print(paste0("For calcCotreatments(): skipping the drug combination consisting of ", drug1, " and ", drug2, "."))
    
    # Write the error to the appropriate file. 
    sample_error_file = paste0(results_dir, "NetworkScore/Logs/", sample, "_log.txt")
    fileConn = file(sample_error_file)
    writeLines(paste0("For calcCotreatments(), there was a problem calculating the drug-drug interactions for the drug combination consisting of ", paste(drugs, sep="", collapse=" "), "."), fileConn)
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
    regex1 = paste0("^\\[", drug1, " co-treated with ", drug2, "] results in (de|in)creased (activity|expression) of") # Regular expression for co-treatments. This excludes the interactions of the form "[drug1 co-treated with drug2] promotes|inhibits the reaction." Maybe in a future version. 
    
    regex2 = paste0("(", drug2, "|", drug1, ")", " (?!\\banalog\\b)") # This will be set to negative, because we do NOT want the analogs. Just the drug itself. 
    regex3 = "(susceptibility to|modified form)" # This will be set to negative, because we do NOT want modified forms of proteins. 
    regex4 = paste0("\\[", drug2, " (promotes|inhibits)") # This will be set to negative, because we do NOT want something like the situation "Dactinomycin inhibits the reaction [Tretinoin promotes the reaction [IFNG protein results in increased expression of IRF1 protein]]" which is one of the results that comes up when drug1 is dactinomycin and drug2 is tretinoin.
    # We will use the regexPipes library to do this. 
    interactions = ctd_subset$Interaction %>% regexPipes::grep(regex1, perl = T, value = T) %>% regexPipes::grep(regex2, perl = T, value = T) %>% regexPipes::grep(regex3, value = T, invert = T) %>% regexPipes::grep(regex4, value = T, invert = T)
    
    # interactions_table = data.frame(Target=character(), InhibitionOrPromotion=character()) # To be filled with results of 2) and 3).
    # I think this is already accomplished in the overall_reaction_scores table.
    
    if(length(interactions)==0) {
      # No co-treatment information on drug1 and drug2. 
      return(0)
      
    } else {
      
      # 2) Now take all those reactions generated in step 1) and determine whether they are good or bad.
      # Because of the possibility of duplicates, we will set this up as a data frame and then remove duplicate rows.
      reactions = mclapply(interactions, FUN = determineRxnGoodOrBad, map = map, mc.cores = num_cores)
      reactions = as.data.frame(do.call(rbind, reactions))
      colnames(reactions) = c("Target", "DrugDirection", "TargetRankChange")  
      reactions$DrugDirection = as.integer(as.character(reactions$DrugDirection))
      reactions$TargetRankChange = as.integer(as.character(reactions$TargetRankChange))
      
      # Remove the duplicate rows from the reactions data frame.
      reactions = reactions %>% distinct()
      # Remove targets not on the network. 
      reactions = reactions[reactions$Target %in% map$GeneSymbol,,drop=F]
      
      if(nrow(reactions)==0) {
        return(0)
      } else {
        # Now calculate the overall reaction scores, i.e. the final assessment of the reaction.
        reaction_scores = jubilee.mcsapply(1:nrow(reactions), FUN = calcOverallRxnScore, reactions = reactions, map = map, mc.cores = num_cores)
        
        # 4) The final score will be the sum of the reaction_scores.
        return(sum(reaction_scores))
      } # End else there are still targets on the network. 
      
    } # End else there are interactions between the drugs.
    
  } # End else loading the CTD chem-gene ixns file was successful.
  
}

