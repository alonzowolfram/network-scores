getTargets = function(drug, disease_map, sample, drug_target_data) {
  # For drug "drug," get the targets from ctd_all_drugs. (In TCGA_data_analysis.Rmd, ctd_all_drugs is named ctd_drug_targets. Oops! We'll leave it as ctd_all_drugs for calcNetworkScore() and all the functions it calls so we don't break anything.)
  # Never mind, let's change it to all_drug_effects, because we'll need ctd_all_drugs for the actual CTD output, which will be used by calcCotreatments and calcDrugDrugChemIxns. 
  # 15 May 2021: never mind, let's change it to drug_target_data. 
  print(paste0("Getting targets for drug ", drug, "."))
  
  possibleError = tryCatch(
    {
  targets = drug_target_data[[drug]]
  # Remove targets not in the disease map generated for the sample.
  # For example, if the disease map was generated for a sample measured using a platform
  # that only measured ~17000 genes, remove the targets not covered in those 17000 genes
  # even if the reference DB (e.g. STRING-DB) has data for 19000 genes. 
  # This is because if it wasn't covered, we can't assign a weight to it. 2022-06-26. 
  targets = targets[targets$Target %in% disease_map$GeneSymbol,,drop=F]
    }, error=function(cond) {
      message(paste0("For getTargets(), there was a problem retrieving the chemical-gene interactions for the drug ", drug, ":"))
      message(cond)
      cond
    }
  )
  
  if(inherits(possibleError, "error") | !exists("targets")) {
    # Write the error to the appropriate file. 
    sample_error_file = paste0(results_dir, "Logs/", sample, "_log.txt")
    fileConn = file(sample_error_file)
    writeLines(paste0("For getTargets(), there was a problem retrieving the chemical-gene interactions for the drug ", drug, "."), fileConn)
    close(fileConn)
    
    # Return an empty data.frame. 
    return(data.frame())
  } else {
    return(targets)
  }

}

# --- OLD FUNCTIONS, PREVIOUSLY CALLED BY calcNetworkScore(). NO LONGER IN USE. ---
getTargetsOld = function(drug, disease_map, sample) {
  # For drug "drug," get the targets. 
  
  # Subset the CTD data frame to include only the rows corresponding to the current drug. 
  # ctd_subset = ctd[ctd$ChemicalName==drug,,drop=F]
  # ABOVE IS THE OLD WAY. WE WILL USE AN HTTP BATCH QUERY FOR THIS ONE.
  # https://www.r-bloggers.com/getting-data-from-an-online-source/
  #query = paste0("http://ctdbase.org/tools/batchQuery.go?inputType=chem&inputTerms=", drug, "&report=cgixns&actionTypes=activity&actionTypes=expression&format=csv")
  #ctd_subset = data.table::fread(query)
  
  print(paste0("Getting targets for drug ", drug, "."))
  
  # Load the chemical-gene interactions from CTD (saved as a CSV file in /Data/CTD/ALMANAC.)
  possibleError = tryCatch(
    {
      # ctd_all_drugs should already be loaded in the global environment. 
      # Keep only the interactions for the drug "drug" and found in humans.
      regex = paste0("^", drug, "$")
      ctd_subset = ctd_all_drugs[ctd_all_drugs$OrganismID=="9606" & (ctd_all_drugs$ChemicalName==drug | ctd_all_drugs$Input==drug | base::grepl(regex, ctd_all_drugs$Input, ignore.case = T) | base::grepl(regex, ctd_all_drugs$ChemicalName, ignore.case = T)),,drop=F]
    }, error=function(cond) {
      message(paste0("For getTargets(), there was a problem retrieving the chemical-gene interactions for the drug ", drug, ":"))
      message(cond)
    }
  )
  
  if(inherits(possibleError, "error") | !exists("ctd_subset")) {
    # Write the error to the appropriate file. 
    sample_error_file = paste0(results_dir, "Logs/", sample, "_log.txt")
    fileConn = file(sample_error_file)
    writeLines(paste0("For getTargets(), there was a problem retrieving the chemical-gene interactions for the drug ", drug, "."), fileConn)
    close(fileConn)
    
    return(NA)
  }
  
  else {
  
  # We want only the interactions that include "[chemical name] (and not an analog) ... results in decreased/increased action/expression of [X protein/mRNA]."
  # regex = paste0(drug, " (?!\\banalog\\b)(binds to and results in |results in(\\s| (de|in)creased [^\\s]+ of and results in ))(de|in)creased (activity|expression) of") # Yes, there should be no space between "(?!\\banalog\\b)" and "results."
    # Standardize the drug name. 
    drug = ctd_subset$ChemicalName
    drug = drug[complete.cases(drug)][1]
  
  regex1 = paste0(" (?!\\banalog\\b)") # This will be set to negative, because we do NOT want the analogs. Just the drug itself. 
  regex2 = "(the reaction \\[|susceptibility to|co-treated|modified form)" # This will be set to negative, because we do NOT want the other drug-chemical interactions that the drug in question affects OR direct gene effects affected by other drugs. Just the direct gene effects. We also don't want co-treatments (for right now.) Finally, we don't want the modified forms of proteins. 
  regex3 = "results in (de|in)creased (activity|expression) of" # This will be set to positive, because we want to make sure that this phrase is present. 
  regex4 = "\\[" # We do NOT want something like "Carboplatin results in increased activity of [NFKB1 protein binds to RELA protein]"
  # We will use the regexPipes library to do this. 
  interactions = ctd_subset$Interaction %>% regexPipes::grep(drug, fixed = T, value = T) %>% regexPipes::grep(regex1, perl = T, value = T) %>% regexPipes::grep(regex2, value = T, invert = T) %>% regexPipes::grep(regex3, value = T) %>% regexPipes::grep(regex4, value = T, invert = T)
  # ctd_subset = ctd_subset[grep(regex, ctd_subset$Interaction, perl = T),,drop=F]
  
  # Build a data frame to hold the targets and the direction the drug affects the target in.
  # We won't distinguish between mRNA and protein here. 
  # Populate the data frame with the targets.
  targets = mclapply(interactions, FUN = extractTarget, drug = drug, mc.cores = num_cores)
  targets = as.data.frame(do.call(rbind, targets))
  if(nrow(targets) > 0) {
    colnames(targets) = c("Drug", "Target", "Effect")
    targets$Effect = as.integer(as.character(targets$Effect))
    
    # Remove duplicates
    targets = unique(targets)
    
    # Remove targets not in the DEGs. 
    targets = targets[targets$Target %in% disease_map$GeneSymbol,,drop=F]
  }
  
  rm(ctd_subset)
  
  return(targets)
  }
}
# --- END OLD FUNCTIONS. ---

extractTarget = function(interaction, drug) {
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
  
  return(c(drug, target, direction))
}