generateRandomNetworkProximity = function(i, disease_module_dist, drug_targets_j_dist, disease_module_ensembl, drug_targets_j_ensembl, ref_network, degree_all_nodes, order_ref_network, order_disease_module) {
  # This is the function that has to be repeated 10 000 times, so we'd better streamline it. ... 
  
  # Ensure reproducibility.
  set.seed(i)
  
  # Randomly generate a group of proteins with the same size and degree distribution as disease_proteins.
  # https://stats.stackexchange.com/questions/286062/distribution-matching-by-subsampling
  group1 = degree_all_nodes # All proteins.
  group2 = disease_module_dist # Disease proteins. 
  label = factor(c(rep("group1", order_ref_network), rep("group2", order_disease_module)))
  var = c(group1, group2)
  
  var_cut = cut(var, breaks = 100)
  
  # Simple density estimators.
  p_density = table(var_cut[label == "group1"])/length(group1)
  q_density = table(var_cut[label == "group2"])/length(group2)
  
  # Evaluate density of both distributions at each point of group 1 .
  p = p_density[var_cut[label == "group1"]]
  q = q_density[var_cut[label == "group1"]]
  
  # Graphing for sanity check.
  #mydata = data.frame(var,label)
  #mydata$matching_vec = rep(NA, nrow(mydata))
  #mydata$matching_vec[matching_idx] = "group1 (matching)"
  #mydata$matching_vec[label == "group2"] = "group2"
  #mydata$matching_vec = factor(mydata$matching_vec)
  #ggplot(subset(mydata, !is.na(matching_vec))) +  aes(x = var , fill=matching_vec)  + geom_histogram( position="identity", alpha=0.5 , bins = 100 )
  
  # Randomly sample from group 1 (the distributions of all proteins in the reference STRING network.)
  matching_idx = sample(which(label == "group1"), sum(label=="group2"), FALSE, q/p)
  random_disease_proteins = group1[matching_idx] %>% names
  
  # Randomly generate a group of proteins with the same size and degree distribution as drug_targets.
  group3 = drug_targets_j_dist # Drug targets. 
  label = factor(c(rep("group1", order_ref_network), rep("group3", length(drug_targets_j_dist))))
  var = c(group1, group3)
  
  var_cut = cut(var, breaks = 100)
  
  # Simple density estimators.
  p_density = table(var_cut[label == "group1"])/length(group1)
  q_density = table(var_cut[label == "group3"])/length(group3)
  
  # Evaluate density of both distributions at each point of group 1 .
  p = p_density[var_cut[label == "group1"]]
  q = q_density[var_cut[label == "group1"]]
  
  # Randomly sample from group 1 (the distributions of all proteins in the reference STRING network.)
  matching_idx = sample(which(label == "group1"), sum(label=="group3"), FALSE, q/p)
  random_drug_targets = group1[matching_idx] %>% names 
  
  # Trying something out from journal entry 05/16/2021.
  d_x_y = distances(ref_network, v = V(ref_network)[name %in% random_drug_targets], to = V(ref_network)[name %in% random_disease_proteins])
  # OK, looks like the "from" vertices are the rows, and the "to" vertices are the columns. 
  # Since we want the shortest distance to each disease node, we'll take the minimum of each column. 
  min_d_x_y = apply(d_x_y, 2, min, na.rm = T)
  d_X_Y = sum(min_d_x_y) * (1/order_disease_module)
  
  #d_X_Y = jubilee.mcsapply(random_disease_proteins, FUN = calculateShortestPathLength, drug_targets = random_drug_targets, ref_network = ref_network, mc.cores = num_cores) %>% # 
  #  sum(na.rm = T) %>% 
  #  c(., (1/size_disease_module)) %>% 
  #  prod
  
  return(d_X_Y)
}

# NOT USED ANYMORE.
#calculateShortestPathLength = function(disease_protein, drug_targets, ref_network) {
#  return(distances(ref_network, v = V(ref_network)[name==disease_protein], to = V(ref_network)[name %in% drug_targets]) %>% min)
#}