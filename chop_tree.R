require(phytools)
require(phangorn)

chop_tree <- function(tree, scale=T) {
  #### Cut up phylogenetic tree into segments #################
  times <- node.depth.edgelength(tree) # date of each node
  
  # Line up date of each node with the split points in the tree
  split_points <- sort(unique(times))
  node_time <- match(times, split_points)
  
  # Create a sequence of nodes, respecting temporal order
  node_seq <- seq(from = 1, to = length(node_time))
  node_seq <- node_seq[order(node_time)]
  
  # Find the "parent" node for each node and tip of the tree
  parent <- Ancestors(tree, node_seq, type = "parent")
  
  # Parent time indicates amount of time since the parent node, scaled by the total depth of the tree
  parent_time <- rep(NA, length(node_seq))
  #parent_time[1] <- -99 # placeholder for ancestral state
  
  for (i in 2:length(parent_time)) {
    if (scale == T) parent_time[i] <- (node.depth.edgelength(tree)[node_seq[i]] - node.depth.edgelength(tree)[parent[i]]) / max(node.depth.edgelength(tree))
    
    if (scale == F)parent_time[i] <- (node.depth.edgelength(tree)[node_seq[i]] - node.depth.edgelength(tree)[parent[i]])
  }
  
  N_seg <- length(node_seq) # total num segments in the tree
  
  return(
    list(
      N_seg = N_seg,
      node_seq = node_seq,
      parent = parent,
      ts = parent_time
      )
  )
}

