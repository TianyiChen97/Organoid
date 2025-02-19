detect_communities <- function(graph) {
  
  # Convert the adjacency matrix to an igraph object
  if (class(graph) != "igraph") {
    graph <- graph_from_adjacency_matrix(A, weighted = FALSE, mode = "undirected", diag = FALSE)
  }
  
  # Apply different community detection algorithms
  # 0. Leiden Method
  leiden_communities <- cluster_leiden(graph)
  num_leiden <- length(sizes(leiden_communities))
  
  # 1. Louvain Method
  louvain_communities <- cluster_louvain(graph)
  num_louvain <- length(sizes(louvain_communities))
  
  # # 2. Edge Betweenness Method
  # NB: This fails for a trivial graph!
  # edge_betweenness_communities <- cluster_edge_betweenness(graph)
  # num_edge_betweenness <- length(sizes(edge_betweenness_communities))
  
  # 3. Infomap Method
  infomap_communities <- cluster_infomap(graph)
  num_infomap <- length(sizes(infomap_communities))
  
  # 4. Walktrap Method
  walktrap_communities <- cluster_walktrap(graph)
  num_walktrap <- length(sizes(walktrap_communities))
  
  # 5. Label Propagation Method
  label_propagation_communities <- cluster_label_prop(graph)
  num_label_propagation <- length(sizes(label_propagation_communities))
  
  # 6. Leading Eigenvector Method
  leading_eigenvector_communities <- cluster_leading_eigen(graph)
  num_leading_eigenvector <- length(sizes(leading_eigenvector_communities))
  
  # Return the number of communities detected by each method as a named vector
  return(c(Leiden = num_leiden,
           Louvain = num_louvain,
           #Edge_Betweenness = num_edge_betweenness,
           Infomap = num_infomap,
           Walktrap = num_walktrap,
           Label_Propagation = num_label_propagation,
           Leading_Eigenvector = num_leading_eigenvector))
}

