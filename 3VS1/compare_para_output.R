edge_list_1 = read.csv("~/Desktop/Research /PhDresearch/Hopkins_Organoid/M07357/adjacency_edges_Aug_29_2025_M07357_try2.csv")
edge_list_2 = read.csv("~/Desktop/Research /PhDresearch/Hopkins_Organoid/M07357/adjacency_edges_Sep_1_2025_M07357_para2.csv")

filenames = unique(edge_list_2$File)

for (file in filenames) {

  # Define the list of wells to iterate over
  well_list <- paste0("well", sprintf("%03d", 0:5))
  
  # --- Main Loop ---
  for (well in well_list) {
    
    # 1. Subset the data for the current file and well
    subset_data <- subset(edge_list_2, File == file & Well == well)

    n_rows <- n_cols <- unique(subset_data$dim)
    
    if (nrow(subset_data) == 0) {
      next # Move to the next well
    }
    
    adj <- sparseMatrix(
        i    = subset_data$Row + 1,
        j    = subset_data$Column + 1,
        dims = c(n_rows, n_cols)
      )
      
      # Create the graph and find the largest connected component (LCC)
      g <- graph_from_adjacency_matrix(adj, mode = "undirected", diag = FALSE)
      comps <- components(g, mode = "weak")
      big <- which.max(comps$csize)
      g_lcc <- induced_subgraph(g, V(g)[comps$membership == big])
      print(paste(extract_network_path(file),well,vcount(g_lcc)))

  }
}
