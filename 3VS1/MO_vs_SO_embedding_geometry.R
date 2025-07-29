edge_list = read_csv("~/Desktop/Research /PhDresearch/Hopkins_Organoid/MO VS SO_2025_May_graphs/adjacency_edges_May_31_2025_ecr_results_no_window.csv")
filenames = unique(edge_list$File)
filenames = filenames[-2]

for (file in filenames) {
  par(mfrow = c(2, 3), mar = c(3, 3, 3, 1), oma = c(4, 0, 3, 0))
  
  # Define the list of wells to iterate over
  well_list <- paste0("well", sprintf("%03d", 0:5))
  
  # --- Main Loop ---
  for (well in well_list) {
    
    # 1. Subset the data for the current file and well
    subset_data <- subset(edge_list, File == file & Well == well)
    
    # Initial check for any data at all
    if (nrow(subset_data) == 0) {
      plot(1, type = "n", axes = FALSE, xlab = "", ylab = "")
      title(main = well)
      text(1, 1, "No data available", col = "red")
      next # Move to the next well
    }
    
    # --- Graph Construction and ASE ---
    # This block is wrapped in tryCatch to handle potential errors
    tryCatch({
      # Create the adjacency matrix
      n_rows <- n_cols <- unique(subset_data$dim)
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
      
      # 2. Check if the size of the LCC is less than 11
      if (vcount(g_lcc) < 11) {
        # If so, plot the specified message and skip the analysis
        plot(1, type = "n", axes = FALSE, xlab = "", ylab = "")
        title(main = well)
        text(1, 1, "Size of largest\nconnected component < 11", col = "red", cex = 0.9)
      } else {
        # --- Proceed with ASE only if LCC size is 11 or greater ---
        
        # Perform Adjacency Spectral Embedding (ASE) on the LCC
        ase <- full.ase(g_lcc, 10)
        
        # Find the elbow to determine the embedding dimension
        elb <- getElbows(ase$eval, plot = F)
        
        if (length(elb) < 2) {
          stop("Embedding dimension is less than 2.")
        }
        
        embedding_dim <- 2
        
        # --- Conditional Plotting ---
        
        # 3. If the embedding dimension is 2, create a 2D plot
        if (embedding_dim == 2) {
          plot(ase$Xhat[, 1], ase$Xhat[, 2], 
               main = well, 
               xlab = "Dimension 1", 
               ylab = "Dimension 2",
               pch = 16,
               col = "blue")
        } 
        # 4. If the embedding dimension is 3 or more, create a 3D plot
        else if (embedding_dim >= 3) {
          scatterplot3d(
            x = ase$Xhat[, 1], 
            y = ase$Xhat[, 2], 
            z = ase$Xhat[, 3],
            main = well,
            xlab = "Dim 1",
            ylab = "Dim 2",
            zlab = "Dim 3",
            color = "blue",
            pch = 16,
            type = "p",
            grid = TRUE,
            box = TRUE
          )
        } 
        # Handle cases where the dimension is less than 2
        else {
          plot(1, type = "n", axes = FALSE, xlab = "", ylab = "")
          title(main = well)
          text(1, 1, paste("Embedding dim < 2\n(dim =", embedding_dim, ")"), col = "orange")
        }
      }
      
    }, error = function(e) {
      # If any other error occurs, plot a placeholder with the error message
      plot(1, type = "n", axes = FALSE, xlab = "", ylab = "")
      title(main = well)
      text(1, 1, "An error occurred", col = "red", cex = 0.8)
      message(paste("Error processing", well, ":", e$message))
    })
  }
  
  # --- Add the Main Title to the Entire Figure ---
  # Use mtext() to write the extracted path part in the outer margin (top)
  mtext(extract_path_part(file), 
        side = 3,       # 3 = top
        line = 1,       # Position in the margin
        outer = TRUE,   # Use the outer margin
        cex = 1.5,      # Character expansion (font size)
        font = 2)       # Font style (2 = bold)
  
  # Optional: Reset the plotting device layout to default
  # par(mfrow = c(1, 1), oma = c(0, 0, 0, 0))
  
  
}
