edge_list = read_csv("~/Desktop/Research /PhDresearch/Hopkins_Organoid/MO VS SO_2025_May_graphs/adjacency_edges_May_31_2025_ecr_results_no_window.csv")
filenames = unique(edge_list$File)
filenames = filenames[-2]

filenames
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
        
        # --- CALCULATE DENSITY MATRIX ---
        # Get dense matrix of LCC for the function
        A_lcc <- as.matrix(as_adjacency_matrix(g_lcc))
        dens_mat <- get_clique_density_matrix(A_lcc)
        
        # Format the matrix text for the legend
        mat_text <- c(
          paste("CC:", sprintf("%.2f", dens_mat[1,1])),
          paste("CP:", sprintf("%.2f", dens_mat[1,2])),
          paste("PP:", sprintf("%.2f", dens_mat[2,2]))
        )
        
        
        # --- Proceed with ASE only if LCC size is 11 or greater ---
        
        # Perform Adjacency Spectral Embedding (ASE) on the LCC
        ase <- full.ase(g_lcc, 10)
        
        # Find the elbow to determine the embedding dimension
        elb <- getElbows(ase$eval, plot = F)
        
        if (length(elb) < 2) {
          stop("Embedding dimension is less than 2.")
        }
        
        embedding_dim <- elb[2]
        # --- Conditional Plotting ---
        
        # 3. If the embedding dimension is 2, create a 2D plot
        if (embedding_dim == 2) {
          plot(ase$Xhat[, 1], ase$Xhat[, 2], 
               main = well, 
               xlab = "Dimension 1", 
               ylab = "Dimension 2",
               pch = 16,
               col = "blue")
          
          legend("topright", legend = mat_text, bty = "n", cex = 2, title="Density")
          
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
          
          legend("topright", 
                 legend = mat_text, bty = "n", 
                 cex = 2, 
                 title="Density",
                 inset=c(0.05, 0))
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


# --- PART 1: Calculate Average per File ---

# Initialize list to store results per file
file_averages <- list()

for (file in filenames) {
  
  file_key <- extract_path_part(file)
  
  # List to store density matrices for valid wells in this specific file
  valid_matrices <- list()
  well_list <- paste0("well", sprintf("%03d", 0:5))
  
  for (well in well_list) {
    subset_data <- subset(edge_list, File == file & Well == well)
    
    if (nrow(subset_data) == 0) next
    
    tryCatch({
      # Construct Graph
      n_rows <- n_cols <- unique(subset_data$dim)
      adj <- sparseMatrix(
        i    = subset_data$Row + 1,
        j    = subset_data$Column + 1,
        dims = c(n_rows, n_cols)
      )
      
      g <- graph_from_adjacency_matrix(adj, mode = "undirected", diag = FALSE)
      
      # Get LCC
      comps <- components(g, mode = "weak")
      g_lcc <- induced_subgraph(g, V(g)[comps$membership == which.max(comps$csize)])
      
      # Filter by size >= 11
      if (vcount(g_lcc) >= 11) {
        A_lcc <- as.matrix(as_adjacency_matrix(g_lcc))
        dens_mat <- get_clique_density_matrix(A_lcc)
        valid_matrices[[length(valid_matrices) + 1]] <- dens_mat
      }
      
    }, error = function(e) {
      # Silent error handling for cleanliness
    })
  }
  
  # Calculate Average for this File
  if (length(valid_matrices) > 0) {
    avg_matrix <- Reduce("+", valid_matrices) / length(valid_matrices)
    file_averages[[file_key]] <- avg_matrix
  } else {
    file_averages[[file_key]] <- NA
  }
}


# --- PART 2: Aggregate by Experiment Group ---

target_groups <- c("M07357", "M08438", "M07359")
group_averages <- list()

for (group in target_groups) {
  
  # 1. Find all file keys that start with this group ID (e.g., "M07357/...")
  matching_keys <- grep(paste0("^", group), names(file_averages), value = TRUE)
  
  # 2. Extract the matrices, filtering out any NA entries
  group_matrices <- file_averages[matching_keys]
  group_matrices <- group_matrices[ !sapply(group_matrices, function(x) all(is.na(x))) ]
  
  # 3. Calculate the group average
  if (length(group_matrices) > 0) {
    group_mean <- Reduce("+", group_matrices) / length(group_matrices)
    group_averages[[group]] <- group_mean
  } else {
    group_averages[[group]] <- NA
  }
}

# --- Output Results ---
print(names(file_averages)) # Printing names to verify keys

print(group_averages)








