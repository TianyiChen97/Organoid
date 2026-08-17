# --- Updated Helper Function ---
get_clique_density_stats <- function(A) {
  # 1. Convert to Graph
  A_mat <- as.matrix(A)
  g <- graph_from_adjacency_matrix(A_mat, mode = "undirected", diag = FALSE)
  
  # 2. Identify Core Nodes 
  kc <- coreness(g)
  core_index <- which(kc == max(kc))
  
  # 3. Identify Periphery Nodes
  all_indices <- 1:nrow(A_mat)
  periph_index <- setdiff(all_indices, core_index)
  
  # 4. Get Node Counts
  n_c <- length(core_index)
  n_p <- length(periph_index)
  
  # 5. Calculate Edges
  edges_CC <- if(n_c > 1) sum(A_mat[core_index, core_index]) / 2 else 0
  edges_PP <- if(n_p > 1) sum(A_mat[periph_index, periph_index]) / 2 else 0
  edges_CP <- if(n_c > 0 && n_p > 0) sum(A_mat[core_index, periph_index]) else 0
  
  # 6. Calculate Max Possible Edges
  max_CC <- n_c * (n_c - 1) / 2
  max_PP <- n_p * (n_p - 1) / 2
  max_CP <- n_c * n_p
  
  # 7. Calculate Densities
  dens_CC <- ifelse(max_CC > 0, edges_CC / max_CC, 0)
  dens_PP <- ifelse(max_PP > 0, edges_PP / max_PP, 0)
  dens_CP <- ifelse(max_CP > 0, edges_CP / max_CP, 0)
  
  # 8. Create Matrix
  M <- matrix(c(dens_CC, dens_CP, 
                dens_CP, dens_PP), 
              nrow=2, byrow=TRUE)
  colnames(M) <- c("C", "P")
  rownames(M) <- c("C", "P")
  
  # Return a list with all stats
  return(list(matrix = M, n_core = n_c, n_periph = n_p))
}


for (file in filenames) {
  # Setup grid of plots
  par(mfrow = c(2, 3), mar = c(1, 1, 3, 1), oma = c(2, 0, 3, 0))
  
  well_list <- paste0("well", sprintf("%03d", 0:5))
  
  for (well in well_list) {
    
    subset_data <- subset(edge_list, File == file & Well == well)
    
    # 1. Setup Blank Plot Canvas
    plot(1, type = "n", axes = FALSE, xlab = "", ylab = "", xlim=c(0, 10), ylim=c(0, 10))
    title(main = well, cex.main = 1.5)
    box() # Optional: puts a box around the stats
    
    if (nrow(subset_data) == 0) {
      text(5, 5, "No data available", col = "red", cex = 1.2)
      next 
    }
    
    tryCatch({
      # Create adjacency matrix and LCC
      n_rows <- n_cols <- unique(subset_data$dim)
      adj <- sparseMatrix(
        i    = subset_data$Row + 1,
        j    = subset_data$Column + 1,
        dims = c(n_rows, n_cols)
      )
      
      # Original Graph
      g <- graph_from_adjacency_matrix(adj, mode = "undirected", diag = FALSE)
      n_original_nodes <- vcount(g) 
      
      comps <- components(g, mode = "weak")
      big <- which.max(comps$csize)
      g_lcc <- induced_subgraph(g, V(g)[comps$membership == big])
      n_lcc_nodes <- vcount(g_lcc) 
      
      if (n_lcc_nodes < 11) {
        text(5, 5, "LCC size < 11", col = "red", cex = 1.2)
      } else {
        
        # --- CALCULATE STATISTICS ---
        A_lcc <- as.matrix(as_adjacency_matrix(g_lcc))
        stats_res <- get_clique_density_stats(A_lcc)
        
        dens_mat <- stats_res$matrix
        n_core   <- stats_res$n_core
        n_periph <- stats_res$n_periph
        
        # Ratios
        ratio_core <- n_core / n_lcc_nodes
        ratio_peri <- n_periph / n_lcc_nodes
        
        # --- DISPLAY TEXT ---
        # We format the numbers into a clean block of text
        # x=5, y=5 centers it in the 0-10 coordinate system
        
        lines_to_print <- c(
          "--- Densities ---",
          paste0("CC: ", sprintf("%.3f", dens_mat[1,1])),
          paste0("CP: ", sprintf("%.3f", dens_mat[1,2])),
          paste0("PP: ", sprintf("%.3f", dens_mat[2,2])),
          "",
          "--- Counts ---",
          paste0("Core: ", n_core),
          paste0("Peri: ", n_periph),
          paste0("LCC Size: ", n_lcc_nodes),
          paste0("Orig Size: ", n_original_nodes),
          "",
          "--- Ratios ---",
          paste0("Core / LCC: ", sprintf("%.3f", ratio_core)),
          paste0("Peri / LCC: ", sprintf("%.3f", ratio_peri))
        )
        
        text(5, 5, labels = paste(lines_to_print, collapse = "\n"), 
             cex = 1.3, family = "mono") 
      }
      
    }, error = function(e) {
      text(5, 5, "Error", col = "red")
      message(paste("Error processing", well, ":", e$message))
    })
  }
  
  mtext(extract_network_path(file), 
        side = 3, 
        line = 1, 
        outer = TRUE, 
        cex = 1.5, 
        font = 2)
}
