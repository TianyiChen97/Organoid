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


# --- Load Data ---
edge_list = read.csv("~/Desktop/Research /PhDresearch/Hopkins_Organoid/M07357/adjacency_edges_Aug_29_2025_M07357_para2.csv")
filenames = unique(edge_list$File)

# --- Main Loop ---
for (file in filenames) {
  par(mfrow = c(2, 3), mar = c(3, 3, 3, 1), oma = c(4, 0, 3, 0))
  
  well_list <- paste0("well", sprintf("%03d", 0:5))
  
  for (well in well_list) {
    
    subset_data <- subset(edge_list, File == file & Well == well)
    
    if (nrow(subset_data) == 0) {
      plot(1, type = "n", axes = FALSE, xlab = "", ylab = "")
      title(main = well)
      text(1, 1, "No data available", col = "red")
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
        plot(1, type = "n", axes = FALSE, xlab = "", ylab = "")
        title(main = well)
        text(1, 1, "Size of largest\nconnected component < 11", col = "red", cex = 0.9)
      } else {
        
        # --- CALCULATE DENSITY & COUNTS ---
        A_lcc <- as.matrix(as_adjacency_matrix(g_lcc))
        
        # Get stats
        stats_res <- get_clique_density_stats(A_lcc)
        dens_mat <- stats_res$matrix
        n_core   <- stats_res$n_core
        n_periph <- stats_res$n_periph
        
        # Calculate Ratios
        ratio_core <- n_core / n_lcc_nodes
        ratio_peri <- n_periph / n_lcc_nodes
        
        # Updated Legend Text
        mat_text <- c(
          paste("CC:", sprintf("%.2f", dens_mat[1,1])),
          paste("CP:", sprintf("%.2f", dens_mat[1,2])),
          paste("PP:", sprintf("%.2f", dens_mat[2,2])),
          paste("Core:", n_core),
          paste("Peri:", n_periph),
          paste("LCC:", n_lcc_nodes),
          paste("Orig:", n_original_nodes),
          paste("C/LCC:", sprintf("%.2f", ratio_core)), # Added Ratio
          paste("P/LCC:", sprintf("%.2f", ratio_peri))  # Added Ratio
        )
        
        # --- Perform ASE ---
        ase <- full.ase(g_lcc, 10)
        elb <- getElbows(ase$eval, plot = F)
        
        if (length(elb) < 2) {
          stop("Embedding dimension is less than 2.")
        }
        
        embedding_dim <- elb[2]
        
        # --- Plotting ---
        if (embedding_dim == 2) {
          plot(ase$Xhat[, 1], ase$Xhat[, 2], 
               main = well, 
               xlab = "Dimension 1", 
               ylab = "Dimension 2",
               pch = 16,
               col = "blue")
          
          # Legend with extra stats
          legend("topright", legend = mat_text, bty = "n", cex = 1.1, title="Stats")
          
        } else if (embedding_dim >= 3) {
          sp <- scatterplot3d(
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
                 cex = 1.1, 
                 title="Stats",
                 inset=c(0.05, 0))
          
        } else {
          plot(1, type = "n", axes = FALSE, xlab = "", ylab = "")
          title(main = well)
          text(1, 1, paste("Embedding dim < 2\n(dim =", embedding_dim, ")"), col = "orange")
        }
      }
      
    }, error = function(e) {
      plot(1, type = "n", axes = FALSE, xlab = "", ylab = "")
      title(main = well)
      text(1, 1, "An error occurred", col = "red", cex = 0.8)
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
