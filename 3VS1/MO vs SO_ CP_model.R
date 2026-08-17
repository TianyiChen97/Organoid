edge_list = read_csv("~/Desktop/Research /PhDresearch/Hopkins_Organoid/MO VS SO_2025_May_graphs/adjacency_edges_May_31_2025_ecr_results_no_window.csv")
filenames = unique(edge_list$File)
filenames = filenames[-c(1,2)]
filenames

# --- 1. Helper Function (Unchanged) ---
get_clique_density_stats <- function(A) {
  A_mat <- as.matrix(A)
  g <- graph_from_adjacency_matrix(A_mat, mode = "undirected", diag = FALSE)
  
  # Identify Core Nodes 
  kc <- coreness(g)
  core_index <- which(kc == max(kc))
  
  # Identify Periphery Nodes
  all_indices <- 1:nrow(A_mat)
  periph_index <- setdiff(all_indices, core_index)
  
  # Node Counts
  n_c <- length(core_index)
  n_p <- length(periph_index)
  
  # Calculate Edges
  edges_CC <- if(n_c > 1) sum(A_mat[core_index, core_index]) / 2 else 0
  edges_PP <- if(n_p > 1) sum(A_mat[periph_index, periph_index]) / 2 else 0
  edges_CP <- if(n_c > 0 && n_p > 0) sum(A_mat[core_index, periph_index]) else 0
  
  # Calculate Densities
  max_CC <- n_c * (n_c - 1) / 2
  max_PP <- n_p * (n_p - 1) / 2
  max_CP <- n_c * n_p
  
  dens_CC <- ifelse(max_CC > 0, edges_CC / max_CC, 0)
  dens_PP <- ifelse(max_PP > 0, edges_PP / max_PP, 0)
  dens_CP <- ifelse(max_CP > 0, edges_CP / max_CP, 0)
  
  # Create Matrix
  M <- matrix(c(dens_CC, dens_CP, 
                dens_CP, dens_PP), 
              nrow=2, byrow=TRUE)
  
  return(list(matrix = M, n_core = n_c, n_periph = n_p))
}

# --- 2. Main Processing Loop ---

# List to store the averaged stats for each file
file_stats_storage <- list()

for (file in filenames) {
  par(mfrow = c(2, 3), mar = c(1, 1, 3, 1), oma = c(2, 0, 3, 0))
  
  # Generate the clean key for this file (e.g., "M07357/Network/000389")
  # We use a regex to extract the part starting with Mxxxxx up to the parent folder of the .pkl
  # This replicates the behavior of your 'extract_path_part'
  file_key <- sub(".*(M[0-9]{5}/.*/[0-9]+)/.*", "\\1", file)
  
  # Lists to collect well-level stats for THIS file
  dens_matrices <- list()
  counts_list   <- list() # Stores c(n_core, n_peri, n_lcc, prop_core, prop_peri)
  
  well_list <- paste0("well", sprintf("%03d", 0:5))
  
  for (well in well_list) {
    subset_data <- subset(edge_list, File == file & Well == well)
    
    # Setup Blank Plot Canvas
    plot(1, type = "n", axes = FALSE, xlab = "", ylab = "", xlim=c(0, 10), ylim=c(0, 10))
    title(main = well, cex.main = 1.5)
    box()
    
    if (nrow(subset_data) == 0) {
      text(5, 5, "No data available", col = "red", cex = 1.2)
      next 
    }
    
    tryCatch({
      # Graph Construction
      n_rows <- n_cols <- unique(subset_data$dim)
      adj <- sparseMatrix(
        i    = subset_data$Row + 1,
        j    = subset_data$Column + 1,
        dims = c(n_rows, n_cols)
      )
      
      g <- graph_from_adjacency_matrix(adj, mode = "undirected", diag = FALSE)
      n_original_nodes <- vcount(g) 
      
      comps <- components(g, mode = "weak")
      g_lcc <- induced_subgraph(g, V(g)[comps$membership == which.max(comps$csize)])
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
        
        ratio_core <- n_core / n_lcc_nodes
        ratio_peri <- n_periph / n_lcc_nodes
        
        # --- STORE FOR FILE AVERAGE ---
        dens_matrices[[length(dens_matrices) + 1]] <- dens_mat
        
        # We store counts and ratios as a vector
        # 1:Core, 2:Peri, 3:LCC, 4:PropCore, 5:PropPeri
        counts_list[[length(counts_list) + 1]] <- c(n_core, n_periph, n_lcc_nodes, ratio_core, ratio_peri)
        
        # --- PLOT TEXT ---
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
    })
  }
  
  # --- CALCULATE FILE AVERAGE ---
  if (length(dens_matrices) > 0) {
    # Average Matrix
    avg_mat <- Reduce("+", dens_matrices) / length(dens_matrices)
    
    # Average Counts/Ratios (Reduce works on list of vectors too)
    avg_counts <- Reduce("+", counts_list) / length(counts_list)
    
    # Store in the master list using the clean key
    file_stats_storage[[file_key]] <- list(
      mat = avg_mat,
      stats = avg_counts # [Core, Peri, LCC, PropCore, PropPeri]
    )
  }
  
  mtext(file_key, side = 3, line = 1, outer = TRUE, cex = 1.5, font = 2)
}


# A. Print Individual File Averages
cat("\n=======================================================\n")
cat("   INDIVIDUAL FILE AVERAGES (Averaged over Wells)\n")
cat("=======================================================\n")

for (key in names(file_stats_storage)) {
  data <- file_stats_storage[[key]]
  mat <- data$mat
  stats <- data$stats # [Core, Peri, LCC, PropCore, PropPeri]
  
  cat(sprintf("\nFile: %s\n", key))
  cat("  CC: ", sprintf("%.3f", mat[1,1]), "  CP: ", sprintf("%.3f", mat[1,2]), "  PP: ", sprintf("%.3f", mat[2,2]), "\n")
  cat("  Core: ", sprintf("%.1f", stats[1]), "  Peri: ", sprintf("%.1f", stats[2]), "  LCC: ", sprintf("%.1f", stats[3]), "\n")
  cat("  PropCore: ", sprintf("%.1f%%", stats[4]*100), "  PropPeri: ", sprintf("%.1f%%", stats[5]*100), "\n")
}


# --- 3. Aggregate by Group (The "Simple" Way) ---

target_groups <- c("M07357", "M08438", "M07359")

cat("\n=======================================================\n")
cat("   FINAL GROUP AVERAGES \n")
cat("=======================================================\n")

for (group in target_groups) {
  
  # 1. Find all file keys that start with this group ID
  # This relies on file_key being "M07357/..." 
  matching_keys <- grep(paste0("^", group), names(file_stats_storage), value = TRUE)
  
  if (length(matching_keys) > 0) {
    # 2. Extract data for these files
    group_data <- file_stats_storage[matching_keys]
    
    # Extract matrices and stats vectors into separate lists for averaging
    list_of_mats  <- lapply(group_data, function(x) x$mat)
    list_of_stats <- lapply(group_data, function(x) x$stats)
    
    # 3. Compute Group Averages
    group_avg_mat <- Reduce("+", list_of_mats) / length(list_of_mats)
    group_avg_stats <- Reduce("+", list_of_stats) / length(list_of_stats)
    
    # Unpack for printing
    n_core_avg <- group_avg_stats[1]
    n_peri_avg <- group_avg_stats[2]
    n_lcc_avg  <- group_avg_stats[3]
    prop_core_avg <- group_avg_stats[4]
    prop_peri_avg <- group_avg_stats[5]
    
    # 4. Print
    cat(sprintf("\nGroup: %s (n=%d files)\n", group, length(matching_keys)))
    cat("-------------------------\n")
    cat("  Density Matrix:\n")
    cat(sprintf("    CC: %.4f\n", group_avg_mat[1,1]))
    cat(sprintf("    CP: %.4f\n", group_avg_mat[1,2]))
    cat(sprintf("    PP: %.4f\n", group_avg_mat[2,2]))
    cat("\n")
    cat("  Average Node Counts:\n")
    cat(sprintf("    Core:       %.2f\n", n_core_avg))
    cat(sprintf("    Periphery:  %.2f\n", n_peri_avg))
    cat(sprintf("    Total LCC:  %.2f\n", n_lcc_avg))
    cat("\n")
    cat("  Average Proportions:\n")
    cat(sprintf("    Core %%:     %.2f%%\n", prop_core_avg * 100))
    cat(sprintf("    Peri %%:     %.2f%%\n", prop_peri_avg * 100))
    
  } else {
    cat(sprintf("\nGroup: %s - No data found.\n", group))
  }
}
