edge_list = read.csv('/Users/tianyichen/Desktop/Research /PhDresearch/Hopkins_Organoid/Codes/R/adjacency_edges.csv')
filenames = unique(edge_list$File)


library(igraph)
library(Matrix)
library(stringr)

# ===================================================================
# 1. HELPER FUNCTIONS
# ===================================================================

# Function to extract "Mxxxxx/Type/xxxxxx" from the path
extract_network_path <- function(input_string) {
  # Regex explanation:
  # M\d+        -> Matches "M" followed by digits (e.g., M07914)
  # /           -> Literal slash
  # [A-Za-z]+   -> Matches letters (e.g., Network, Stimulation, ActivityScan)
  # /           -> Literal slash
  # \d+         -> Matches digits (e.g., 000297)
  pattern <- "M\\d+/[A-Za-z]+/\\d+"
  
  result <- str_extract(input_string, pattern)
  return(result)
}

get_directed_clique_density_stats <- function(A) {
  # [PASTE THE UPDATED FUNCTION FROM PREVIOUS TURN HERE]
  # (I will assume the function definition is already loaded in your environment 
  #  as defined in our previous step. If not, copy-paste it here.)
  
  # Ensure standard matrix format
  A_mat <- as.matrix(A)
  g <- graph_from_adjacency_matrix(A_mat, mode = "directed", diag = FALSE)
  kc <- coreness(g, mode = "all") 
  core_index <- which(kc == max(kc))
  all_indices <- 1:nrow(A_mat)
  periph_index <- setdiff(all_indices, core_index)
  n_c <- length(core_index)
  n_p <- length(periph_index)
  
  Acc <- if(n_c > 0) A_mat[core_index, core_index, drop=FALSE] else matrix(0,0,0)
  App <- if(n_p > 0) A_mat[periph_index, periph_index, drop=FALSE] else matrix(0,0,0)
  Acp <- if(n_c > 0 && n_p > 0) A_mat[core_index, periph_index, drop=FALSE] else matrix(0,0,0)
  Apc <- if(n_c > 0 && n_p > 0) A_mat[periph_index, core_index, drop=FALSE] else matrix(0,0,0)
  
  edges_CC <- sum(Acc); edges_PP <- sum(App); edges_CP <- sum(Acp); edges_PC <- sum(Apc)
  total_edges <- edges_CC + edges_PP + edges_CP + edges_PC
  
  max_CC <- n_c * (n_c - 1); max_PP <- n_p * (n_p - 1)
  max_CP <- n_c * n_p; max_PC <- n_p * n_c 
  
  dens_CC <- ifelse(max_CC > 0, edges_CC / max_CC, 0)
  dens_PP <- ifelse(max_PP > 0, edges_PP / max_PP, 0)
  dens_CP <- ifelse(max_CP > 0, edges_CP / max_CP, 0)
  dens_PC <- ifelse(max_PC > 0, edges_PC / max_PC, 0)
  
  rho_global <- 0
  if (total_edges > 0) rho_global <- sum(A_mat * t(A_mat)) / total_edges
  
  rho_CC <- 0; if (edges_CC > 0) rho_CC <- sum(Acc * t(Acc)) / edges_CC
  rho_PP <- 0; if (edges_PP > 0) rho_PP <- sum(App * t(App)) / edges_PP
  
  rho_inter <- 0
  edges_inter <- edges_CP + edges_PC
  if (edges_inter > 0 && n_c > 0 && n_p > 0) {
    matches_upper_right <- sum(Acp * t(Apc))
    matches_lower_left  <- sum(Apc * t(Acp))
    rho_inter <- (matches_upper_right + matches_lower_left) / edges_inter
  }
  
  M <- matrix(c(dens_CC, dens_CP, dens_PC, dens_PP), nrow=2, byrow=TRUE)
  colnames(M) <- c("To_Core", "To_Periph"); rownames(M) <- c("From_Core", "From_Periph")
  
  return(list(matrix = M, n_core = n_c, n_periph = n_p, rho_global = rho_global, 
              rho_CC = rho_CC, rho_PP = rho_PP, rho_CP = rho_inter, rho_PC = rho_inter))
}

# ===================================================================
# 2. CONFIGURATION & FILE SETUP
# ===================================================================

# Assuming 'edge_list' is already loaded with the new 3Vs1 data
# If you need to reload it:
# edge_list <- read.csv("path/to/your/new_3Vs1_edge_list.csv")

filenames <- unique(edge_list$File)
MAIN_IDENTIFIER <- "3Vs1"
output_dir <- paste0(MAIN_IDENTIFIER, "_DIRECTED_CP_STATS_outputs")

if (!dir.exists(output_dir)) {
  dir.create(output_dir)
}

# ===================================================================
# 3. MAIN LOOP
# ===================================================================

for (file in filenames) {
  
  # --- Title and Filename Logic ---
  
  # Extracts: "M07914/ActivityScan/000297"
  path_part <- extract_network_path(file)
  
  # Check if extraction worked, fallback if NULL
  if (is.na(path_part)) {
    path_part <- basename(dirname(file)) # Fallback to folder name
  }
  
  full_title_string <- path_part
  
  # Create safe filename
  pdf_filename_base <- full_title_string
  pdf_filename_base <- gsub("/", "_", pdf_filename_base) 
  pdf_filename_base <- gsub(" ", "_", pdf_filename_base)
  pdf_filename <- paste0(pdf_filename_base, '_DIRECTED_CP.pdf')
  
  final_pdf_path <- file.path(output_dir, pdf_filename)
  
  # --- Plotting ---
  
  pdf(file = final_pdf_path, width = 14, height = 12)
  
  # Grid: 2 rows, 3 columns (for Wells 0-5)
  par(mfrow = c(2, 3), mar = c(1, 1, 3, 1), oma = c(2, 0, 3, 0))
  
  well_list <- paste0("well", sprintf("%03d", 0:5))
  
  for (well in well_list) {
    
    subset_data <- subset(edge_list, File == file & Well == well)
    
    plot(1, type = "n", axes = FALSE, xlab = "", ylab = "", xlim=c(0, 10), ylim=c(0, 10))
    title(main = well, cex.main = 1.5)
    box() 
    
    if (nrow(subset_data) == 0) {
      text(5, 5, "No data available", col = "red", cex = 1.2)
      next 
    }
    
    tryCatch({
      n_rows <- n_cols <- unique(subset_data$dim)
      adj <- sparseMatrix(
        i    = subset_data$Row + 1,
        j    = subset_data$Column + 1,
        dims = c(n_rows, n_cols)
      )
      
      # Directed Graph & LCC
      g <- graph_from_adjacency_matrix(adj, mode = "directed", diag = FALSE)
      n_original_nodes <- vcount(g) 
      
      comps <- components(g, mode = "weak")
      big <- which.max(comps$csize)
      g_lcc <- induced_subgraph(g, V(g)[comps$membership == big])
      n_lcc_nodes <- vcount(g_lcc) 
      
      if (n_lcc_nodes < 11) {
        text(5, 5, "LCC size < 11", col = "red", cex = 1.2)
      } else {
        
        # Calculate Stats
        A_lcc <- as.matrix(as_adjacency_matrix(g_lcc))
        stats_res <- get_directed_clique_density_stats(A_lcc)
        
        dens_mat <- stats_res$matrix
        n_core   <- stats_res$n_core
        n_periph <- stats_res$n_periph
        ratio_core <- n_core / n_lcc_nodes
        ratio_peri <- n_periph / n_lcc_nodes
        
        # Text Output
        lines_to_print <- c(
          "--- Densities (Dir) ---",
          paste0("C -> C: ", sprintf("%.3f", dens_mat["From_Core","To_Core"])),
          paste0("C -> P: ", sprintf("%.3f", dens_mat["From_Core","To_Periph"])),
          paste0("P -> C: ", sprintf("%.3f", dens_mat["From_Periph","To_Core"])), 
          paste0("P -> P: ", sprintf("%.3f", dens_mat["From_Periph","To_Periph"])),
          "",
          "--- Block Symmetry ---",
          paste0("Rho CC: ", sprintf("%.3f", stats_res$rho_CC)),
          paste0("Rho CP: ", sprintf("%.3f", stats_res$rho_CP)),
          paste0("Rho PC: ", sprintf("%.3f", stats_res$rho_PC)),
          paste0("Rho PP: ", sprintf("%.3f", stats_res$rho_PP)),
          "",
          "--- Counts ---",
          paste0("Core: ", n_core),
          paste0("Peri: ", n_periph),
          paste0("LCC Size: ", n_lcc_nodes),
          paste0("Orig Size: ", n_original_nodes),
          "",
          "--- Global Symmetry ---",
          paste0("Mutual/All: ", sprintf("%.3f", stats_res$rho_global)),
          "",
          "--- Ratios ---",
          paste0("Core / LCC: ", sprintf("%.3f", ratio_core)),
          paste0("Peri / LCC: ", sprintf("%.3f", ratio_peri))
        )
        
        text(5, 5, labels = paste(lines_to_print, collapse = "\n"), 
             cex = 1.1, family = "mono") 
      }
      
    }, error = function(e) {
      text(5, 5, "Error", col = "red")
      message(paste("Error processing", well, ":", e$message))
    })
  }
  
  # Main Title: Uses the extracted Mxxxxx/Type/xxxxxx path
  mtext(full_title_string, side = 3, line = 1, outer = TRUE, cex = 1.5, font = 2)
  dev.off()
  print(paste("Processed:", full_title_string))
}
