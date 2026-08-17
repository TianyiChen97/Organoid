library(igraph)
library(Matrix)

get_directed_clique_density_stats <- function(A) {
  # Ensure standard matrix format
  A_mat <- as.matrix(A)
  
  # Create directed graph for coreness calculation
  g <- graph_from_adjacency_matrix(A_mat, mode = "directed", diag = FALSE)
  
  # Identify Core Nodes 
  # mode="all" considers both in and out edges for structural importance
  kc <- coreness(g, mode = "all") 
  core_index <- which(kc == max(kc))
  
  # Identify Periphery Nodes
  all_indices <- 1:nrow(A_mat)
  periph_index <- setdiff(all_indices, core_index)
  
  # Node Counts
  n_c <- length(core_index)
  n_p <- length(periph_index)
  
  # --- EXTRACT SUB-MATRICES (drop=FALSE preserves dimensions) ---
  
  # Core -> Core
  Acc <- if(n_c > 0) A_mat[core_index, core_index, drop=FALSE] else matrix(0,0,0)
  
  # Periph -> Periph
  App <- if(n_p > 0) A_mat[periph_index, periph_index, drop=FALSE] else matrix(0,0,0)
  
  # Core -> Periph
  Acp <- if(n_c > 0 && n_p > 0) A_mat[core_index, periph_index, drop=FALSE] else matrix(0,0,0)
  
  # Periph -> Core
  Apc <- if(n_c > 0 && n_p > 0) A_mat[periph_index, core_index, drop=FALSE] else matrix(0,0,0)
  
  # --- 1. EDGE COUNTS (Sums) ---
  edges_CC <- sum(Acc)
  edges_PP <- sum(App)
  edges_CP <- sum(Acp)
  edges_PC <- sum(Apc)
  
  total_edges <- edges_CC + edges_PP + edges_CP + edges_PC
  
  # --- 2. DENSITIES ---
  max_CC <- n_c * (n_c - 1)
  max_PP <- n_p * (n_p - 1)
  max_CP <- n_c * n_p
  max_PC <- n_p * n_c 
  
  dens_CC <- ifelse(max_CC > 0, edges_CC / max_CC, 0)
  dens_PP <- ifelse(max_PP > 0, edges_PP / max_PP, 0)
  dens_CP <- ifelse(max_CP > 0, edges_CP / max_CP, 0)
  dens_PC <- ifelse(max_PC > 0, edges_PC / max_PC, 0)
  
  # --- 3. RECIPROCITY (RHO) CALCULATIONS ---
  
  # A. Global Reciprocity: sum(A * A') / sum(A)
  # Uses element-wise multiplication (*) with transpose
  rho_global <- 0
  if (total_edges > 0) {
    rho_global <- sum(A_mat * t(A_mat)) / total_edges
  }
  
  # B. Rho CC (Core Internal)
  rho_CC <- 0
  if (edges_CC > 0) {
    # Matches within core / edges within core
    rho_CC <- sum(Acc * t(Acc)) / edges_CC
  }
  
  # C. Rho PP (Periph Internal)
  rho_PP <- 0
  if (edges_PP > 0) {
    rho_PP <- sum(App * t(App)) / edges_PP
  }
  
  # D. Rho Inter-Block (Interaction between CP and PC)
  # Numerator: Count of mutual pairs * 2 (or explicitly sum both match directions)
  # Denominator: Total edges in both CP and PC blocks
  rho_inter <- 0
  edges_inter <- edges_CP + edges_PC
  
  if (edges_inter > 0 && n_c > 0 && n_p > 0) {
    # Check matches from C->P perspective
    matches_upper_right <- sum(Acp * t(Apc))
    
    # Check matches from P->C perspective (should be identical sum)
    matches_lower_left  <- sum(Apc * t(Acp))
    
    # Sum of mutual endpoints / Total endpoints
    rho_inter <- (matches_upper_right + matches_lower_left) / edges_inter
  }
  
  # --- FORMAT OUTPUT ---
  M <- matrix(c(dens_CC, dens_CP, 
                dens_PC, dens_PP), 
              nrow=2, byrow=TRUE)
  
  colnames(M) <- c("To_Core", "To_Periph")
  rownames(M) <- c("From_Core", "From_Periph")
  
  return(list(
    matrix = M, 
    n_core = n_c, 
    n_periph = n_p,
    rho_global = rho_global,
    rho_CC = rho_CC,
    rho_PP = rho_PP,
    rho_CP = rho_inter, # CP and PC share the same reciprocity
    rho_PC = rho_inter
  ))
}



get_outdegree_core_stats <- function(A, top_n = 10) {
  A_mat <- as.matrix(A)
  n_total <- nrow(A_mat)
  
  # --- 1. Identify Core by Top Out-Degree ---
  # Assuming A[i, j] = 1 means edge from i -> j, rowSums is out-degree
  out_degree <- rowSums(A_mat)
  
  # Ensure we don't request more nodes than exist
  k <- min(top_n, n_total)
  
  # Get indices of the top k nodes with highest out-degree
  # order() puts the indices of the largest values first
  core_index <- order(out_degree, decreasing = TRUE)[1:k]
  
  # --- 2. Identify Periphery Nodes ---
  all_indices <- 1:n_total
  periph_index <- setdiff(all_indices, core_index)
  
  # Node Counts
  n_c <- length(core_index)
  n_p <- length(periph_index)
  
  # --- 3. Calculate Directed Edges ---
  # Core -> Core
  edges_CC <- if(n_c > 0) sum(A_mat[core_index, core_index]) else 0
  
  # Core -> Periphery (Out-flow from Core)
  edges_CP <- if(n_c > 0 && n_p > 0) sum(A_mat[core_index, periph_index]) else 0
  
  # Periphery -> Core (In-flow to Core)
  edges_PC <- if(n_c > 0 && n_p > 0) sum(A_mat[periph_index, core_index]) else 0
  
  # Periphery -> Periphery
  edges_PP <- if(n_p > 0) sum(A_mat[periph_index, periph_index]) else 0
  
  # --- 4. Calculate Densities ---
  
  # Max possible edges (Directed: n*(n-1) for diagonal blocks)
  max_CC <- n_c * (n_c - 1)
  max_PP <- n_p * (n_p - 1)
  
  # Max possible edges for off-diagonal blocks
  max_CP <- n_c * n_p
  max_PC <- n_p * n_c 
  
  dens_CC <- ifelse(max_CC > 0, edges_CC / max_CC, 0)
  dens_PP <- ifelse(max_PP > 0, edges_PP / max_PP, 0)
  dens_CP <- ifelse(max_CP > 0, edges_CP / max_CP, 0)
  dens_PC <- ifelse(max_PC > 0, edges_PC / max_PC, 0)
  
  # Create Matrix (Rows = Source, Cols = Target)
  M <- matrix(c(dens_CC, dens_CP, 
                dens_PC, dens_PP), 
              nrow=2, byrow=TRUE)
  
  # Labeling for clarity
  colnames(M) <- c("To_Core", "To_Periph")
  rownames(M) <- c("From_Core", "From_Periph")
  
  return(list(matrix = M, n_core = n_c, n_periph = n_p, core_indices = core_index))
}


extract_network_path <- function(input_string) {
  # This regular expression looks for the literal string "Network/"
  # followed by one or more digits (\\d+).
  pattern <- "Network/\\d+"
  
  # str_extract will find and return the first match of the pattern.
  result <- str_extract(input_string, pattern)
  
  return(result)
}
extract_path_part <- function(filepath) {
  # Use sub() to find the pattern and return only that part
  sub(".*/(M\\d+/[A-Za-z]+?/\\d+)/.*", "\\1", filepath)
}


# ------------------------------------------------------------------
# CONFIGURATION
# ------------------------------------------------------------------
active_file_path <- "~/Desktop/Research /PhDresearch/Hopkins_Organoid/M07357/adjacency_edges_Aug_29_2025_M07357_para2.csv"
edge_list <- read.csv(active_file_path)

filenames = unique(edge_list$File)


parameter_choice <- sub(".*M07357_([a-zA-Z0-9]+)\\.csv", "\\1", active_file_path)
if (nchar(parameter_choice) > 10 | parameter_choice == active_file_path) {
  parameter_choice <- "Default" 
}
MAIN_IDENTIFIER <- "M07357"

# New Suffix for the DIRECTED Statistics Plot
# Added "DIRECTED" to ensure filename is distinct
TITLE_SUFFIX <- paste0("PCinSS:", parameter_choice, " + DIRECTED_CP_STATS")

# ------------------------------------------------------------------
# 2. DIRECTED CP STATISTICS PLOTTING LOOP
# ------------------------------------------------------------------
# Assumes 'filenames', 'edge_list', 'extract_network_path', etc. are defined

for (file in filenames) {
  
  # ----------------------------------------------------
  # 2. Filename and Directory Logic
  # ----------------------------------------------------
  
  path_part <- extract_network_path(file)
  
  full_title_string <- paste(MAIN_IDENTIFIER, path_part, sep = '/')
  full_title_string <- paste(full_title_string, TITLE_SUFFIX, sep = " | ")
  
  # Create unique filename
  pdf_filename_base <- full_title_string
  pdf_filename_base <- gsub("/", "_", pdf_filename_base) 
  pdf_filename_base <- gsub(" \\| ", "_", pdf_filename_base) 
  pdf_filename_base <- gsub(": ", "_", pdf_filename_base) 
  pdf_filename_base <- gsub(" \\+ ", "_", pdf_filename_base) 
  pdf_filename <- paste0(pdf_filename_base, '.pdf')
  
  output_dir <- "M07357_DIRECTED_CP_STATS_outputs" 
  
  if (!dir.exists(output_dir)) {
    dir.create(output_dir)
  }
  
  final_pdf_path <- file.path(output_dir, pdf_filename)
  
  # Increase Height to fit new stats
  pdf(file = final_pdf_path, width = 14, height = 12)
  
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
      # Create adjacency matrix
      n_rows <- n_cols <- unique(subset_data$dim)
      adj <- sparseMatrix(
        i    = subset_data$Row + 1,
        j    = subset_data$Column + 1,
        dims = c(n_rows, n_cols)
      )
      
      # Directed Graph
      g <- graph_from_adjacency_matrix(adj, mode = "directed", diag = FALSE)
      n_original_nodes <- vcount(g) 
      
      # LCC (Weakly Connected)
      comps <- components(g, mode = "weak")
      big <- which.max(comps$csize)
      g_lcc <- induced_subgraph(g, V(g)[comps$membership == big])
      n_lcc_nodes <- vcount(g_lcc) 
      
      if (n_lcc_nodes < 11) {
        text(5, 5, "LCC size < 11", col = "red", cex = 1.2)
      } else {
        
        # Calculate Stats using NEW function
        A_lcc <- as.matrix(as_adjacency_matrix(g_lcc))
        stats_res <- get_directed_clique_density_stats(A_lcc)
        
        dens_mat <- stats_res$matrix
        n_core   <- stats_res$n_core
        n_periph <- stats_res$n_periph
        
        ratio_core <- n_core / n_lcc_nodes
        ratio_peri <- n_periph / n_lcc_nodes
        
        # Build Text Output
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
  
  mtext(full_title_string, side = 3, line = 1, outer = TRUE, cex = 1.5, font = 2)
  dev.off()
}



file = filenames [1]
print(file)
well = "well000"

subset_data <- subset(edge_list, File == file & Well == well)
if (nrow(subset_data) == 0) {
  next
}
n_rows <-   n_cols   <- unique(subset_data$dim)

sparse_matrix <- sparseMatrix(
  i = subset_data$Row + 1,  
  j = subset_data$Column + 1,
  dims = c(n_rows, n_cols)
)

g <- graph_from_adjacency_matrix(sparse_matrix, weighted = NULL, mode = "directed", diag = FALSE)

components <- igraph::clusters(g, mode="weak")
biggest_cluster_id <- which.max(components$csize)
vert_ids <- V(g)[components$membership == biggest_cluster_id]
lcc_2 <- igraph::induced_subgraph(g, vert_ids)

A2 <- lcc_2[]

diag(A2)

sum(A2 * t(A2)) / sum(A2)



get_clique_density_stats(lcc)
get_directed_clique_density_stats(lcc)

get_outdegree_core_stats(lcc,30)

A = lcc[]
in_degree = colSums(A)
out_degree = rowSums(A) 
#par(mfrow=c(1,2))
#hist(in_degree, ylim = c(0,800), breaks = 100)
#hist(out_degree, ylim = c(0,800), breaks = 100)
#length(which(in_degree == 0))
#out_degree[which(in_degree == 0)]
#in_degree[which(out_degree > max(out_degree)-100)]
#length(which(out_degree == 0))
#dim(A)
degree_matrix <- cbind(in_degree, out_degree)
rownames(degree_matrix) <- paste0("Node_", seq_len(nrow(degree_matrix)))
colnames(degree_matrix) <- c("In-degree", "Out-degree")

library(pheatmap)

# Only show row labels for high-out-degree nodes 
show_labels <- function(degree_matrix, top_n = 10) {
  top_nodes <- order(degree_matrix[,'Out-degree'], decreasing = TRUE)[1:top_n]
  labels <- rep("", nrow(degree_matrix))
  labels[top_nodes] <- paste0("Node_", top_nodes)
  labels
}
pheatmap(
  degree_matrix,
  cluster_rows = FALSE,
  cluster_cols = FALSE,
  scale = "none",
  color = colorRampPalette(c("white", "blue"))(100),
  show_rownames = TRUE,
  labels_row = show_labels(degree_matrix, top_n = 20),  # only label top 10 by out degree
  angle_col = 0,  
  fontsize_row = 6,
  fontsize_col = 12,
  main = "In-Degree and Out-Degree per Node",
  legend = TRUE,
  border_color = NA
)



