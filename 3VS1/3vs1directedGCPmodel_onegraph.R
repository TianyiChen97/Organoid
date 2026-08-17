# 3vs1
edge_list = read.csv('/Users/tianyichen/Desktop/Research /PhDresearch/Hopkins_Organoid/Codes/R/adjacency_edges.csv')
filenames = unique(edge_list$File)
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

A = g[]
components <- igraph::clusters(g, mode="weak")
biggest_cluster_id <- which.max(components$csize)
vert_ids <- V(g)[components$membership == biggest_cluster_id]
lcc_1 <- igraph::induced_subgraph(g, vert_ids)

A1 <- lcc_1[]
sum(A1 * t(A1)) / sum(A1)


get_clique_density_stats(lcc)
get_directed_clique_density_stats(lcc)

get_outdegree_core_stats(lcc,400)

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





library(igraph)

get_directed_clique_density_stats <- function(A) {
  A_mat <- as.matrix(A)
  
  # Create directed graph
  g <- graph_from_adjacency_matrix(A_mat, mode = "directed", diag = FALSE)
  
  # Identify Core Nodes 
  # mode="all" considers both in and out edges for structural importance, 
  # keeping the logic consistent with the undirected version.
  kc <- coreness(g, mode = "all") 
  core_index <- which(kc == max(kc))
  
  # Identify Periphery Nodes
  all_indices <- 1:nrow(A_mat)
  periph_index <- setdiff(all_indices, core_index)
  
  # Node Counts
  n_c <- length(core_index)
  n_p <- length(periph_index)
  
  # --- Calculate Directed Edges (Sums) ---
  # Note: We do NOT divide by 2 for diagonal blocks in directed graphs
  
  # Core -> Core
  edges_CC <- if(n_c > 0) sum(A_mat[core_index, core_index]) else 0
  
  # Core -> Periphery
  edges_CP <- if(n_c > 0 && n_p > 0) sum(A_mat[core_index, periph_index]) else 0
  
  # Periphery -> Core (Distinct from CP now)
  edges_PC <- if(n_c > 0 && n_p > 0) sum(A_mat[periph_index, core_index]) else 0
  
  # Periphery -> Periphery
  edges_PP <- if(n_p > 0) sum(A_mat[periph_index, periph_index]) else 0
  
  # --- Calculate Densities ---
  
  # Max possible edges in directed graph is n*(n-1) for diagonal blocks
  max_CC <- n_c * (n_c - 1)
  max_PP <- n_p * (n_p - 1)
  
  # Max possible edges for off-diagonal is just n_row * n_col
  max_CP <- n_c * n_p
  max_PC <- n_p * n_c # Mathematically same as above, written for clarity
  
  dens_CC <- ifelse(max_CC > 0, edges_CC / max_CC, 0)
  dens_PP <- ifelse(max_PP > 0, edges_PP / max_PP, 0)
  dens_CP <- ifelse(max_CP > 0, edges_CP / max_CP, 0)
  dens_PC <- ifelse(max_PC > 0, edges_PC / max_PC, 0)
  
  # Create Matrix (Rows = Source, Cols = Target)
  # Row 1: Core -> Core, Core -> Periph
  # Row 2: Periph -> Core, Periph -> Periph
  M <- matrix(c(dens_CC, dens_CP, 
                dens_PC, dens_PP), 
              nrow=2, byrow=TRUE)
  
  colnames(M) <- c("To_Core", "To_Periph")
  rownames(M) <- c("From_Core", "From_Periph")
  
  return(list(matrix = M, n_core = n_c, n_periph = n_p))
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
