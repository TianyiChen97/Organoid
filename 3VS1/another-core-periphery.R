library(igraph)

# 1. Define the graph structure
n_core <- 80
n_periphery <- 61 # Total vertices will be 141 as in your data
n_total <- n_core + n_periphery

# Start with an empty adjacency matrix
adj_matrix <- matrix(0, nrow = n_total, ncol = n_total)

# --- REVISION 1: Make the core very dense, but not a perfect clique ---
# This lowers the minimum degree of core nodes to help prevent a jump.
# We use a high probability (p=0.98) to make it "clique-like".
core_subgraph <- sample_gnp(n = n_core, p = 1, directed = FALSE)
adj_matrix[1:n_core, 1:n_core] <- as_adjacency_matrix(core_subgraph, sparse = FALSE)

# --- REVISION 2: Create a smoother, randomized gradient of connections ---
# This vector defines how many core nodes each periphery node will connect to.
# It ramps up smoothly from a low number to a high number.
num_connections_vector <- round(seq(from = 1, to = n_core ,  by = 2))
num_connections_vector <- sort(rep(num_connections_vector,round(n_periphery/length(num_connections_vector))))
length(num_connections_vector)

# Create the connections for the periphery
for (i in 1:n_periphery) {
  periphery_node_idx <- n_core + i
  num_connections <- num_connections_vector[i]
  
  # Randomly sample which core nodes to connect to.
  # This ensures the connections are spread evenly across the core.
  target_core_nodes <- sample(1:n_core, size = num_connections, replace = FALSE)
  
  adj_matrix[periphery_node_idx, target_core_nodes] <- 1
}

# Make the final matrix symmetric
adj_matrix <- adj_matrix + t(adj_matrix)
adj_matrix[adj_matrix > 0] <- 1

sort(apply(adj_matrix, 1, sum))
# You can now proceed with this revised 'adj_matrix'
# For example, check the degree distribution:
final_graph <- graph_from_adjacency_matrix(adj_matrix, mode = "undirected", diag = FALSE,weighted = NULL)
plot(sort(degree(final_graph)), main = "Degree Distribution of Revised Model", ylab = "Degree")

eigen_decomp <- eigen(adj_matrix)
eigenvectors <- eigen_decomp$vectors

# 4. Create the Plot
# Color the nodes by their degree to visualize the smooth transition
color_palette <- colorRampPalette(c("SkyBlue2", "red"))
# order(target_degrees) ensures the colors match the sorted degrees
node_colors <- color_palette(length(target_degrees))[rank(degree(g))]

# Plot the most informative eigenvectors
plot(eigenvectors[, 1], eigenvectors[, 2],
     col = node_colors,
     pch = 19,
     xlab = "Eigenvector 1",
     ylab = "Eigenvector 2",
     main = "Spectral Embedding from Target Degree Sequence"
)


