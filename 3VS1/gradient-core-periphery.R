

library(igraph)

# 1. Define the graph structure with a gradient
n_core <- 80
n_periphery <- 60
n_total <- n_core + n_periphery

adj_matrix <- matrix(0, nrow = n_total, ncol = n_total)

adj_matrix[1:n_core, 1:n_core] <- 1
diag(adj_matrix) <- 0 # Ensure no self-loops

for (i in 1:n_periphery) {
  # This makes periphery nodes connect to an increasing number of core nodes
  # Node 11 connects to 1 core node, node 12 to 2, ..., up to node 20 connecting to all 10.
  # We repeat this pattern for nodes 21-30 and 31-40
  
  num_connections <- (i - 1) %% n_core + 20
  print(num_connections)
  # The current periphery node's index
  periphery_node_idx <- n_core + i
  
  # Connect this periphery node to the first 'num_connections' core nodes
  adj_matrix[periphery_node_idx, 1:num_connections] <- 1
  
}

# Make the matrix symmetric
adj_matrix <- adj_matrix + t(adj_matrix)
adj_matrix[adj_matrix > 0] <- 1

sort(apply(adj_matrix, 2, sum))
plot(sort(apply(adj_matrix, 2, sum)))


# 2. Get the adjacency matrix eigenvectors (the correct method from before)
eigen_decomp <- eigen(adj_matrix)
eigenvectors <- eigen_decomp$vectors

sort(abs(eigen_decomp$values))

node_colors <- c(rep('red',n_core),rep('green',n_periphery))


g <- graph_from_adjacency_matrix(adj_matrix)
plot(eigenvectors[, 1], eigenvectors[, 140],
     col = node_colors,
     pch = 19,
     xlab = "Eigenvector 1",
     ylab = "Eigenvector 140",
     main = "Simulated"
)

