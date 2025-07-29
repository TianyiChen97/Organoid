library(igraph)

# 1. Define the Target Degree Sequence (from your data)
target_degrees <- c(
  1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   2,   3,   3,   9,
  9,  10,  11,  11,  12,  13,  14,  15,  16,  17,  18,  21,  21,  22,  22,  22,
  24,  24,  28,  28,  29,  29,  30,  33,  34,  38,  40,  42,  45,  49,  50,  51,
  52,  55,  56,  56,  58,  65,  65,  66,  74,  75,  78,  79,  80,  80,  82,  82,
  82,  82,  83,  83,  83,  83,  84,  84,  84,  85,  85,  85,  86,  87,  87,  87,
  88,  88,  89,  89,  90,  90,  90,  91,  91,  92,  92,  93,  94,  94,  95,  96,
  96,  96,  96,  97,  97,  97,  97,  98,  98, 100, 101, 101, 102, 102, 103, 103,
  103, 103, 104, 105, 106, 107, 108, 112, 113, 113, 114, 115, 115, 116, 117, 118,
  119, 120, 120, 122, 124, 125, 125, 125, 125, 125, 127, 128, 140
)

# 2. Generate a graph with this exact degree sequence
# Note: The sum of degrees must be even, which it is in this case.
# The 'vl' method is for constructing a graph from a valid degree sequence.
g <- sample_degseq(target_degrees, method = "vl")

# 3. Perform Adjacency Spectral Embedding
adj_matrix <- as_adjacency_matrix(g, sparse = FALSE)
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
