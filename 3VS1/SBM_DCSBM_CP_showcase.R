library(igraph)
sample_bernoulli <- function(P) {
  n <- nrow(P)
  # Generate 0s and 1s based on the probability matrix P
  # We generate for the whole matrix to keep it simple, then symmetrize
  adj_vec <- rbinom(n * n, size = 1, prob = P)
  adj <- matrix(adj_vec, nrow = n, ncol = n)
  
  # Make it symmetric (undirected graph)
  adj[lower.tri(adj)] <- t(adj)[lower.tri(adj)]
  diag(adj) <- 0 # No self-loops
  
  return(graph_from_adjacency_matrix(adj, mode = "undirected"))
}
# --- Global Settings ---
n <- 200
# Define colors: Red for Core, Blue for Periphery
# We generate a vector of colors for plotting later
get_colors <- function(n, type) {
  if (type == "block") {
    # Discrete colors for SBM/DC-SBM
    c(rep("SkyBlue2", n/2), rep("red", n/2))
  } else {
    # Gradient colors for GCP
    colorRampPalette(c("SkyBlue2", "red"))(n)
  }
}

# ==========================================
# Model 1: Standard Stochastic Block Model (SBM)
# Signature: Tight Clusters (Blobs)
# ==========================================
# Two distinct blocks with homogeneous probabilities
n_sbm <- n
p_matrix <- matrix(c(0.05, 0.05,  # Periphery-Periphery, Periphery-Core
                     0.05, 0.80), # Core-Periphery, Core-Core
                   nrow = 2)
g_sbm <- sample_sbm(n_sbm, pref.matrix = p_matrix, block.sizes = c(n/2, n/2))

# Spectral Embedding
ase_sbm <- eigen(as_adjacency_matrix(g_sbm, sparse=FALSE))$vectors[, 1:2]


# ==========================================
# Model 2: Degree-Corrected SBM (DC-SBM)
# Signature: Linear Rays (Lines)
# ==========================================
# Two blocks, but nodes have varying 'popularity' (theta)
theta <- runif(n, min = 0.2, max = 1.0) # Random degree corrections
# Assign block vectors: Block 1 is [1, 0], Block 2 is [0, 1] (simplified structure)
# We simulate this by constructing the P matrix directly
B_dcsbm <- matrix(c(0.1, 0.1, 
                    0.1, 0.9), nrow=2)
# Assign half nodes to block 1, half to block 2
block_assign <- c(rep(1, n/2), rep(2, n/2))

# Construct Prob Matrix: P_ij = theta_i * theta_j * B_{g_i, g_j}
P_dcsbm <- matrix(0, n, n)
for(i in 1:n) {
  for(j in 1:n) {
    P_dcsbm[i,j] <- theta[i] * theta[j] * B_dcsbm[block_assign[i], block_assign[j]]
  }
}
P_dcsbm[P_dcsbm > 1] <- 1
g_dcsbm <- sample_bernoulli(P_dcsbm)

# Spectral Embedding
ase_dcsbm <- eigen(as_adjacency_matrix(g_dcsbm, sparse=FALSE))$vectors


# ==========================================
# Model 3: Gradient Core-Periphery (GCP)
# Signature: Parabolic Curve
# ==========================================
# Latent positions on a manifold curve X(t) = [t, gamma*t^2]

t_seq <- seq(0, 1, length.out = n)
gamma <- 10.0
# Latent positions X
X_gcp <- cbind(t_seq, gamma * t_seq^2)
# Generate probabilities P = X * X^T
P_gcp <- X_gcp %*% t(X_gcp)
P_gcp[P_gcp > 1] <- 1
g_gcp <- sample_bernoulli(P_gcp)

# Spectral Embedding
ase_gcp <- svd(as_adjacency_matrix(g_gcp, sparse=FALSE))
plot(ase_gcp$u[,1], ase_gcp$u[,2], col = get_colors(n, "gradient"), pch=19,
     main = "3. Gradient Core-Periphery", xlab = "Dim 1", ylab = "Dim 2",
     sub = "Signature: Parabolic Curve")

# ==========================================
# PLOTTING THE COMPARISON
# ==========================================
par(mfrow = c(1, 3)) # 3 plots in one row

# Plot 1: SBM
plot(ase_sbm[,1], ase_sbm[,2], col = get_colors(n, "block"), pch=19,
     main = "1. Standard SBM", xlab = "Dim 1", ylab = "Dim 2",
     sub = "Signature: Discrete Blobs")

# Plot 2: DC-SBM
plot(ase_dcsbm[,1], ase_dcsbm[,2], col = get_colors(n, "block"), pch=19,
     main = "2. Degree-Corrected SBM", xlab = "Dim 1", ylab = "Dim 2",
     sub = "Signature: Linear Rays")

# Plot 3: Gradient Core-Periphery
plot(ase_gcp$u[,1], ase_gcp$u[,2], col = get_colors(n, "gradient"), pch=19,
     main = "3. Gradient Core-Periphery", xlab = "Dim 1", ylab = "Dim 2",
     sub = "Signature: Parabolic Curve")

par(mfrow = c(1, 1)) # Reset layout
