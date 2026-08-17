library(Matrix)
library(igraph)
library(mclust)
library(randnet)
library(irlba)
library(ggplot2)
library(MASS)
library(mvtnorm)
library(dplyr)
library(RColorBrewer)
library(readr)
library(patchwork)  
library(iGraphMatch) 
library(scatterplot3d)
library(stringr)

# ===================================================================
# 1. HELPER FUNCTIONS
# ===================================================================

# User-defined regex for path extraction
extract_network_path <- function(input_string) {
  pattern <- "M\\d+/[A-Za-z]+/\\d+"
  result <- str_extract(input_string, pattern)
  return(result)
}

# --- Embedding Functions ---

full.ase <- function(A, d, diagaug=TRUE, doptr=FALSE) {
  require(irlba)
  if (doptr) { g <- ptr(A); A <- g[] } else { A <- A[] }
  A = as.matrix(A)
  if (diagaug) { diag(A) <- rowSums(A) / (nrow(A)-1) }
  A.svd <- irlba(A,d)
  Xhat <- A.svd$u %*% diag(sqrt(A.svd$d))
  Xhat.R <- NULL
  if (!isSymmetric(A)) { Xhat.R <- A.svd$v %*% diag(sqrt(A.svd$d)) }
  return(list(eval=A.svd$d, Xhat=Matrix(Xhat), Xhat.R=Xhat.R))
}

full.lse <- function(A, d) {
  A <- A[]
  A = as.matrix(A)
  deg.seq <- rowSums(A)
  # Prevent division by zero
  deg.seq[deg.seq == 0] <- 1 
  L.svd <- irlba(A/sqrt(outer(deg.seq,deg.seq)),d)
  Xhat <- L.svd$u %*% diag(sqrt(L.svd$d))
  Xhat.R <- NULL
  if (!isSymmetric(A)) { Xhat.R <- L.svd$v %*% diag(sqrt(L.svd$d)) }
  return(list(eval=L.svd$d, Xhat=(Xhat), Xhat.R=(Xhat.R)))
}

getElbows <- function(dat, n = 3, threshold = FALSE, plot = TRUE, main="") {
  if (is.matrix(dat)) { d <- sort(apply(dat,2,sd), decreasing=TRUE) } else { d <- sort(dat,decreasing=TRUE) }
  if (!is.logical(threshold)) d <- d[d > threshold]
  p <- length(d)
  if (p == 0) stop("d must have elements larger than threshold")
  lq <- rep(0.0, p)
  for (q in 1:p) {
    mu1 <- mean(d[1:q])
    mu2 <- mean(d[-(1:q)])
    sigma2 <- (sum((d[1:q] - mu1)^2) + sum((d[-(1:q)] - mu2)^2)) / (p - 1 - (q < p))
    lq[q] <- sum( dnorm(  d[1:q ], mu1, sqrt(sigma2), log=TRUE) ) + sum( dnorm(d[-(1:q)], mu2, sqrt(sigma2), log=TRUE) )
  }
  q <- which.max(lq)
  if (n > 1 && q < (p-1)) { q <- c(q, q + getElbows(d[(q+1):p], n-1, plot=FALSE)) }
  if (plot==TRUE) {
    if (is.matrix(dat)) { sdv <- d; plot(sdv,type="b",xlab="dim",ylab="stdev",main=main); points(q,sdv[q],col=2,pch=19) } 
    else { plot(dat, type="b",main=main); points(q,dat[q],col=2,pch=19) }
  }
  return(q)
}

# --- Undirected CP Stats Function ---

get_clique_density_stats <- function(A) {
  A_mat <- as.matrix(A)
  # Mode is UNDIRECTED here
  g <- graph_from_adjacency_matrix(A_mat, mode = "undirected", diag = FALSE)
  
  kc <- coreness(g)
  core_index <- which(kc == max(kc))
  periph_index <- setdiff(1:nrow(A_mat), core_index)
  
  n_c <- length(core_index); n_p <- length(periph_index)
  
  # Edges (divided by 2 for undirected)
  edges_CC <- if(n_c > 1) sum(A_mat[core_index, core_index]) / 2 else 0
  edges_PP <- if(n_p > 1) sum(A_mat[periph_index, periph_index]) / 2 else 0
  edges_CP <- if(n_c > 0 && n_p > 0) sum(A_mat[core_index, periph_index]) else 0
  
  max_CC <- n_c * (n_c - 1) / 2
  max_PP <- n_p * (n_p - 1) / 2
  max_CP <- n_c * n_p
  
  dens_CC <- ifelse(max_CC > 0, edges_CC / max_CC, 0)
  dens_PP <- ifelse(max_PP > 0, edges_PP / max_PP, 0)
  dens_CP <- ifelse(max_CP > 0, edges_CP / max_CP, 0)
  
  M <- matrix(c(dens_CC, dens_CP, dens_CP, dens_PP), nrow=2, byrow=TRUE)
  colnames(M) <- c("C", "P"); rownames(M) <- c("C", "P")
  
  return(list(matrix = M, n_core = n_c, n_periph = n_p))
}

# ===================================================================
# 2. CONFIGURATION & DATA LOADING
# ===================================================================

# Load the MOvsSO dataset
edge_list <- read_csv("~/Desktop/Research /PhDresearch/Hopkins_Organoid/MO VS SO_2025_May_graphs/adjacency_edges_May_31_2025_ecr_results_no_window.csv")

# Filter filenames as requested
filenames = unique(edge_list$File)
filenames = filenames[-c(1,2)] # Removing first two files as per your logic


file = filenames[9]
well = 'well001'

EMBEDDING_METHOD = "ASE"

if (EMBEDDING_METHOD == "LSE") {
  embedding_function <- full.lse 
} else {
  embedding_function <- full.ase
}

subset_data <- subset(edge_list, File == file & Well == well)


n_rows <- n_cols <- unique(subset_data$dim)
adj <- sparseMatrix(i=subset_data$Row+1, j=subset_data$Column+1, dims=c(n_rows, n_cols))
  
  # Undirected Graph & LCC
g <- graph_from_adjacency_matrix(adj, mode="undirected", diag=FALSE)
comps <- components(g, mode="weak")
g_lcc <- induced_subgraph(g, V(g)[comps$membership == which.max(comps$csize)])
  

get_clique_density_stats(g_lcc)

g_lcc_embed_result <- embedding_function(g_lcc, 10)

hist(apply(g_lcc[], 2, sum))
sort(apply(g_lcc[], 2, sum))

plot(g_lcc_embed_result$Xhat[,1], g_lcc_embed_result$Xhat[,2], 
     main=well, xlab="Dim 1", ylab="Dim 2", pch=16, col="blue")




library(igraph)

# 1. Define the parameters from your empirical graph
n_core <- 72
n_periph <- 164
total_nodes <- n_core + n_periph

# The block connectivity (preference) matrix
# Based on your C/P density matrix
pref_matrix <- matrix(
  c(0.9510955, 0.2653286, 
    0.2653286, 0.03052521), 
  nrow = 2, 
  byrow = TRUE
)

# Block sizes
block_sizes <- c(n_core, n_periph)

# 2. Generate the standard SBM random graph
# We use a seed so you can reproduce these exact results
set.seed(123) 
g_sbm_standard <- sample_sbm(
  n = total_nodes, 
  pref.matrix = pref_matrix, 
  block.sizes = block_sizes, 
  directed = FALSE
)

# 3. Extract the degree distribution
sbm_degrees <- degree(g_sbm_standard)

# Sort the degrees so you can easily see the clustering 
# and compare it to the smooth spread of your real data
sorted_sbm_degrees <- sort(sbm_degrees, decreasing = TRUE)

# Print the result
print(sorted_sbm_degrees)

hist(sorted_sbm_degrees)



A_mat <- as.matrix(as_adjacency_matrix(g_lcc))
n <- nrow(A_mat)
degrees <- rowSums(A_mat)

# 1. Get Memberships (X matrix)
kc <- coreness(g_lcc)
max_kc <- max(kc)
core_index <- which(kc == max_kc)
periph_index <- setdiff(1:n, core_index)

X <- matrix(0, nrow = n, ncol = 2)
X[core_index, 1] <- 1
X[periph_index, 2] <- 1

# 2. The Connectivity Matrix (B matrix)
# Plugging in the results from your get_clique_density_stats function
B <- matrix(
  c(0.9510955, 0.2653286, 
    0.2653286, 0.03052521), 
  nrow = 2, 
  byrow = TRUE
)

# 3. Estimate Degree Correction Parameters (Theta matrix)
# Calculate the mean degree within each block
mean_deg_core <- mean(degrees[core_index])
mean_deg_periph <- mean(degrees[periph_index])

# Initialize and assign theta values
theta_vec <- numeric(n)
theta_vec[core_index] <- degrees[core_index] / mean_deg_core
theta_vec[periph_index] <- degrees[periph_index] / mean_deg_periph

Theta <- diag(theta_vec)

# 4. Generate the Expected Probability Matrix (P)
# P = Theta * X * B * X^T * Theta
P_matrix <- Theta %*% X %*% B %*% t(X) %*% Theta

# 5. Clean up the P matrix for graph generation
# Cap probabilities at 1 (extreme hubs might slightly exceed 1)
P_matrix[P_matrix > 1] <- 1
# Remove self-loops
diag(P_matrix) <- 0


# Assuming B, theta_vec, core_index, and periph_index are already in your environment

# 1. Eigendecomposition of the B matrix to find the base latent positions
eig_B <- eigen(B)

# Note: Your B matrix has a negative determinant, meaning it has one negative eigenvalue. 
# Technically, this forms a Generalized Random Dot Product Graph (GRDPG). 
# We use the absolute values of the eigenvalues to construct the 2D embedding space coordinates.
Z_base <- eig_B$vectors %*% diag(sqrt(abs(eig_B$values)))

# 2. Assign base latent positions to all nodes
n <- length(theta_vec)
Y_base <- matrix(0, nrow = n, ncol = 2)

# Z_base[1,] is the base position for Core; Z_base[2,] is for Periphery
for (i in core_index) { Y_base[i, ] <- Z_base[1, ] }
for (i in periph_index) { Y_base[i, ] <- Z_base[2, ] }

# 3. Apply the degree corrections (Theta) to scale the positions
# Y = Theta * Y_base
Y_scaled <- diag(theta_vec) %*% Y_base

# 4. Plotting the results
plot(Y_scaled[, 1], Y_scaled[, 2], 
     pch = 16, col = rgb(0.5, 0.5, 0.5, alpha = 0.5), # semi-transparent gray dots
     xlab = "Latent Dimension 1", ylab = "Latent Dimension 2",
     main = "Theoretical ASE of DC-2-SBM",
     xlim = range(0, Y_scaled[, 1]) * 1.1, 
     ylim = range(0, Y_scaled[, 2]) * 1.1)

# 5. Draw the rays from the origin
max_theta_core <- max(theta_vec[core_index])
max_theta_periph <- max(theta_vec[periph_index])

segments(x0 = 0, y0 = 0, 
         x1 = Z_base[1, 1] * max_theta_core, y1 = Z_base[1, 2] * max_theta_core, 
         col = "blue", lty = 2, lwd = 2) # Core Ray
segments(x0 = 0, y0 = 0, 
         x1 = Z_base[2, 1] * max_theta_periph, y1 = Z_base[2, 2] * max_theta_periph, 
         col = "green", lty = 2, lwd = 2) # Periphery Ray

# 6. Add the two base latent positions (where theta = 1) as red circles
points(Z_base[, 1], Z_base[, 2], col = "red", pch = 1, cex = 2.5, lwd = 3)

#points(-P_embed_result$Xhat[,1], P_embed_result$Xhat[,2]-0.6)


eig_P <- eigen(P_matrix, symmetric = TRUE)

# We take the top 2 eigenvalues by absolute magnitude 
# (This is standard practice since B might have a negative eigenvalue, creating a GRDPG)
top_indices <- order(abs(eig_P$values), decreasing = TRUE)[1:2]
evals_P <- eig_P$values[top_indices]
evecs_P <- eig_P$vectors[, top_indices]

# 2. Construct the ASE for P
Y_ASE <- evecs_P %*% diag(sqrt(abs(evals_P)))

# 3. Orthogonal Procrustes Alignment
# Align Y_ASE to our theoretical Y_scaled to handle the rotational ambiguity
svd_align <- svd(t(Y_ASE) %*% Y_scaled)
W <- svd_align$u %*% t(svd_align$v)
Y_ASE_aligned <- Y_ASE %*% W

# 4. Add the aligned ASE of P to the existing plot
# We use black crosses (pch = 4) to distinguish them
points(Y_ASE_aligned[, 1], Y_ASE_aligned[, 2], 
       col = "black", pch = 4, cex = 0.8, lwd = 1.5)



# 1 the exact rank-2 P matrix
P_pure <- Theta %*% X %*% B %*% t(X) %*% Theta

# 2. Perform the eigendecomposition on P_pure instead of the hollowed P_matrix
eig_P_pure <- eigen(P_pure)

# 3. Extract top 2 eigenvectors/eigenvalues
top_indices <- order(abs(eig_P_pure$values), decreasing = TRUE)[1:2]
evals_P <- eig_P_pure$values[top_indices]
evecs_P <- eig_P_pure$vectors[, top_indices]

# 4. Construct the pure ASE
Y_ASE_pure <- evecs_P %*% diag(sqrt(abs(evals_P)))

# Instead of Orthogonal Procrustes via SVD, we use the General Linear Solution
# We want W such that Y_ASE_pure %*% W = Y_scaled
# The exact solution is: W = (Y_ASE^T * Y_ASE)^(-1) * Y_ASE^T * Y_scaled

W_exact <- solve(t(Y_ASE_pure) %*% Y_ASE_pure) %*% t(Y_ASE_pure) %*% Y_scaled

# Apply the exact transformation
Y_ASE_aligned_exact <- Y_ASE_pure %*% W_exact

# Plot these as purple diamonds
points(Y_ASE_aligned_exact[, 1], Y_ASE_aligned_exact[, 2], 
       col = "purple", pch = 18, cex = 1.5)



n <-  nrow(P_matrix)
U <- matrix(0, nrow = n, ncol = n)
U[col(U) > row(U)] <- runif(n*(n-1)/2)
U <- (U + t(U))
diag(U) <- runif(n)
A <- (U < P_matrix) + 0 ;
diag(A) <- 0

A_embed_result <- embedding_function(A, 10)
A_embed_result$Xhat

Y_emp <- as.matrix(A_embed_result$Xhat)[,1:2]
Y_emp_aligned <- align_grdpg_o11(Y_emp, Y_scaled)$Y_aligned
points(Y_emp_aligned[,1],Y_emp_aligned[,2])

g_lcc_embed_result <- embedding_function(g_lcc, 10)
Y_emp <- as.matrix(g_lcc_embed_result$Xhat)[,1:2]
original_Y_aligned = align_grdpg_o11(Y_emp, Y_scaled)$Y_aligned

points(original_Y_aligned[,1],original_Y_aligned[,2],col = 'blue')


