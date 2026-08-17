library(igraph)
library(irlba)
library(Matrix)

# =====================================================================
# 0. HELPER FUNCTIONS
# =====================================================================

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

full.ase <- function(A, d, diagaug=TRUE, doptr=FALSE) {
  require(irlba)
  if (doptr) { g <- ptr(A); A <- g[] } else { A <- A[] }
  A = as.matrix(A)
  if (diagaug) { diag(A) <- rowSums(A) / (nrow(A)-1) }
  A.svd <- irlba(A, d)
  Xhat <- A.svd$u %*% diag(sqrt(A.svd$d))
  Xhat.R <- NULL
  if (!isSymmetric(A)) { Xhat.R <- A.svd$v %*% diag(sqrt(A.svd$d)) }
  return(list(eval=A.svd$d, Xhat=Matrix(Xhat), Xhat.R=Xhat.R))
}

align_grdpg_o11 <- function(Y_source, Y_target) {
  loss_func <- function(x, sign_type, Y_s, Y_t) {
    cx <- cosh(x); sx <- sinh(x)
    if (sign_type == 1) W <- matrix(c(cx, sx, sx, cx), 2, 2, byrow=TRUE)
    else if (sign_type == 2) W <- matrix(c(-cx, -sx, sx, cx), 2, 2, byrow=TRUE)
    else if (sign_type == 3) W <- matrix(c(cx, sx, -sx, -cx), 2, 2, byrow=TRUE)
    else W <- matrix(c(-cx, -sx, -sx, -cx), 2, 2, byrow=TRUE)
    return(sum((Y_s %*% W - Y_t)^2))
  }
  
  best_loss <- Inf; best_W <- NULL
  for (st in 1:4) {
    opt <- optimize(loss_func, interval = c(-10, 10), sign_type = st, Y_s = Y_source, Y_t = Y_target)
    if (opt$objective < best_loss) {
      best_loss <- opt$objective
      cx <- cosh(opt$minimum); sx <- sinh(opt$minimum)
      if (st == 1) best_W <- matrix(c(cx, sx, sx, cx), 2, 2, byrow=TRUE)
      if (st == 2) best_W <- matrix(c(-cx, -sx, sx, cx), 2, 2, byrow=TRUE)
      if (st == 3) best_W <- matrix(c(cx, sx, -sx, -cx), 2, 2, byrow=TRUE)
      if (st == 4) best_W <- matrix(c(-cx, -sx, -sx, -cx), 2, 2, byrow=TRUE)
    }
  }
  return(list(W = best_W, Y_aligned = as.matrix(Y_source %*% best_W)))
}

# =====================================================================
# 1. PARAMETER ESTIMATION FROM EMPIRICAL GRAPH
# =====================================================================

A_mat <- as.matrix(as_adjacency_matrix(g_lcc))
n <- nrow(A_mat)
degrees <- rowSums(A_mat)

# Memberships (X matrix)
kc <- coreness(g_lcc)
core_index <- which(kc == max(kc))
periph_index <- setdiff(1:n, core_index)

X <- matrix(0, nrow = n, ncol = 2)
X[core_index, 1] <- 1
X[periph_index, 2] <- 1

# Connectivity Matrix (B) - Dynamically calculated!
B <- get_clique_density_stats(g_lcc)$matrix

# Degree Corrections (Theta)
theta_vec <- numeric(n)
theta_vec[core_index] <- degrees[core_index] / mean(degrees[core_index])
theta_vec[periph_index] <- degrees[periph_index] / mean(degrees[periph_index])
Theta <- diag(theta_vec)

# Construct Pure and Hollow P matrices
P_pure <- Theta %*% X %*% B %*% t(X) %*% Theta

P_matrix <- P_pure

summary(as.vector(P_matrix))
P_matrix[P_matrix > 1] <- 1
diag(P_matrix) <- 0

# =====================================================================
# 2. GENERATE LATENT POSITIONS & EMBEDDINGS
# =====================================================================

# A. Theoretical Geometry
eig_B <- eigen(B)
Z_base <- eig_B$vectors %*% diag(sqrt(abs(eig_B$values)))

Y_scaled <- matrix(0, nrow = n, ncol = 2)
for (i in core_index) Y_scaled[i, ] <- Z_base[1, ] * theta_vec[i]
for (i in periph_index) Y_scaled[i, ] <- Z_base[2, ] * theta_vec[i]

# B. ASE of P_pure (Exact Rank-2, No Diagonal Augmentation)
ase_pure <- full.ase(P_pure, d=2, diagaug=FALSE)
Y_pure <- as.matrix(ase_pure$Xhat)
# Use exact general linear mapping for the noiseless matrix to avoid optimization artifacts
W_exact <- solve(t(Y_pure) %*% Y_pure) %*% t(Y_pure) %*% Y_scaled
Y_pure_aligned <- Y_pure %*% W_exact

# C. ASE of P_hollow (Diagonal Augmentation applied)
ase_hollow <- full.ase(P_matrix, d=2, diagaug=TRUE)
Y_hollow_aligned <- align_grdpg_o11(as.matrix(ase_hollow$Xhat), Y_scaled)$Y_aligned

# D. ASE of Simulated Graph (A_sim)
U <- matrix(0, nrow = n, ncol = n)
U[col(U) > row(U)] <- runif(n*(n-1)/2)
U <- (U + t(U))
diag(U) <- runif(n)
A_sim <- (U < P_matrix) + 0 
diag(A_sim) <- 0

ase_sim <- full.ase(A_sim, d=2, diagaug=TRUE)
Y_sim_aligned <- align_grdpg_o11(as.matrix(ase_sim$Xhat), Y_scaled)$Y_aligned

# E. ASE of Original Empirical Graph (g_lcc)
ase_emp <- full.ase(A_mat, d=2, diagaug=TRUE)
Y_emp_aligned <- align_grdpg_o11(as.matrix(ase_emp$Xhat), Y_scaled)$Y_aligned

# =====================================================================
# 3. PLOTTING
# =====================================================================

# Calculate dynamic bounding box to ensure all points fit
all_x <- c(0, Y_scaled[,1], Y_pure_aligned[,1], Y_hollow_aligned[,1], Y_sim_aligned[,1], Y_emp_aligned[,1], Z_base[,1])
all_y <- c(0, Y_scaled[,2], Y_pure_aligned[,2], Y_hollow_aligned[,2], Y_sim_aligned[,2], Y_emp_aligned[,2], Z_base[,2])

xlims <- range(all_x) + c(-0.1, 0.1) * diff(range(all_x))
ylims <- range(all_y) + c(-0.1, 0.1) * diff(range(all_y))

# 1. Initialize Plot with Theoretical Dots
plot(Y_scaled[, 1], Y_scaled[, 2], 
     pch = 16, col = "black", cex = 1,
     xlab = "ASE Dimension 1", ylab = "ASE Dimension 2",
     main = "DC-2-SBM: ASE Embeddings",
     xlim = xlims, ylim = ylims)

# 2. Draw Theoretical Rays
max_t_core <- max(theta_vec[core_index])
max_t_periph <- max(theta_vec[periph_index])
segments(0, 0, Z_base[1, 1] * max_t_core, Z_base[1, 2] * max_t_core, col = "black", lty = 2, lwd = 2)
segments(0, 0, Z_base[2, 1] * max_t_periph, Z_base[2, 2] * max_t_periph, col = "black", lty = 2, lwd = 2)

# 3. Add Base Latent Positions (Red Circles)
points(Z_base[, 1], Z_base[, 2], col = "red", pch = 1, cex = 2.5, lwd = 3)

# 4. Add ASE of P_pure (Open Magenta Squares - Will perfectly frame the grey dots)
#points(Y_pure_aligned[, 1], Y_pure_aligned[, 2], col = "magenta", pch = 0, cex = 1.3, lwd = 1.5)

# 5. Add ASE of P_hollow (Black Crosses - Shows the diagonal bending artifact)
points(Y_hollow_aligned[, 1], Y_hollow_aligned[, 2], col = "black", pch = 4, cex = 0.8, lwd = 1.5)

# 6. Add ASE of Simulated Graph (Orange Triangles, Semi-transparent)
#points(Y_sim_aligned[, 1], Y_sim_aligned[, 2], col = rgb(1, 0.5, 0, alpha=0.5), pch = 17, cex = 1)

# 7. Add ASE of Original Empirical Graph (Blue Dots, Semi-transparent)
points(Y_emp_aligned[, 1], Y_emp_aligned[, 2], col = rgb(0, 0, 1, alpha=0.6), pch = 16, cex = 1)

# 8. Legend
legend("bottomleft", 
       legend = c("Theta_hat * X_hat * |B_hat|^{1/2}", 
                  "ASE of normalized Phat (cap and hollow)", 
                  "ASE of original Graph", 
                  "Base latent positions"),
       col = c("black", "black", rgb(0, 0, 1, alpha=0.6), "red"),
       pch = c(16, 4, 16, 1),
       pt.cex = c(1, 0.8, 1, 2.5),
       pt.lwd = c(1, 1.5, 1, 3),
       bg = "white", cex=0.85)

