# 1. Subset the data for the current file and well
subset_data <- subset(edge_list, File == filenames[1] & Well == 'well001')

n_rows <- n_cols <- unique(subset_data$dim)
adj <- sparseMatrix(
    i    = subset_data$Row + 1,
    j    = subset_data$Column + 1,
    dims = c(n_rows, n_cols)
  )
  
  # Create the graph and find the largest connected component (LCC)
  g <- graph_from_adjacency_matrix(adj, mode = "undirected", diag = FALSE)
  comps <- components(g, mode = "weak")
  big <- which.max(comps$csize)
  g_lcc <- induced_subgraph(g, V(g)[comps$membership == big])
    # --- Proceed with ASE only if LCC size is 11 or greater ---
    
    # Perform Adjacency Spectral Embedding (ASE) on the LCC
ase <- full.ase(g_lcc, 10)
get_clique_density_stats(g_lcc)

lse <- full.lse(g_lcc, 10)


scatterplot3d(
  x = ase$Xhat[, 1], 
  y = ase$Xhat[, 2], 
  z = ase$Xhat[, 3],
  main = well,
  xlab = "Dim 1",
  ylab = "Dim 2",
  zlab = "Dim 3",
  color = 'blue',
  pch = 16,
  type = "p",
  grid = TRUE,
  box = TRUE
)


scatterplot3d(
  x = lse$Xhat[, 1], 
  y = lse$Xhat[, 2], 
  z = lse$Xhat[, 3],
  main = well,
  xlab = "Dim 1",
  ylab = "Dim 2",
  zlab = "Dim 3",
  color = 'blue',
  pch = 16,
  type = "p",
  grid = TRUE,
  box = TRUE
)


# 1. Create a color function (Blue = Low Degree, Red = High Degree)
rbPal <- colorRampPalette(c('blue', 'red'))

# 2. Map degrees to colors (Using log scale to handle the heavy tail)
# This divides the log-degrees into 100 bins and assigns a color to each
deg_colors <- rbPal(100)[as.numeric(cut(log(deg.seq), breaks = 100))]

# 3. Plot
scatterplot3d(
  x = lse$Xhat[, 1], 
  y = lse$Xhat[, 2], 
  z = lse$Xhat[, 3],
  main = well,
  xlab = "Dim 1",
  ylab = "Dim 2",
  zlab = "Dim 3",
  color = deg_colors,  # <--- Apply the color vector here
  pch = 16,
  type = "p",
  grid = TRUE,
  box = TRUE
)

# Optional: Add a legend to interpret the colors
legend("topleft", 
       legend = c("Low Degree", "High Degree"),
       col = c("blue", "red"), 
       pch = 16, 
       bty = "n")

library(igraph)

kc <- coreness(g_lcc)
core_idx <- which(kc == max(kc))
periph_idx <- which(kc < max(kc))

ord <- c(core_idx, periph_idx)
A_sorted <- as.matrix(as_adjacency_matrix(g_lcc))[ord, ord]

par(mfrow = c(1, 1), mar = c(2, 2, 2, 2))
image(t(A_sorted)[, nrow(A_sorted):1], 
      axes = FALSE, 
      col = c("white", "black"), 
      main = "Adjacency Matrix: Core vs. Periphery")

cut_frac <- length(core_idx) / nrow(A_sorted)
abline(v = cut_frac, col = "red", lwd = 2)
abline(h = 1 - cut_frac, col = "red", lwd = 2)

deg_sorted <- rowSums(A_sorted)
deg_root <- sqrt(deg_sorted)
deg_root[deg_root == 0] <- 1 

B_sorted <- A_sorted / outer(deg_root, deg_root)
ss = irlba(B_sorted,10)
ss$d
non_zero_values <- B_sorted[B_sorted > 0]
saturation_point <- quantile(non_zero_values, 0.95) # Cap at the 95th percentile of actual edges

# Create the plot with the breaks
par(mfrow = c(1, 1), mar = c(2, 2, 2, 2))

# We use 'breaks' to force the color scale to saturate early
image(t(B_sorted)[, nrow(B_sorted):1], 
      axes = FALSE, 
      col = grey(seq(1, 0, length.out = 100)), 
      breaks = c(seq(0, saturation_point, length.out = 100), max(B_sorted)),
      main = "Normalized B (Contrast Adjusted)")

# Add the red structure lines
cut_frac <- length(which(kc == max(kc))) / nrow(B_sorted)
abline(v = cut_frac, col = "red", lwd = 2)
abline(h = 1 - cut_frac, col = "red", lwd = 2)


A <- g_lcc[]
A = as.matrix(A)
deg.seq <- rowSums(A)
sort(deg.seq)
B <- A/sqrt(outer(deg.seq,deg.seq))

L.svd <- irlba(B,10)

Rank1_B = L.svd$u[,1] %*% t(L.svd$u[,1])



(B - Rank1_B)/B

Xhat <- L.svd$u %*% diag(sqrt(L.svd$d))    


aase 


plot(L.svd$d)
ase$eval
plot(ase$Xhat[,1],ase$Xhat[,2])
hist(Xhat[,3])

plot(1:nrow(ase$Xhat),sort(ase$Xhat[,1]))

# 1. Plot Dim 2 vs Dim 3 directly to see the hidden structure
# (Ignoring the dominant Dim 1 "degree" effect)
plot(Xhat[,2], Xhat[,3], main="Latent Positions: Dim 2 vs Dim 3")

# 2. Check localization: Do a few nodes hog all the mass?
# If the top few values are huge and the rest are 0, it's localization.
plot(sort(L.svd$u[,2]^2), main="Squared entries of Eigenvector 2")



# 1. Identify the culprits (the indices of the two nodes)
bad_nodes <- order(abs(L.svd$u[,2]), decreasing = TRUE)[1:2]
# Or check both tails if one is negative:


bad_nodes <- which(L.svd$u[,2]^2 > 0.1)

print(bad_nodes)

deg.seq[bad_nodes]

# 2. Remove them from your adjacency matrix
A_clean <- A[-bad_nodes, -bad_nodes]

# 3. Re-run your SVD on the clean matrix
deg.seq.clean <- rowSums(A_clean)
B_clean <- A_clean / sqrt(outer(deg.seq.clean, deg.seq.clean))
L.svd.clean <- irlba(B_clean, 10)
L.svd.clean$d
Xhat.clean <- L.svd.clean$u %*% diag(sqrt(L.svd.clean$d))

# 4. Plot again
scatterplot3d(Xhat.clean[,1], Xhat.clean[,2], Xhat.clean[,3], main="Cleaned Graph")

    