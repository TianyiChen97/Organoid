# Organoid Graph Analysis — Complete Reference

This document contains all reusable helper functions and the analysis pipeline for organoid MEA graph data. It is self-contained: any agent or collaborator can reproduce the full analysis on a new dataset using only this file.

---

## Required Libraries

```r
library(Matrix)
library(igraph)
library(irlba)
library(ggplot2)
library(scatterplot3d)
library(stringr)
library(readr)
library(patchwork)
```

---

## 1. Helper Functions

### 1.1 Adjacency Spectral Embedding (ASE)

Computes the truncated SVD of the (optionally diagonal-augmented) adjacency matrix and returns the scaled left singular vectors as the embedding.

- `A`: adjacency matrix or igraph object (via `A[]`)
- `d`: number of dimensions to embed into
- `diagaug`: if `TRUE`, replaces diagonal with `rowSums(A) / (n-1)` (default)
- Returns: `eval` (singular values), `Xhat` (left embedding), `Xhat.R` (right embedding, non-NULL only for asymmetric/directed graphs)

```r
full.ase <- function(A, d, diagaug=TRUE, doptr=FALSE) {
  require(irlba)
  if (doptr) { g <- ptr(A); A <- g[] } else { A <- A[] }
  A <- as.matrix(A)
  if (diagaug) { diag(A) <- rowSums(A) / (nrow(A) - 1) }
  A.svd <- irlba(A, d)
  Xhat <- A.svd$u %*% diag(sqrt(A.svd$d))
  Xhat.R <- NULL
  if (!isSymmetric(A)) { Xhat.R <- A.svd$v %*% diag(sqrt(A.svd$d)) }
  return(list(eval = A.svd$d, Xhat = Matrix(Xhat), Xhat.R = Xhat.R))
}
```

### 1.2 Laplacian Spectral Embedding (LSE)

Computes the truncated SVD of the normalized Laplacian `D^{-1/2} A D^{-1/2}`.

- `A`: adjacency matrix or igraph object
- `d`: number of dimensions
- Isolated nodes (degree 0) have their degree set to 1 to avoid division by zero.

```r
full.lse <- function(A, d) {
  A <- A[]
  A <- as.matrix(A)
  deg.seq <- rowSums(A)
  deg.seq[deg.seq == 0] <- 1
  L.svd <- irlba(A / sqrt(outer(deg.seq, deg.seq)), d)
  Xhat <- L.svd$u %*% diag(sqrt(L.svd$d))
  Xhat.R <- NULL
  if (!isSymmetric(A)) { Xhat.R <- L.svd$v %*% diag(sqrt(L.svd$d)) }
  return(list(eval = L.svd$d, Xhat = Xhat, Xhat.R = Xhat.R))
}
```

### 1.3 Elbow Selection (`getElbows`)

Selects up to `n` elbow points from a scree plot (eigenvalues or singular values) using a Gaussian mixture likelihood criterion.

- `dat`: vector of eigenvalues/singular values, or a matrix (uses column SDs)
- `n`: number of elbows to find (default 3)
- Returns: vector of elbow indices (e.g., `c(2, 5, 8)`)

**Usage convention in this project:** use `getElbows(embed_result$eval, plot=FALSE)[2]` (the second elbow) as the embedding dimension.

```r
getElbows <- function(dat, n = 3, threshold = FALSE, plot = TRUE, main = "") {
  if (is.matrix(dat)) {
    d <- sort(apply(dat, 2, sd), decreasing = TRUE)
  } else {
    d <- sort(dat, decreasing = TRUE)
  }
  if (!is.logical(threshold)) d <- d[d > threshold]
  p <- length(d)
  if (p == 0) stop("d must have elements larger than threshold")
  lq <- rep(0.0, p)
  for (q in 1:p) {
    mu1 <- mean(d[1:q])
    mu2 <- mean(d[-(1:q)])
    sigma2 <- (sum((d[1:q] - mu1)^2) + sum((d[-(1:q)] - mu2)^2)) / (p - 1 - (q < p))
    lq[q] <- sum(dnorm(d[1:q], mu1, sqrt(sigma2), log = TRUE)) +
             sum(dnorm(d[-(1:q)], mu2, sqrt(sigma2), log = TRUE))
  }
  q <- which.max(lq)
  if (n > 1 && q < (p - 1)) {
    q <- c(q, q + getElbows(d[(q + 1):p], n - 1, plot = FALSE))
  }
  if (plot == TRUE) {
    if (is.matrix(dat)) {
      sdv <- d
      plot(sdv, type = "b", xlab = "dim", ylab = "stdev", main = main)
      points(q, sdv[q], col = 2, pch = 19)
    } else {
      plot(dat, type = "b", main = main)
      points(q, dat[q], col = 2, pch = 19)
    }
  }
  return(q)
}
```

### 1.4 Undirected Core–Periphery Statistics (`get_clique_density_stats`)

Identifies core nodes via **k-core decomposition** on the undirected graph (nodes with maximum coreness = core). Computes the 2×2 block density matrix.

- Input `A`: adjacency matrix of the LCC (undirected, binary)
- Core identification: `coreness(g)` → nodes where `kc == max(kc)`
- Edge counts divided by 2 for undirected graphs
- Returns:
  - `matrix`: 2×2 symmetric density matrix `[[CC, CP], [CP, PP]]`
  - `n_core`, `n_periph`: block sizes

```r
get_clique_density_stats <- function(A) {
  A_mat <- as.matrix(A)
  g <- graph_from_adjacency_matrix(A_mat, mode = "undirected", diag = FALSE)

  kc <- coreness(g)
  core_index <- which(kc == max(kc))
  periph_index <- setdiff(1:nrow(A_mat), core_index)

  n_c <- length(core_index)
  n_p <- length(periph_index)

  edges_CC <- if (n_c > 1) sum(A_mat[core_index, core_index]) / 2 else 0
  edges_PP <- if (n_p > 1) sum(A_mat[periph_index, periph_index]) / 2 else 0
  edges_CP <- if (n_c > 0 && n_p > 0) sum(A_mat[core_index, periph_index]) else 0

  max_CC <- n_c * (n_c - 1) / 2
  max_PP <- n_p * (n_p - 1) / 2
  max_CP <- n_c * n_p

  dens_CC <- ifelse(max_CC > 0, edges_CC / max_CC, 0)
  dens_PP <- ifelse(max_PP > 0, edges_PP / max_PP, 0)
  dens_CP <- ifelse(max_CP > 0, edges_CP / max_CP, 0)

  M <- matrix(c(dens_CC, dens_CP, dens_CP, dens_PP), nrow = 2, byrow = TRUE)
  colnames(M) <- c("C", "P"); rownames(M) <- c("C", "P")

  return(list(matrix = M, n_core = n_c, n_periph = n_p))
}
```

### 1.5 Directed Core–Periphery Statistics (`get_directed_clique_density_stats`)

Identifies core nodes via **k-core decomposition with `mode="all"`** (considering both in-degree and out-degree). Computes the 2×2 asymmetric block density matrix and block-level reciprocity statistics.

- Input `A`: adjacency matrix of the LCC (directed, binary)
- Core identification: `coreness(g, mode = "all")` → nodes where `kc == max(kc)`
- Edge counts are **not** divided by 2 (directed edges)
- Reciprocity `rho = sum(A * t(A)) / sum(A)` measures the fraction of edges that are reciprocated
- Returns:
  - `matrix`: 2×2 asymmetric density matrix `[[C→C, C→P], [P→C, P→P]]`
  - `n_core`, `n_periph`: block sizes
  - `rho_global`: global reciprocity
  - `rho_CC`, `rho_PP`: within-block reciprocity
  - `rho_CP`, `rho_PC`: inter-block reciprocity (same value, computed from both directions)

```r
get_directed_clique_density_stats <- function(A) {
  A_mat <- as.matrix(A)
  g <- graph_from_adjacency_matrix(A_mat, mode = "directed", diag = FALSE)

  kc <- coreness(g, mode = "all")
  core_index <- which(kc == max(kc))
  periph_index <- setdiff(1:nrow(A_mat), core_index)

  n_c <- length(core_index)
  n_p <- length(periph_index)

  Acc <- if (n_c > 0) A_mat[core_index, core_index, drop = FALSE] else matrix(0, 0, 0)
  App <- if (n_p > 0) A_mat[periph_index, periph_index, drop = FALSE] else matrix(0, 0, 0)
  Acp <- if (n_c > 0 && n_p > 0) A_mat[core_index, periph_index, drop = FALSE] else matrix(0, 0, 0)
  Apc <- if (n_c > 0 && n_p > 0) A_mat[periph_index, core_index, drop = FALSE] else matrix(0, 0, 0)

  edges_CC <- sum(Acc); edges_PP <- sum(App)
  edges_CP <- sum(Acp); edges_PC <- sum(Apc)
  total_edges <- edges_CC + edges_PP + edges_CP + edges_PC

  max_CC <- n_c * (n_c - 1); max_PP <- n_p * (n_p - 1)
  max_CP <- n_c * n_p; max_PC <- n_p * n_c

  dens_CC <- ifelse(max_CC > 0, edges_CC / max_CC, 0)
  dens_PP <- ifelse(max_PP > 0, edges_PP / max_PP, 0)
  dens_CP <- ifelse(max_CP > 0, edges_CP / max_CP, 0)
  dens_PC <- ifelse(max_PC > 0, edges_PC / max_PC, 0)

  # Global reciprocity
  rho_global <- 0
  if (total_edges > 0) rho_global <- sum(A_mat * t(A_mat)) / total_edges

  # Block reciprocities
  rho_CC <- 0; if (edges_CC > 0) rho_CC <- sum(Acc * t(Acc)) / edges_CC
  rho_PP <- 0; if (edges_PP > 0) rho_PP <- sum(App * t(App)) / edges_PP

  # Inter-block reciprocity
  rho_inter <- 0
  edges_inter <- edges_CP + edges_PC
  if (edges_inter > 0 && n_c > 0 && n_p > 0) {
    matches_upper_right <- sum(Acp * t(Apc))
    matches_lower_left  <- sum(Apc * t(Acp))
    rho_inter <- (matches_upper_right + matches_lower_left) / edges_inter
  }

  M <- matrix(c(dens_CC, dens_CP, dens_PC, dens_PP), nrow = 2, byrow = TRUE)
  colnames(M) <- c("To_Core", "To_Periph"); rownames(M) <- c("From_Core", "From_Periph")

  return(list(
    matrix = M, n_core = n_c, n_periph = n_p,
    rho_global = rho_global,
    rho_CC = rho_CC, rho_PP = rho_PP,
    rho_CP = rho_inter, rho_PC = rho_inter
  ))
}
```

### 1.6 Out-Degree Core Statistics (`get_outdegree_core_stats`)

Alternative core identification: instead of k-core decomposition, defines the core as the top `top_n` nodes ranked by out-degree. Useful for comparing against the k-core definition.

```r
get_outdegree_core_stats <- function(A, top_n = 10) {
  A_mat <- as.matrix(A)
  n_total <- nrow(A_mat)
  out_degree <- rowSums(A_mat)
  k <- min(top_n, n_total)
  core_index <- order(out_degree, decreasing = TRUE)[1:k]
  periph_index <- setdiff(1:n_total, core_index)
  n_c <- length(core_index); n_p <- length(periph_index)

  edges_CC <- if (n_c > 0) sum(A_mat[core_index, core_index]) else 0
  edges_CP <- if (n_c > 0 && n_p > 0) sum(A_mat[core_index, periph_index]) else 0
  edges_PC <- if (n_c > 0 && n_p > 0) sum(A_mat[periph_index, core_index]) else 0
  edges_PP <- if (n_p > 0) sum(A_mat[periph_index, periph_index]) else 0

  max_CC <- n_c * (n_c - 1); max_PP <- n_p * (n_p - 1)
  max_CP <- n_c * n_p; max_PC <- n_p * n_c

  dens_CC <- ifelse(max_CC > 0, edges_CC / max_CC, 0)
  dens_PP <- ifelse(max_PP > 0, edges_PP / max_PP, 0)
  dens_CP <- ifelse(max_CP > 0, edges_CP / max_CP, 0)
  dens_PC <- ifelse(max_PC > 0, edges_PC / max_PC, 0)

  M <- matrix(c(dens_CC, dens_CP, dens_PC, dens_PP), nrow = 2, byrow = TRUE)
  colnames(M) <- c("To_Core", "To_Periph"); rownames(M) <- c("From_Core", "From_Periph")

  return(list(matrix = M, n_core = n_c, n_periph = n_p, core_indices = core_index))
}
```

### 1.7 Path Extraction Helper

Extracts `Mxxxxx/Type/xxxxxx` (e.g., `M07914/Network/000297`) from full file paths in the edge list.

```r
extract_network_path <- function(input_string) {
  pattern <- "M\\d+/[A-Za-z]+/\\d+"
  result <- str_extract(input_string, pattern)
  return(result)
}
```

---

## 2. Analysis Pipeline

### 2.1 Input Data Format

The edge list CSV must have these columns:

| Column | Description |
|--------|-------------|
| `File` | Full path identifying the recording session (contains `Mxxxxx/Type/xxxxxx`) |
| `Well` | Well identifier, e.g. `well000` through `well005` |
| `Row` | 0-indexed row of edge source |
| `Column` | 0-indexed column of edge target |
| `dim` | Dimension of the adjacency matrix (number of electrodes) |

### 2.2 Graph Construction & LCC

For each `(File, Well)` pair:

```r
subset_data <- subset(edge_list, File == file & Well == well)
n_rows <- n_cols <- unique(subset_data$dim)

# Build sparse adjacency (note: 0-indexed → 1-indexed)
adj <- sparseMatrix(
  i = subset_data$Row + 1,
  j = subset_data$Column + 1,
  dims = c(n_rows, n_cols)
)

# Undirected graph + LCC
g <- graph_from_adjacency_matrix(adj, mode = "undirected", diag = FALSE)
comps <- components(g, mode = "weak")
g_lcc <- induced_subgraph(g, V(g)[comps$membership == which.max(comps$csize)])

# Directed graph + LCC (same LCC logic, just directed construction)
g_dir <- graph_from_adjacency_matrix(adj, mode = "directed", diag = FALSE)
comps_dir <- components(g_dir, mode = "weak")
g_lcc_dir <- induced_subgraph(g_dir, V(g_dir)[comps_dir$membership == which.max(comps_dir$csize)])

# Skip analysis if LCC too small
if (vcount(g_lcc) < 11) next
```

### 2.3 Embedding Analysis

```r
# Choose method
EMBEDDING_METHOD <- "ASE"  # or "LSE"
embedding_function <- if (EMBEDDING_METHOD == "LSE") full.lse else full.ase

# Embed into 10 dimensions, then select via elbow
embed_result <- embedding_function(g_lcc, 10)
elb <- getElbows(embed_result$eval, plot = FALSE)
embedding_dim <- elb[2]  # second elbow

# Plot
if (embedding_dim == 2) {
  plot(embed_result$Xhat[, 1], embed_result$Xhat[, 2],
       pch = 16, col = "blue", xlab = "Dim 1", ylab = "Dim 2")
} else if (embedding_dim >= 3) {
  scatterplot3d(embed_result$Xhat[, 1], embed_result$Xhat[, 2], embed_result$Xhat[, 3],
                color = "blue", pch = 16, xlab = "D1", ylab = "D2", zlab = "D3")
}
```

### 2.4 Undirected CP Analysis

```r
A_lcc <- as.matrix(as_adjacency_matrix(g_lcc))
stats <- get_clique_density_stats(A_lcc)

# stats$matrix    — 2×2 density matrix (CC, CP, PP)
# stats$n_core    — number of core nodes
# stats$n_periph  — number of periphery nodes
```

### 2.5 Directed CP Analysis

```r
A_lcc_dir <- as.matrix(as_adjacency_matrix(g_lcc_dir))
stats_dir <- get_directed_clique_density_stats(A_lcc_dir)

# stats_dir$matrix     — 2×2 asymmetric density matrix (C→C, C→P, P→C, P→P)
# stats_dir$n_core     — number of core nodes
# stats_dir$n_periph   — number of periphery nodes
# stats_dir$rho_global — global reciprocity
# stats_dir$rho_CC     — core internal reciprocity
# stats_dir$rho_PP     — periphery internal reciprocity
# stats_dir$rho_CP     — inter-block reciprocity
```

### 2.6 Output Directory Convention

For a dataset with identifier `XXXX`:

| Analysis | Output Folder |
|----------|---------------|
| Embedding (ASE/LSE) plots | `XXXX_Embed_outputs/` |
| Undirected CP stats | `XXXX_UNDIRECTED_CP_STATS_outputs/` |
| Directed CP stats | `XXXX_DIRECTED_CP_STATS_outputs/` |

PDF naming: `{path_part}_{method}.pdf` where `path_part` is `Mxxxxx_Type_xxxxxx` extracted from the file path.

---

## 3. Datasets Analyzed

| Identifier | Description | Edge List Path |
|------------|-------------|----------------|
| `3Vs1` | 3 vs 1 organoid dataset | `adjacency_edges.csv` (local R directory) |
| `MOvsSO` | Mature vs Single organoid (May 2025) | `adjacency_edges_May_31_2025_ecr_results_no_window.csv` |
| `M07357` | Individual organoid M07357 | `adjacency_edges_Aug_29_2025_M07357_para2.csv` |
