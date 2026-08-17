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

MAIN_IDENTIFIER <- "MOvsSO"

# Choose your embedding method here ("ASE" or "LSE")
EMBEDDING_METHOD <- "ASE" 

if (EMBEDDING_METHOD == "LSE") {
  embedding_function <- full.lse 
} else {
  embedding_function <- full.ase
}

# Create Output Directories
output_dir_embed <- paste0(MAIN_IDENTIFIER, "_Embed_outputs")
output_dir_stats <- paste0(MAIN_IDENTIFIER, "_UNDIRECTED_CP_STATS_outputs")

if (!dir.exists(output_dir_embed)) dir.create(output_dir_embed)
if (!dir.exists(output_dir_stats)) dir.create(output_dir_stats)

print(paste("Processing", length(filenames), "files for dataset:", MAIN_IDENTIFIER))

# ===================================================================
# 3. ANALYSIS A: EMBEDDING (ASE/LSE)
# ===================================================================

message("Starting Embedding Analysis...")

for (file in filenames) {
  
  # Title Logic
  path_part <- extract_network_path(file)
  if (is.na(path_part)) path_part <- basename(dirname(file))
  
  full_title_string <- paste(path_part, EMBEDDING_METHOD, sep = " | ")
  
  # File Logic
  pdf_filename_base <- gsub("/", "_", path_part)
  pdf_filename_base <- gsub(" ", "_", pdf_filename_base)
  pdf_filename <- paste0(pdf_filename_base, "_", EMBEDDING_METHOD, 'fixdim=2.pdf')
  final_pdf_path <- file.path(output_dir_embed, pdf_filename)
  
  pdf(file = final_pdf_path, width = 14, height = 12)
  par(mfrow = c(2, 3), mar = c(3, 3, 3, 1), oma = c(4, 0, 3, 0))
  
  well_list <- paste0("well", sprintf("%03d", 0:5))
  
  for (well in well_list) {
    subset_data <- subset(edge_list, File == file & Well == well)
    
    if (nrow(subset_data) == 0) {
      plot(1, type="n", axes=F, xlab="", ylab=""); title(main=well); text(1,1,"No data", col="red"); next
    }
    
    tryCatch({
      n_rows <- n_cols <- unique(subset_data$dim)
      adj <- sparseMatrix(i=subset_data$Row+1, j=subset_data$Column+1, dims=c(n_rows, n_cols))
      
      # Undirected Graph & LCC
      g <- graph_from_adjacency_matrix(adj, mode="undirected", diag=FALSE)
      comps <- components(g, mode="weak")
      g_lcc <- induced_subgraph(g, V(g)[comps$membership == which.max(comps$csize)])
      
      if (vcount(g_lcc) < 11) {
        plot(1, type="n", axes=F, xlab="", ylab=""); title(main=well)
        text(1, 1, "LCC < 11", col = "red")
      } else {
        # Run Embedding
        embed_result <- embedding_function(g_lcc, 10)
        elb <- getElbows(embed_result$eval, plot=F)
        
        if (length(elb) < 2) stop("Embedding dim < 2")
        #embedding_dim <- elb[2]
        embedding_dim <- 2
        if (embedding_dim == 2) {
          plot(embed_result$Xhat[,1], embed_result$Xhat[,2], 
               main=well, xlab="Dim 1", ylab="Dim 2", pch=16, col="blue")
        } else if (embedding_dim >= 3) {
          scatterplot3d(embed_result$Xhat[,1], embed_result$Xhat[,2], embed_result$Xhat[,3],
                        main=well, xlab="D1", ylab="D2", zlab="D3", color="blue", pch=16)
        } else {
          plot(1, type="n", axes=F); title(main=well); text(1,1,"Dim Issue", col="orange")
        }
      }
    }, error = function(e) {
      plot(1, type="n", axes=F); title(main=well); text(1,1,"Error", col="red")
    })
  }
  mtext(full_title_string, side=3, line=1, outer=TRUE, cex=1.5, font=2)
  dev.off()
}

# ===================================================================
# 4. ANALYSIS B: UNDIRECTED CP STATISTICS
# ===================================================================

message("Starting Undirected GCP Analysis...")

for (file in filenames) {
  
  # Title Logic
  path_part <- extract_network_path(file)
  if (is.na(path_part)) path_part <- basename(dirname(file))
  
  full_title_string <- paste(path_part, "Undirected_CP", sep = " | ")
  
  # File Logic
  pdf_filename_base <- gsub("/", "_", path_part)
  pdf_filename_base <- gsub(" ", "_", pdf_filename_base)
  pdf_filename <- paste0(pdf_filename_base, "_Undir_CP.pdf")
  final_pdf_path <- file.path(output_dir_stats, pdf_filename)
  
  pdf(file = final_pdf_path, width = 14, height = 12)
  par(mfrow = c(2, 3), mar = c(1, 1, 3, 1), oma = c(2, 0, 3, 0))
  
  for (well in well_list) {
    subset_data <- subset(edge_list, File == file & Well == well)
    
    # Blank Plot
    plot(1, type = "n", axes = FALSE, xlab = "", ylab = "", xlim=c(0, 10), ylim=c(0, 10))
    title(main = well, cex.main = 1.5); box()
    
    if (nrow(subset_data) == 0) {
      text(5, 5, "No data available", col = "red"); next
    }
    
    tryCatch({
      n_rows <- n_cols <- unique(subset_data$dim)
      adj <- sparseMatrix(i=subset_data$Row+1, j=subset_data$Column+1, dims=c(n_rows, n_cols))
      
      # Undirected Graph & LCC
      g <- graph_from_adjacency_matrix(adj, mode = "undirected", diag = FALSE)
      n_orig <- vcount(g)
      comps <- components(g, mode = "weak")
      g_lcc <- induced_subgraph(g, V(g)[comps$membership == which.max(comps$csize)])
      n_lcc <- vcount(g_lcc)
      
      if (n_lcc < 11) {
        text(5, 5, "LCC size < 11", col = "red", cex = 1.2)
      } else {
        # Calculate Undirected Stats
        stats_res <- get_clique_density_stats(as.matrix(as_adjacency_matrix(g_lcc)))
        dens_mat <- stats_res$matrix
        n_core <- stats_res$n_core
        n_periph <- stats_res$n_periph
        
        lines_to_print <- c(
          "--- Densities (Undir) ---",
          paste0("CC: ", sprintf("%.3f", dens_mat[1,1])),
          paste0("CP: ", sprintf("%.3f", dens_mat[1,2])),
          paste0("PP: ", sprintf("%.3f", dens_mat[2,2])),
          "",
          "--- Counts ---",
          paste0("Core: ", n_core),
          paste0("Peri: ", n_periph),
          paste0("LCC Size: ", n_lcc),
          paste0("Orig Size: ", n_orig),
          "",
          "--- Ratios ---",
          paste0("Core / LCC: ", sprintf("%.3f", n_core/n_lcc)),
          paste0("Peri / LCC: ", sprintf("%.3f", n_periph/n_lcc))
        )
        
        text(5, 5, labels = paste(lines_to_print, collapse = "\n"), cex = 1.3, family = "mono")
      }
    }, error = function(e) {
      text(5, 5, "Error", col = "red")
    })
  }
  mtext(full_title_string, side=3, line=1, outer=TRUE, cex=1.5, font=2)
  dev.off()
}

print("Processing Complete.")
