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
edge_list = read_csv("~/Desktop/Research /PhDresearch/Hopkins_Organoid/MO VS SO_2025_May_graphs/adjacency_edges_May_31_2025_ecr_results_no_window.csv")
filenames = unique(edge_list$File)
filenames = filenames[-2]
#Single organoid plate – M08438 
#Multi-organoid plate – M07359 


full.ase <- function(A, d, diagaug=TRUE, doptr=FALSE) {
  require(irlba)
  
  # doptr
  if (doptr) {
    g <- ptr(A)
    A <- g[]
  } else {
    A <- A[]
  }
  A = as.matrix(A)
  
  # diagaug
  if (diagaug) {
    diag(A) <- rowSums(A) / (nrow(A)-1)
  }
  
  A.svd <- irlba(A,d)
  Xhat <- A.svd$u %*% diag(sqrt(A.svd$d))
  Xhat.R <- NULL
  
  if (!isSymmetric(A)) {
    Xhat.R <- A.svd$v %*% diag(sqrt(A.svd$d))
  }
  
  return(list(eval=A.svd$d, Xhat=Matrix(Xhat), Xhat.R=Xhat.R))
}



full.lse <- function(A, d) {
  deg.seq <- rowSums(A)
  L.svd <- irlba(A/sqrt(outer(deg.seq,deg.seq)),d)
  Xhat <- L.svd$u %*% diag(sqrt(L.svd$d))
  Xhat.R <- NULL
  
  if (!isSymmetric(A)) {
    Xhat.R <- L.svd$v %*% diag(sqrt(L.svd$d))
  }
  
  return(list(eval=L.svd$d, Xhat=(Xhat), Xhat.R=(Xhat.R)))
}

getElbows <- function(dat, n = 3, threshold = FALSE, plot = TRUE, main="") {
  if (is.matrix(dat)) {
    d <- sort(apply(dat,2,sd), decreasing=TRUE)
  } else {
    d <- sort(dat,decreasing=TRUE)
  }
  
  if (!is.logical(threshold))
    d <- d[d > threshold]
  
  p <- length(d)
  if (p == 0)
    stop(paste("d must have elements that are larger than the threshold ",
               threshold), "!", sep="")
  
  lq <- rep(0.0, p)                     # log likelihood, function of q
  for (q in 1:p) {
    mu1 <- mean(d[1:q])
    mu2 <- mean(d[-(1:q)])              # = NaN when q = p
    sigma2 <- (sum((d[1:q] - mu1)^2) + sum((d[-(1:q)] - mu2)^2)) /
      (p - 1 - (q < p))
    lq[q] <- sum( dnorm(  d[1:q ], mu1, sqrt(sigma2), log=TRUE) ) +
      sum( dnorm(d[-(1:q)], mu2, sqrt(sigma2), log=TRUE) )
  }
  
  q <- which.max(lq)
  if (n > 1 && q < (p-1)) {
    q <- c(q, q + getElbows(d[(q+1):p], n-1, plot=FALSE))
  }
  
  if (plot==TRUE) {
    if (is.matrix(dat)) {
      sdv <- d # apply(dat,2,sd)
      plot(sdv,type="b",xlab="dim",ylab="stdev",main=main)
      points(q,sdv[q],col=2,pch=19)
    } else {
      plot(dat, type="b",main=main)
      points(q,dat[q],col=2,pch=19)
    }
  }
  
  return(q)
}
plot_GMM_contours <- function(mc, contour_or_not = TRUE) {
  
  
  # Use the data directly if it has 2 columns, otherwise apply PCA to extract the first two components
  if(ncol(mc$data) == 2) {
    Xhat <- mc$data
  } else {
    pca_result <- prcomp(mc$data, scale. = TRUE)
    Xhat <- pca_result$x[, 1:2]
  }
  
  clusters <- mc$classification  # Cluster labels
  
  # Create a data frame for plotting
  plot_data <- data.frame(
    X1 = Xhat[, 1],
    X2 = Xhat[, 2],
    Cluster = as.factor(clusters)
  )
  
  num_clusters <- length(unique(plot_data$Cluster))
  
  # Base plot with scatter points colored by cluster
  p <- ggplot(plot_data, aes(x = X1, y = X2, color = Cluster)) +
    geom_point(size = 3, alpha = 0.7) +
    labs(title = "GMM Clustering", x = "Dimension 1", y = "Dimension 2") +
    theme_minimal()
  
  # Use an appropriate palette: if more than 9 clusters, generate a custom palette
  if (num_clusters > 9) {
    p <- p + scale_color_manual(values = colorRampPalette(brewer.pal(9, "Set1"))(num_clusters))
  } else {
    p <- p + scale_color_brewer(palette = "Set1")
  }
  
  # If contour_or_not is TRUE, compute and add contour lines
  if (contour_or_not) {
    contour_data <- data.frame()
    
    for (k in 1:mc$G) {
      mean_k <- mc$parameters$mean[, k]  # Mean for cluster k
      sigma_k <- mc$parameters$variance$sigma[,,k]  # Covariance matrix for cluster k
      
      # Create a grid for density estimation
      x_seq <- seq(min(Xhat[, 1]), max(Xhat[, 1]), length.out = 100)
      y_seq <- seq(min(Xhat[, 2]), max(Xhat[, 2]), length.out = 100)
      grid <- expand.grid(X1 = x_seq, X2 = y_seq)
      
      # Compute the multivariate density on the grid
      density_values <- dmvnorm(grid, mean = mean_k, sigma = sigma_k)
      grid$Density <- density_values
      grid$Cluster <- as.factor(k)
      
      contour_data <- rbind(contour_data, grid)
    }
    
    p <- p + geom_contour(data = contour_data,
                          aes(x = X1, y = X2, z = Density, color = Cluster),
                          bins = 5)
  }
  
  return(p)
}


filenames

stage0 <- c('000386', '000384')  # ActivityScan
stage1 <- c('000387', '000385')  # Baseline #1 (pre‑stim)
stage2 <- c('000395', '000390')  # Baseline #2 (pre‑stim)
stage3 <- c('000396', '000391')  # During stimulation #1 (random‑site)
stage4 <- c('000397', '000392')  # After stimulation #1
stage5 <- c('000398', '000393')  # During stimulation #2 (single‑region)
stage6 <- c('000399', '000394')  # After stimulation #2

get_stage <- function(file_path) {
  if (grepl(paste(stage0, collapse = "|"), file_path)) {
    return(0)
  } else if (grepl(paste(stage1, collapse = "|"), file_path)) {
    return(1)
  } else if (grepl(paste(stage2, collapse = "|"), file_path)) {
    return(2)
  } else if (grepl(paste(stage3, collapse = "|"), file_path)) {
    return(3)
  } else if (grepl(paste(stage4, collapse = "|"), file_path)) {
    return(4)
  } else if (grepl(paste(stage5, collapse = "|"), file_path)) {
    return(5)
  } else if (grepl(paste(stage6, collapse = "|"), file_path)) {
    return(6)
  } else {
    return(NA) # Assign NA if no stage matches
  }
}


extract_chip_id <- function(path) {
  sub(
    ".*?/([^/]+)/(?:Network|Stimulation|ActivityScan)/.*",
    "\\1",
    path,
    perl = TRUE
  )
}



# 0) Choose your penalty multiplier
kappa <- 3
Kmax <- 50
# 1) Initialize summary_table with two extra integer columns:
summary_table <- data.frame(
  File                        = character(),
  Well                        = character(),
  Num_edges                   = integer(),
  Num_nodes                   = integer(),
  Num_Communities_Leiden_0.8  = integer(),
  Num_Communities_Leiden_1    = integer(),
  Num_Communities_Leiden_2    = integer(),
  Num_Communities_BHMC        = integer(),
  Num_Communities_ASE_GMM     = integer(),
  Num_Communities_LSE_GMM     = integer(),
  Num_Communities_ASE_BIC_kappa = integer(),
  Num_Communities_LSE_BIC_kappa = integer(),
  GMM_ASE                     = I(list()),
  GMM_LSE                     = I(list()),
  stringsAsFactors            = FALSE
)

total_files <- length(filenames)

for (i in seq_along(filenames))  {
  file <-  filenames[i]
  for (well in paste0("well", sprintf("%03d", 0:5))) {
    
    cat(sprintf("Processing file %d/%d: %s, %s...\n", i, total_files, file, well))
    
    
    subset_data <- subset(edge_list, File==file & Well==well)
    if (nrow(subset_data)==0){
      cat("  -> Skipped (no data)\n")
      next
    }

    # build graph
    n_rows <- n_cols <- unique(subset_data$dim)
    adj    <- sparseMatrix(
      i    = subset_data$Row+1,
      j    = subset_data$Column+1,
      dims = c(n_rows,n_cols)
    )
    g     <- graph_from_adjacency_matrix(adj, mode="undirected", diag=FALSE)
    comps <- clusters(g, mode="weak")
    big   <- which.max(comps$csize)
    g     <- induced_subgraph(g, V(g)[comps$membership==big])
    
    n_nodes <- vcount(g)
    n_edges <- ecount(g)
    
    if (n_nodes < 51 ) {
      summary_table <- rbind(summary_table, data.frame(
        File                         = file,
        Well                         = well,
        Num_edges                    = n_edges,
        Num_nodes                    = n_nodes,
        Num_Communities_Leiden_0.8   = NA,
        Num_Communities_Leiden_1     = NA,
        Num_Communities_Leiden_2     = NA,
        Num_Communities_BHMC         = NA,
        Num_Communities_ASE_GMM      = NA,
        Num_Communities_LSE_GMM      = NA,
        Num_Communities_ASE_BIC_kappa = NA,
        Num_Communities_LSE_BIC_kappa = NA,
        GMM_ASE                      = I(list(NA)),
        GMM_LSE                      = I(list(NA))
      ))
      cat("  -> Skipped (too few nodes)\n")
      next
    }
    
    # Leiden communities
    leiden_k <- lapply(c(0.8,1,2), function(res) {
      comm <- cluster_leiden(
        graph                = g,
        resolution_parameter = res,
        objective_function   = "modularity"
      )
      length(unique(membership(comm)))
    })
    
    # BHMC
    bhmc_k <- tryCatch(BHMC.estimate(g[], Kmax)$K[1], error=function(e) NA)
    
    # ASE–GMM
    ase_GMM_result <- tryCatch({
      ase   <- full.ase(g[], 10)
      elb   <- getElbows(ase$eval, plot = F)
      Mclust(ase$Xhat[,1:elb[2]], G=2:Kmax, verbose= T)
    }, error=function(e) NA)
    ase_k <-  ase_GMM_result$G
    
    # LSE–GMM
    lse_GMM_result <- tryCatch({
      lse   <- full.lse(g[], 10)
      elb   <- getElbows(lse$eval, plot = F)
      Mclust(lse$Xhat[,1:elb[2]], G=2:Kmax, verbose=T)
    }, error=function(e) NA)
    lse_k <-  lse_GMM_result$G
    
    # ───────────────────────────────────────────────────────────
    # 2) Custom‐BIC with extra penalty for ASE
    if (!is.na(ase_GMM_result$G)) {
      BIC_mat_ase    <- ase_GMM_result$BIC
      modelNames_ase <- colnames(BIC_mat_ase)
      Gvals_ase      <- as.numeric(rownames(BIC_mat_ase))
      n_ase          <- nrow(ase_GMM_result$data)  # or nrow(ase$Xhat)
      d_ase          <- ncol(ase_GMM_result$data)  # or length of ASE dims
      # build df_mat for ASE
      df_mat_ase <- outer(Gvals_ase, modelNames_ase,
                          Vectorize(function(G, mod) {
                            var_pars  <- nMclustParams(mod, d_ase, G)
                            mean_pars <- G * d_ase
                            prop_pars <- G - 1
                            var_pars + mean_pars + prop_pars
                          })
      )
      dimnames(df_mat_ase) <- list(rownames(BIC_mat_ase), modelNames_ase)
      extra_penalty_ase <- (kappa-1) * df_mat_ase * log(n_ase)
      customBIC_ase     <- BIC_mat_ase - extra_penalty_ase
      best_idx_ase      <- which(customBIC_ase == max(customBIC_ase, na.rm=TRUE), arr.ind=TRUE)
      num_of_communities_ase_BIC_kappa <- 
        as.numeric(rownames(customBIC_ase)[best_idx_ase[1]])
    } else {
      num_of_communities_ase_BIC_kappa <- NA
    }
    
    # 3) Custom‐BIC with extra penalty for LSE
    if (!is.na(lse_GMM_result$G)) {
      BIC_mat_lse    <- lse_GMM_result$BIC
      modelNames_lse <- colnames(BIC_mat_lse)
      Gvals_lse      <- as.numeric(rownames(BIC_mat_lse))
      n_lse          <- nrow(lse_GMM_result$data)
      d_lse          <- ncol(lse_GMM_result$data)
      df_mat_lse <- outer(Gvals_lse, modelNames_lse,
                          Vectorize(function(G, mod) {
                            var_pars  <- nMclustParams(mod, d_lse, G)
                            mean_pars <- G * d_lse
                            prop_pars <- G - 1
                            var_pars + mean_pars + prop_pars
                          })
      )
      dimnames(df_mat_lse) <- list(rownames(BIC_mat_lse), modelNames_lse)
      extra_penalty_lse <- (kappa-1) * df_mat_lse * log(n_lse)
      customBIC_lse     <- BIC_mat_lse - extra_penalty_lse
      best_idx_lse      <- which(customBIC_lse == max(customBIC_lse, na.rm=TRUE), arr.ind=TRUE)
      num_of_communities_lse_BIC_kappa <- 
        as.numeric(rownames(customBIC_lse)[best_idx_lse[1]])
    } else {
      num_of_communities_lse_BIC_kappa <- NA
    }
    # ───────────────────────────────────────────────────────────
    
    # 4) Finally rbind everything—atomic + list‑columns + new BIC columns
    summary_table <- rbind(summary_table, data.frame(
      File                           = file,
      Well                           = well,
      Num_edges                      = n_edges,
      Num_nodes                      = n_nodes,
      Num_Communities_Leiden_0.8     = leiden_k[[1]],
      Num_Communities_Leiden_1       = leiden_k[[2]],
      Num_Communities_Leiden_2       = leiden_k[[3]],
      Num_Communities_BHMC           = bhmc_k,
      Num_Communities_ASE_GMM        = ase_k,
      Num_Communities_LSE_GMM        = lse_k,
      Num_Communities_ASE_BIC_kappa  = num_of_communities_ase_BIC_kappa,
      Num_Communities_LSE_BIC_kappa  = num_of_communities_lse_BIC_kappa,
      GMM_ASE                        = I(list(ase_GMM_result)),
      GMM_LSE                        = I(list(lse_GMM_result))
    ))
    
    cat(sprintf("  -> Finished processing file %d/%d: %s, %s\n", i, total_files, file, well))
  }
}

example_GMM <- summary_table$GMM_ASE[[4]]
example_GMM <- summary_table$GMM_LSE[[4]]

plot_GMM_contours(example_GMM, contour_or_not = F)

BIC_mat    <- example_GMM$BIC
modelNames <- colnames(BIC_mat)
Gvals      <- as.numeric(rownames(BIC_mat))
n          <- nrow(example_GMM$data)
d          <- ncol(example_GMM$data)
kappa      <- 3   # how harsh you want the penalty

df_mat <- matrix(
  nrow = length(Gvals),
  ncol = length(modelNames),
  dimnames = list(rownames(BIC_mat), colnames(BIC_mat))
)


for (i in seq_along(Gvals)) {
  G <- Gvals[i]
  for (j in seq_along(modelNames)) {
    mod      <- modelNames[j]
    var_pars <- nMclustParams(mod, d, G)      # # covariance parameters
    mean_pars<- G * d                          # # mean parameters
    prop_pars<- G - 1                          # # mixing proportions
    df_mat[i,j] <- var_pars + mean_pars + prop_pars
  }
}


extra_penalty <- (kappa - 1) * df_mat * log(n)
customBIC     <- BIC_mat - extra_penalty

best_idx <- which(customBIC == max(customBIC, na.rm=TRUE), arr.ind=TRUE)
best_G    <- as.numeric(rownames(customBIC)[ best_idx[1] ])
best_mod  <- colnames(customBIC)[                best_idx[2] ]
#cat("With κ =", kappa, "→ choose model", best_mod, "with G =", best_G, "\n")

gmm_custom <- Mclust(
  example_GMM$data,
  G          = best_G,
  modelNames = best_mod,
  verbose    = F
)

plot_GMM_contours(gmm_custom, contour_or_not = F)
