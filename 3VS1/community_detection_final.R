library('randnet')
library(irlba)
library(mclust)
library(igraph)
library(ggplot2)
library(MASS)
library(mvtnorm)
library(dplyr)
library(mvtnorm)
library(RColorBrewer)
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


stage1 <- c('000297', '000289')
stage2 <- c('000298', '000290')
stage3 <- c('000300', '000293')
stage4 <- c('000302', '000296')
stage5 <- c('000299', '000291')
stage6 <- c('000301', '000295')

get_stage <- function(file_path) {
  if (grepl(paste(stage1, collapse = "|"), file_path)) {
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

edge_list = read.csv('/Users/tianyichen/Desktop/Research /PhDresearch/Hopkins_Organoid/Codes/R/adjacency_edges.csv')
filenames = unique(edge_list$File)

# Initialize a data frame to store the summary
summary_table <- data.frame(
  File = character(),
  Well = character(),
  Num_edges = character(),
  Num_nodes = character(),
  Num_Communities_BHMC = character(),
  Num_Communities_ASE_GMM = character(),
  GMM = I(list()),
  stringsAsFactors = FALSE
)


file = filenames[1]

well = "well000"

subset_data <- subset(edge_list, File == file & Well == well)

if (nrow(subset_data) == 0) {
  next
}

n_rows <-   n_cols   <- unique(subset_data$dim)

# Create a sparse adjacency matrix
sparse_matrix <- sparseMatrix(
  i = subset_data$Row + 1,  
  j = subset_data$Column + 1,
  dims = c(n_rows, n_cols)
)


g <- graph_from_adjacency_matrix(sparse_matrix, weighted = NULL, mode = "undirected", diag = FALSE)

components <- igraph::clusters(g, mode="weak")
biggest_cluster_id <- which.max(components$csize)
vert_ids <- V(g)[components$membership == biggest_cluster_id]
g <- igraph::induced_subgraph(g, vert_ids)

n = vcount(g)
num_edges = ecount(g)

A = g[]

ed = 10

ase = full.ase(g, ed)
ase$eval
elbow = getElbows(ase$eval)
plot(ase$Xhat[,1], ase$Xhat[,2], main = 'ASE embedding')
Kmax = 30
ase_GMM =  Mclust(ase$Xhat[,1:elbow[2]], G=2:Kmax, verbose=T)
p2 = plot_GMM_contours(ase_GMM, contour_or_not = F)
print(p2)

lse = full.lse( A ,d=10)
elbow_lse = getElbows(lse$eval)
lse_GMM = Mclust(lse$Xhat[,1:elbow_lse[2]], G=2:Kmax, verbose=T)
p1 = plot_GMM_contours(lse_GMM, contour_or_not = F)
print(p1)



BIC_mat    <- lse_GMM$BIC
modelNames <- colnames(BIC_mat)
Gvals      <- as.numeric(rownames(BIC_mat))
n          <- nrow(lse$Xhat)
d          <- ncol(lse$Xhat[,1:elbow_lse[2]])
kappa      <- 3   # how harsh you want the penalty

# 2. Build df_mat with rows=G, cols=modelNames
df_mat <- matrix(
  nrow = length(Gvals),
  ncol = length(modelNames),
  dimnames = list(rownames(BIC_mat), colnames(BIC_mat))
)

df_mat

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


df_mat

# 3. Compute custom–BIC with penalty κ·p·log(n)
extra_penalty <- (kappa - 1) * df_mat * log(n)
customBIC     <- BIC_mat - extra_penalty

# 4. Find the best (G, model)
best_idx <- which(customBIC == max(customBIC, na.rm=TRUE), arr.ind=TRUE)
best_G    <- as.numeric(rownames(customBIC)[ best_idx[1] ])
best_mod  <- colnames(customBIC)[                best_idx[2] ]
cat("With κ =", kappa, "→ choose model", best_mod, "with G =", best_G, "\n")


lse_custom <- Mclust(
  lse$Xhat[, 1:elbow_lse[2]],
  G          = best_G,
  modelNames = best_mod,
  verbose    = TRUE
)


p3 = plot_GMM_contours(lse_custom, contour_or_not = F)
print(p3)



head(lse_GMM$BIC)
lse_GMM$BIC
plot(lse_GMM, what = "BIC")



comm <- cluster_leiden(
  graph              = g,
  resolution_parameter = 0.5,      # try e.g. 0.5, 1.0, 2.0
  objective_function  = "modularity",
  n_iterations        = -1         # use default full refinement
)
comm$nb_clusters



comm <- cluster_leiden(
  graph              = g,
  resolution_parameter = 0.5,      # try e.g. 0.5, 1.0, 2.0
  objective_function  = "modularity")
comm$nb_clusters



comm <- cluster_leiden(
  graph              = g,
  resolution_parameter = 1,      # try e.g. 0.5, 1.0, 2.0
  objective_function  = "modularity")
comm$nb_clusters


comm <- cluster_leiden(
  graph              = g,
  resolution_parameter = 2,      # try e.g. 0.5, 1.0, 2.0
  objective_function  = "modularity")
comm$nb_clusters
