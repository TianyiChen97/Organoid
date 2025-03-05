library('randnet')
library(irlba)
library(mclust)
library(igraph)
library(ggplot2)
library(MASS)
library(mvtnorm)
library(dplyr)



plot_GMM_contours <- function(mc) {
  library(mvtnorm)
  # Extract GMM model results
  clusters <- mc$classification  # Cluster labels
  Xhat <- mc$data[,1:2]  # Use first two dimensions for visualization
  
  # Create a data frame for plotting
  plot_data <- data.frame(
    X1 = Xhat[,1],
    X2 = Xhat[,2],
    Cluster = as.factor(clusters)
  )
  
  # Generate contour data for each cluster using estimated parameters
  contour_data <- data.frame()
  
  for (k in 1:mc$G) {
    mean_k <- mc$parameters$mean[, k]  # Mean of cluster k
    sigma_k <- mc$parameters$variance$sigma[,,k]  # Covariance matrix of cluster k
    
    # Create a grid for density estimation
    x_seq <- seq(min(Xhat[,1]), max(Xhat[,1]), length.out = 100)
    y_seq <- seq(min(Xhat[,2]), max(Xhat[,2]), length.out = 100)
    grid <- expand.grid(X1 = x_seq, X2 = y_seq)
    
    # Compute multivariate density on the grid
    density_values <- dmvnorm(grid, mean = mean_k, sigma = sigma_k)
    grid$Density <- density_values
    grid$Cluster <- as.factor(k)
    
    contour_data <- rbind(contour_data, grid)
  }
  
  # Plot scatter with GMM contours
  ggplot(plot_data, aes(x = X1, y = X2, color = Cluster)) +
    geom_point(size = 3, alpha = 0.7) +  # Scatter plot of ASE embeddings
    geom_contour(data = contour_data, aes(x = X1, y = X2, z = Density, color = Cluster), bins = 5) +  # GMM contour lines
    labs(title = "GMM Clustering with Contours on ASE Embedding",
         x = "ASE Dimension 1",
         y = "ASE Dimension 2") +
    theme_minimal() +
    scale_color_brewer(palette = "Set1")
}
full.ase <- function(A, d, diagaug=TRUE, doptr=FALSE) {
  require(RSpectra)
  
  # doptr
  if (doptr) {
    g <- ptr(A)
    A <- g[]
  } else {
    A <- A[]
  }
  
  # diagaug
  if (diagaug) {
    diag(A) <- rowSums(A) / (nrow(A)-1)
  }
  
  A.svd <- svds(A,k=d)
  Xhat <- A.svd$u %*% diag(sqrt(A.svd$d))
  Xhat.R <- NULL
  
  if (!isSymmetric(A)) {
    Xhat.R <- A.svd$v %*% diag(sqrt(A.svd$d))
  }
  
  return(list(eval=A.svd$d, Xhat=Matrix(Xhat), Xhat.R=Xhat.R))
}
scree_plot_tc <- function(kk){
  
  file = summary_table$File[kk]
  well = summary_table$Well[kk]
  
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
  
  g <- graph_from_adjacency_matrix(sparse_matrix, weighted = NULL, mode = "directed", diag = FALSE)
  
  components <- igraph::clusters(g, mode="weak")
  biggest_cluster_id <- which.max(components$csize)
  vert_ids <- V(g)[components$membership == biggest_cluster_id]
  g <- igraph::induced_subgraph(g, vert_ids)
  
  dmax = 100 
  
  ase = full.ase(g, d=dmax, diagaug=TRUE, doptr=FALSE)
  
  elb = getElbows(ase$eval, plot=T, main = paste('Scree plot for', summary_table$Group[kk] , 'Stage' , summary_table$Stage[kk]   ) )
  
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
doMclustASE_directed = function(g, dmax=100) {
  #pacman::p_load(irlba, gmmase, mclust)
  n = vcount(g)
  Kmax = n/5
  ase = full.ase(g, d=dmax, diagaug=TRUE, doptr=FALSE)
  elb = getElbows(ase$eval, plot=F)
  #dhat = max(elb[1],2)
  dhat = elb[1]
  Xhat = cbind(ase$Xhat[,1:dhat],ase$Xhat.R[,1:dhat])
  #Xhat = ase$Xhat[,1:dhat]
  mc = Mclust(Xhat, G=2:Kmax, verbose=T)
  Khat = mc$G
  return( list(mc=mc,Khat=Khat) )
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


for (file in filenames) {
  for (well in c("well000", "well001", "well002", "well003", "well004", "well005")) {
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
    
    g <- graph_from_adjacency_matrix(sparse_matrix, weighted = NULL, mode = "directed", diag = FALSE)
    
    components <- igraph::clusters(g, mode="weak")
    biggest_cluster_id <- which.max(components$csize)
    vert_ids <- V(g)[components$membership == biggest_cluster_id]
    g <- igraph::induced_subgraph(g, vert_ids)
    
    n = vcount(g)
    num_edges = ecount(g)
    
    if ( n  < 100 ) {
      # Add a row to the summary table with NULL communities
      summary_table <- rbind(summary_table, data.frame(
        File = file,
        Well = well,
        Num_edges = num_edges,
        Num_nodes = n,
        Num_Communities_BHMC = NA,
        Num_Communities_ASE_GMM = NA,
        GMM = I(list(NA))
      ))
      next
    }
    

    bhmc <- tryCatch({
      #print(c(file,well))
      as.numeric(BHMC.estimate(g[], n/5 )$K)[1]
    }, error = function(e) {
      NA  
    })
    
    CEPYP <- tryCatch({
      GMM_result = doMclustASE_directed(g, dmax = 10)
      
      # Extract Khat safely
      if (is.list(GMM_result) && "Khat" %in% names(GMM_result)) {
        as.numeric(GMM_result$Khat)[1] 
      } else {
        stop("GMM_result does not contain Khat")
      }
    }, error = function(e) {
      cat("Error in GMM clustering:", e$message, "\n")
      return(NA)
    })
    
    
    summary_table <- rbind(summary_table, data.frame(
      File = file,
      Well = well,
      Num_edges = num_edges,
      Num_nodes = n,
      Num_Communities_BHMC = bhmc,
      Num_Communities_ASE_GMM = CEPYP,  
      GMM = I(list(GMM_result)))
    )
  }
}

summary_table$Group <- ifelse(grepl("M07914", summary_table$File), "M07914",
                              ifelse(grepl("M07915", summary_table$File), "M07915", NA))



summary_file <- "summary_table.RData"

# Save the dataframe as an .RData file
save(summary_table, file = summary_file)

filtered_table <- subset(summary_table, !is.na(Num_Communities_BHMC) & !is.na(Group))

# Split into two groups
group_1 <- subset(filtered_table, Group == "M07914")
group_2 <- subset(filtered_table, Group == "M07915")


#Single organoid plate – M07914 
#Multi-organoid plate – M07915 

# Compute summary statistics for each group
group_1_mean <- mean(group_1$Num_Communities_BHMC)
group_2_mean <- mean(group_2$Num_Communities_BHMC)

group_1_median <- median(group_1$Num_Communities_BHMC)
group_2_median <- median(group_2$Num_Communities_BHMC)


# Print summary statistics
cat("Group M07914: Mean =", group_1_mean, ", Median =", group_1_median, "\n")
cat("Group M07915: Mean =", group_2_mean, ", Median =", group_2_median, "\n")


wilcox.test(group_1$Num_Communities_BHMC, group_2$Num_Communities_BHMC)

wilcox.test(group_1$Num_Communities_ASE_GMM, group_2$Num_Communities_ASE_GMM)


hist(group_1$Num_Communities_BHMC)
hist(group_2$Num_Communities_BHMC)

group_1_mean <- mean(na.omit(group_1$Num_Communities_ASE_GMM))
group_2_mean <- mean(na.omit(group_2$Num_Communities_ASE_GMM))

group_1_median <- median(na.omit(group_1$Num_Communities_ASE_GMM))
group_2_median <- median(na.omit(group_2$Num_Communities_ASE_GMM))


# Print summary statistics
cat("Group M07914: Mean =", group_1_mean, ", Median =", group_1_median, "\n")
cat("Group M07915: Mean =", group_2_mean, ", Median =", group_2_median, "\n")


#33 38 7 2 




mc = summary_table$GMM[2][[1]]$mc
plot_GMM_contours(mc)


mc = summary_table$GMM[33][[1]]$mc
plot_GMM_contours(mc)




# Apply function to create new Stage column

# View the updated dataframe
head(summary_table)



library(dplyr)

# Group by 'Group' and 'Stage' and compute mean and count of non-NA values
summary_table_grouped_BHMC <- summary_table %>%
  group_by(Group, Stage) %>%
  summarise(
    max_NUM_communities = max(Num_Communities_BHMC, na.rm = TRUE),
    median_NUM_communities = median(Num_Communities_BHMC, na.rm = TRUE),
    mean_NUM_communities = mean(Num_Communities_BHMC, na.rm = TRUE),
    num_nonNA = sum(!is.na(Num_Communities_BHMC))
    , .groups = "drop")  # Removes unnecessary grouping

# View the result
print(summary_table_grouped_BHMC)



plot_GMM_contours(summary_table$GMM)




library(dplyr)

# Group by 'Group' and 'Stage' and compute mean and count of non-NA values
summary_table_grouped_ASE_GMM <- summary_table %>%
  group_by(Group, Stage) %>%
  summarise(
    max_NUM_communities = max(Num_Communities_ASE_GMM, na.rm = TRUE),
    median_NUM_communities = median(Num_Communities_ASE_GMM, na.rm = TRUE),
    mean_NUM_communities = mean(Num_Communities_ASE_GMM, na.rm = TRUE),
    num_nonNA = sum(!is.na(Num_Communities_ASE_GMM))
    , .groups = "drop")  # Removes unnecessary grouping

# View the result
print(summary_table_grouped_ASE_GMM)


summary_table_grouped <- summary_table %>%
  group_by(Group, Stage) %>%
  summarise(mean_NUM_communities = mean(Num_Communities_ASE_GMM, na.rm = TRUE)) %>%
  ungroup()

# View the result
print(summary_table_grouped)


