library('randnet')
library('igraph')
library(irlba)
library(mclust)
library(irlba)

##ASE for a network A with embedding dimension d
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


## Time-stamp: <Sun Mar 20, 2016 11:49:22 YP>

getElbows <- function(dat, n = 3, threshold = FALSE, plot = TRUE, main="") {
  ## Given a decreasingly sorted vector, return the given number of elbows
  ##
  ## Args:
  ##   dat: a input vector (e.g. a vector of standard deviations), or a input feature matrix.
  ##   n: the number of returned elbows.
  ##   threshold: either FALSE or a number. If threshold is a number, then all
  ##   the elements in d that are not larger than the threshold will be ignored.
  ##   plot: logical. When T, it depicts a scree plot with highlighted elbows.
  ##
  ## Return:
  ##   q: a vector of length n.
  ##
  ## Reference:
  ##   Zhu, Mu and Ghodsi, Ali (2006), "Automatic dimensionality selection from
  ##   the scree plot via the use of profile likelihood", Computational
  ##   Statistics & Data Analysis, Volume 51 Issue 2, pp 918-930, November, 2006. 
  
  #  if (is.unsorted(-d))
  
  
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

doMclustASE = function(g, dmax=100) {
  #pacman::p_load(irlba, gmmase, mclust)
  n = vcount(g)
  Kmax = n/2
  ase = full.ase(g, d=dmax, diagaug=TRUE, doptr=FALSE)
  elb = getElbows(ase$eval, plot=F)
  #dhat = max(elb[1],2)
  dhat = elb[1]
  Xhat = ase$Xhat[,1:dhat]
  mc = Mclust(Xhat, G=2:Kmax, verbose=FALSE)
  Khat = mc$G
  return(Khat)
}


edge_list = read.csv('/Users/tianyichen/Desktop/Research /PhDresearch/Hopkins_Organoid/Codes/R/adjacency_edges.csv')
filenames = unique(edge_list$File)

# Define the threshold for density
density_threshold <- 0.01  # Adjust as needed

# Initialize a data frame to store the summary
summary_table <- data.frame(
  File = character(),
  Well = character(),
  Num_edges = character(),
  Num_nodes = character(),
  Density = character(),
  Num_Communities_BHMC = character(),
  Num_Communities_ASE_GMM = character(),
  stringsAsFactors = FALSE
)


for (file in filenames) {
  for (well in c("well000", "well001", "well002", "well003", "well004", "well005")) {
    # Subset the data for the current file and well
    subset_data <- subset(edge_list, File == file & Well == well)
    
    if (nrow(subset_data) == 0) {
      # If no data is found for the well, skip it
      next
    }
    
    # Extract dimensions
    n_rows <-   n_cols   <- unique(subset_data$dim)
    
    # Create a sparse adjacency matrix
    sparse_matrix <- sparseMatrix(
      i = subset_data$Row + 1,  # Convert to 1-based indexing
      j = subset_data$Column + 1,
      dims = c(n_rows, n_cols)
    )
    
    # Calculate the density of the matrix
    num_edges <- sum(sparse_matrix)
    density <- num_edges / (n_rows * n_cols)
    # Check if density meets the threshold
    if (density < density_threshold) {
      # Add a row to the summary table with NULL communities
      summary_table <- rbind(summary_table, data.frame(
        File = file,
        Well = well,
        Num_edges = num_edges,
        Num_nodes = n_rows,
        Density = density,
        Num_Communities_BHMC = NA,
        Num_Communities_ASE_GMM = NA
      ))
      next
    }
    
    sample_graph = sparse_matrix
    # Estimate the number of communities using BHMC
    g <- graph_from_adjacency_matrix(sparse_matrix, weighted = NULL, mode = "directed", diag = FALSE)
    
    components <- igraph::clusters(g, mode="strong")
    biggest_cluster_id <- which.max(components$csize)
    vert_ids <- V(g)[components$membership == biggest_cluster_id]

    g <- igraph::induced_subgraph(g, vert_ids)
    n = vcount(g)
    
    #g <- graph_from_adjacency_matrix(sparse_matrix, weighted = NULL, mode = "undirected", diag = FALSE)
    
    bhmc <- tryCatch({
      #print(c(file,well))
      as.numeric(BHMC.estimate(g[], n/2 )$K)[1]
    }, error = function(e) {
      NA  # If an error occurs, set bhmc to NA
    })
    
    # Estimate number of communities using ASE GMM with error handling
    CEPYP <- tryCatch({
      #print(c(file,well))
      as.numeric(doMclustASE(g, dmax = 100))[1]
    }, error = function(e) {
      NA  # If an error occurs, set CEPYP to NA
    })
    
    print(c( n, density ,bhmc, CEPYP))

    # Add the result to the summary table
    summary_table <- rbind(summary_table, data.frame(
      File = file,
      Well = well,
      Num_edges = num_edges,
      Num_nodes = n_rows,
      Density = density,
      Num_Communities_BHMC = bhmc,
      Num_Communities_ASE_GMM = CEPYP
    ))
  }
}



# Add a new column for group based on whether the filename contains "M07914" or "M07915"
summary_table$Group <- ifelse(grepl("M07914", summary_table$File), "M07914",
                              ifelse(grepl("M07915", summary_table$File), "M07915", NA))

# Remove rows with NA in Num_Communities or Group
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


wilcox.test(group_1$Num_Communities, group_2$Num_Communities)


group_1_mean <- mean(group_1$Num_Communities_ASE_GMM)
group_2_mean <- mean(group_2$Num_Communities_ASE_GMM)

group_1_median <- median(group_1$Num_Communities_ASE_GMM)
group_2_median <- median(group_2$Num_Communities_ASE_GMM)


# Print summary statistics
cat("Group M07914: Mean =", group_1_mean, ", Median =", group_1_median, "\n")
cat("Group M07915: Mean =", group_2_mean, ", Median =", group_2_median, "\n")



hist(group_1$Num_Communities_ASE_GMM, breaks = 30)

bhmc

num_communities <-  bhmc$K


