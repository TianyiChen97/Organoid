library(Rcpp)
library(doParallel)
library(foreach)

# ---------------------------------------------------------------------------
# Part 1: C++ Greedy Optimization Function (Unchanged)
# ---------------------------------------------------------------------------
cppFunction('IntegerVector borgattiCpp(IntegerMatrix A){

  double n = A.nrow();
  double m = n*(n-1)/2;
  int k = 0;
  
  IntegerVector C(n);
  for(int i = 0; i < n; i++){
    if(rand()%10 > 5){ C(i) = 1; k++; } else { C(i) = 0; }
  }
  
  int kk = k;
  double sum_A = 0; 
  double sum_CP = 0.5*k*(k-1)+(n-k)*k; 
  double sum_ACP = 0;
  
  for (int i = 1; i < n; i++){
      for(int j = 0; j < i; j++){
        sum_A   += A(i,j);
        sum_ACP += A(i,j) * (C(i) + C(j) - C(i)*C(j)); 
      }
    }
  
  double obj_last = (sum_ACP - sum_A / m * sum_CP) / 
    (sqrt(sum_A - sum_A/m*sum_A) * sqrt(sum_CP - sum_CP/m*sum_CP) );
                    
  int ind = 1;
  int n_iter = 0;
  int max_iter = 100;
  
  IntegerVector Ctest = clone(C);
  IntegerVector v(n);
  for(int i=0; i<n;i++){v(i) = i;}
  
  double obj_new = 0;
  
  while(ind > 0 & n_iter < max_iter){
    ind = 0;
    int N = n;
    for(int a=0;a<n;a++){
      int b = round(rand()%N);
      int c = v(a); v(a)=v(b); v(b)=c;
    }
    
    for(int i = 0; i < n; i++){
      Ctest = clone(C);
      Ctest(v(i)) = 1 - Ctest(v(i));
      
      if(Ctest(v(i))==0){ kk = k - 1; } else { kk = k + 1; }

      double sum_CP_test = 0.5*kk*(kk-1) + (n-kk)*kk; 
      double sum_ACP_test = sum_ACP;
    
      for(int j = 0; j < n; j++){
         sum_ACP_test += A(v(i),j) * (Ctest(v(i)) + Ctest(j) - Ctest(v(i))*Ctest(j)); 
         sum_ACP_test -= A(v(i),j) * (C(v(i)) + C(j) - C(v(i))*C(j)); 
      }

      obj_new = (sum_ACP_test - sum_A / m * sum_CP_test) / 
       (sqrt(sum_A - sum_A/m*sum_A) * sqrt(sum_CP_test - sum_CP_test/m*sum_CP_test) );

      if(obj_new > obj_last){
        C(v(i)) = Ctest(v(i));
        sum_ACP = sum_ACP_test;
        ind++;
        obj_last = obj_new;
        k = kk;
      }
    }
    n_iter++;
  }
  return( C );
}')

# ---------------------------------------------------------------------------
# Part 2: Main Divide-and-Conquer Function (Updated with Metric)
# ---------------------------------------------------------------------------
dacCP <- function(A, m=0.10, iters=1000, no_cores = detectCores()-1){
  
  A <- as.matrix(A)
  n <- nrow(A)
  
  # --- Step 1: Divide and Conquer ---
  if(.Platform$OS.type == "unix") {
    cl <- makeCluster(no_cores, type="FORK")
  } else {
    cl <- makeCluster(no_cores, type="PSOCK")
    clusterExport(cl, "borgattiCpp")
  }
  registerDoParallel(cl)
  
  out <- foreach(i=1:iters, .combine="+", .packages='Rcpp') %dopar% {
    idx <- sample(1:n, round(n*m))
    B   <- A[idx, idx]
    while(sum(B) == 0){
      idx <- sample(1:n, round(n*m))
      B   <- A[idx, idx]
    }
    C_sub <- as.logical( borgattiCpp(B) )
    samp <- core <- rep(0, n)
    samp[idx] <- 1
    core[idx[C_sub]] <- 1
    return(cbind(core, samp))
  }
  stopCluster(cl)
  
  # --- Step 2: Thresholding ---
  props <- out[,1] / pmax(out[,2], 1)
  
  M      <- matrix(0, nrow=n, ncol=3)
  M[,1]  <- 1:n
  M[,2]  <- props
  M <- M[order(M[,2], decreasing=TRUE),] 
  
  d   <- -diff(M[,2])  
  idx_cut <- which.max(d) 
  M[1:idx_cut, 3] <- 1  
  M <- M[order(M[,1], decreasing=FALSE),]
  
  final_core_labels <- M[,3]
  
  # --- Step 3: Calculate Final Metric rho(A, c) ---
  # Implements Eq(1) efficiently [cite: 68]
  
  # Total pairs (m in paper notation)
  total_pairs <- n * (n - 1) / 2
  
  # 1. Properties of Ideal Matrix (Delta)
  # k = number of core nodes
  k <- sum(final_core_labels)
  # sum_CP = number of 1s in ideal matrix (Core-Core edges + Core-Periph edges)
  sum_CP <- 0.5 * k * (k - 1) + k * (n - k)
  
  # 2. Properties of Observed Matrix (A)
  # sum_A = number of edges in A (assuming unweighted/binary as per paper)
  sum_A <- sum(A) / 2
  
  # 3. Covariance Term (sum_ACP)
  # Sum of A where Ideal is 1. This equals Total Edges - Edges strictly inside Periphery.
  # We find edges in periphery to subtract (faster for sparse A)
  periph_indices <- which(final_core_labels == 0)
  if(length(periph_indices) > 1){
    edges_in_periph <- sum(A[periph_indices, periph_indices]) / 2
  } else {
    edges_in_periph <- 0
  }
  sum_ACP <- sum_A - edges_in_periph
  
  # 4. Pearson Correlation Formula
  # Numerator: E[XY] - E[X]E[Y]
  num <- sum_ACP - (sum_A * sum_CP / total_pairs)
  
  # Denominator: sigma_X * sigma_Y
  # Note: For binary matrices, sum(X^2) = sum(X)
  den_A  <- sqrt(sum_A  - (sum_A  * sum_A  / total_pairs))
  den_CP <- sqrt(sum_CP - (sum_CP * sum_CP / total_pairs))
  
  rho <- num / (den_A * den_CP)
  
  return(list(
    core_nodes = which(final_core_labels == 1),
    periphery_nodes = which(final_core_labels == 0),
    coreness_scores = props,
    rho = rho
  ))
}


get_rho_fast <- function(A, core_idx) {
  # A: Adjacency matrix (sparse or dense)
  # core_idx: Vector of indices belonging to the core
  
  n <- nrow(A)
  total_pairs <- n * (n - 1) / 2
  
  # 1. Calculate properties of the Ideal Model (Y)
  k <- length(core_idx)
  # Count of 1s in ideal model: All pairs except Periphery-Periphery
  # Formula: Total Pairs - Periphery Pairs
  n_periph <- n - k
  pairs_periph <- n_periph * (n_periph - 1) / 2
  sum_CP <- total_pairs - pairs_periph
  
  # 2. Calculate properties of the Data (X)
  sum_A <- sum(A) / 2 # Total edges
  
  # 3. Calculate Cross Product (XY)
  # This is the number of edges that fall in the "Core" region (Core-Core or Core-Periph)
  # Easier to calc: Total Edges - Edges inside Periphery
  periph_idx <- setdiff(1:n, core_idx)
  if(length(periph_idx) > 1){
    edges_periph <- sum(A[periph_idx, periph_idx]) / 2
  } else {
    edges_periph <- 0
  }
  sum_ACP <- sum_A - edges_periph
  
  # 4. Pearson Correlation
  # Covariance numerator
  num <- sum_ACP - (sum_A * sum_CP / total_pairs)
  
  # Standard Deviations
  # Variance of binary vector X = p(1-p) * N, but using raw sums is safer for consistency
  sd_A  <- sqrt(sum_A  - (sum_A^2  / total_pairs))
  sd_CP <- sqrt(sum_CP - (sum_CP^2 / total_pairs))
  
  return(num / (sd_A * sd_CP))
}


calc_rho_simple <- function(A, core_index) {
  # 0. Get network size and create binary labels from indices
  n <- nrow(A)
  core_labels <- rep(0, n)       # Initialize all as Periphery (0)
  core_labels[core_index] <- 1   # Set Core nodes to 1
  
  # 1. Create the perfect Core-Periphery matrix Delta
  # Delta_ij = 1 if i or j is Core, 0 if both are Periphery
  Delta <- matrix(0, n, n)
  
  # Outer product logic: (c_i + c_j - c_i*c_j) is 1 if either is 1
  Delta <- outer(core_labels, core_labels, function(x, y) x + y - x*y)
  diag(Delta) <- 0 # Ignore self-loops
  
  # 2. Extract upper triangles
  vec_A <- A[upper.tri(A)]
  vec_Delta <- Delta[upper.tri(Delta)]
  
  # 3. Calculate Correlation
  return(cor(vec_A, vec_Delta))
}

# 1. Create a random symmetric adjacency matrix (example)
n <- 100
A <- matrix(rbinom(n*n, 1, 0.1), n, n)
A[lower.tri(A)] <- t(A)[lower.tri(A)]
diag(A) <- 0

# 2. Run the function
result <- dacCP(A, m=0.1, iters=500)
result

