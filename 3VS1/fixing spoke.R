library('randnet')
library(irlba)
library(mclust)
library(igraph)
library(ggplot2)
library(MASS)
library(mvtnorm)
library(dplyr)
## Fix the eigenspoke

file = filenames [1]

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


## directed version:
g <- graph_from_adjacency_matrix(sparse_matrix, weighted = NULL, mode = "directed", diag = FALSE)
components <- igraph::clusters(g, mode="weak")
biggest_cluster_id <- which.max(components$csize)
vert_ids <- V(g)[components$membership == biggest_cluster_id]
g <- igraph::induced_subgraph(g, vert_ids)

plot(g, 
     vertex.size = 1,           # Adjust vertex size
     vertex.color = "lightblue", # Set vertex color
     vertex.label.cex = 0.1,     # Adjust label size
     edge.arrow.size = 0.5) 

n = vcount(g)
num_edges = ecount(g)

A = g[]



in_degrees <- degree(g, mode = "in")

# For out-degrees (number of outgoing edges)
out_degrees <- degree(g, mode = "out")

table(in_degrees)
table(out_degrees)




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




ase = full.ase(g, d = 20, diagaug = F, spoke = F , spoke_coef = 0.01)


elb = getElbows( ase$eval )

dhat = max(elb[1],2)
dhat
Xhat = cbind(ase$Xhat[,1:dhat],ase$Xhat.R[,1:dhat])

pairs(as.matrix(Xhat))

Kmax = 200
mc = Mclust(Xhat, G=2:Kmax, verbose=T)

mc$G

if(ncol(mc$data) == 2) {
  Xhat <- mc$data
} else {
  pca_result <- prcomp(mc$data, scale. = F)
  Xhat <- pca_result$x[, 1:2]
}

plot(Xhat[,1],Xhat[,2])

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

print(p)

p1 = plot_GMM_contours(mc, contour_or_not = F)

#e = getElbows(sort(ase$eval,decreasing=TRUE ))
#ase$eval

dhat = e[1]

## for directed ASE using both left and right singular vectors, we see eigenspoke.
plot(ase$Xhat[,1:e[1]], ase$Xhat.R[,1:e[1]], main = paste( 'original directed graph ',
                                                           'added constant = ' , cc ))

hist(ase$Xhat[,1], breaks = 40)
plot(density(ase$Xhat[,1]))
plot(density(ase$Xhat.R[,1]))


ase = full.ase(g, 20, spoke = T , spoke_coef = 2)
e = getElbows(sort(ase$eval,decreasing=TRUE ))
plot(ase$Xhat[,1:e[1]], ase$Xhat.R[,1:e[1]])




##### undirected version 
g <- graph_from_adjacency_matrix(sparse_matrix, weighted = NULL, mode = "undirected", diag = FALSE)

components <- igraph::clusters(g, mode="weak")
biggest_cluster_id <- which.max(components$csize)
vert_ids <- V(g)[components$membership == biggest_cluster_id]
g <- igraph::induced_subgraph(g, vert_ids)

n = vcount(g)
num_edges = ecount(g)

A = g[]




### try Lisa's regularization

# 1. Compute the degree for each vertex.
# You may use igraph's degree() function, which is equivalent to the rowSums of the adjacency matrix.
deg_seq <- degree(g)

# 2. Estimate d.
# Here d is theoretically defined as: d = max_ij (n * p_ij).
# One simple estimator is to let d be the maximum observed degree.
d_est <- max(deg_seq)*0.00001
# Note: For real data you might want to use a model-based estimator if you have estimates of the p_ij.

# 3. Compute the lambda values for each vertex.
# For each vertex i, lambda_i = min(2*d_est / d_i, 1)
lambda <- pmin(2 * d_est / deg_seq, 1)

lambda

W <- sqrt(outer(lambda, lambda)) * A

ase = full.ase( W , 20, spoke = F , spoke_coef = 0.01)
#ase$eval
#getElbows(ase$eval)
plot(ase$Xhat[,1], ase$Xhat[,2], main = 'ASE embedding')












deg.seq <- rowSums(A)
L.svd <- irlba(A/sqrt(outer(deg.seq,deg.seq)),10)
Xhat <- L.svd$u %*% diag(sqrt(L.svd$d))



deg.seq <- rowSums(A)
L.svd <- irlba(A/sqrt(outer(deg.seq,deg.seq)),10)
Xhat <- L.svd$u %*% diag(sqrt(L.svd$d))
Xhat.R <- NULL

if (!isSymmetric(A)) {
  Xhat.R <- L.svd$v %*% diag(sqrt(L.svd$d))
}


plot(Xhat[,1],Xhat[,2])


par(mfrow=c(2,2))
lse = full.lse(A,d=10)
lse$eval
getElbows(lse$eval)
plot(lse$Xhat[,1],lse$Xhat[,2],main='LSE embedding')

Kmax = 200
lse_GMM = Mclust(lse$Xhat[,1:2], G=2:Kmax, verbose=T)

p1 = plot_GMM_contours(lse_GMM)


ase = full.ase(g, 20, spoke = F , spoke_coef = 0.01)
ase$eval
getElbows(ase$eval)
plot(ase$Xhat[,1], ase$Xhat[,2], main = 'ASE embedding')

ase_GMM =  Mclust(ase$Xhat[,1:2], G=2:Kmax, verbose=T)

p2 = plot_GMM_contours(ase_GMM)

library(gridExtra)
grid.arrange(p1, p2, ncol = 2)



require(RSpectra)
eigs(A, k =10, which = "LM")$values
svds(A, k = 10)$d
irlba(A, nv = 10)$d
#eigs(A, k =10, which = "LM")$vectors[,10] - svds(A, k=10)$u[,10]

#eA$values

e = getElbows(ase$eval)
#e = getElbows(sort(ase$eval,decreasing=TRUE ))

#ase$eval
#dhat = e[1]


## for directed ASE using both left and right singular vectors, we see eigenspoke.
plot(ase$Xhat[,1], ase$Xhat[,2])

hist(ase$Xhat[,1], breaks = 40)
plot(density(ase$Xhat[,1]))
plot(density(ase$Xhat[,20]))

cc=0

ase = full.ase(g, 20, spoke = T , spoke_coef = cc)
#ase$eval
#e = getElbows(ase$eval)
#e = getElbows(sort(ase$eval,decreasing=TRUE ))
plot(ase$Xhat[,1], ase$Xhat[,20], main = paste( 'symmetrized graph ',  'added constant = ' , cc ) )




GMM_result = doMclustASE_directed(g, dmax = 10 ,  sk = T , s_c = 0.1 )
GMM_result$Khat

plot_GMM_contours(GMM_result$mc)


