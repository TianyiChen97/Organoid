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

n = vcount(g)
num_edges = ecount(g)



ase = full.ase(g, 20, spoke = F , spoke_coef = 0.01)

#A = g[]
#eA = eigs(A, k =10)
#eA$values


e = getElbows(sort(ase$eval,decreasing=TRUE ))
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
## directed version:
g <- graph_from_adjacency_matrix(sparse_matrix, weighted = NULL, mode = "undirected", diag = FALSE)

components <- igraph::clusters(g, mode="weak")
biggest_cluster_id <- which.max(components$csize)
vert_ids <- V(g)[components$membership == biggest_cluster_id]
g <- igraph::induced_subgraph(g, vert_ids)

n = vcount(g)
num_edges = ecount(g)



ase = full.ase(g, 20, spoke = F , spoke_coef = 0.01)

#ase$eval
A = g[]
eigs(A, k =10, which = "LM")$values
svds(A, k = 10)$d

eigs(A, k =10, which = "LM")$vectors[,10] - svds(A, k=10)$u[,10]


#eA$values

#e = getElbows(ase$eval)
#e = getElbows(sort(ase$eval,decreasing=TRUE ))

#ase$eval
#dhat = e[1]


## for directed ASE using both left and right singular vectors, we see eigenspoke.
plot(ase$Xhat[,1], ase$Xhat[,20])

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


