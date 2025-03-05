
library('randnet')
library('igraph')
library(irlba)
library(mclust)

#file = "/cis/project/organoid/Dec_29_2024_ecr_results/M07915/Stimulation/000291/data.raw_20241213_18h15m.pkl"

file = "/cis/project/organoid/Dec_29_2024_ecr_results/M07914/Stimulation/000299/data.raw_20241213_18h15m.pkl"
#well =  "well003"   
well = "well001"

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

# Estimate the number of communities using BHMC
g <- graph_from_adjacency_matrix(sparse_matrix, weighted = NULL, mode = "directed", diag = FALSE)



components <- igraph::clusters(g, mode="weak")
biggest_cluster_id <- which.max(components$csize)

# ids
vert_ids <- V(g)[components$membership == biggest_cluster_id]
vcount(g)
# subgraph
g <- igraph::induced_subgraph(g, vert_ids)

dmax= 100
n = vcount(g)
doMclustASE(g, dmax = 100)

svds(g[], k = 71)

#Kmax = n/2
ase = full.ase(g, d=dmax, diagaug=TRUE, doptr=FALSE)
#ase$Xhat[,1]
#ase$Xhat.R[,1]
elb = getElbows(ase$eval, plot=F)
#dhat = max(elb[1],2)
dhat=elb[1]
Xhat = cbind(ase$Xhat[,1:dhat],ase$Xhat.R[,1:dhat])
#Xhat = ase$Xhat[,1:dhat]
plot(Xhat[,1],Xhat[,2],xlab = 'Left singular vector', ylab = 'right singular vector' )
hist(Xhat[,2],breaks = 40)

sum(g[])

#BIC <- mclustBIC(Xhat[,1])
#plot(BIC)

Kmax = n/2
mc = Mclust(Xhat, G=2:Kmax, verbose=TRUE)
plot(mc)
mc$G
0








