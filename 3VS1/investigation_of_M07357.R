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
well_list <- paste0("well", sprintf("%03d", 0:5))


file = filenames[1]
well=well_list[1]

subset_data <- subset(edge_list, File == file & Well == well)

n_rows <- n_cols <- unique(subset_data$dim)
adj <- sparseMatrix(
    i    = subset_data$Row + 1,
    j    = subset_data$Column + 1,
    dims = c(n_rows, n_cols)
)
  
  # Create the graph and find the largest connected component (LCC)
g <- graph_from_adjacency_matrix(adj, mode = "undirected", diag = FALSE)
plot(g,vertex.label.cex = 0.5,vertex.size = 0.5)

comps <- components(g, mode = "weak")
big <- which.max(comps$csize)
g_lcc <- induced_subgraph(g, V(g)[comps$membership == big])
plot(g_lcc, 
     vertex.size = 0.5,           # Adjust vertex size
     vertex.color = "lightblue", # Set vertex color
     vertex.label.cex = 0.1,     # Adjust label size
     edge.arrow.size = 0.25) 
plot(g_lcc,vertex.label.cex = 0.5,vertex.size = 0.5)

largest_clqs <- largest_cliques(g_lcc)
clique_v = union(largest_clqs[[1]],largest_clqs[[2]])
all_v <- V(g_lcc)
periphery_v <- setdiff(all_v, clique_v)


degree(g_lcc, v = clique_v)
periphery_degrees <- degree(g_lcc, v = periphery_v)
periphery_degrees 

# 3. Count connections from each periphery node to the core
# This is done efficiently by subsetting the adjacency matrix and summing the rows.
adj_matrix_subset <- g_lcc[periphery_v, clique_v]
core_connection_counts <- rowSums(as.matrix(adj_matrix_subset))


# 4. Create the final data frame
periphery_analysis_df <- data.frame(
  Vertex_ID = periphery_v,
  Total_Degree = periphery_degrees,
  Core_Connections = core_connection_counts
)


# 5. Sort the table by the number of core connections to see the gradient
periphery_analysis_df <- periphery_analysis_df[order(periphery_analysis_df$Core_Connections), ]




v_colors <- rep("black", vcount(g_lcc)) # Default color for all vertices
e_colors <- rep("grey", ecount(g_lcc)) # Default color for all edges

clique_edges <- E(g_lcc)[clique_v %--% clique_v]
clique_peri_edges <- E(g_lcc)[clique_v %--% periphery_v]

e_colors[clique_edges] <- "red"
e_colors[clique_peri_edges] <- "green"
v_colors[ periphery_v] <- "green"             # Set clique vertices to red


# 3. Plot the graph with the new color attributes
plot(g_lcc,
     vertex.label.cex = 0.5,
     vertex.size = .5,          # Increased size slightly for visibility
     vertex.label.color = v_colors,  # Use your custom vertex colors
     edge.color = e_colors     # Use your custom edge colors
)

ase <- full.ase(g_lcc[], 10)
s_m= as.matrix(g_lcc[])

nrow(s_m)
sort(apply(s_m,2,sum))
plot(sort(apply(s_m,2,sum)), main = "Degree Distribution of Revised Model", ylab = "Degree")


largest_clqs[[3]]

ev <- eigen(s_m)
ev$values
plot(1:141,ev$values)

all_vertices <- 1:nrow(s_m)

other_vertices <- setdiff(all_vertices, clique_vertices)

# 3. Create the new order: clique vertices first, then the rest
new_order <- c(clique_vertices, other_vertices)

# 4. Reorder the matrix using the new index vector
perm_sm<- s_m[new_order, new_order]

heatmap(s_m, scale = "none",Rowv = NA, Colv = NA )
heatmap(perm_sm, scale = "none",Rowv = NA, Colv = NA )

heatmap(perm_sm, scale = "none")

par(mfrow=c(1,2))
plot(ev$vectors[,1],ev$vectors[,141], col = v_colors,main = 'M07357 LCC spectral embedding')
plot(ev$vectors[,1],ev$vectors[,2], col = v_colors, main = '')


plot(ase$Xhat[,1],ase$Xhat[,3])
ase$Xhat[,c(1,2)]


plot(1:n, ase$Xhat[,1])
plot(1:n, eigen(perm_sm)$vectors[,1])

