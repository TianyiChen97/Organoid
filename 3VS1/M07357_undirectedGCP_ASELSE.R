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
library(igraph)
library(Matrix)
library(iGraphMatch) # Assuming full.ase and getElbows are from this package
library(scatterplot3d)
library(stringr)

extract_network_path <- function(input_string) {
  # This regular expression looks for the literal string "Network/"
  # followed by one or more digits (\\d+).
  pattern <- "Network/\\d+"
  
  # str_extract will find and return the first match of the pattern.
  result <- str_extract(input_string, pattern)
  
  return(result)
}
extract_path_part <- function(filepath) {
  # Use sub() to find the pattern and return only that part
  sub(".*/(M\\d+/[A-Za-z]+?/\\d+)/.*", "\\1", filepath)
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



full.lse <- function(A, d) {
  A <- A[]
  A = as.matrix(A)
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

get_clique_density_stats <- function(A) {
  # 1. Convert to Graph
  A_mat <- as.matrix(A)
  g <- graph_from_adjacency_matrix(A_mat, mode = "undirected", diag = FALSE)
  
  # 2. Identify Core Nodes 
  kc <- coreness(g)
  core_index <- which(kc == max(kc))
  
  # 3. Identify Periphery Nodes
  all_indices <- 1:nrow(A_mat)
  periph_index <- setdiff(all_indices, core_index)
  
  # 4. Get Node Counts
  n_c <- length(core_index)
  n_p <- length(periph_index)
  
  # 5. Calculate Edges
  edges_CC <- if(n_c > 1) sum(A_mat[core_index, core_index]) / 2 else 0
  edges_PP <- if(n_p > 1) sum(A_mat[periph_index, periph_index]) / 2 else 0
  edges_CP <- if(n_c > 0 && n_p > 0) sum(A_mat[core_index, periph_index]) else 0
  
  # 6. Calculate Max Possible Edges
  max_CC <- n_c * (n_c - 1) / 2
  max_PP <- n_p * (n_p - 1) / 2
  max_CP <- n_c * n_p
  
  # 7. Calculate Densities
  dens_CC <- ifelse(max_CC > 0, edges_CC / max_CC, 0)
  dens_PP <- ifelse(max_PP > 0, edges_PP / max_PP, 0)
  dens_CP <- ifelse(max_CP > 0, edges_CP / max_CP, 0)
  
  # 8. Create Matrix
  M <- matrix(c(dens_CC, dens_CP, 
                dens_CP, dens_PP), 
              nrow=2, byrow=TRUE)
  colnames(M) <- c("C", "P")
  rownames(M) <- c("C", "P")
  
  # Return a list with all stats
  return(list(matrix = M, n_core = n_c, n_periph = n_p))
}

#edge_list = read.csv("~/Desktop/Research /PhDresearch/Hopkins_Organoid/M07357/adjacency_edges_Aug_28_2025_M07357.csv")
#edge_list = read.csv("~/Desktop/Research /PhDresearch/Hopkins_Organoid/M07357/adjacency_edges_Aug_29_2025_M07357_try2.csv")
#edge_list = read.csv("~/Desktop/Research /PhDresearch/Hopkins_Organoid/M07357/adjacency_edges_Aug_29_2025_M07357_para2.csv")
#edge_list = read.csv("~/Desktop/Research /PhDresearch/Hopkins_Organoid/M07357/adjacency_edges_Sep_1_2025_M07357_para3.csv")


# --- FILE LOADING ---
active_file_path <- "~/Desktop/Research /PhDresearch/Hopkins_Organoid/M07357/adjacency_edges_Aug_29_2025_M07357_para2.csv"
edge_list <- read.csv(active_file_path)

# --- PARAMETER EXTRACTION ---
parameter_choice <- sub(".*M07357_([a-zA-Z0-9]+)\\.csv", "\\1", active_file_path)
if (nchar(parameter_choice) > 10 | parameter_choice == active_file_path) {
  parameter_choice <- "Default" 
}

filenames = unique(edge_list$File)


EMBEDDING_METHOD <- "LSE" # <--- **CHANGE THIS TO "LSE" WHEN NEEDED**

TITLE_SUFFIX <- paste0("PCinSS:", parameter_choice, " + ", EMBEDDING_METHOD) 

MAIN_IDENTIFIER <- "M07357"

if (EMBEDDING_METHOD == "LSE") {
  # Assuming your LSE function is defined as 'full.lse'
  embedding_function <- full.lse 
} else {
  # Default to ASE
  embedding_function <- full.ase
}






## Basic graph statistics here.



## Embedding result here

for (file in filenames) {
  
  # A. Define the part of the title that comes from the file (e.g., "Network/000389")
  path_part <- extract_network_path(file)
  
  # B. Construct the full title string used for the MTEXT and the PDF name
  # ORDER: M07357/Network/000389 | PCinSS:para2 + ASE
  full_title_string <- paste(MAIN_IDENTIFIER, path_part, sep = '/')
  full_title_string <- paste(full_title_string, TITLE_SUFFIX, sep = " | ")
  
  # C. Create the intended PDF file name (replace invalid characters with '_')
  # Example: "M07357_Network_000389_PCinSS_para2_ASE.pdf"
  pdf_filename_base <- full_title_string
  pdf_filename_base <- gsub("/", "_", pdf_filename_base) 
  pdf_filename_base <- gsub(" \\| ", "_", pdf_filename_base) 
  pdf_filename_base <- gsub(": ", "_", pdf_filename_base) 
  pdf_filename_base <- gsub(" \\+ ", "_", pdf_filename_base) 
  pdf_filename <- paste0(pdf_filename_base, '.pdf')
  
  output_dir <- 'M07357_Embed_outputs'
  
  # E. Create the folder if it doesn't exist
  if (!dir.exists(output_dir)) {
    dir.create(output_dir)
  }
  
  # F. Construct the final absolute path for the PDF
  final_pdf_path <- file.path(output_dir, pdf_filename)
  
  # G. Open the PDF device
  pdf(file = final_pdf_path, width = 10, height = 7)
  
  par(mfrow = c(2, 3), mar = c(3, 3, 3, 1), oma = c(4, 0, 3, 0))
  
  # Define the list of wells to iterate over
  well_list <- paste0("well", sprintf("%03d", 0:5))
  
  # --- Main Loop ---
  for (well in well_list) {
    
    # 1. Subset the data for the current file and well
    subset_data <- subset(edge_list, File == file & Well == well)
    
    # Initial check for any data at all
    if (nrow(subset_data) == 0) {
      plot(1, type = "n", axes = FALSE, xlab = "", ylab = "")
      title(main = well)
      text(1, 1, "No data available", col = "red")
      next # Move to the next well
    }
    
    # --- Graph Construction and Embedding ---
    tryCatch({
      # Create the adjacency matrix
      n_rows <- n_cols <- unique(subset_data$dim)
      adj <- sparseMatrix(
        i    = subset_data$Row + 1,
        j    = subset_data$Column + 1,
        dims = c(n_rows, n_cols)
      )
      
      # Create the graph and find the largest connected component (LCC)
      g <- graph_from_adjacency_matrix(adj, mode = "undirected", diag = FALSE)
      comps <- components(g, mode = "weak")
      big <- which.max(comps$csize)
      g_lcc <- induced_subgraph(g, V(g)[comps$membership == big])
      
      # 2. Check if the size of the LCC is less than 11
      if (vcount(g_lcc) < 11) {
        # If so, plot the specified message and skip the analysis
        plot(1, type = "n", axes = FALSE, xlab = "", ylab = "")
        title(main = well)
        text(1, 1, "Size of largest\nconnected component < 11", col = "red", cex = 0.9)
      } else {
        # --- Proceed with Embedding only if LCC size is 11 or greater ---
        
        # **MODIFIED:** Call the selected embedding function (full.ase or full.lse)
        embed_result <- embedding_function(g_lcc, 10)
        
        # Find the elbow to determine the embedding dimension
        elb <- getElbows(embed_result$eval, plot = F)
        
        if (length(elb) < 2) {
          stop("Embedding dimension is less than 2.")
        }
        
        embedding_dim <- elb[2]
        # --- Conditional Plotting (uses Xhat from the result) ---
        
        # 3. If the embedding dimension is 2, create a 2D plot
        if (embedding_dim == 2) {
          plot(embed_result$Xhat[, 1], embed_result$Xhat[, 2], 
               main = well, 
               xlab = "Dimension 1", 
               ylab = "Dimension 2",
               pch = 16,
               col = "blue")
        } 
        # 4. If the embedding dimension is 3 or more, create a 3D plot
        else if (embedding_dim >= 3) {
          scatterplot3d(
            x = embed_result$Xhat[, 1], 
            y = embed_result$Xhat[, 2], 
            z = embed_result$Xhat[, 3],
            main = well,
            xlab = "Dim 1",
            ylab = "Dim 2",
            zlab = "Dim 3",
            color = "blue",
            pch = 16,
            type = "p",
            grid = TRUE,
            box = TRUE
          )
        } 
        # Handle cases where the dimension is less than 2
        else {
          plot(1, type = "n", axes = FALSE, xlab = "", ylab = "")
          title(main = well)
          text(1, 1, paste("Embedding dim < 2\n(dim =", embedding_dim, ")"), col = "orange")
        }
      }
      
    }, error = function(e) {
      # If any other error occurs, plot a placeholder with the error message
      plot(1, type = "n", axes = FALSE, xlab = "", ylab = "")
      title(main = well)
      text(1, 1, "An error occurred", col = "red", cex = 0.8)
      message(paste("Error processing", well, ":", e$message))
    })
  }
  
  # --- Add the Main Title to the Entire Figure ---
  mtext(full_title_string, 
        side = 3, 
        line = 1, 
        outer = TRUE, 
        cex = 1.5, 
        font = 2) 
  
  dev.off()
}



# Extract the parameter choice
parameter_choice <- sub(".*M07357_([a-zA-Z0-9]+)\\.csv", "\\1", active_file_path)
if (nchar(parameter_choice) > 10 | parameter_choice == active_file_path) {
  parameter_choice <- "Default" 
}
MAIN_IDENTIFIER <- "M07357"

# New Suffix for the Statistics Plot
TITLE_SUFFIX <- paste0("PCinSS:", parameter_choice, " + CP_STATS")

# ------------------------------------------------------------------
# 2. CP STATISTICS PLOTTING LOOP (MODIFIED for Saving and Title)
# ------------------------------------------------------------------

for (file in filenames) {
  
  # ----------------------------------------------------
  # 2. Filename and Directory Logic
  # ----------------------------------------------------
  
  # A. Define the part of the title that comes from the file (e.g., "Network/000389")
  path_part <- extract_network_path(file)
  
  # B. Construct the full title string for the figure
  # ORDER: M07357/Network/000389 | PCinSS:para2 + CP_STATS
  full_title_string <- paste(MAIN_IDENTIFIER, path_part, sep = '/')
  full_title_string <- paste(full_title_string, TITLE_SUFFIX, sep = " | ")
  
  # C. Create the intended PDF file name (replace invalid characters with '_')
  # Example: "M07357_Network_000389_PCinSS_para2_CP_STATS.pdf"
  pdf_filename_base <- full_title_string
  pdf_filename_base <- gsub("/", "_", pdf_filename_base) 
  pdf_filename_base <- gsub(" \\| ", "_", pdf_filename_base) 
  pdf_filename_base <- gsub(": ", "_", pdf_filename_base) 
  pdf_filename_base <- gsub(" \\+ ", "_", pdf_filename_base) 
  pdf_filename <- paste0(pdf_filename_base, '.pdf')
  
  # D. Define the output folder (as requested: "M07357_CP_STATS_outputs")
  output_dir <- "M07357_CP_STATS_outputs" 
  
  # E. Create the folder if it doesn't exist
  if (!dir.exists(output_dir)) {
    dir.create(output_dir)
  }
  
  # F. Construct the final absolute path for the PDF
  final_pdf_path <- file.path(output_dir, pdf_filename)
  
  # G. Open the PDF device
  pdf(file = final_pdf_path, width = 10, height = 7)
  
  # Setup grid of plots
  par(mfrow = c(2, 3), mar = c(1, 1, 3, 1), oma = c(2, 0, 3, 0))
  
  well_list <- paste0("well", sprintf("%03d", 0:5))
  
  for (well in well_list) {
    
    subset_data <- subset(edge_list, File == file & Well == well)
    
    # 1. Setup Blank Plot Canvas
    plot(1, type = "n", axes = FALSE, xlab = "", ylab = "", xlim=c(0, 10), ylim=c(0, 10))
    title(main = well, cex.main = 1.5)
    box() # Optional: puts a box around the stats
    
    if (nrow(subset_data) == 0) {
      text(5, 5, "No data available", col = "red", cex = 1.2)
      next 
    }
    
    tryCatch({
      # Create adjacency matrix and LCC
      n_rows <- n_cols <- unique(subset_data$dim)
      adj <- sparseMatrix(
        i    = subset_data$Row + 1,
        j    = subset_data$Column + 1,
        dims = c(n_rows, n_cols)
      )
      
      # Original Graph
      g <- graph_from_adjacency_matrix(adj, mode = "undirected", diag = FALSE)
      n_original_nodes <- vcount(g) 
      
      comps <- components(g, mode = "weak")
      big <- which.max(comps$csize)
      g_lcc <- induced_subgraph(g, V(g)[comps$membership == big])
      n_lcc_nodes <- vcount(g_lcc) 
      
      if (n_lcc_nodes < 11) {
        text(5, 5, "LCC size < 11", col = "red", cex = 1.2)
      } else {
        
        # --- CALCULATE STATISTICS ---
        A_lcc <- as.matrix(as_adjacency_matrix(g_lcc))
        # This requires the function get_clique_density_stats() to be defined
        stats_res <- get_clique_density_stats(A_lcc)
        
        dens_mat <- stats_res$matrix
        n_core   <- stats_res$n_core
        n_periph <- stats_res$n_periph
        
        # Ratios
        ratio_core <- n_core / n_lcc_nodes
        ratio_peri <- n_periph / n_lcc_nodes
        
        # --- DISPLAY TEXT ---
        lines_to_print <- c(
          "--- Densities ---",
          paste0("CC: ", sprintf("%.3f", dens_mat[1,1])),
          paste0("CP: ", sprintf("%.3f", dens_mat[1,2])),
          paste0("PP: ", sprintf("%.3f", dens_mat[2,2])),
          "",
          "--- Counts ---",
          paste0("Core: ", n_core),
          paste0("Peri: ", n_periph),
          paste0("LCC Size: ", n_lcc_nodes),
          paste0("Orig Size: ", n_original_nodes),
          "",
          "--- Ratios ---",
          paste0("Core / LCC: ", sprintf("%.3f", ratio_core)),
          paste0("Peri / LCC: ", sprintf("%.3f", ratio_peri))
        )
        
        text(5, 5, labels = paste(lines_to_print, collapse = "\n"), 
             cex = 1.3, family = "mono") 
      }
      
    }, error = function(e) {
      text(5, 5, "Error", col = "red")
      message(paste("Error processing", well, ":", e$message))
    })
  }
  
  # --- Add the Main Title to the Entire Figure ---
  # Use the full title string including M07357, path, and stats type
  mtext(full_title_string, 
        side = 3, 
        line = 1, 
        outer = TRUE, 
        cex = 1.5, 
        font = 2)
  
  # Close PDF device
  dev.off()
}






####### BElow is ASE with CP stats

for (file in filenames) {
  par(mfrow = c(2, 3), mar = c(3, 3, 3, 1), oma = c(4, 0, 3, 0))
  
  well_list <- paste0("well", sprintf("%03d", 0:5))
  
  for (well in well_list) {
    
    subset_data <- subset(edge_list, File == file & Well == well)
    
    if (nrow(subset_data) == 0) {
      plot(1, type = "n", axes = FALSE, xlab = "", ylab = "")
      title(main = well)
      text(1, 1, "No data available", col = "red")
      next 
    }
    
    tryCatch({
      # Create adjacency matrix and LCC
      n_rows <- n_cols <- unique(subset_data$dim)
      adj <- sparseMatrix(
        i    = subset_data$Row + 1,
        j    = subset_data$Column + 1,
        dims = c(n_rows, n_cols)
      )
      
      # Original Graph
      g <- graph_from_adjacency_matrix(adj, mode = "undirected", diag = FALSE)
      n_original_nodes <- vcount(g) 
      
      comps <- components(g, mode = "weak")
      big <- which.max(comps$csize)
      g_lcc <- induced_subgraph(g, V(g)[comps$membership == big])
      n_lcc_nodes <- vcount(g_lcc) 
      
      if (n_lcc_nodes < 11) {
        plot(1, type = "n", axes = FALSE, xlab = "", ylab = "")
        title(main = well)
        text(1, 1, "Size of largest\nconnected component < 11", col = "red", cex = 0.9)
      } else {
        
        # --- CALCULATE DENSITY & COUNTS ---
        A_lcc <- as.matrix(as_adjacency_matrix(g_lcc))
        
        # Get stats
        stats_res <- get_clique_density_stats(A_lcc)
        dens_mat <- stats_res$matrix
        n_core   <- stats_res$n_core
        n_periph <- stats_res$n_periph
        
        # Calculate Ratios
        ratio_core <- n_core / n_lcc_nodes
        ratio_peri <- n_periph / n_lcc_nodes
        
        # Updated Legend Text
        mat_text <- c(
          paste("CC:", sprintf("%.2f", dens_mat[1,1])),
          paste("CP:", sprintf("%.2f", dens_mat[1,2])),
          paste("PP:", sprintf("%.2f", dens_mat[2,2])),
          paste("Core:", n_core),
          paste("Peri:", n_periph),
          paste("LCC:", n_lcc_nodes),
          paste("Orig:", n_original_nodes),
          paste("C/LCC:", sprintf("%.2f", ratio_core)), # Added Ratio
          paste("P/LCC:", sprintf("%.2f", ratio_peri))  # Added Ratio
        )
        
        # --- Perform ASE ---
        ase <- full.ase(g_lcc, 10)
        elb <- getElbows(ase$eval, plot = F)
        
        if (length(elb) < 2) {
          stop("Embedding dimension is less than 2.")
        }
        
        embedding_dim <- elb[2]
        
        # --- Plotting ---
        if (embedding_dim == 2) {
          plot(ase$Xhat[, 1], ase$Xhat[, 2], 
               main = well, 
               xlab = "Dimension 1", 
               ylab = "Dimension 2",
               pch = 16,
               col = "blue")
          
          # Legend with extra stats
          legend("topright", legend = mat_text, bty = "n", cex = 1.1, title="Stats")
          
        } else if (embedding_dim >= 3) {
          sp <- scatterplot3d(
            x = ase$Xhat[, 1], 
            y = ase$Xhat[, 2], 
            z = ase$Xhat[, 3],
            main = well,
            xlab = "Dim 1",
            ylab = "Dim 2",
            zlab = "Dim 3",
            color = "blue",
            pch = 16,
            type = "p",
            grid = TRUE,
            box = TRUE
          )
          
          legend("topright", 
                 legend = mat_text, bty = "n", 
                 cex = 1.1, 
                 title="Stats",
                 inset=c(0.05, 0))
          
        } else {
          plot(1, type = "n", axes = FALSE, xlab = "", ylab = "")
          title(main = well)
          text(1, 1, paste("Embedding dim < 2\n(dim =", embedding_dim, ")"), col = "orange")
        }
      }
      
    }, error = function(e) {
      plot(1, type = "n", axes = FALSE, xlab = "", ylab = "")
      title(main = well)
      text(1, 1, "An error occurred", col = "red", cex = 0.8)
      message(paste("Error processing", well, ":", e$message))
    })
  }
  
  mtext(extract_network_path(file), 
        side = 3, 
        line = 1, 
        outer = TRUE, 
        cex = 1.5, 
        font = 2)
}
