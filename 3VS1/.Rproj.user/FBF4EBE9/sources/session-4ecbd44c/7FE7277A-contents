library('randnet')
edge_list = read.csv('/Users/tianyichen/Desktop/Research /PhDresearch/Hopkins_Organoid/Codes/R/adjacency_edges.csv')
filenames = unique(edge_list$File)

# Define the threshold for density
density_threshold <- 0.00001  # Adjust as needed

# Initialize a data frame to store the summary
summary_table <- data.frame(
  File = character(),
  Well = character(),
  Num_Communities = character(),
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
        Num_Communities = NA
      ))
      next
    }
    
    # Estimate the number of communities using BHMC
    bhmc <- BHMC.estimate((sparse_matrix + t(sparse_matrix))/2 , 20 )
    num_communities <-  bhmc$K
    
    # Add the result to the summary table
    summary_table <- rbind(summary_table, data.frame(
      File = file,
      Well = well,
      Num_Communities = num_communities
    ))
  }
}




# Add a new column for group based on whether the filename contains "M07914" or "M07915"
summary_table$Group <- ifelse(grepl("M07914", summary_table$File), "M07914",
                              ifelse(grepl("M07915", summary_table$File), "M07915", NA))

# Remove rows with NA in Num_Communities or Group
filtered_table <- subset(summary_table, !is.na(Num_Communities) & !is.na(Group))

# Split into two groups
group_1 <- subset(filtered_table, Group == "M07914")
group_2 <- subset(filtered_table, Group == "M07915")

# Compute summary statistics for each group
group_1_mean <- mean(group_1$Num_Communities)
group_2_mean <- mean(group_2$Num_Communities)

group_1_median <- median(group_1$Num_Communities)
group_2_median <- median(group_2$Num_Communities)

# Print summary statistics
cat("Group M07914: Mean =", group_1_mean, ", Median =", group_1_median, "\n")
cat("Group M07915: Mean =", group_2_mean, ", Median =", group_2_median, "\n")

wilcox.test(group_1$Num_Communities, group_2$Num_Communities)




bhmc

num_communities <-  bhmc$K


