#file <-  '/cis/project/organoid/May_31_2025_ecr_results_no_window/250314/M07359/Network/000394/data.raw_20250601_13h43m.pkl'
well <-  'well003'


par(mfrow = c(2, 3), mar = c(3, 3, 2, 2) + 0.1) # Default is c(5.1, 4.1, 4.1, 2.1)

file <- filenames[1]
for (well in paste0("well", sprintf("%03d", 0:5))) 
well <-  'well005'

subset_data <- subset(edge_list, File==file & Well==well)

n_rows <- n_cols <- unique(subset_data$dim)
adj    <- sparseMatrix(
  i    = subset_data$Row+1,
  j    = subset_data$Column+1,
  dims = c(n_rows,n_cols)
)
g     <- graph_from_adjacency_matrix(adj, mode="undirected", diag=FALSE)
comps <- clusters(g, mode="weak")
big   <- which.max(comps$csize)
g     <- induced_subgraph(g, V(g)[comps$membership==big])

ase   <- full.ase(g[], 10)

A <- as.matrix(g[])
#svd(A)$d

elb   <- getElbows(ase$eval, plot = F)
elb
#Mclust(ase$Xhat[,1:elb[2]], G=2:Kmax, verbose= T)
plot(ase$Xhat[,1],ase$Xhat[,2])
#plot(ase$Xhat[,1],ase$Xhat[,3])
#plot(ase$Xhat[,2],ase$Xhat[,3])
library(scatterplot3d)

scatterplot3d(
  ase$Xhat[,1], ase$Xhat[,2], ase$Xhat[,3],
  main = well,
  color = "blue",        # Color of the points
  pch = 16,              # Point character (filled circle)
  type = "p",            # Type of plot: "p" for points, "l" for lines, "h" for histogram-like
  grid = TRUE,           # Add a grid to the plot
  box = TRUE             # Draw a box around the plot
)




example_GMM <- summary_table$GMM_ASE[[41]]

plot_GMM_contours(example_GMM, contour_or_not = F)

BIC_mat    <- example_GMM$BIC
modelNames <- colnames(BIC_mat)
Gvals      <- as.numeric(rownames(BIC_mat))
n          <- nrow(example_GMM$data)
d          <- ncol(example_GMM$data)
kappa      <- 7   # how harsh you want the penalty

df_mat <- matrix(
  nrow = length(Gvals),
  ncol = length(modelNames),
  dimnames = list(rownames(BIC_mat), colnames(BIC_mat))
)


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


extra_penalty <- (kappa - 1) * df_mat * log(n)
customBIC     <- BIC_mat - extra_penalty

best_idx <- which(customBIC == max(customBIC, na.rm=TRUE), arr.ind=TRUE)
best_G    <- as.numeric(rownames(customBIC)[ best_idx[1] ])
best_mod  <- colnames(customBIC)[                best_idx[2] ]
#cat("With κ =", kappa, "→ choose model", best_mod, "with G =", best_G, "\n")

gmm_custom <- Mclust(
  example_GMM$data,
  G          = best_G,
  modelNames = best_mod,
  verbose    = F
)

plot_GMM_contours(gmm_custom, contour_or_not = F)









library(dplyr)
library(tidyr)

#Single organoid plate – M08438 
#Multi-organoid plate – M07359 


summary_table <- summary_table[-7,]

# Add Group and Stage, filter valid rows
filtered_table <- summary_table %>%
  mutate(
    Group = case_when(
      grepl("M08438", File) ~ "M08438",
      grepl("M07359", File) ~ "M07359",
      TRUE ~ NA_character_
    ),
    Stage = sapply(File, get_stage)
  ) %>%
  filter(!is.na(Num_Communities_BHMC), !is.na(Group))

hist((filtered_table$Num_edges[filtered_table$Group=='M08438']))

hist((filtered_table$Num_edges[filtered_table$Group=='M07359']))

cols_to_test <- names(Filter(is.numeric, filtered_table[filtered_table$Group == "M08438", ]))
cols_to_test <- setdiff(cols_to_test, "Stage")

# Wilcoxon tests
test_results <- lapply(cols_to_test, function(col) {
  wilcox.test(
    filtered_table %>% filter(Group == "M08438") %>% pull(col),
    filtered_table %>% filter(Group == "M07359") %>% pull(col),
    paired = FALSE, exact = F
  )
})

pval_df <- data.frame(
  Column = cols_to_test,
  p_value = sapply(test_results, `[[`, "p.value")
)

# Create summary table with medians and p-values
summary_table_out <- filtered_table %>%
  group_by(Group) %>%
  summarise(across(all_of(cols_to_test), median, na.rm = TRUE)) %>%
  pivot_longer(-Group, names_to = "Column", values_to = "Median") %>%
  pivot_wider(names_from = Group, values_from = Median) %>%
  left_join(pval_df, by = "Column") %>%
  mutate(across(where(is.numeric), round, digits = 2))
summary_table_out <- filtered_table %>%
  group_by(Group) %>%
  summarise(across(all_of(cols_to_test), mean, na.rm = TRUE)) %>%
  pivot_longer(-Group, names_to = "Column", values_to = "Mean") %>%
  pivot_wider(names_from = Group, values_from = Mean) %>%
  left_join(pval_df, by = "Column") %>%
  mutate(across(where(is.numeric), round, digits = 2))

# Display
print(summary_table_out)


summary_table$Group <- ifelse(grepl("M08438", summary_table$File), "M08438",
                              ifelse(grepl("M07359", summary_table$File), "M07359", NA))

summary_table$Stage <- sapply(summary_table$File, get_stage)
filtered_table <- subset(summary_table, !is.na(Num_Communities_BHMC))
Num_community_table <- filtered_table[, c(3:12,15,16)]

grouped_summary_2 <- Num_community_table %>%
  group_by(Group, Stage) %>%
  summarise(
    across(where(is.numeric), mean, na.rm = TRUE),
    Count = n(),
    .groups = 'drop'
  )


# Function to create a table for each stage
create_stage_table <- function(stage_number) {
  data <- grouped_summary_2 %>%
    filter(Stage == stage_number) %>%
    pivot_longer(
      cols = names(.)[!(names(.) %in% c("Group", "Stage"))],
      names_to = "Variable",
      values_to = "Value"
    ) %>%
    pivot_wider(names_from = Group, values_from = Value, values_fill = list(Value = NA))
  return(data)
}


# Create tables for each stage and print them
stage_values <- sort(unique(grouped_summary_2$Stage))
stage_tables <- lapply(stage_values, create_stage_table)
names(stage_tables) <- paste0("Stage_", stage_values)


# Print all tables
for (name in names(stage_tables)) {
  print(stage_tables[[name]] %>% mutate(across(where(is.numeric), ~ round(.x, 0))))
}

summary_table[summary_table$Stage == 4, c(3:12,15,16) ]
  
