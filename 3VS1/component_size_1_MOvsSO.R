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

#Single organoid plate – M08438 
#Multi-organoid plate – M07359 


stage0 <- c('000386', '000384')  # ActivityScan
stage1 <- c('000387', '000385')  # Baseline #1 (pre‑stim)
stage2 <- c('000395', '000390')  # Baseline #2 (pre‑stim)
stage3 <- c('000396', '000391')  # During stimulation #1 (random‑site)
stage4 <- c('000397', '000392')  # After stimulation #1
stage5 <- c('000398', '000393')  # During stimulation #2 (single‑region)
stage6 <- c('000399', '000394')  # After stimulation #2

get_stage <- function(file_path) {
  if (grepl(paste(stage0, collapse = "|"), file_path)) {
    return(0)
  } else if (grepl(paste(stage1, collapse = "|"), file_path)) {
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


extract_chip_id <- function(path) {
  sub(
    ".*?/([^/]+)/(?:Network|Stimulation|ActivityScan)/.*",
    "\\1",
    path,
    perl = TRUE
  )
}

files <- filenames[3:8] ## M08438
files <- filenames[9:13] ## M07359




wells <- paste0("well", sprintf("%03d", 0:5))

library(purrr)

# build & sort the index using base R
index_df <- expand.grid(
  file = files,
  well = wells,
  stringsAsFactors = FALSE
)

# get stage as numeric
index_df$stage <- sapply(index_df$file, get_stage)

# extract the numeric suffix of 'well'
index_df$well_num <- as.integer(sub("well", "", index_df$well))

# reorder rows by stage, then well_num
index_df <- index_df[order(index_df$stage, index_df$well_num), ]
rownames(index_df) <- NULL
# now loop in that order
plots <- vector("list", nrow(index_df))
for (i in seq_len(nrow(index_df))) {
  file  <- index_df$file[i]
  well  <- index_df$well[i]
  stage <- index_df$stage[i]
  Chip  <- extract_chip_id(file)
  
  subset_data <- subset(edge_list, File == file & Well == well)
  if (nrow(subset_data) == 0) {
    message("Skipping empty: Stage=", stage, " Well=", well)
    next
  }
  
  n_rows <- n_cols <- unique(subset_data$dim)
  sparse_matrix <- sparseMatrix(
    i    = subset_data$Row + 1,
    j    = subset_data$Column + 1,
    dims = c(n_rows, n_cols)
  )
  g <- graph_from_adjacency_matrix(sparse_matrix,
                                   weighted = NULL,
                                   mode     = "directed",
                                   diag     = FALSE)
  
  comps  <- clusters(g, mode = "weak")
  df <- map_df(seq_along(comps$csize), function(cid) {
    verts <- V(g)[comps$membership == cid]
    subg  <- induced_subgraph(g, verts)
    tibble(
      component    = cid,
      size         = comps$csize[cid],
      density      = edge_density(subg, loops = FALSE)
    )
  })
  
  df_summary <- df %>%
    group_by(size) %>%
    summarize(
      freq         = n(),
      mean_density = mean(density, na.rm = TRUE),
      .groups      = "drop"
    ) %>%
    mutate(mean_density = ifelse(is.nan(mean_density), NA_real_, mean_density))
  
  p <- ggplot(df_summary, aes(size, freq, fill = mean_density)) +
    geom_col(width = 0.5, colour = "black", show.legend = FALSE) +
    geom_text(
      aes(x = size + .5 ,y = freq / 2 + 10 , label = freq),
      size = 3,
      color = 'red'
    ) +
    scale_fill_gradient(
      low      = "grey",
      high     = "black",
      na.value = "white"
    ) +
    labs(
      x     = "Size",
      y     = "Freq",
      title = sprintf("Stage=%s %s  n=%d", stage, well, vcount(g))
    ) +
    theme_minimal(base_size = 8) +
    theme(
      axis.title     = element_text(size = 8),
      axis.text      = element_text(size = 6),
      plot.title     = element_text(size = 9),
      panel.border   = element_rect(colour = "black", fill = NA),
      panel.background = element_blank()
    )
  plots[[i]] <- p
}

dfp <- tibble(
  stage = index_df$stage,
  plot  = plots
) %>%
  filter(!map_lgl(plot, is.null))

df_rows <- dfp %>%
  group_by(stage) %>%
  summarise(plots_list = list(plot), .groups = "drop")
df_rows %>%
  mutate(n_plots = lengths(plots_list)) %>%      # how many plots per stage?
  arrange(stage)


df_rows <- df_rows %>%
  mutate(
    padded = purrr::map(plots_list, ~{
      pls <- .x
      if (length(pls) < 6) {
        pls <- c(pls, rep(list(plot_spacer()), 6 - length(pls)))
      }
      pls
    })
  )

stage_rows <- purrr::map(df_rows$padded, ~ wrap_plots(.x, ncol = 6))

final_plot <- wrap_plots(stage_rows, ncol = 1) +
  plot_annotation(
    title    = extract_chip_id(file[1])
  ) 

# Draw it
final_plot


