
# define the sets
files <- filenames[1:6] ## M07915
files <- filenames[7:12] ## M07914

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

df_rows <- df_rows %>%
  mutate(
    padded = map(plots_list, ~{
      pls <- .x
      if (length(pls) < 6) {
        pls <- c(pls, rep(list(plot_spacer()), 6 - length(pls)))
      }
      pls
    })
  )

stage_rows <- map(df_rows$padded, ~ wrap_plots(.x, ncol = 6))

final_plot <- wrap_plots(stage_rows, ncol = 1) +
  plot_annotation(
    title    = "M07914"
  ) 

# Draw it
final_plot
files
