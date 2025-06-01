## let 's see whether there is multiple large connected component, 
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
library(dplyr)
library(ggplot2)
library(purrr)
library(patchwork)

stage1 <- c('000297', '000289')
stage2 <- c('000298', '000290')
stage4 <- c('000300', '000293')
stage6 <- c('000302', '000296')
stage3 <- c('000299', '000291')
stage5 <- c('000301', '000295')

get_stage <- function(file_path) {
  if (grepl(paste(stage1, collapse = "|"), file_path)) {
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

file = filenames[1]
well = "well000"

stage = get_stage(file)
Chip = extract_chip_id(file)

subset_data <- subset(edge_list, File == file & Well == well)
if (nrow(subset_data) == 0) {
  next
}
n_rows <-   n_cols   <- unique(subset_data$dim)

sparse_matrix <- sparseMatrix(
  i = subset_data$Row + 1,  
  j = subset_data$Column + 1,
  dims = c(n_rows, n_cols)
)

g <- graph_from_adjacency_matrix(sparse_matrix, weighted = NULL, mode = "directed", diag = FALSE)

n = vcount(g)

comps <- igraph::clusters(g, mode="weak")

# unique component IDs
comp_ids <- seq_along(comps$csize)

# compute size and density for each component
result <- lapply(comp_ids, function(cid) {
  verts    <- V(g)[comps$membership == cid]
  subg     <- induced_subgraph(g, verts)
  list(
    component = cid,
    size      = comps$csize[cid],
    density   = edge_density(subg, loops = FALSE)
  )
})

# turn into a data.frame
df <- do.call(rbind, lapply(result, as.data.frame))
rownames(df) <- NULL

df_summary <- df %>%
  group_by(size) %>%
  summarize(
    freq         = n(),                    # how many components of this size
    mean_density = mean(density, na.rm=TRUE) # average density among them
  )


ggplot(df_summary, aes(x = size, y = freq, fill = mean_density)) +
  geom_col(width = 2,
           colour = "black",    show.legend = FALSE        # <— this also hides the fill legend
    # outline each bar in black
           ) +    # bar‐width = 1 vertex
  scale_fill_gradient(
    low     = "grey",
    high    = "black",
    na.value = "white",        # <--- missing (NA) densities become pure white
    name    = "Avg. density"
  ) +
  labs(
    x = "Component size (number of vertices)",
    y = "Freq",
    title = paste0(
                  well,
                   " Stage = ", stage,
                   " Chip = ", Chip,'\nTotal num of vertices = ',n
                   )
  ) +
  theme_minimal(base_size = 8) +              # smaller base text size
  theme(
    axis.title   = element_text(size = 8),    # axis labels
    axis.text    = element_text(size = 6),    # tick labels
    plot.title   = element_text(size = 9),    # plot title
    plot.caption = element_text(size = 6),    # if you have captions
    panel.border = element_rect(colour = "black", fill = NA),
    panel.background = element_blank()
  ) 



library(Matrix)
library(dplyr)
library(igraph)
library(ggplot2)
library(patchwork)   # for wrap_plots()

# helper fns (assume you have these)
# get_stage <- function(path) { ... }
# extract_chip_id <- function(path) { ... }

plots <- list()
i_plot <- 1

for (file in filenames[1:6]) {
  stage <- get_stage(file)
  Chip  <- extract_chip_id(file)
  for (well in paste0("well", sprintf("%03d", 0:5))) {
    
    # subset & bail if empty
    subset_data <- subset(edge_list, File == file & Well == well)
    if (nrow(subset_data) == 0) {
      
      print(paste('Stage = ',stage,'  ',well))
      next
    }
    
    # build graph
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
    
    # compute components
    n      <- vcount(g)
    comps  <- clusters(g, mode = "weak")
    result <- lapply(seq_along(comps$csize), function(cid) {
      verts <- V(g)[comps$membership == cid]
      subg  <- induced_subgraph(g, verts)
      data.frame(
        component = cid,
        size      = comps$csize[cid],
        density   = edge_density(subg, loops = FALSE)
      )
    })
    df <- bind_rows(result)
    
    # summarize
    df_summary <- df %>%
      group_by(size) %>%
      summarize(
        freq         = n(),
        mean_density = mean(density, na.rm = TRUE),
        .groups      = "drop"
      ) %>%
      mutate(mean_density = ifelse(is.nan(mean_density), NA_real_, mean_density))
    
    # make plot
    p <- ggplot(df_summary, aes(x = size, y = freq, fill = mean_density)) +
      geom_col(width = 0.5 ,
               colour = "black", show.legend = FALSE      # outline each bar in black
      ) +    # bar‐width = 1 vertex
      scale_fill_gradient(
        low     = "grey",
        high    = "black",
        na.value = "white",        # <--- missing (NA) densities become pure white
        name    = "Avg. density"
      ) +
      labs(
        x = "Size",
        y = "Freq",
        title = paste0(
          well,
          " Stage = ", stage,
          " Chip = ", Chip,'n = ',n
        )
      ) +
      theme_minimal(base_size = 8) +              # smaller base text size
      theme(
        axis.title   = element_text(size = 8),    # axis labels
        axis.text    = element_text(size = 6),    # tick labels
        plot.title   = element_text(size = 9),    # plot title
        plot.caption = element_text(size = 6),    # if you have captions
        panel.border = element_rect(colour = "black", fill = NA),
        panel.background = element_blank()
      ) 
    
    
    
    plots[[i_plot]] <- p
    i_plot <- i_plot + 1
  }
}
i_plot

wrap_plots(plots, ncol = 5, nrow = 6)
