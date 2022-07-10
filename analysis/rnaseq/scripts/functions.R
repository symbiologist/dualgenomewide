# Plot gene
sleuth_gene_plot <- function(sleuth_table, 
                             feature_id) {
  
  sleuth_table_subset <- sleuth_table %>% filter(feature_id %in% feature_id)
  
  sleuth_table_means <- sleuth_table_subset %>% group_by(feature_id, condition, time) %>% summarize(tpm = mean(tpm))
  
  sleuth_table_subset %>% 
    ggplot(aes(x = time, y = tpm, color = condition)) + 
    geom_point(size = 2) + 
    geom_line(data = sleuth_table_means) + 
    scale_color_few() + 
    labs(x = 'Time',
         y = 'TPM')
}

# scale matrix (across rows)
scale_matrix_rows <- function(matrix, center = TRUE, scale = TRUE) {
  
  output <- matrix %>% 
    t() %>% # transpose
    scale(center = center, scale = scale) %>% # center and scale
    t() # transpose back
  
  output
}
# convert matrix to (long) table
matrix2table <- function(matrix, rownames = 'feature_id') {
  output <- matrix %>% 
    as.data.frame() %>% 
    rownames_to_column() %>% 
    as_tibble() %>% 
    pivot_longer(cols = -rowname, names_to = 'sample', values_to = 'scaled_expression')
  
  colnames(output)[1] <- rownames
  
  output
}

# Plotting scaled_expression clusters over time
cluster_plots <- function(input_table, 
                          which_algorithm = 'none', 
                          group = 'feature_id',
                          color_palette = 'carto', 
                          colors = 'Prism',
                          breaks = seq(0, max(input_table$time), by = 2), 
                          line_size = 0.1, 
                          line_alpha = 'auto', 
                          alpha_scaling = 1,
                          grid = FALSE,
                          ylim = c(-2, 2),
                          title_prefix = 'Cluster',
                          subtitle = FALSE,
                          exclude = TRUE,
                          hide_axis = FALSE) {
  
  # grouping variable
  input_table$group <- input_table[[group]]
  
  # use only numbered clusters
  if (which_algorithm == 'none') {
    data_subset <- input_table 
  } else {
    data_subset <- input_table %>% filter(algorithm == which_algorithm)
  }
  
  if (exclude) {
    data_subset <- data_subset %>% filter(!is.na(as.integer(cluster)))
  } 
  
  clusters <- data_subset$cluster %>% unique() %>% sort()
  
  if(color_palette == 'carto') {
    color_scheme <- rev(rcartocolor::carto_pal(n = length(clusters), name = colors))
  } else {
    color_scheme <- colors
  }
  
  map(1:length(clusters), function(i) {
    
    plot_subset <- data_subset %>% 
      filter(cluster == clusters[i])
    
    if(line_alpha == 'auto') {
      # adjust opacity to inverse square relationship of data plotted (since plotting takes up area)
      data_points <- plot_subset %>% pull(group) %>% n_distinct()
      line_alpha <- alpha_scaling * sqrt(1/data_points) 
    }
    
    p <- plot_subset %>% 
      ggplot(aes(x = time, 
                 y = scaled_expression, 
                 group = group)) + 
      geom_line(size = line_size, alpha = line_alpha, color = color_scheme[i]) + 
      theme_publication(grid = grid) +
      theme(plot.title = element_text(size = rel(1.2), face = 'plain')) +
      scale_x_continuous(breaks = breaks) + ylim(ylim) +
      ggtitle(paste(title_prefix, clusters[i])) +
      xlab('') +
      ylab('')
    
    if (subtitle) {
      subtitle <- paste0('n = ', n_distinct(plot_subset$group))
      p <- p + 
        labs(subtitle = subtitle) +
        theme(plot.subtitle = element_text(hjust = 0.5))
    } 
    
    if (hide_axis) {
      p <- p + 
        theme(axis.line = element_line(color = 'white'),
              axis.text = element_text(color = 'white'),
              axis.ticks = element_line(color = 'white'))
    }
  
      p
  })
}

# project correlations
project_clusters <- function(expression_table, # table of unassigned features, algorithm, cluster, time, and (scaled) expression
                             reference_table, # table to map onto
                             cor_method = 'pearson', # correlation method
                             threshold = 0.5,
                             unassigned = 0) { # minimum correlation threshold for assignment
  
  features_to_assign <- expression_table$feature_id %>% unique() # features to map
  
  reference_matrix <- reference_table %>% 
    filter(cluster != unassigned) %>% 
    select(cluster, scaled_expression, time) %>% 
    pivot_wider(names_from = cluster, values_from = scaled_expression) %>% 
    column_to_rownames('time') %>% 
    as.matrix()
  
  map(features_to_assign, function(i) {
    feature_subset <- arrange(filter(expression_table, feature_id == i), time)
    original_cluster <- feature_subset %>% pull(cluster) %>% unique()
    expression_vector <- feature_subset$scaled_expression
    
    correlation_results <- cor(expression_vector, reference_matrix, method = cor_method)
    max_cor <- max(correlation_results)
    top_cluster <- colnames(correlation_results)[which(correlation_results == max_cor)]  
    
    if(max_cor >= threshold) {
      cluster_assigned <- TRUE
      cluster_assignment <- top_cluster
    } else {
      cluster_assigned <- FALSE
      cluster_assignment <- unassigned 
    }
    
    tibble('feature_id' = i,
           'old_cluster' = original_cluster,
           'top_cluster' = top_cluster,
           'correlation' = max_cor,
           'assigned' = cluster_assigned,
           'new_cluster' = as.character(cluster_assignment))
  }) %>% bind_rows()
}