library(cowplot)

seurat_statsplot <- function(seurat, x = 'nCount_RNA', y = 'nFeature_RNA', marginal.type = 'density', size = 0.5, alpha = .05,
                             title = seurat@project.name, cells = NULL, logx = T, logy = T, vline = NULL, hline = NULL, limits = c(100, NULL), median_lines = T, ...){
  
  df <- FetchData(seurat, vars = c(x, y), cells = cells, ...)
  
  library(ggExtra)
  p <- df %>% ggplot(aes(x = df[,1], y = df[,2])) + geom_jitter(alpha = alpha, color = 'dodgerblue4', size = size) + xlab(x) + ylab(y) + theme_publication()
  
  x_max <- max(df[,x])
  y_max <- max(df[,y])
  
  x_median <- median(df[,1])
  y_median <- median(df[,2])
  
  if(logx) {
    p <- p + scale_x_continuous(trans = 'log10', limits = c(limits[1], x_max))
  }
  if(logy) {
    p <- p + scale_y_continuous(trans = 'log10', limits = c(limits[1], y_max))
  }
  if(!is.null(vline)) {
    p <- p + geom_vline(xintercept = vline, linetype = 'dashed')
  }
  if(!is.null(hline)) {
    p <- p + geom_hline(yintercept = hline, linetype = 'dashed')
  }
  if(median_lines) {
    p <- p + geom_vline(xintercept = x_median, linetype = 'dashed', alpha = 0.5) +
      geom_hline(yintercept = y_median, linetype = 'dashed', alpha = 0.5) +
      annotate('text', x = x_median, y = y_max, label = x_median, hjust = 1.5) +
      annotate('text', x = x_max, y = y_median, label = y_median, vjust = -0.5)
  }
  if(!is.null(title)) {
    p <- p + ggtitle(title)
  }
  
  ggMarginal(p, type = marginal.type, margins = "both", size = 4, fill = c('#D55E00'))
  
}

seurat_embedmeta <- function(seurat, 
                             feature, 
                             reduction = 'umap', 
                             dims = 1:2, 
                             slot = 'data') {
  embedding <- Seurat::Embeddings(seurat, reduction = reduction)[ , dims]
  df <- as.data.frame(cbind(embedding, Seurat::FetchData(seurat, feature, slot = slot)))
  df
}

seurat_plot <- function(seurat, 
                        feature = 'ident', 
                        size = 0.1, 
                        alpha = 0.5, 
                        show.legend = F,
                        regress = F, 
                        cells.use = NULL, 
                        do.label = T, 
                        reduction = 'umap', 
                        dims = 1:2,
                        color_package = 'hues', 
                        color_scale = 'pimp', 
                        title_params = 'italics',
                        label_color = c('white', 'black'),
                        facet = NULL,
                        ncol = NULL,
                        ...){
  
  if(!is.null(facet)) {
    df <- seurat_embedmeta(seurat, c(feature, facet), reduction = reduction, dims = dims)
    colnames(df) <- c('dim1', 'dim2', 'feature', 'facet')  
  } else {
    df <- seurat_embedmeta(seurat, feature, reduction = reduction, dims = dims)
    colnames(df) <- c('dim1', 'dim2', 'feature')  
  }
  
  if (reduction == 'pca') {
    xlabel <- paste0('PC', dims[1])
    ylabel <- paste0('PC', dims[2])
  } else if (reduction == 'tsne') {
    xlabel <- bquote(italic('t')*'-SNE'*.(dims[1]))
    ylabel <- bquote(italic('t')*'-SNE'*.(dims[2]))
  } else {
    xlabel <- paste0(toupper(reduction), dims[1])
    ylabel <- paste0(toupper(reduction), dims[2])
  }
  
  if(title_params == 'italics') {
    title_params <- element_text(face = 'italic', size = rel(1))
  } else if(title_params == 'bold') {
    title_params <- element_text(face = "bold", size = rel(1.2), hjust = 0.5)
  }
  
  plot <- ggplot2::ggplot(df, aes(x = dim1, y = dim2, color = feature)) +
    geom_point(size = size, alpha = alpha) + xlab(xlabel) + ylab(ylabel) +
    theme_publication(axis = FALSE, grid = FALSE) +
    theme(axis.text = element_blank(),
          axis.ticks = element_blank())
  
  if (color_package == 'carto') {
    plot <- plot + rcartocolor::scale_color_carto_d(palette = color_scale)
  } else if (color_package == 'hues') {
    num_colors <- dplyr::n_distinct(df$feature)
    cat(crayon::red('Generating', num_colors, 'colors using hues package\n'))
    plot <- plot + scale_color_manual(values = generate_colors(num_colors, palette = color_scale))
  } else if (color_package == 'custom') {
    plot <- plot + scale_color_manual(values = color_scale)
  }
  
  if (!show.legend) {
    plot <- plot + theme(legend.position = 'none')
    } else {
    plot <- plot + theme(legend.position = 'right')}
  
  
  if(do.label) {
    centers <- df %>% dplyr::group_by(feature) %>% dplyr::summarize(x = median(dim1), y = median(dim2))
    plot <- plot + shadowtext::geom_shadowtext(data = centers, aes(x = x, y = y, label = feature), color = label_color[1], size = 5, bg.color = label_color[2])
  }
  
  if(!is.null(facet)) {
    plot <- plot + 
      facet_wrap(~facet)
  }
    
  return(plot)
}

### 3d plot
seurat_3d <- function(seurat, 
                      feature = 'ident', 
                      reduction = 'umap3d', 
                      dims = 1:3,
                      size = 2, 
                      alpha = 1, 
                      showlegend = T, 
                      continuous = 'auto',
                      plot_colors = 'auto', 
                      na_break = T,
                      log_transform = F, 
                      title = feature) {
  
  # get embedding and feature to color by
  df <- seurat_embedmeta(seurat, feature, dims = dims, reduction = reduction)
  colnames(df) <- c('dim1', 'dim2', 'dim3', 'feature')
  
  if ((log_transform) & is.numeric(df$feature)) {
    df$feature <- log10(df$feature + 1)
  }
  
  # determine continous
  if (continuous == 'auto') {
    if (is.numeric(df$feature) & !is.integer(df$feature)) {
      continuous <- TRUE
    } else {
      continuous <- FALSE
    }
  }
  
  if (plot_colors == 'auto') {
    if(continuous) {
      plot_colors <- c('grey80', 'dodgerblue4')
    } else {
      plot_colors <- generate_colors(dplyr::n_distinct(df$feature), palette = 'pimp')
    }
  }
  
  # plot
  xlabel <- paste0(toupper(reduction), dims[1])
  ylabel <- paste0(toupper(reduction), dims[2])
  zlabel <- paste0(toupper(reduction), dims[3])
  
  if (continuous) {
    cat(red('Continuous scale\n'))
    p <- plot_ly(df,
                 x = ~dim1,
                 y = ~dim2,
                 z = ~dim3,
                 color = ~feature,
                 marker = list(colorscale = plot_colors,
                               showscale = T,
                               size = size,
                               opacity = alpha))
    
  } else {
    cat(red('Discrete scale\n'))
    
    df$feature <- as.factor(df$feature)
    
    p <- plot_ly(df,
                 x = ~dim1,
                 y = ~dim2,
                 z = ~dim3,
                 color = ~feature,
                 colors = plot_colors,
                 marker = list(size = size,
                               opacity = alpha),
                 text = ~paste('Cluster: ', feature))
  }
  
  p <- p %>%
    layout(title = title,
           scene = list(xaxis = list(title = xlabel),
                        yaxis = list(title = ylabel),
                        zaxis = list(title = zlabel)),
           showlegend = showlegend)
  p
}

seurat_feature <- function(seurat, gene, show = T, alpha = 0.9, size = 0.7, na.break = T, reduction = 'umap', dims = 1:2,
                           title = gene, remove_axis = F, title_params = 'italics', colors = 1, scaled = FALSE, guide = F) {
  
  if (length(colors) == 1) {
    if(colors == 1) {
      colors = c('grey90', 'grey90', 'darkorchid4')
    } else if (colors == 2) {
      colors <- c('grey90','wheat3', 'darkorchid4')
    } else if (colors == 3) {
      colors <- c('grey90', 'dodgerblue4', 'orangered')
    }
  }
  
  reduction.matrix <- Seurat::Embeddings(seurat, reduction = reduction)[,dims]
  
  if (reduction == 'pca') {
    xlabel <- paste0('PC', dims[1])
    ylabel <- paste0('PC', dims[2])
  } else if (reduction == 'tsne') {
    xlabel <- bquote(italic('t')*'-SNE'*.(dims[1]))
    ylabel <- bquote(italic('t')*'-SNE'*.(dims[2]))
  } else {
    xlabel <- paste0(toupper(reduction), dims[1])
    ylabel <- paste0(toupper(reduction), dims[2])
  }
  
  if(title_params == 'italics') {
    title_params <- element_text(face = 'italic', size = rel(1), hjust = 0.5)
  } else if(title_params == 'bold') {
    title_params <- element_text(face = "bold", size = rel(1.2), hjust = 0.5)
  }
  
  expression_matrix <- as.matrix(FetchData(seurat, gene))
  
  if(scaled) {
    expression_matrix[,1] <- as.numeric(scale(expression_matrix[,1]))
  }
  
  # Set 0's to a separate color
  if(na.break) {
    expression_matrix[expression_matrix == 0] <- NA
  }
  
  gene_frame <- data.frame(reduction.matrix, 'expression' = expression_matrix[,1]) %>% dplyr::arrange(!is.na(expression), expression)
  
  plot <- ggplot(gene_frame, aes(x = gene_frame[,1],
                                 y = gene_frame[,2],
                                 color = expression)) +
    geom_point(size = size, alpha = alpha) +
    xlab(xlabel) + ylab(ylabel) +
    theme_minimal() +
    ggtitle(title) +
    theme(plot.title = title_params,
          panel.grid = element_blank(),
          axis.text = element_blank())
  
  if (remove_axis) {
    plot <- plot + theme(axis.title = element_blank())
  }
  
  if (scaled) {
    plot <- plot +
      scale_color_gradient2(low = colors[2], mid = colors[1], high = colors[3], labels = NULL)
  } else {
    plot <- plot +
      scale_color_gradient(na.value = colors[1], low = colors[2], high = colors[3], labels = NULL)
  }
  
  if(guide) {
    plot <- plot + theme(legend.position = 'right')
  } else {
    plot <- plot + theme(legend.position = 'none')
  }
  
  if(show == T) {print(plot)}
  return(plot)
}

seurat_features <- function(seurat, gene.list, alpha = NULL, size = NULL, na.break = T, reduction = 'umap',
                            dims = 1:2, label_size = 10, nrow = NULL, ncol = NULL, titles = gene.list, remove_axis = ifelse(length(gene.list > 3), T, F),
                            title_params = 'italics', guide = F, colors = 1, scaled = FALSE, ...) {
  
  if(is.null(alpha)) {
    alpha <- ifelse(ncol(seurat) < 2000, 0.5, 0.25)
    print(c('alpha ', alpha))
  }
  if(is.null(size)) {
    size <- ifelse(ncol(seurat) < 2000, 1.5, 7500/ncol(seurat))
    print(c('size ', size))
  }
  
  plots <- lapply(1:length(gene.list), function(i){
    seurat_feature(seurat, gene.list[i], size = size, show = F, alpha = alpha, na.break = na.break, reduction = reduction,
                   dims = dims, title = titles[i], remove_axis = remove_axis, title_params = title_params, guide = guide, colors = colors, scaled = scaled, ...)
  })
  plot_grid(plotlist = plots, label_size = label_size, nrow = nrow, ncol = ncol)
}

umap_density <- function(coordinates_df, 
                         x_coord = 'UMAP_1',
                         y_coord = 'UMAP_2',
                         n = 100, 
                         expansion = 0.2) {
  
  coordinates_df <- coordinates_df %>% select(one_of(c(x_coord, y_coord)))
  colnames(coordinates_df) <- c('x_coord', 'y_coord')
  
  kdeout <- coordinates_df %>% 
    with( 
      MASS::kde2d(x_coord, 
                  y_coord, 
                  n = n,
                  lims = c(
                    scales::expand_range(range(x_coord), expansion),
                    scales::expand_range(range(y_coord), expansion)
                  )
      )
    )
  
  kde_df <- kdeout %>% 
    .[c("x", "y")] %>% 
    cross_df() %>% 
    #rename("UMAP_1" = "x", "UMAP_2" = "y") %>% 
    mutate(density = as.vector(kdeout$z)) 
  
  colnames(kde_df) <- c(x, y, 'density')
  
  kde_df
  
}