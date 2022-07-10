library(ggthemes) # themes
library(patchwork) # compose ggplots
library(rcartocolor) # color palette
library(RColorBrewer) # color palette
library(ggforce) # plotting

theme_publication <- function(base_size = 14, 
                              base_family = 'Arial', 
                              axis = T, 
                              grid = F, 
                              legend = F) {
  
  if(axis) {
    axis_element <- element_line(color = 'black')
  } else {
    axis_element <- element_blank()
  }
  
  if(grid) {
    grid_element <- element_line(color = 'grey90')
  } else {
    grid_element <- element_blank()
  }
  
  if(legend) {
    legend_position = 'right'
  } else {
    legend_position = 'none'
  }
  
  theme_foundation(base_size = base_size, base_family = base_family) + 
    theme(plot.title = element_text(face = "bold", size = 18, hjust = 0.5),
          text = element_text(),
          panel.background = element_rect(color = NA),
          plot.background = element_rect(color = NA),
          panel.border = element_blank(),
          axis.title = element_text(face = "plain", size = rel(1)),
          axis.title.y = element_text(angle=90, vjust = 0.5),
          axis.title.x = element_text(vjust = 0),
          axis.text = element_text(color = 'black'),
          axis.line = axis_element,
          panel.grid.major = grid_element,
          panel.grid.minor = element_blank(),
          legend.key = element_rect(color = NA),
          legend.key.size = unit(0.4, "cm"),
          legend.background = element_rect(fill = NA),
          legend.spacing = unit(0, "cm"),
          legend.title = element_blank(), #element_text(face="italic"),
          legend.position = legend_position,
          strip.background = element_rect(color="#f0f0f0",fill="#f0f0f0"),
          strip.text = element_text(face="bold"))
}

theme_publication2 <- function(base_size = 6, 
                              base_family = 'Arial', 
                              axis = T, 
                              grid = F, 
                              legend.position = 'top') {
  
  if(axis) {
    axis_element <- element_line(color = 'black')
  } else {
    axis_element <- element_blank()
  }
  
  if(grid) {
    grid_element <- element_line(color = 'grey90')
  } else {
    grid_element <- element_blank()
  }

  
  theme_foundation(base_size = base_size, base_family = base_family) + 
    theme(plot.title = element_text(face = 'plain', size = 8, hjust = 0.5),
          text = element_text(),
          panel.background = element_rect(color = NA),
          plot.background = element_rect(color = NA),
          panel.border = element_blank(),
          axis.title = element_text(size = 7),
          axis.title.y = element_text(angle=90, vjust = 0.5),
          axis.title.x = element_text(vjust = 0),
          axis.text = element_text(color = 'black', size = 6),
          axis.line = element_line(size = 0.5),
          panel.grid.major = grid_element,
          panel.grid.minor = element_blank(),
          legend.key = element_rect(color = NA),
          legend.key.size = unit(0.4, "cm"),
          legend.background = element_rect(fill = NA),
          legend.spacing = unit(0, "cm"),
          legend.title = element_blank(),
          legend.position = legend.position,
          strip.background = element_rect(color="#f0f0f0",fill="#f0f0f0"),
          strip.text = element_text(face="bold"))
}

save_figure <- function(plot = last_plot(), 
                        filename, 
                        device = c('png', 'pdf'), 
                        dpi = 600, 
                        w = 5, 
                        h = 5, 
                        units = 'in', 
                        directory = 'figures', 
                        size = 'custom', 
                        gg = TRUE,
                        overwrite = TRUE) {
  
  # Some default sizes
  if(size == 'small') {
    w = 5
    h = 5
  } else if(size == 'medium') {
    w = 7
    h = 7
  } else if(size == 'large'){
    w = 10
    h = 10
  } else if(size == 'wide') {
    h = 5
    w = 8
  } else if(size == 'xwide') {
    h = 5
    w = 12
  } else if(size == 'long') {
    h = 12
    w = 8
  } else if(size == 'xlong') {
    h = 12
    w = 5
  }
  
  # prevent overwrites
  i = 0
  prefix <- filename
  map(device, function(x) {
    if(!overwrite){
      while(file.exists(paste0(directory, filename,'.', x)))
      {
        i = i + 1
        ifilename = paste0(prefix,i)
      }
    }
    cat('Saving as', paste0(filename, '.', x), '\n')
    
    filename = paste0(filename, '.', x)
    if(gg) { # for ggplot objects
      if(x == 'pdf') {
        x <- cairo_pdf
      }
      ggsave(filename = filename, plot = plot, dpi = dpi, path = directory, width = w, height = h, units = units)}
    else { # for others
      dev.copy(png, file = paste0(directory, filename), width = w, height = h, units = units, res = dpi)
      dev.off() }
    
  })
}