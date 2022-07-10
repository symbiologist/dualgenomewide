# reference tables to use 
sgrna_table <- read_tsv('analysis/screen/input/ScreenProcessing/library_tables/CRISPRi_v2_human_lincRNA_unique_merged_librarytable.txt')
lncrna_annotations <- read_csv('analysis/screen/input/crincl_reference/csv_tables/aah7111-TableS1.csv')
unified_metadata <- read_tsv('analysis/reference/output/02_unified_reference/unified_metadata.tsv.gz')
lhid_table <- unified_metadata %>% select(feature_id, gene_name) %>% unique()


# Runs readtable, calculate_score, fdr, mark_hits, and volcano_plot
screen_analysis <- function(folder,
                            compact = T, 
                            drop_na = T,
                            library = 'auto',
                            experiment = 'auto',
                            desired_fdr = 0.05, 
                            thresholds = seq(0, 30, by = 0.05)) {

  screen_table <- screen_readtable(folder = folder,
                                   compact = compact,
                                   drop_na = drop_na,
                                   library = library,
                                   experiment = experiment) %>% 
    screen_calculate_score() %>% 
    screen_mark_hits(desired_fdr = desired_fdr)
  
  fdr <- screen_table %>% 
    screen_fdr(desired_fdr = desired_fdr,
               thresholds = thresholds)
  
  list('table' = screen_table,
       'fdr' = fdr)
  
  
}

# Check replicates
screen_replicates <- function(folder,
                              exclude_control = TRUE,
                              correlation = TRUE,
                              trendline = TRUE) {
  
  # raw counts (convert to cpm)
  counts_table <- read_tsv(file.path(folder, 'screen_mergedcountstable.txt'),
                           skip = 3,
                           col_names = c('feature_id', 'Sample1_Rep1', 'Sample1_Rep2', 'Sample2_Rep1', 'Sample2_Rep2')) %>%
      mutate_if(is.numeric, function(.) (1e6*./sum(.)))
  
  experiment_name <- paste(basename(folder), basename(dirname(folder)))
  
  if(exclude_control) {
    counts_table <- counts_table %>% filter(str_detect(feature_id, 'non-targeting', negate = TRUE))
  }
  
  sample1_cpm <- counts_table %>% 
    scatter_stats(x = 'Sample1_Rep1',
                  y = 'Sample1_Rep2',
                  xlab = 'log10 CPM + 1 (Replicate 1)',
                  ylab = 'log10 CPM + 1 (Replicate 2)',
                  title = paste0(experiment_name, ': Sample 1'),
                  correlation = correlation,
                  trendline = trendline, 
                  log_transform = TRUE,
                  pseudocount = 1)
  
  sample2_cpm <- counts_table %>% 
    scatter_stats(x = 'Sample2_Rep1',
                  y = 'Sample2_Rep2',
                  xlab = 'log10 CPM + 1 (Replicate 1)',
                  ylab = 'log10 CPM + 1 (Replicate 2)',
                  title = paste0(experiment_name, ': Sample 2'),
                  correlation = correlation,
                  trendline = trendline, 
                  log_transform = TRUE,
                  pseudocount = 1)
  
  # phenotype 
  screen_table <- screen_readtable(folder, compact = F) %>% screen_calculate_score()
  
  if(exclude_control) {
    screen_table <- screen_table %>% filter(group == 'Gene')
  }
  
  pheno <- screen_table %>% 
    scatter_stats(x = 'rep1_pheno',
                  y = 'rep2_pheno',
                  xlab = 'Phenotype (Replicate 1)',
                  ylab = 'Phenotype (Replicate 2)',
                  title = paste0(experiment_name, ': Sample 1'),
                  correlation = correlation,
                  trendline = trendline, 
                  log_transform = FALSE)
  
  # score
  score <- screen_table %>% 
    scatter_stats(x = 'rep1_score',
                  y = 'rep2_score',
                  xlab = 'Score (Replicate 1)',
                  ylab = 'Score (Replicate 2)',
                  title = paste0(experiment_name, ': Sample 1'),
                  correlation = correlation,
                  trendline = trendline, 
                  log_transform = FALSE)
  
  list('counts' = counts_table,
       'screen_table' = screen_table,
       'Sample1' = sample1_cpm,
       'Sample2' = sample2_cpm,
       'Pheno' = pheno,
       'Score' = score)
  
} 

# Reads in table from results directory
screen_readtable <- function(folder, 
                             compact = T, 
                             drop_na = T,
                             assay = 'auto',
                             library = 'auto',
                             experiment = 'auto') {
  
  # read table and rename columns
  screen_table <- read_tsv(file.path(folder, 'screen_genetable_collapsed.txt'), 
                     skip = 4, 
                     col_names = c('feature_id', 
                                   'transcripts', 
                                   'rep1_mw', 
                                   'rep1_pheno', 
                                   'rep1_grna_count_mw', 
                                   'rep1_grna_count_pheno',
                                   'transcripts_2_dup',
                                   'rep2_mw', 
                                   'rep2_pheno', 
                                   'rep2_grna_count_mw', 
                                   'rep2_grna_count_pheno',
                                   'transcripts_3_dup',
                                   'avg_mw', 
                                   'avg_pheno', 
                                   'avg_grna_count_mw', 
                                   'avg_grna_count_pheno')) %>% 
    select(-ends_with('dup'))
  
  # compact table
  if(compact) {
    screen_table <- screen_table %>% select(feature_id, transcripts, starts_with('avg')) 
  }
  
  # remove dropouts
  if(drop_na) {
    screen_table <- screen_table %>% drop_na()
  }

  #add column to identify negative controls
  screen_table <- screen_table %>% mutate(group = ifelse(str_detect(feature_id, 'pseudo_'), 'Control', 'Gene'))
  
  # Categorize phenotypes
  screen_table <- screen_table %>% mutate(direction = ifelse(avg_pheno >= 0, 'Positive', 'Negative'))
  
  # Add experiment info
  # Automatically detect assay direction from results folder
  if (assay == 'auto') {
    assay <- basename(dirname(folder))
  }
  screen_table <- screen_table %>% mutate(assay = assay)
  
  # if lncRNA library, add information from annotation
  if(library == 'auto') {
    if(any(str_detect(screen_table$feature_id, 'LH0'))) {
      library <- 'lncRNA'
    } else {
      library <- 'Coding'
    }
  }
  
  #convert LHID to gene names from CRiNCL supp table
  if(library == 'lncRNA') {
    screen_table <- right_join(lhid_table, screen_table) %>% mutate(library = library) 
  } else {
    screen_table <- screen_table %>% mutate(gene_name = feature_id, library = library) 
  }
  
  # Automatically detect experiment name from results folder
  if (experiment == 'auto') {
    experiment <- basename(folder)
  }
  screen_table <- screen_table %>% mutate(experiment = experiment)
  
  
  screen_table %>% select(feature_id, gene_name, everything())
}

screen_calculate_score <- function(screen_table) {
  
  # detect if replicates are present
  replicates <- any(str_detect(colnames(screen_table), 'rep'))
  
  if(replicates) {
    rep1_sd <- screen_table %>% filter(group == 'Control') %>% pull(rep1_pheno) %>% sd()
    rep2_sd <- screen_table %>% filter(group == 'Control') %>% pull(rep2_pheno) %>% sd()
    
    screen_table <- screen_table %>% mutate(rep1_score = abs(rep1_pheno) * -log10(rep1_mw)/rep1_sd,
                                            rep2_score = abs(rep2_pheno) * -log10(rep2_mw)/rep2_sd)
  }
  
  # calculate average distribution from negative controls
  control_sd <- screen_table %>% filter(group == 'Control') %>% pull(avg_pheno) %>% sd()
  
  # add average score column
  screen_table <- screen_table %>% 
    mutate(control_sd = control_sd,
           score = abs(avg_pheno) * -log10(avg_mw)/control_sd) %>% 
    arrange(-score)
  
  screen_table
}

# calculates fdr from negative control distribution
screen_fdr <- function(screen_table, 
                        desired_fdr = 0.05, 
                        thresholds = seq(0, 30, by = 0.05)) {
  
  fdr_table <- map_df(thresholds, function(threshold) {
    total_hits <- screen_table %>% filter(score >= threshold) %>% nrow()
    gene_hits <- screen_table %>% filter(score >= threshold, group == 'Gene') %>% nrow()
    control_hits <- screen_table %>% filter(score >= threshold, group == 'Control') %>% nrow()
    
    list(threshold = threshold, 
         total_hits = total_hits, 
         gene_hits = gene_hits, 
         control_hits = control_hits,
         fdr = control_hits/total_hits)
  }) %>% mutate(significant = ifelse(fdr <= desired_fdr, TRUE, FALSE))
  
  score_threshold <- fdr_table %>% filter(significant) %>% pull(threshold) %>% min()
  
  list('fdr_table' = fdr_table,
       'score_threshold' = score_threshold)
}

screen_mark_hits <- function(screen_table, desired_fdr = 0.05) {
  
  threshold <- screen_fdr(screen_table, desired_fdr = desired_fdr)$score_threshold
  
  # mark hits
  screen_table <- screen_table %>% 
    mutate(threshold = threshold,
           status = ifelse(score >= threshold, 'Hit', 'Non-hit'),
           legend = paste(status, group),
           direction = ifelse(status != 'Hit', 'None', direction),
           hit_direction = ifelse(status == 'Hit', paste(direction, status), status))
  
}

screen_fdr_plot <- function(screen_table, 
                             fdr = 0.05,
                             thresholds = seq(0, 30, by = 0.05)) {
  
  fdr_table <- screen_fdr(screen_table)[['fdr_table']]
  
  fdr_table %>% 
    ggplot(aes(x = threshold, y = fdr, color = experiment)) + 
    geom_line(size = 1, alpha = 0.5) + 
    labs(x = 'Score threshold',
         y = 'FDR') +
    geom_hline(yintercept = 0.05, linetype = 'dashed', color = 'grey50')
}

screen_volcano <- function(screen_table, 
                           top_genes = 25, 
                           extra_genes = NULL, 
                           size = 1, 
                           alpha = 0.5, 
                           xlab = 'log2 enrichment', 
                           ylab = '-log10 p-value', 
                           ylim = 'auto',
                           title = '', 
                           color_by = 'legend', 
                           colors = c('darkred', 'blue', 'grey75', 'grey75'), 
                           label_size = 4, 
                           fontface = 'plain') {
  
  # Get experiment info
  library <- screen_table %>% pull(library) %>% unique()
  assay <- screen_table %>% pull(assay) %>% unique()
  experiment <- screen_table %>% pull(experiment) %>% unique()
  threshold <- screen_table %>% pull(threshold) %>% unique()
  #
  title <- paste(assay, experiment) 
  
  if(ylim == 'auto') {
    ylim <- c(0, ceiling(max(-log10(screen_table$avg_mw))))
  }
  
  color_vector <- unlist(screen_table[,color_by])
  
  plot <- screen_table %>% ggplot() + geom_point(aes(x = avg_pheno, y = -log10(avg_mw), color = color_vector), size = size, alpha = alpha) +
    xlab(xlab) + ylab(ylab) + ylim(ylim) + ggtitle(title) + 
    theme_publication() + scale_color_manual(values = colors) + 
    theme(plot.title = element_text(hjust = 0.5, face = 'plain'), 
          axis.title = element_text(hjust = 0.5), legend.title = element_blank(),
          text = element_text(size = 14)) +
    ggtitle(title)
  
  if(top_genes > 0 | !is.null(extra_genes)) {
    
    interesting_genes <- c(unique(c(screen_table %>% group_by(direction) %>% top_n(top_genes, score) %>% pull(feature_id), extra_genes)))
    if(library == 'lncRNA') {
      plot <- plot + 
        geom_text_repel(data = screen_table %>% 
                          filter(feature_id %in% interesting_genes, legend == 'hit gene'), aes(x = avg_pheno, y = -log10(avg_mw), label = gene_name), size = label_size)
    }
    else {
      plot <- plot + 
        geom_text_repel(data = screen_table %>% filter(feature_id %in% interesting_genes, legend == 'hit gene'), 
                        aes(x = avg_pheno, y = -log10(avg_mw), label = feature_id), size = label_size,
                                     fontface = fontface)
    }
    
    #calculate distribution from pseudo genes
    control_sd <- screen_table %>% filter(group == 'Control') %>% pull(avg_pheno) %>% sd()
    
    # mark hits
    score_threshold <- control_sd * threshold
    plot + stat_function(fun = function(x) {abs(score_threshold/x)}, geom = 'line', linetype = 'dashed')
  }
}


screen_origin <- function(screen_table) {
  
  sgrna_table_subset <- sgrna_table_subset %>% select(Sublibrary, feature_id = `Gene ID`) %>% unique()
  
  origin_of_hits <- inner_join(table_s2_subset, table)
  origin_of_hits %>% filter(legend == 'hit gene') %>% group_by(Sublibrary) %>% count()
}
