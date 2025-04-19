#' Create a Chromosome Ideogram
#'
#' @param genome_data GRanges object with chromosome information
#' @param highlight_regions Optional GRanges object with regions to highlight
#' @return A Gviz ideogram plot
#' @examples
#' # ideogram <- plot_chromosome_ideogram(genome_data)
plot_chromosome_ideogram <- function(genome_data, highlight_regions = NULL) {
  require(Gviz)
  
  # Create chromosome axis track
  axisTrack <- GenomeAxisTrack()
  
  # Create ideogram track
  ideoTrack <- IdeogramTrack(genome = "custom", chromosome = unique(seqnames(genome_data)))
  
  # Create data track if highlight regions are provided
  if (!is.null(highlight_regions)) {
    dataTrack <- AnnotationTrack(highlight_regions, 
                                 name = "Regions of Interest", 
                                 fill = "red", 
                                 col = "black")
    
    # Plot tracks together
    plotTracks(list(ideoTrack, axisTrack, dataTrack))
  } else {
    # Plot without highlight regions
    plotTracks(list(ideoTrack, axisTrack))
  }
}

#' Create Expression Heatmap
#'
#' @param expression_data Matrix or data frame of expression values (genes in rows, samples in columns)
#' @param cluster_rows Logical, whether to cluster rows
#' @param cluster_columns Logical, whether to cluster columns
#' @param annotation_data Optional data frame with sample annotations
#' @return A ComplexHeatmap object
#' @examples
#' # heatmap <- plot_expression_heatmap(expression_data, TRUE, TRUE)
plot_expression_heatmap <- function(expression_data, 
                                   cluster_rows = TRUE, 
                                   cluster_columns = TRUE,
                                   annotation_data = NULL) {
  require(ComplexHeatmap)
  require(viridis)
  
  # Convert to matrix if data frame
  if (is.data.frame(expression_data)) {
    row_names <- expression_data[[1]]  # Assuming first column has row names
    expression_matrix <- as.matrix(expression_data[, -1])
    rownames(expression_matrix) <- row_names
  } else {
    expression_matrix <- expression_data
  }
  
  # Create color mapping for heatmap
  col_fun <- viridis(100)
  
  # Create column annotation if provided
  if (!is.null(annotation_data)) {
    column_annotation <- HeatmapAnnotation(df = annotation_data)
    
    # Create heatmap with annotation
    hm <- Heatmap(expression_matrix,
                 name = "Expression",
                 col = col_fun,
                 cluster_rows = cluster_rows,
                 cluster_columns = cluster_columns,
                 top_annotation = column_annotation,
                 show_row_names = TRUE,
                 show_column_names = TRUE)
  } else {
    # Create heatmap without annotation
    hm <- Heatmap(expression_matrix,
                 name = "Expression",
                 col = col_fun,
                 cluster_rows = cluster_rows,
                 cluster_columns = cluster_columns,
                 show_row_names = TRUE,
                 show_column_names = TRUE)
  }
  
  return(hm)
}

#' Create Genome Coverage Plot
#'
#' @param coverage_data Data frame with chromosome, position, and coverage columns
#' @param chromosomes Vector of chromosomes to plot (default: all)
#' @param bin_size Size of bins for aggregating coverage (default: 1000)
#' @return A ggplot object
#' @examples
#' # coverage_plot <- plot_genome_coverage(coverage_data, c("chr1", "chr2"))
plot_genome_coverage <- function(coverage_data, chromosomes = NULL, bin_size = 1000) {
  require(tidyverse)
  
  # Filter chromosomes if specified
  if (!is.null(chromosomes)) {
    coverage_data <- coverage_data %>%
      filter(chromosome %in% chromosomes)
  }
  
  # Bin data for better visualization
  binned_data <- coverage_data %>%
    mutate(bin = floor(position / bin_size) * bin_size) %>%
    group_by(chromosome, bin) %>%
    summarize(mean_coverage = mean(coverage, na.rm = TRUE),
              .groups = "drop")
  
  # Create plot
  p <- ggplot(binned_data, aes(x = bin, y = mean_coverage)) +
    geom_line() +
    facet_wrap(~chromosome, scales = "free_x") +
    theme_minimal() +
    labs(x = "Genomic Position", 
         y = "Mean Coverage", 
         title = "Genome Coverage Profile")
  
  return(p)
}

#' Create Synteny Plot
#'
#' @param synteny_data Data frame with synteny blocks information
#' @param genome1_name Name of the first genome
#' @param genome2_name Name of the second genome
#' @return A ggplot object
#' @examples
#' # synteny_plot <- plot_synteny(synteny_data, "Oscheius_sp1", "Oscheius_sp2")
plot_synteny <- function(synteny_data, genome1_name, genome2_name) {
  require(ggplot2)
  
  p <- ggplot() +
    geom_segment(data = synteny_data,
                aes(x = genome1_start, xend = genome1_end,
                    y = genome1_chr, yend = genome1_chr,
                    color = factor(block_id)),
                size = 3) +
    geom_segment(data = synteny_data,
                aes(x = genome2_start, xend = genome2_end,
                    y = genome2_chr, yend = genome2_chr,
                    color = factor(block_id)),
                size = 3) +
    geom_segment(data = synteny_data,
                aes(x = genome1_start + (genome1_end - genome1_start)/2,
                    xend = genome2_start + (genome2_end - genome2_start)/2,
                    y = genome1_chr, 
                    yend = genome2_chr,
                    color = factor(block_id)),
                alpha = 0.3) +
    facet_grid(rows = vars(genome_name), scales = "free", space = "free") +
    theme_minimal() +
    theme(legend.position = "none") +
    labs(x = "Genomic Position",
         y = "Chromosome",
         title = paste("Synteny between", genome1_name, "and", genome2_name))
  
  return(p)
}
