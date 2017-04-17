suppressMessages(library(argparser))
suppressMessages(library(ggplot2))
suppressMessages(library(ggthemes))
suppressMessages(library(RColorBrewer))
suppressMessages(library(scales))
suppressMessages(library(reshape2))
suppressMessages(library(plyr))
suppressMessages(library(dplyr))
options(stringsAsFactors = F)
source('/home/public/scripts/RNAseq/R/quantification/quant_plot.R')

p <- arg_parser("rseqc plot")
p <- add_argument(p, "--sample_inf", help = "sample information")
p <- add_argument(p, '--read_distribution_dir', help = 'read distribution directory')
p <- add_argument(p, '--genebody_cov_dir', help = 'genebody coverage directory')
p <- add_argument(p, '--inner_distance_dir', help = 'inner distance directory')
p <- add_argument(p, '--reads_duplication_dir', help = 'reads duplication directory')
argv <- parse_args(p)

## for test
#source('../../atom/R/quantification/quant_plot.R')
# setwd('C:\\work\\scripts\\test\\rseqc_test')
# read_distribution_df <- read.delim('read_distribution/read_distribution.summary.txt', header = T)
# group_inf_df <- read.delim('group_sample', header = F)
# out_dir <- './'


## read arguments
sample_inf_file <- argv$sample_inf
read_distribution_dir <- argv$read_distribution_dir
group_inf_df <- read.delim(sample_inf_file, header = F)
sample_number <- length(group_inf_df$V2)

## read distribution plot
read_distribution_file <- file.path(read_distribution_dir, 'read_distribution.summary.txt')
read_distribution_df <- read.delim(read_distribution_file, header = T)
read_distribution_df$Sample <- factor(read_distribution_df$Sample, levels = group_inf_df$V2)
read_distribution_filter_df <- dplyr::filter(read_distribution_df, Total_bases > 0)

### plot
theme_set(theme_onmath() + theme(panel.grid.minor = element_blank(), panel.grid.major = element_blank()))
onmath_color <- colorRampPalette(onmath_color_basel)(length(unique(read_distribution_df$Group)))
p <- ggplot(read_distribution_filter_df, aes(Sample, Tag_count, fill = Group)) + 
  geom_bar(colour = 'black', position = "fill", stat = 'identity') + 
  scale_fill_manual('Genomic features', values = onmath_color) +
  scale_y_continuous(labels = percent_format(), expand = c(0, 0)) + 
  labs(x = NULL, y = NULL) 

### output
sample_number <- length(unique(read_distribution_df$Sample))
plot_width <- 2 + sample_number

ggsave(filename = paste(read_distribution_dir,'read_distribution.png', sep = '/'), plot = p, type="cairo-png", width = plot_width, height = 6, dpi = 300)
ggsave(filename = paste(read_distribution_dir,'read_distribution.pdf', sep = '/'), plot = p, width = plot_width, height = 6)

## genebody coverage

### read files
#gene_cov_files <- file.path('./genebody_coverage/', paste(group_inf_df$V2, "geneBodyCoverage.txt", sep = '.'))
#file.exists(gene_cov_files)

genebody_cov_dir <- argv$genebody_cov_dir
gene_cov_files <- file.path(genebody_cov_dir, paste(group_inf_df$V2, "geneBodyCoverage.txt", sep = '.'))
gene_cov_list = list()
for (i in seq(length(gene_cov_files))){
  each_file_df <- read.delim(gene_cov_files[i])
  gene_cov_list[[i]] = each_file_df
}
gene_cov_df <- ldply(gene_cov_list, data.frame)
gene_cov_num_df <- gene_cov_df[,2:dim(gene_cov_df)[2]]
gene_cov_per_df <- gene_cov_num_df/rowSums(gene_cov_num_df)
colnames(gene_cov_per_df) <- seq(1,100,1)
gene_cov_per_df$samples <- gene_cov_df$Percentile
gene_cov_mdf <- melt(gene_cov_per_df, id = 'samples')

### plot
theme_set(theme_onmath())
onmath_color <- colorRampPalette(onmath_color_basel)(length(group_inf_df$V2))

gene_body_cov_plot <- ggplot(gene_cov_mdf, aes(variable, value, colour = samples, group = samples)) + 
  geom_point(size = 1) + geom_line() +
  scale_y_continuous(labels = percent_format()) +
  scale_x_discrete(breaks = seq(0,100,5)) +
  scale_color_manual(values = onmath_color) +
  labs(x = "Gene body percentile (5'->3')", y = 'Coverage')

ggsave(filename = paste(genebody_cov_dir,'genebody_coverage.png', sep = '/'), plot = gene_body_cov_plot, type="cairo-png", width = 8, height = 6, dpi = 300)
ggsave(filename = paste(genebody_cov_dir,'genebody_coverage.pdf', sep = '/'), plot = gene_body_cov_plot, width = 8, height = 6)

## inner distance
# inner_distance_files <- file.path('./inner_distance/', paste(group_inf_df$V2, "inner_distance_freq.txt", sep = '.'))
inner_distance_dir <- argv$inner_distance_dir
inner_distance_files <- file.path(inner_distance_dir, paste(group_inf_df$V2, "inner_distance_freq.txt", sep = '.'))
inner_distance_list = list()
for (i in seq(length(inner_distance_files))){
  each_file_df <- read.delim(inner_distance_files[i], header = F)
  each_file_df$sample = group_inf_df$V2[i]
  each_file_df$proportion = each_file_df$V3/sum(each_file_df$V3)
  inner_distance_list[[i]] = each_file_df
}
inner_distance_df <- ldply(inner_distance_list, data.frame)
inner_distance_df$sample <- factor(inner_distance_df$sample, levels = group_inf_df$V2)

theme_set(theme_onmath_border() + 
            theme(panel.grid.minor = element_blank(), 
                  panel.grid.major = element_blank(), 
                  strip.background = element_blank()
                  )
)

inner_distance_plot <- ggplot(inner_distance_df, aes(x = V2, y = proportion, fill = sample )) + 
  geom_bar(stat = 'identity', width = 5, colour = 'black') + 
  scale_x_continuous(breaks = seq(-200, 200, 50)) +
  scale_fill_manual(values = onmath_color) +
  scale_y_continuous(labels = percent_format(), expand = c(0, 0)) +
  labs(x = 'Inner distance (bp)', y = 'Percent of reads') +
  facet_wrap(~sample , ncol = 6)

inner_distance_plot_height <- ceiling(sample_number/6)*4
if (sample_number < 6) {
  inner_distance_plot_width <- sample_number*4 + 1
} else {
  inner_distance_plot_width = 26
}

ggsave(filename = paste(inner_distance_dir,'inner_distance.png', sep = '/'), plot = inner_distance_plot, type="cairo-png", width = inner_distance_plot_width, height = inner_distance_plot_height, dpi = 300)
ggsave(filename = paste(inner_distance_dir,'inner_distance.pdf', sep = '/'), plot = inner_distance_plot, width = inner_distance_plot_width, height = inner_distance_plot_height)


## reads duplication
# reads_duplication_dir <- './read_duplication/'
reads_duplication_dir <- argv$reads_duplication_dir
reads_duplication_files <- file.path(reads_duplication_dir, paste(group_inf_df$V2, "pos.DupRate.xls", sep = '.'))
reads_duplication_list = list()
for (i in seq(length(reads_duplication_files))){
  each_file_df <- read.delim(reads_duplication_files[i])
  each_file_df$sample = group_inf_df$V2[i]
  each_file_df$proportion = each_file_df$UniqReadNumber/sum(each_file_df$UniqReadNumber)
  each_file_filter_df <- dplyr::filter(each_file_df, Occurrence <= 10)
  add_row <- c('>10', 0, as.character(group_inf_df$V2[i]), 0)
  add_row[2] <- sum(each_file_df$UniqReadNumber)-sum(each_file_filter_df$UniqReadNumber)
  add_row[4] <- (1-sum(each_file_filter_df$proportion))
  each_file_plot_df <- rbind(each_file_filter_df, add_row)
  each_file_plot_df$proportion <- as.numeric(each_file_plot_df$proportion)
  each_file_plot_df$Occurrence <- factor(each_file_plot_df$Occurrence, levels = each_file_plot_df$Occurrence)
  reads_duplication_list[[i]] = each_file_plot_df
}

reads_duplication_df <- ldply(reads_duplication_list, data.frame)
reads_duplication_df$sample <- factor(reads_duplication_df$sample, levels = group_inf_df$V2)
max_proportion <- round(max(reads_duplication_df$proportion) + 0.1, 1)
onmath_color <- colorRampPalette(onmath_color_basel)(length(group_inf_df$V2))

theme_set(theme_onmath())
onmath_color <- colorRampPalette(onmath_color_basel)(length(group_inf_df$V2))
reads_duplication_plot <- ggplot(reads_duplication_df, aes(Occurrence, proportion, group = sample, colour = sample)) + 
  geom_point() + geom_line() +
  scale_y_continuous(breaks = seq(0,max_proportion,0.1), limits = c(0,max_proportion), 
                     labels = percent_format(), expand = c(0, 0)) + 
  scale_color_manual(values = onmath_color) + 
  labs(x = 'Occurrence of reads', y = 'Percent of reads') 

ggsave(filename = paste(reads_duplication_dir,'reads_duplication.png', sep = '/'), plot = reads_duplication_plot, type="cairo-png", width = 6, height = 6, dpi = 300)
ggsave(filename = paste(reads_duplication_dir,'reads_duplication.pdf', sep = '/'), plot = reads_duplication_plot, width = 6, height = 6)
