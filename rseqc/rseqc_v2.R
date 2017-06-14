suppressMessages(library(argparser))
suppressMessages(library(ggplot2))
suppressMessages(library(ggthemes))
suppressMessages(library(RColorBrewer))
suppressMessages(library(scales))
suppressMessages(library(reshape2))
suppressMessages(library(plyr))
suppressMessages(library(dplyr))
suppressMessages(library(ggthemr))
options(stringsAsFactors = F)
source('/public/scripts/RNAseq/R/rseqc/RNAseq_plot_lib.R')

p <- arg_parser("rseqc plot")
p <- add_argument(p, "--sample_inf", help = "sample information")
p <- add_argument(p, '--read_distribution_dir', help = 'read distribution directory')
p <- add_argument(p, '--genebody_cov_dir', help = 'genebody coverage directory')
p <- add_argument(p, '--inner_distance_dir', help = 'inner distance directory')
p <- add_argument(p, '--reads_duplication_dir', help = 'reads duplication directory')
argv <- parse_args(p)

## for test
# source('C:\\work\\scripts\\atom\\R\\quantification\\quant_plot.R')
# setwd('C:\\work\\pipe\\rseqc')
# source('./RNAseq_plot_lib.R')
#read_distribution_df <- read.delim('read_distribution/read_distribution.summary.txt', header = T)
#group_inf_df <- read.delim('group_sample', header = F)
#out_dir <- './'
# sample_inf_file <- './group_sample'
# read_distribution_dir <- './rseqc/read_distribution/'
# genebody_cov_dir <- './rseqc/genebody_coverage/'
# inner_distance_dir <- './rseqc/inner_distance/'
# reads_duplication_dir <- './rseqc/read_duplication/'

## read arguments
sample_inf_file <- argv$sample_inf
read_distribution_dir <- argv$read_distribution_dir
genebody_cov_dir <- argv$genebody_cov_dir
inner_distance_dir <- argv$inner_distance_dir
reads_duplication_dir <- argv$reads_duplication_dir

group_inf_df <- read.delim(sample_inf_file, header = F)
sample_number <- length(group_inf_df$V2)

## read distribution plot
read_distribution_file <- file.path(read_distribution_dir, 'read_distribution.summary.txt')
read_distribution_df <- read.delim(read_distribution_file, header = T)
read_distribution_df$Sample <- factor(read_distribution_df$Sample, levels = group_inf_df$V2)
read_distribution_filter_df <- dplyr::filter(read_distribution_df, Total_bases > 0)

read_distribution_filter_df$Group <- factor(read_distribution_filter_df$Group, levels = rev(unique(read_distribution_filter_df$Group)))

for (each_sample in group_inf_df$V2) {
  each_sample_pie_plot_data <- dplyr::filter(read_distribution_filter_df, Sample == each_sample)
  each_sample_pie_plot_data$portion <- each_sample_pie_plot_data$Tag_count/sum(each_sample_pie_plot_data$Tag_count)
  each_sample_pie_plot_data$label<-paste(each_sample_pie_plot_data$Group," (",percent(each_sample_pie_plot_data$portion),")",sep="")
  each_sample_pie_plot_data$label <- factor(each_sample_pie_plot_data$label, levels = rev(each_sample_pie_plot_data$label))
  each_sample_out_name <- paste(each_sample, 'read_distribution.pie', sep = '.')
  each_sample_out_prefix <- file.path(read_distribution_dir, each_sample_out_name)
  read_distribution_pie_plot(each_sample_pie_plot_data, each_sample_out_prefix)
}

read_distribution_bar_plot_out_prefix <- file.path(read_distribution_dir, 'read_distribution.bar')
read_distribution_bar_plot(read_distribution_filter_df, read_distribution_bar_plot_out_prefix)

## genebody coverage

### read files
#gene_cov_files <- file.path('./genebody_coverage/', paste(group_inf_df$V2, "geneBodyCoverage.txt", sep = '.'))
#file.exists(gene_cov_files)

gene_cov_files <- file.path(genebody_cov_dir, paste(group_inf_df$V2, "geneBodyCoverage.txt", sep = '.'))
gene_cov_list = list()
for (i in seq(length(gene_cov_files))){
  each_file_df <- read.delim(gene_cov_files[i])
  # each_file_df[1,1] <- gsub('^V([0-9])','\\1',each_file_df[1,1])
  each_file_df[1,1] <- group_inf_df$V2[i]
  gene_cov_list[[i]] = each_file_df
}
gene_cov_df <- ldply(gene_cov_list, data.frame)
gene_cov_num_df <- gene_cov_df[,2:dim(gene_cov_df)[2]]
gene_cov_per_df <- gene_cov_num_df/rowSums(gene_cov_num_df)
colnames(gene_cov_per_df) <- seq(1,100,1)
gene_cov_per_df$sample <- gene_cov_df$Percentile
gene_cov_mdf <- melt(gene_cov_per_df, id = 'sample')

for (each_sample in unique(gene_cov_mdf$sample)) {
  each_sample_cov_plot_data <- dplyr::filter(gene_cov_mdf, sample == each_sample)
  each_sample_out_name <- paste(each_sample, 'genebody_coverage.point', sep = '.')
  each_sample_out_prefix <- file.path(genebody_cov_dir, each_sample_out_name)
  gene_body_cov_plot(each_sample_cov_plot_data, each_sample_out_prefix)
}

gene_body_merged_plot_prefix <- file.path(genebody_cov_dir, 'genebody_coverage.point')
gene_body_cov_plot(gene_cov_mdf, gene_body_merged_plot_prefix)

## inner distance

inner_distance_files <- file.path(inner_distance_dir, paste(group_inf_df$V2, "inner_distance_freq.txt", sep = '.'))
inner_distance_list = list()
for (i in seq(length(inner_distance_files))){
  each_file_df <- read.delim(inner_distance_files[i], header = F)
  each_file_df$sample = group_inf_df$V2[i]
  each_file_df$proportion = each_file_df$V3/sum(each_file_df$V3)
  inner_distance_list[[i]] = each_file_df
}
inner_distance_df <- ldply(inner_distance_list, data.frame)

for (each_sample in unique(inner_distance_df$sample)) {
  each_sample_inner_distance_plot_data <- dplyr::filter(inner_distance_df, sample == each_sample)
  each_sample_out_name <- paste(each_sample, 'inner_distance.bar', sep = '.')
  each_sample_out_prefix <- file.path(inner_distance_dir, each_sample_out_name)
  inner_distance_plot(each_sample_inner_distance_plot_data, each_sample_out_prefix)
}

inner_distance_df$sample <- factor(inner_distance_df$sample, levels = group_inf_df$V2)
inner_distance_out_prefix <- file.path(inner_distance_dir, 'inner_distance.bar')
inner_distance_plot(inner_distance_df, inner_distance_out_prefix)

## reads duplication
# reads_duplication_dir <- './read_duplication/'

reads_duplication_pos_files <- file.path(reads_duplication_dir, paste(group_inf_df$V2, "pos.DupRate.xls", sep = '.'))
reads_duplication_seq_files <- file.path(reads_duplication_dir, paste(group_inf_df$V2, "seq.DupRate.xls", sep = '.'))
reads_duplication_pos_list = list()
reads_duplication_seq_list = list()

get_dup_df <- function(dup_files, df_list) {
  for (i in seq(length(dup_files))){
    each_file_df <- read.delim(dup_files[i])
    each_file_df$sample = group_inf_df$V2[i]
    each_file_df$proportion = each_file_df$UniqReadNumber/sum(each_file_df$UniqReadNumber)
    each_file_filter_df <- dplyr::filter(each_file_df, Occurrence <= 10)
    add_row <- c('>10', 0, as.character(group_inf_df$V2[i]), 0)
    add_row[2] <- sum(each_file_df$UniqReadNumber)-sum(each_file_filter_df$UniqReadNumber)
    add_row[4] <- (1-sum(each_file_filter_df$proportion))
    each_file_plot_df <- rbind(each_file_filter_df, add_row)
    each_file_plot_df$proportion <- as.numeric(each_file_plot_df$proportion)
    each_file_plot_df$Occurrence <- factor(each_file_plot_df$Occurrence, levels = each_file_plot_df$Occurrence)
    df_list[[i]] = each_file_plot_df
  }
  df_list
}

reads_duplication_pos_list <- get_dup_df(reads_duplication_pos_files, reads_duplication_pos_list)
reads_duplication_pos_df <- ldply(reads_duplication_pos_list, data.frame)
reads_duplication_seq_list <- get_dup_df(reads_duplication_seq_files, reads_duplication_seq_list)
reads_duplication_seq_df <- ldply(reads_duplication_seq_list, data.frame)
reads_duplication_pos_df$method <- 'position'
reads_duplication_seq_df$method <- 'sequence'
reads_duplication_df <- rbind(reads_duplication_pos_df, reads_duplication_seq_df)

for (each_sample in unique(reads_duplication_df$sample)) {
  each_sample_duplication_plot_data <- dplyr::filter(reads_duplication_df, sample == each_sample)
  each_sample_out_name <- paste(each_sample, 'reads_duplication.point', sep = '.')
  each_sample_out_prefix <- file.path(reads_duplication_dir, each_sample_out_name)
  reads_duplication_plot(each_sample_duplication_plot_data, each_sample_out_prefix)
}

reads_duplication_df$sample <- factor(reads_duplication_df$sample, levels = group_inf_df$V2)
reads_duplication_out_prefix <- file.path(reads_duplication_dir, 'reads_duplication.point')
reads_duplication_plot(reads_duplication_df, reads_duplication_out_prefix)
