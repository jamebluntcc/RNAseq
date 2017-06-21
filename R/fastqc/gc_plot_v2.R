#2017-03-01
options(warn = -1)
options(stringsAsFactors = F)
suppressMessages(require('tidyverse',quietly = T))
suppressMessages(require('reshape2',quietly = T))
suppressMessages(require('scales',quietly = T))
suppressMessages(require('argparser',quietly = T))
source('/public/scripts/RNAseq/R/rseqc/RNAseq_plot_lib.R')

#----plot----
p <- arg_parser('gc plot')
p <- add_argument(p,'--gc_dir',help = 'gc stats')
p <- add_argument(p, '--sample_inf', help = 'sample info file')
p <- add_argument(p,'--out_dir',help = 'path to save plot')
argv <- parse_args(parser = p)

# for test
# source('../../scripts/atom/R/rseqc/RNAseq_plot_lib.R')
# setwd('C:\\work\\pipe\\fastqc')
# gc_dir <- './gc_plot/'
# sample_inf <- './group_sample'

# read parameters
gc_dir <- argv$gc_dir
sample_inf <- argv$sample_inf
out_dir <- argv$out_dir


group_inf_df <- read.delim(sample_inf, header = F)
sample_number <- length(group_inf_df$V2)

gc_files <- file.path(gc_dir, paste(group_inf_df$V2, "gc.txt", sep = '.'))
#file.exists(gc_files)
gc_file_list <- list()
for (i in seq(sample_number)) {
  each_sample_gc_df <- read.delim(gc_files[i])
  each_sample_gc_df[,2:dim(each_sample_gc_df)[2]] <- each_sample_gc_df[,2:dim(each_sample_gc_df)[2]] / 100
  each_sample_gc_df$sample <- group_inf_df$V2[i]
  each_sample_gc_df[is.na(each_sample_gc_df)] <- 0
  gc_file_list[[i]] <- each_sample_gc_df
  each_sample_out_name <- paste(group_inf_df$V2[i], 'gc_distribution.line', sep = '.')
  each_sample_out_path <- file.path(gc_dir, each_sample_out_name)
  rs_each_sample_gc_df <- melt(each_sample_gc_df,id=c('X.Base', 'sample'))
  gc_line_plot(rs_each_sample_gc_df, each_sample_out_path)
}
gc_file_df <- ldply(gc_file_list, data.frame)
rs_gc_file_df <- melt(gc_file_df,id=c('X.Base', 'sample'))

gc_plot_out <- file.path(gc_dir, 'gc_distribution.line')
rs_gc_file_df$sample <- factor(rs_gc_file_df$sample, levels = group_inf_df$V2)
gc_line_plot(rs_gc_file_df, gc_plot_out)
