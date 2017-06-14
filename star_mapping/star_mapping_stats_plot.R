suppressMessages(library(ggplot2))
suppressMessages(library(reshape2))
suppressMessages(library(scales))
suppressMessages(library(argparser))
options(stringsAsFactors = F)
source('/public/scripts/RNAseq/R/quantification/quant_plot.R')
# source('C:\\work\\scripts\\atom\\R\\quantification\\quant_plot.R')
#
# setwd('C:\\work\\project\\mRNA\\2017\\OM-mRNA-20-Wheat-P20170502\\mapping')

p <- arg_parser("star mapping stats plot")
p <- add_argument(p, "--mapping_stats", help = "star mapping plot data")
p <- add_argument(p, "--sample_inf", help = "group sample file")
p <- add_argument(p, '--out_dir', help = 'output directory')
argv <- parse_args(p)

plot_data <- read.delim(argv$mapping_stats)
sample_df <- read.delim(argv$sample_inf, header = F)
plot_data$ummapped_reads = plot_data$total_reads - plot_data$unique_mapped_reads - plot_data$multiple_mapped_reads
plot_data2 <- plot_data[, c(1,3:5)]
melt_plot_data2 <- melt(plot_data2,id = c('Sample'))
melt_plot_data2$variable <- factor(melt_plot_data2$variable, levels = c('ummapped_reads', 'multiple_mapped_reads', 'unique_mapped_reads'))
melt_plot_data2$Sample <- factor(melt_plot_data2$Sample, levels = sample_df$V2)

theme_set(theme_onmath() +
            theme(panel.grid.minor = element_blank(),
                  panel.grid.major = element_blank(),
                  axis.text.x = element_text(color = "black",face = "bold",angle = 90,hjust = 1,vjust = 0.5, size = rel(0.9)))
          )

p <- ggplot(melt_plot_data2, aes(x =Sample, y = value, fill = variable)) +
  geom_bar(colour = 'black', position = "fill", stat = 'identity', width = 0.75) +
  scale_y_continuous(labels = percent_format(), expand = c(0, 0), breaks = seq(0,1,0.1)) +
  labs(x = NULL, y = NULL) +
  scale_fill_brewer('', palette = 'RdYlBu')

sample_number = dim(plot_data)[1]
plot_width = 6 + sample_number/10
plot_height = 6 + sample_number/20
ggsave(paste(argv$out_dir, 'mapping_stats_plot.png', sep = '/'), plot = p, type = 'cairo-png', width = plot_width, height = plot_height, dpi = 300)
ggsave(paste(argv$out_dir, 'mapping_stats_plot.pdf', sep = '/'), plot = p, width = plot_width, height = plot_height)
