suppressMessages(library(ggplot2))
suppressMessages(library(reshape2))
suppressMessages(library(dplyr))
suppressMessages(library(RColorBrewer))
suppressMessages(library(argparser))
suppressMessages(library(xlsx))
options(stringsAsFactors = F)
source('/public/scripts/RNAseq/R/quantification/quant_plot.R')


p <- arg_parser("snp stats plot")
p <- add_argument(p, "--snp_stats", help = "star mapping plot data")
p <- add_argument(p, "--sample_inf", help = "group sample file")
p <- add_argument(p, '--out_dir', help = 'output directory')
argv <- parse_args(p)

snp_stats <- argv$snp_stats
sample_inf <- argv$sample_inf
out_dir <- argv$out_dir

## for test

#source('C:\\work\\scripts\\atom\\R\\quantification\\quant_plot.R')
# out_dir <- './'
# snp_stats <- './tstv_by_sample.0.dat'

samples <- read.delim(sample_inf, header = F)
colnames(samples) <- c('condition', 'sample')

snp_summary <- read.delim(snp_stats)
snp_summary <- snp_summary[,c(8,2:7)]
colnames(snp_summary) <- c('Sample_id', 'ts/tv', 'het/hom', 'nSNPs', 'nIndels', 'Avarage_depth', 'nSingletons')
write.table(snp_summary, file = paste(out_dir, 'snp_stats.txt', sep = '/'), quote=F, row.names = F, sep = '\t')
write.xlsx(snp_summary, file = paste(out_dir, 'snp_stats.xlsx', sep = '/'), sheetName = "snp_stats", append = FALSE, row.names = F)

plot_data <- melt(snp_summary)
plot_data <- filter(plot_data, variable == "nSNPs" | variable == "nIndels" | variable == "nSingletons")
plot_data$Sample_id <- factor(plot_data$Sample_id, levels = samples$sample)

pal <- c(brewer.pal(3,'Set1'))

theme_set(theme_onmath()+theme(axis.text.x = element_text(angle = -90,color = 'black',
                                                          vjust = 0.5,hjust = 0,
                                                          size = rel(1.2)),
                               axis.text.y = element_text(size = rel(1.2)),
                               legend.text = element_text(size = rel(0.8)),
                               legend.key = element_blank(),
                               legend.title = element_blank()))

p <- ggplot(plot_data, aes(x = Sample_id, y = value, fill = variable)) +
  geom_bar(stat = 'identity', position = 'dodge', color = 'black') +
  guides(title = '')+
  scale_fill_manual(values = pal) +
  xlab("") + ylab('Number')

sample_number = dim(plot_data)[1]
plot_width = 6 + sample_number/10
plot_height = 6 + sample_number/20

ggsave(paste(out_dir, 'snp_stats_plot.png', sep = '/'), plot = p, type = 'cairo-png', width = plot_width, height = plot_height, dpi = 300)
ggsave(paste(out_dir, 'snp_stats_plot.pdf', sep = '/'), plot = p, width = plot_width, height = plot_height)
