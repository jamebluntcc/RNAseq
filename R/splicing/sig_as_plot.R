suppressMessages(library(ggplot2))
suppressMessages(library(reshape2))
suppressMessages(library(dplyr))
suppressMessages(library(RColorBrewer))
suppressMessages(library(argparser))
suppressMessages(library(ggthemr))
suppressMessages(library(scales))
suppressMessages(library(argparser))
source('/public/scripts/RNAseq/R/quantification/quant_plot.R')
options(stringsAsFactors = F)

p <- arg_parser("as summary plot")
p <- add_argument(p, "--as_summary", help = "as summary plot data")
p <- add_argument(p, '--out_dir', help = 'output directory')
argv <- parse_args(p)

# for test
# source('C:\\work\\scripts\\atom\\R\\quantification\\quant_plot.R')
# out_dir = './'
# as_summary_file <- './B_vs_L.plot.data'

as_summary_file <- argv$as_summary
out_dir <- argv$out_dir

as_summary_name <- basename(as_summary_file)
compare_name <- unlist(strsplit(as_summary_name,split = "\\."))[1]
group_vector <- unlist(strsplit(compare_name,split = "_vs_"))


plot_data <- read.delim(as_summary_file, header = F)
plot_data_1 <- plot_data[,1:3]
plot_data_2 <- plot_data[,c(1,4:5)]
colnames(plot_data_1) <- c('EventType', group_vector[1], group_vector[2])
colnames(plot_data_2) <- c('EventType', group_vector[1], group_vector[2])
plot_data_1$reads_type <- 'Junction_Counts'
plot_data_2$reads_type <- 'ALL_Counts'
plot_data_merge <- rbind(plot_data_1, plot_data_2)
plot_data_merge_reshape <- melt(plot_data_merge, id = c('reads_type', 'EventType'))

ggthemr('fresh')
pal <- swatch()[c(2,4)]
ggthemr_reset()

enrich_theme <- theme_bw()+theme(
  legend.key = element_blank(),
  axis.text.x = element_text(color = "black",face = "bold", size = rel(1.2)),
  axis.text.y = element_text(color = "black",face = "bold", size = rel(1.2)),
  axis.title.y = element_text(color = "black",face = "bold",size = rel(1.2)),
  axis.title.x = element_text(color = "black",face = "bold",size = rel(1.2)),
  panel.grid.minor.x  = element_line(colour = "black"),
  panel.grid.major = element_blank(),
  strip.background = element_blank(),
  strip.text.x = element_text(color = "black",face = "bold",size = rel(1.2)),
  legend.text = element_text(size = rel(0.8))
)
theme_set(enrich_theme)

as_plot <- ggplot(plot_data_merge_reshape, aes(x = EventType, y = value, fill = variable)) + 
  geom_bar(stat = 'identity', position = position_dodge(), color = 'black') +
  facet_grid(. ~ reads_type) +
  scale_fill_manual(values = pal, guide = guide_legend(title = "Group")) +
  ylab('Significant Event Number') +
  xlab('Event Type')

ggsave(paste(out_dir, 'as_summary_plot.png', sep = '/'), plot = as_plot, type = 'cairo-png', width = 8, height = 8, dpi = 300)
ggsave(paste(out_dir, 'as_summary_plot.pdf', sep = '/'), plot = as_plot, width = 8, height = 8)
