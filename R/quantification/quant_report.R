suppressMessages(library(argparser))
suppressMessages(library(plyr))
suppressMessages(library(dplyr))
source('/public/scripts/RNAseq/R/quantification/quant_plot.R')

p <- arg_parser("for expression analysis report plot")
p <- add_argument(p, '--quant_dir',  help = 'quantification directory')
p <- add_argument(p, '--sample_inf', help = 'sample information with sample names and group names')
p <- add_argument(p, '--qvalue',     help = 'diff gene qvalue cutoff', default = 0.05)
p <- add_argument(p, '--logfc',      help = 'diff gene logfc cutoff',  default = 1)
argv <- parse_args(p)

# for test
# setwd('C:\\work\\pipe\\quantification')
# quant_dir <- './'
# sample_inf <- '../fastqc/group_sample'
# source('C:\\work\\scripts\\atom\\R\\quantification\\quant_plot.R')

# read parameters
quant_dir <- argv$quant_dir
sample_inf <- argv$sample_inf
qvalue <- argv$qvalue
logfc <- argv$logfc

samples <- read.delim(sample_inf, stringsAsFactors = F, header = F)
colnames(samples) <- c('condition', 'sample')
all_combine = combn(unique(samples$condition), 2)
diff_dir <- file.path(quant_dir, 'differential_analysis')
exp_dir <- file.path(quant_dir, 'expression_summary/')


diff_df_list = list()
for (i in seq(dim(all_combine)[2])) {

  each_pair <- all_combine[,i]
  each_compare <- paste(each_pair[1], '_vs_', each_pair[2], sep = '')
  each_compare_name = paste(each_compare, 'edgeR.DE_results.txt', sep = '.')
  each_compare_file = file.path(diff_dir, each_compare, each_compare_name)
  each_compare_df <- read.delim(each_compare_file, check.names = F)
  each_compare_plot <- each_compare_df[, c('Gene_ID', 'logFC', 'FDR')]
  each_compare_plot$compare <- each_compare
  diff_df_list[[i]] <- each_compare_plot

}

diff_df <- ldply(diff_df_list, data.frame)
om_volcano_plot(diff_df, 'ALL', logfc, qvalue, exp_dir)
pdf_example_diff_out <- file.path(exp_dir, 'pdf.example.diff.table.txt')
html_example_diff_out <- file.path(exp_dir, 'html.example.diff.table.txt')
pdf_example_diff_df <- each_compare_df[1:100, c('Gene_ID', 'logFC', 'PValue', 'FDR')]

write.table(pdf_example_diff_df, file=pdf_example_diff_out, sep='\t', quote=F, row.names=F)
write.table(each_compare_df, file=html_example_diff_out, sep='\t', quote=F, row.names=F)

diff_exp_df <- filter(diff_df, abs(logFC) >= logfc, FDR <= qvalue)
diff_gene_id <- unique(diff_exp_df$Gene_ID)
gene_count_matrix <- file.path(exp_dir, 'Gene.count.txt')
gene_tpm_matrix <- file.path(exp_dir, 'Gene.tpm.txt')
gene_count_matrix_df <- read.delim(gene_count_matrix, check.names = F)
gene_tpm_matrix_df <- read.delim(gene_tpm_matrix, check.names = F)
out_diff_gene_count_matrix <- filter(gene_count_matrix_df, Gene_ID %in% diff_gene_id)
out_diff_gene_tpm_matrix <- filter(gene_tpm_matrix_df, Gene_ID %in% diff_gene_id)

if (dim(out_diff_gene_count_matrix)[1] > 0) {
  write.table(out_diff_gene_count_matrix, file = paste(exp_dir, 'Diff.gene.count.txt', sep = '/'), quote=F, row.names = F, sep = '\t')
  #write.xlsx(out_diff_gene_count_matrix, file = paste(expression_stat_dir, 'Diff.gene.count.xlsx', sep = '/'), sheetName = "gene.count", append = FALSE, row.names = F)
  write.table(out_diff_gene_tpm_matrix, file = paste(exp_dir, 'Diff.gene.tpm.txt', sep = '/'), quote = F, row.names = F, sep = '\t')
  #write.xlsx(out_diff_gene_tpm_matrix, file = paste(expression_stat_dir, 'Diff.gene.tpm.xlsx', sep = '/'), sheetName = "gene.tpm", append = FALSE, row.names = F)
  ## heatmap
  diff_gene_tpm_matrix <- out_diff_gene_tpm_matrix[,-1]
  rownames(diff_gene_tpm_matrix) <- out_diff_gene_tpm_matrix[,1]
  om_heatmap(plot_data = diff_gene_tpm_matrix, samples = samples, outdir = exp_dir)
}
