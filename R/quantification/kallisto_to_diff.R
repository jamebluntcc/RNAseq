suppressMessages(library(argparser))
suppressMessages(library(tximport))
suppressMessages(library(readr))
suppressMessages(library(edgeR))
suppressMessages(library(dplyr))
suppressMessages(library(RColorBrewer))
suppressMessages(library(gplots))
suppressMessages(library(xlsx))
suppressMessages(library(tibble))
suppressMessages(library(rhdf5))
source('/public/scripts/RNAseq/R/quantification/quant_plot.R')


p <- arg_parser("read kallisto quant files and perform differential analysis")
p <- add_argument(p, '--quant_dir',  help = 'kallisto quantification directory')
p <- add_argument(p, '--sample_inf', help = 'sample information with sample names and group names')
p <- add_argument(p, '--gene2tr',    help = 'gene id and transcript id mapping file')
p <- add_argument(p, '--out_dir',    help = 'diff analyssi output directory')
p <- add_argument(p, '--qvalue',     help = 'diff gene qvalue cutoff', default = 0.05)
p <- add_argument(p, '--logfc',      help = 'diff gene logfc cutoff',  default = 1)
argv <- parse_args(p)

# ## for test
# sample_inf <- 'group_sample2'
# quant_dir <- './kallisto2/'
# gene2tr_file <- 'gene_trans.map'
# outdir <- './'

## read parameters
sample_inf <- argv$sample_inf
quant_dir <- argv$quant_dir
gene2tr_file <- argv$gene2tr
outdir <- argv$out_dir
qvalue <- argv$qvalue
logfc <- argv$logfc

## directory prepare
expression_stat_dir <- file.path(outdir, 'expression_summary')
diff_dir <- file.path(outdir, 'differential_analysis')
dir.create(expression_stat_dir, showWarnings = FALSE)
dir.create(diff_dir, showWarnings = FALSE)

samples <- read.delim(sample_inf, stringsAsFactors = F, header = F)
colnames(samples) <- c('condition', 'sample')
files <- file.path(quant_dir, samples$sample, "abundance.tsv")
names(files) <- samples$sample
gene2tr <- read.delim(gene2tr_file, header = FALSE)
colnames(gene2tr) <- c('gene_id', 'transcript_id')
tx2gene <- gene2tr[,c('transcript_id', 'gene_id')]


## normalized expression matrix
txi <- tximport(files, type = "kallisto", tx2gene = tx2gene)
cts <- txi$counts
y <- DGEList(cts)
normfactors <- calcNormFactors(y)
gene_tpm_matrix <- (txi$abundance)/(normfactors$samples$norm.factors)

## output quant table
out_gene_tpm_matrix <- as.data.frame(gene_tpm_matrix)
out_gene_tpm_matrix <- round(out_gene_tpm_matrix, 3)
out_gene_tpm_matrix <- rownames_to_column(out_gene_tpm_matrix, var="Gene_ID")

out_cts <- as.data.frame(cts)
out_cts <- round(out_cts, 3)
out_cts <- rownames_to_column(out_cts, var = 'Gene_ID')
write.table(out_cts, file = paste(expression_stat_dir, 'Gene.count.txt', sep = '/'), quote=F, row.names = F, sep = '\t')
#write.xlsx(out_cts, file = paste(expression_stat_dir, 'Gene.count.xlsx', sep = '/'), sheetName = "gene.count", append = FALSE, row.names = F)
write.table(out_gene_tpm_matrix, file = paste(expression_stat_dir, 'Gene.tpm.txt', sep = '/'), quote = F, row.names = F, sep = '\t')
#write.xlsx(out_gene_tpm_matrix, file = paste(expression_stat_dir, 'Gene.tpm.xlsx', sep = '/'), sheetName = "gene.tpm", append = FALSE, row.names = F)

## TODO: speed with apply
## TODO: parallel
diff_genes <- c()
all_combine = combn(unique(samples$condition), 2)
for (i in seq(dim(all_combine)[2])) {
  ## get compare samples
  each_pair <- all_combine[,i]
  con1_sample <- samples[samples$condition == each_pair[2],'sample']
  con2_sample <- samples[samples$condition == each_pair[1],'sample']
  each_compare_name = paste(each_pair[2], 'vs', each_pair[1], sep = '_')
  ## get sample data
  each_pair_samples <- samples[samples$sample %in% c(con1_sample, con2_sample),]
  each_pair_files <- file.path(quant_dir, each_pair_samples$sample, "abundance.tsv")
  names(each_pair_files) <- each_pair_samples$sample
  each_pair_txi <- tximport(each_pair_files, type = "kallisto", tx2gene = tx2gene)
  each_pair_cts <- each_pair_txi$counts
  each_pair_cts <- each_pair_cts[rowSums(each_pair_cts)>=2,]
  each_pair_genes <- row.names(each_pair_cts)
  ## normalization
  normMat <- each_pair_txi$length
  normMat <- normMat[each_pair_genes,]
  normMat <- normMat/exp(rowMeans(log(normMat)))
  o <- log(calcNormFactors(each_pair_cts/normMat)) + log(colSums(each_pair_cts/normMat))
  conditions = factor(c(rep(each_pair[1], length(con2_sample)), rep(each_pair[2], length(con1_sample))))
  y <- DGEList(each_pair_cts, group=conditions)
  y$offset <- t(t(log(normMat)) + o)
  y <- estimateDisp(y)
  ## diff
  et <- exactTest(y, pair=each_pair)
  tTags <- topTags(et,n=NULL)
  new_tTags <- tTags$table
  new_tTags <- new_tTags[, !(names(new_tTags) %in% c("logCPM"))]
  each_pair_matrix <- gene_tpm_matrix[,c(con1_sample, con2_sample)]
  merged_df <- merge(each_pair_matrix, new_tTags, by.x = 0, by.y = 0, all.y = T)
  sorted_merged_df <- arrange(merged_df, FDR)
  colnames(sorted_merged_df)[1] <- 'Gene_ID'
  each_pair_outdir <- file.path(diff_dir, each_compare_name)
  dir.create(each_pair_outdir, showWarnings = FALSE)
  out_file_name_prefix <- paste(each_pair_outdir, '/', each_compare_name, sep = '')
  up_regulate_name_prefix <- paste(each_pair_outdir, '/', each_compare_name, '.',each_pair[2], '-UP',  sep = '')
  down_regulate_name_prefix <- paste(each_pair_outdir, '/', each_compare_name, '.',each_pair[1], '-UP', sep = '')
  diff_genes <- c()
  up_regulate_df <- filter(sorted_merged_df, logFC >= logfc, FDR <= qvalue)
  down_regulate_df <- filter(sorted_merged_df, logFC <= -(logfc), FDR <= qvalue)
  diff_genes <- c(diff_genes, up_regulate_df$Gene_ID, down_regulate_df$Gene_ID)
  write.table(sorted_merged_df, file=paste(out_file_name_prefix, 'edgeR.DE_results', 'txt', sep = '.'), sep='\t', quote=F, row.names=F)
  write.table(up_regulate_df, file=paste(up_regulate_name_prefix, 'edgeR.DE_results', 'txt', sep = '.'), sep='\t', quote=F, row.names=F)
  write.table(down_regulate_df, file=paste(down_regulate_name_prefix, 'edgeR.DE_results', 'txt', sep = '.'), sep='\t', quote=F, row.names=F)
#  write.xlsx(sorted_merged_df, file = paste(out_file_name_prefix, 'edgeR.DE_results', 'xlsx', sep = '.'), sheetName = each_compare_name, append = F, row.names = F)
#  write.xlsx(up_regulate_df, file = paste(up_regulate_name_prefix, 'edgeR.DE_results', 'xlsx', sep = '.'), sheetName = each_compare_name, append = F, row.names = F)
#  write.xlsx(down_regulate_df, file = paste(down_regulate_name_prefix, 'edgeR.DE_results', 'xlsx', sep = '.'), sheetName = each_compare_name, append = F, row.names = F)
  ## write diff gene list
  write(as.character(up_regulate_df$Gene_ID), file = paste(up_regulate_name_prefix, 'edgeR.DE_results.diffgenes', 'txt', sep = '.'), sep = '\n')
  write(as.character(down_regulate_df$Gene_ID), file = paste(down_regulate_name_prefix, 'edgeR.DE_results.diffgenes', 'txt', sep = '.'), sep = '\n')
  write(as.character(diff_genes), file = paste(out_file_name_prefix, 'ALL', 'edgeR.DE_results.diffgenes', 'txt', sep = '.'), sep = '\n')
  ## volcano plot
  om_volcano_plot(sorted_merged_df, each_compare_name, logfc, qvalue, each_pair_outdir)
}

## diff gene matrix
diff_genes <- unique(diff_genes)
diff_gene_tpm_matrix <- gene_tpm_matrix[diff_genes,]

## output diff gene quant table
out_diff_gene_count_matrix <- filter(out_cts, Gene_ID %in% diff_genes)
out_diff_gene_tpm_matrix <- filter(out_gene_tpm_matrix, Gene_ID %in% diff_genes)

write.table(out_diff_gene_count_matrix, file = paste(expression_stat_dir, 'Diff.gene.count.txt', sep = '/'), quote=F, row.names = F, sep = '\t')
#write.xlsx(out_diff_gene_count_matrix, file = paste(expression_stat_dir, 'Diff.gene.count.xlsx', sep = '/'), sheetName = "gene.count", append = FALSE, row.names = F)
write.table(out_diff_gene_tpm_matrix, file = paste(expression_stat_dir, 'Diff.gene.tpm.txt', sep = '/'), quote = F, row.names = F, sep = '\t')
#write.xlsx(out_diff_gene_tpm_matrix, file = paste(expression_stat_dir, 'Diff.gene.tpm.xlsx', sep = '/'), sheetName = "gene.tpm", append = FALSE, row.names = F)

# setwd('C:\\work\\project\\mRNA\\2017\\OM-mRNA-20-Wheat-P20170502\\expression')
# diff_gene_tpm_matrix <- read.delim('Diff.gene.tpm.txt', row.names = 1)
# samples <- read.delim('group_sample', stringsAsFactors = F, header = F)
# colnames(samples) <- c('condition', 'sample')
# expression_stat_dir <- './'
# diff_gene_tpm_matrix <- diff_gene_tpm_matrix[1:1000, ]
# colnames(diff_gene_tpm_matrix) <- samples$sample
# source('C:\\work\\scripts\\atom\\R\\quantification\\quant_plot.R')

# diff_gene_tpm_matrix1 <- diff_gene_tpm_matrix[, 1:5]
# diff_gene_tpm_matrix1 <- diff_gene_tpm_matrix1[rowSums(diff_gene_tpm_matrix1) >= 2, ]
# diff_gene_tpm_matrix2 <- diff_gene_tpm_matrix[, 1:10]
# diff_gene_tpm_matrix2 <- diff_gene_tpm_matrix2[rowSums(diff_gene_tpm_matrix2) >= 2, ]
# diff_gene_tpm_matrix3 <- diff_gene_tpm_matrix
# diff_gene_tpm_matrix4 <- cbind(diff_gene_tpm_matrix, diff_gene_tpm_matrix)
# diff_gene_tpm_matrix5 <- cbind(diff_gene_tpm_matrix, diff_gene_tpm_matrix, diff_gene_tpm_matrix)

# diff_gene_tpm_matrix = read.delim(paste(expression_stat_dir, 'Diff.gene.tpm.txt', sep = '/'), row.names = 1)
# colnames(diff_gene_tpm_matrix) <- samples$sample

## heatmap
om_heatmap(plot_data = diff_gene_tpm_matrix, samples = samples, outdir = expression_stat_dir)

## boxplot
om_boxplot(plot_data = gene_tpm_matrix, samples = samples, outdir = expression_stat_dir)

## PCA plot
om_pca_plot(plot_data = gene_tpm_matrix, samples = samples, outdir = expression_stat_dir)

## sample correlation
om_correlation_plot(plot_data = cts, samples = samples, outdir = expression_stat_dir)
