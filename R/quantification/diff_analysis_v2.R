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


p <- arg_parser("perform differential analysis")
p <- add_argument(p, '--kallisto_dir',  help = 'kallisto quantification directory')
p <- add_argument(p, '--tpm_table',  help = 'quantification normalized tpm table')
p <- add_argument(p, '--compare',  help = 'compare name')
p <- add_argument(p, '--sample_inf', help = 'sample information with sample names and group names')
p <- add_argument(p, '--gene2tr',    help = 'gene id and transcript id mapping file')
p <- add_argument(p, '--out_dir', help = 'output directory')
p <- add_argument(p, '--qvalue',     help = 'diff gene qvalue cutoff', default = 0.05)
p <- add_argument(p, '--logfc',      help = 'diff gene logfc cutoff',  default = 1)
argv <- parse_args(p)

# test
# kallisto_dir <- './kallisto/'
# tpm_table <- './expression_summary/Gene.tpm.txt'
# compare <- '3BZH_vs_2BZH'
# outdir <- './differential_analysis/3BZH_vs_2BZH'
# logfc <- 1
# qvalue <- 0.05

# read aguments
kallisto_dir <- argv$kallisto_dir
tpm_table <- argv$tpm_table
compare <- argv$compare
sample_inf <- argv$sample_inf
gene2tr_file <- argv$gene2tr
outdir <- argv$out_dir
qvalue <- argv$qvalue
logfc <- argv$logfc

dir.create(outdir, showWarnings = FALSE)
samples <- read.delim(sample_inf, stringsAsFactors = F, header = F)
colnames(samples) <- c('condition', 'sample')
gene2tr <- read.delim(gene2tr_file, header = FALSE)
colnames(gene2tr) <- c('gene_id', 'transcript_id')
tx2gene <- gene2tr[,c('transcript_id', 'gene_id')]

gene_tpm_matrix <- read.delim(tpm_table, check.names = F,row.names = 1)

each_pair <- unlist(strsplit(compare, split = "_vs_"))
con1_sample <- samples[samples$condition == each_pair[1],'sample']
con2_sample <- samples[samples$condition == each_pair[2],'sample']
each_pair_samples <- samples[samples$sample %in% c(con1_sample, con2_sample),]
each_pair_files <- file.path(kallisto_dir, each_pair_samples$sample, "abundance.tsv")
names(each_pair_files) <- each_pair_samples$sample
each_pair_txi <- tximport(each_pair_files, type = "kallisto", tx2gene = tx2gene)
each_pair_cts <- each_pair_txi$counts
each_pair_cts <- each_pair_cts[, c(con1_sample, con2_sample)]
each_pair_cts <- each_pair_cts[rowSums(each_pair_cts)>=2,]
each_pair_genes <- row.names(each_pair_cts)

## normalization
normMat <- each_pair_txi$length
normMat <- normMat[each_pair_genes,]
normMat <- normMat/exp(rowMeans(log(normMat)))
o <- log(calcNormFactors(each_pair_cts/normMat)) + log(colSums(each_pair_cts/normMat))
conditions = factor(c(rep(each_pair[1], length(con1_sample)), rep(each_pair[2], length(con2_sample))))
y <- DGEList(each_pair_cts, group=conditions)
y$offset <- t(t(log(normMat)) + o)
y <- estimateDisp(y)
## diff
et <- exactTest(y, pair=rev(each_pair))
tTags <- topTags(et,n=NULL)
new_tTags <- tTags$table
new_tTags <- new_tTags[, !(names(new_tTags) %in% c("logCPM"))]
each_pair_matrix <- gene_tpm_matrix[,c(con1_sample, con2_sample)]
merged_df <- merge(each_pair_matrix, new_tTags, by.x = 0, by.y = 0, all.y = T)
sorted_merged_df <- arrange(merged_df, FDR)
colnames(sorted_merged_df)[1] <- 'Gene_ID'

out_file_name_prefix <- paste(outdir, '/', compare, sep = '')
up_regulate_name_prefix <- paste(outdir, '/', compare, '.',each_pair[1], '-UP',  sep = '')
down_regulate_name_prefix <- paste(outdir, '/', compare, '.',each_pair[2], '-UP', sep = '')

diff_genes <- c()
up_regulate_df <- filter(sorted_merged_df, logFC >= logfc, FDR <= qvalue)
down_regulate_df <- filter(sorted_merged_df, logFC <= -(logfc), FDR <= qvalue)
diff_genes <- c(diff_genes, up_regulate_df$Gene_ID, down_regulate_df$Gene_ID)
write.table(sorted_merged_df, file=paste(out_file_name_prefix, 'edgeR.DE_results.txt', sep = '.'), sep='\t', quote=F, row.names=F)
#write.xlsx(sorted_merged_df, file = paste(out_file_name_prefix, 'edgeR.DE_results', 'xlsx', sep = '.'), sheetName = each_compare_name, append = F, row.names = F)
if (dim(up_regulate_df)[1] > 0) {
  write.table(up_regulate_df, file=paste(up_regulate_name_prefix, 'edgeR.DE_results.txt', sep = '.'), sep='\t', quote=F, row.names=F)
  #write.xlsx(up_regulate_df, file = paste(up_regulate_name_prefix, 'edgeR.DE_results', 'xlsx', sep = '.'), sheetName = each_compare_name, append = F, row.names = F)
  write(as.character(up_regulate_df$Gene_ID), file = paste(up_regulate_name_prefix, 'edgeR.DE_results.diffgenes.txt', sep = '.'), sep = '\n')
}
if (dim(down_regulate_df)[1] > 0) {
  write.table(down_regulate_df, file=paste(down_regulate_name_prefix, 'edgeR.DE_results.txt', sep = '.'), sep='\t', quote=F, row.names=F)
  #write.xlsx(down_regulate_df, file = paste(down_regulate_name_prefix, 'edgeR.DE_results', 'xlsx', sep = '.'), sheetName = each_compare_name, append = F, row.names = F)
  write(as.character(down_regulate_df$Gene_ID), file = paste(down_regulate_name_prefix, 'edgeR.DE_results.diffgenes.txt', sep = '.'), sep = '\n')
}
## write diff gene list
if (length(diff_genes) > 0) {
  write(as.character(diff_genes), file = paste(out_file_name_prefix, 'ALL.edgeR.DE_results.diffgenes.txt', sep = '.'), sep = '\n')
}
## volcano plot
om_volcano_plot(sorted_merged_df, compare, logfc, qvalue, outdir)
