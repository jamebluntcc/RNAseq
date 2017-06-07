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

p <- arg_parser("read kallisto quant files generate expression matrix")
p <- add_argument(p, '--kallisto_dir',  help = 'kallisto quantification directory')
p <- add_argument(p, '--sample_inf', help = 'sample information with sample names and group names')
p <- add_argument(p, '--gene2tr',    help = 'gene id and transcript id mapping file')
p <- add_argument(p, '--out_dir',    help = 'diff analyssi output directory')
argv <- parse_args(p)

# ## for test
# source('C:\\work\\scripts\\atom\\R\\quantification\\quant_plot.R')
# sample_inf <- '../fastqc/group_sample'
# quant_dir <- './kallisto/'
# gene2tr_file <- 'gene_trans.map'
# outdir <- './'

## read parameters
sample_inf <- argv$sample_inf
kallisto_dir <- argv$kallisto_dir
gene2tr_file <- argv$gene2tr
outdir <- argv$out_dir

## directory prepare
expression_stat_dir <- file.path(outdir, 'expression_summary')
diff_dir <- file.path(outdir, 'differential_analysis')
dir.create(expression_stat_dir, showWarnings = FALSE)
dir.create(diff_dir, showWarnings = FALSE)

samples <- read.delim(sample_inf, stringsAsFactors = F, header = F)
colnames(samples) <- c('condition', 'sample')
files <- file.path(kallisto_dir, samples$sample, "abundance.tsv")
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

## boxplot
om_boxplot(plot_data = gene_tpm_matrix, samples = samples, outdir = expression_stat_dir)

## PCA plot
om_pca_plot(plot_data = gene_tpm_matrix, samples = samples, outdir = expression_stat_dir)

## sample correlation
om_correlation_plot(plot_data = cts, samples = samples, outdir = expression_stat_dir)

pdf_example_tpm_df <- out_gene_tpm_matrix[1:100, 1:5]
write.table(pdf_example_tpm_df, file = paste(expression_stat_dir, 'pdf.example.Gene.tpm.txt', sep = '/'), quote = F, row.names = F, sep = '\t')
