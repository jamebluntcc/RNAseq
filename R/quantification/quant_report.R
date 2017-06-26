suppressMessages(library(argparser))
suppressMessages(library(plyr))
suppressMessages(library(dplyr))
suppressMessages(library(data.table))
source("/public/scripts/RNAseq/R/quantification/quant_plot.R")

p <- arg_parser("for expression analysis report plot")
p <- add_argument(p, "--quant_dir", help = "quantification directory")
p <- add_argument(p, "--sample_inf", help = "sample information with sample names and group names")
p <- add_argument(p, "--qvalue", help = "diff gene qvalue cutoff", default = 0.05)
p <- add_argument(p, "--logfc", help = "diff gene logfc cutoff", default = 1)
argv <- parse_args(p)

MERGED_VOL_PLOT_NUM = 6
DIFF_HEATMAP_GENE = 40000
CUT_TREE_PER = 20
MIN_CLUSTER_NUM = 10
MIN_CLUSTER_POR = 0.005
# for test setwd('C:\\work\\pipe\\quantification') quant_dir <- './'
# sample_inf <- '../fastqc/group_sample'
# source('C:\\work\\scripts\\atom\\R\\quantification\\quant_plot.R')

# read parameters
quant_dir <- argv$quant_dir
sample_inf <- argv$sample_inf
qvalue <- argv$qvalue
logfc <- argv$logfc

samples <- read.delim(sample_inf, stringsAsFactors = F, header = F)
colnames(samples) <- c("condition", "sample")
all_combine = combn(unique(samples$condition), 2)
diff_dir <- file.path(quant_dir, "differential_analysis")
exp_dir <- file.path(quant_dir, "expression_summary/")


diff_df_list = list()
plot_number = min(c(16, dim(all_combine)[2]))

diff_genes <- c()
for (i in seq(dim(all_combine)[2])) {

  each_pair <- all_combine[, i]
  each_compare <- paste(each_pair[1], "_vs_", each_pair[2], sep = "")
  each_compare_name = paste(each_compare, "edgeR.DE_results.txt", sep = ".")
  each_compare_file = file.path(diff_dir, each_compare, each_compare_name)
  each_compare_df <- fread(each_compare_file)
  each_compare_plot <- each_compare_df[, c("Gene_ID", "logFC", "FDR")]
  each_compare_plot$compare <- each_compare
  each_compare_diff_df <- filter(each_compare_df, abs(logFC) >= logfc, FDR <= qvalue)
  each_diff_genes <- each_compare_diff_df$Gene_ID
  diff_genes <- c(diff_genes, each_diff_genes)
  if (i <= MERGED_VOL_PLOT_NUM) {
    diff_df_list[[i]] <- each_compare_plot
  }
}

diff_df <- ldply(diff_df_list, data.frame)
om_volcano_plot(diff_df, "ALL", logfc, qvalue, exp_dir)
pdf_example_diff_out <- file.path(exp_dir, "pdf.example.diff.table.txt")
html_example_diff_out <- file.path(exp_dir, "html.example.diff.table.txt")
pdf_example_diff_df <- each_compare_df[1:100, c("Gene_ID", "logFC", "PValue", "FDR")]

write.table(pdf_example_diff_df, file = pdf_example_diff_out, sep = "\t", quote = F,
  row.names = F)
write.table(each_compare_df, file = html_example_diff_out, sep = "\t", quote = F,
  row.names = F)

diff_gene_id <- unique(diff_genes)
gene_count_matrix <- file.path(exp_dir, "Gene.count.txt")
gene_tpm_matrix <- file.path(exp_dir, "Gene.tpm.txt")
gene_count_matrix_df <- read.delim(gene_count_matrix, check.names = F)
gene_tpm_matrix_df <- read.delim(gene_tpm_matrix, check.names = F)
out_diff_gene_count_matrix <- filter(gene_count_matrix_df, Gene_ID %in% diff_gene_id)
out_diff_gene_tpm_matrix <- filter(gene_tpm_matrix_df, Gene_ID %in% diff_gene_id)

if (dim(out_diff_gene_count_matrix)[1] > 0) {
  write.table(out_diff_gene_count_matrix, file = paste(exp_dir, "Diff.gene.count.txt",
    sep = "/"), quote = F, row.names = F, sep = "\t")
  write.table(out_diff_gene_tpm_matrix, file = paste(exp_dir, "Diff.gene.tpm.txt",
    sep = "/"), quote = F, row.names = F, sep = "\t")
  diff_gene_tpm_matrix <- out_diff_gene_tpm_matrix[, -1]
  rownames(diff_gene_tpm_matrix) <- out_diff_gene_tpm_matrix[, 1]

  # plot heatmap
  diff_gene_count = dim(diff_gene_tpm_matrix)[1]
  if (diff_gene_count <= DIFF_HEATMAP_GENE) {
    om_heatmap(plot_data=diff_gene_tpm_matrix, samples=samples, outdir=exp_dir)
  } else {
    gene_mean_exp <- sort(rowMeans(diff_gene_tpm_matrix), decreasing = T)
    top_genes <- names(gene_mean_exp[1:DIFF_HEATMAP_GENE])
    diff_gene_tpm_matrix <- diff_gene_tpm_matrix[top_genes, ]
    om_heatmap(plot_data=diff_gene_tpm_matrix, samples=samples, outdir=exp_dir)
    diff_gene_count <- dim(diff_gene_tpm_matrix)[1]
  }

  # cluster plot
  cluster_data_dir <- file.path(exp_dir, "cluster_data")
  dir.create(cluster_data_dir, showWarnings = F)
  diff_matrix <- as.matrix(diff_gene_tpm_matrix)
  log_diff_matrix <- log10(diff_matrix + 1)
  # center rows, mean substracted
  scale_log_diff_matrix = t(scale(t(log_diff_matrix), scale = F))

  # gene clustering according to centered distance values.
  gene_dist = dist(scale_log_diff_matrix, method = "euclidean")
  hc_genes = hclust(gene_dist, method = "complete")
  gene_partition_assignments <- cutree(as.hclust(hc_genes), h = CUT_TREE_PER/100 * max(hc_genes$height))

  max_cluster_count = max(gene_partition_assignments)
  cluster_num_cutoff = max(c(MIN_CLUSTER_NUM, diff_gene_count * MIN_CLUSTER_POR))

  all_partition_list <- list()
  m = 1
  for (i in 1:max_cluster_count) {
    partition_i = (gene_partition_assignments == i)
    partition_data = scale_log_diff_matrix[partition_i, , drop = F]
    cluster_name <- paste("cluster", i, sep = "_")
    partition_data_df <- as.data.frame(partition_data)
    partition_data_df <- cbind(Gene_id = rownames(partition_data_df), partition_data_df)
    write.table(partition_data_df, file = paste(cluster_data_dir, "/", cluster_name,
      ".txt", sep = ""), quote = F, row.names = F, sep = "\t")
    partition_data_df$cluster <- cluster_name
    melt_partition_data_df <- melt(partition_data_df, id = c("cluster", "Gene_id"))
    out_prefix <- file.path(cluster_data_dir, cluster_name)
    cluster_plot(plot_data = melt_partition_data_df, out_prefix = out_prefix)
    if (dim(partition_data)[1] > cluster_num_cutoff) {
      all_partition_list[[m]] <- partition_data_df
      m <- m + 1
    }
  }

  all_cluster_df <- ldply(all_partition_list, data.frame)
  colnames(all_cluster_df) <- colnames(partition_data_df)
  melt_all_cluster_df <- melt(all_cluster_df, id = c("cluster", "Gene_id"))
  melt_all_cluster_df$variable <- factor(melt_all_cluster_df$variable, levels = samples$sample)
  melt_all_cluster_df$cluster <- factor(melt_all_cluster_df$cluster, levels = unique(melt_all_cluster_df$cluster))
  all_cluster_prefix = file.path(exp_dir, 'Diff.genes.cluster')
  cluster_plot(plot_data = melt_all_cluster_df, out_prefix = all_cluster_prefix)
}
