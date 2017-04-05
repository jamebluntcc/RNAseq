suppressMessages(library(dplyr))
suppressMessages(library(argparser))

source('/home/public/scripts/RNAseq/R/enrichment/go_lib.R')
# source('../atom/R/enrichment/enrich_lib.R')

p <- arg_parser("perform go analysis")
p <- add_argument(p, '--quant_dir',  help = 'quantification directory')
p <- add_argument(p, '--go_anno',    help = 'gene id and transcript id mapping file')
p <- add_argument(p, '--topgo_anno', help = 'top go annotation file')
p <- add_argument(p, '--gene_length',help = 'gene length file')
p <- add_argument(p, '--out_dir',    help = 'diff analyssi output directory')
argv <- parse_args(p)

## for test
# quant_dir <- './'
# diff_dir <- file.path(quant_dir, 'differential_analysis')
# go_anno_file <- './gga.txt'
# gene_length_file <- './Gallus_gallus.Galgal4.85.gene.length'
# out_dir <- './enrich/'
#topgo_anno_file <- './gga_gene_go.txt'

## get the parameter
diff_dir <- file.path(argv$quant_dir, 'differential_analysis')
go_anno_file <- argv$go_anno
topgo_anno_file <- argv$topgo_anno
gene_length_file <- argv$gene_length
out_dir <- argv$out_dir

## output directories
go_out_dir <- file.path(out_dir, 'go')
dir.create(go_out_dir, recursive = 1, showWarnings = F)
gene_length_df <- read.delim(gene_length_file, header = F)
go_anno_df <- read.delim(go_anno_file, sep = ',')

## get each compare diff file
all_compares <- basename(list.dirs(diff_dir, recursive = 0))

for (i in seq(length(all_compares))) {
  each_compare <- all_compares[i]
  each_compare_go_out_dir <- file.path(go_out_dir, each_compare)
  dir.create(each_compare_go_out_dir, showWarnings = FALSE)
  each_top_go_out_dir <- file.path(each_compare_go_out_dir, 'DAG')
  dir.create(each_top_go_out_dir, recursive = 1, showWarnings = F)
  up_down_name <- unlist(strsplit(each_compare,split = "_vs_"))
  all_diff_genes <- c()
  for (j in seq(length(up_down_name))) {
    each_reg <- up_down_name[j]
    each_diff_file_name <- paste(each_compare, '.', each_reg, '-UP.edgeR.DE_results.txt', sep = '')
    each_up_down_file <- file.path(diff_dir, each_compare, each_diff_file_name)
    each_up_down_df <- read.delim(each_up_down_file, stringsAsFactors = F)
    each_reg_diff_genes <- each_up_down_df[,'Gene_ID']
    all_diff_genes <- c(all_diff_genes, each_reg_diff_genes)
    each_out_preifx_name <- paste(each_compare, '.', each_reg, '-UP.go.enrichment', sep = '')
    each_reg_out_prefix_path <- file.path(each_compare_go_out_dir, each_out_preifx_name)
    each_reg_enrich_result <- run_goseq(each_reg_diff_genes, gene_length_df, go_anno_df, each_reg_out_prefix_path)
    ## run topgo
    each_top_go_name <- paste(each_reg, 'UP', sep = '-')
    run_topgo(topgo_anno_file, each_reg_diff_genes, each_reg_enrich_result, each_top_go_name,each_top_go_out_dir)
  }
  all_diff_out_preifx_name <- paste(each_compare, '.', 'ALL', '.go.enrichment', sep = '')
  all_diff_out_preifx_path <- file.path(each_compare_go_out_dir, all_diff_out_preifx_name)
  all_diff_enrich_result <- run_goseq(all_diff_genes, gene_length_df, go_anno_df, all_diff_out_preifx_path)
  ## run topgo
  all_top_go_name <- 'ALL'
  run_topgo(topgo_anno_file, all_diff_genes, all_diff_enrich_result, all_top_go_name, each_top_go_out_dir)
}

# library(topGO)
# 
# diff_genes <- all_diff_genes
# 
run_topgo <- function(gene_go_map, diff_genes, enrich_result, name, out_dir) {
  geneID2GO <- readMappings(file = topgo_anno_file)
  geneNames <- names(geneID2GO)
  enrich_result_df <- out_go
  geneList <- factor(as.integer(geneNames %in% diff_genes))
  names(geneList) <- geneNames
  go_catogary_vector <- c('MF','CC','BP')
  enrich_result_df <- filter(enrich_result_df, numInCat >= 5)

  for (i in 1:length(go_catogary_vector)) {
    go_catogary <- go_catogary_vector[i]
    each_enrich_result <- filter(enrich_result_df, ontology == go_catogary)
    each_go_qvalue <- each_enrich_result$qvalue
    each_go_qvalue[which(each_go_qvalue == 0)] <- 1e-100
    names(each_go_qvalue)<-each_enrich_result[,1]
    if (dim(each_enrich_result)[1] < 2) {
      out_info <- paste('Too little gene annotated to ',go_catogary,sep="")
      print(out_info)
    } else {
      GOdata <- new("topGOdata", ontology = go_catogary, allGenes = geneList,
                    annot = annFUN.gene2GO, gene2GO = geneID2GO,nodeSize = 5)
      if (dim(each_enrich_result)[1] <= 10) {
        pdf(file = paste(out_dir,'/',name,'.',go_catogary,'.GO.DAG.pdf',sep = ''),
            width = 8, height = 8)
        showSigOfNodes(GOdata, each_go_qvalue, firstSigNodes = 1, useInfo = 'all')
        dev.off()
        Cairo(file=paste(out_dir,'/',name,'.',go_catogary,'.GO.DAG.png',sep = ''),
              type="png",
              units="in",
              width=8,
              height=8,
              pointsize=12,
              dpi=300,
              bg = "white")
        showSigOfNodes(GOdata, each_go_qvalue, firstSigNodes = 1, useInfo = 'all')
        dev.off()
      } else {
        pdf(file = paste(out_dir,'/',name,'.',go_catogary,'.GO.DAG.pdf',sep = ''),
            width = 8, height = 8)
        showSigOfNodes(GOdata, each_go_qvalue, firstSigNodes = 5, useInfo = 'all')
        dev.off()
        Cairo(file=paste(out_dir,'/',name,'.',go_catogary,'.GO.DAG.png',sep = ''),
              type="png",
              units="in",
              width=8,
              height=8,
              pointsize=12,
              dpi=300,
              bg = "white")
        showSigOfNodes(GOdata, each_go_qvalue, firstSigNodes = 5, useInfo = 'all')
        dev.off()
      }
    }
  }
}
