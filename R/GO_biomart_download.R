suppressMessages(library("biomaRt"))
suppressMessages(library(argparser))

p <- arg_parser("|-- download ensembl go annotation --|")
p <- add_argument(p, "--gene_tr_file", help = "gene transcript map file")
p <- add_argument(p, "--species", help = 'biomart species database name')
p <- add_argument(p, "--output", help = 'output file')
argv <- parse_args(p)

gene_tr_file = argv$gene_tr_file
data_name = argv$species
out_name = argv$output

geneid_info <- read.table(gene_tr_file, header = F, sep = '\t')
gene_id <- unique(geneid_info$V1)
ensembl = useMart("ENSEMBL_MART_ENSEMBL")
ensembl = useDataset(data_name, mart = ensembl)
goids = getBM(attributes=c('ensembl_gene_id','go_id'), filters='ensembl_gene_id', values=gene_id, mart=ensembl)
write.table(goids, file = paste(out_name, ".go.txt", sep = ""), quote = F, sep = ',' , row.names = F )
