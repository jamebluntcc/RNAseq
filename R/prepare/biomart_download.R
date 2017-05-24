##########################################################
###  download ensembl go and gene interpro annotation  ###
##########################################################

suppressMessages(library(biomaRt))
suppressMessages(library(argparser))
suppressMessages(library(GO.db))

p <- arg_parser("|-- download ensembl go and interpro annotation --|")
p <- add_argument(p, "--gene_tr_file", help = "gene transcript map file")
p <- add_argument(p, "--species", help = 'biomart species database name')
p <- add_argument(p, "--output", help = 'output file')
argv <- parse_args(p)

gene_tr_file = argv$gene_tr_file
data_name = argv$species
out_name = argv$output

geneid_info <- read.table(gene_tr_file, header = F, sep = '\t')
gene_id <- unique(geneid_info$V1)

ensembl_animal <- useMart("ENSEMBL_MART_ENSEMBL")
ensembl_animal_database <- listDatasets(ensembl_animal)
animal_name <- paste(data_name, 'gene_ensembl', sep = '_')

ensembl_plant <- useMart("plants_mart", host="plants.ensembl.org")
ensembl_plant_database <- listDatasets(ensembl_plant)
plant_name <- paste(data_name, 'eg_gene', sep = '_')

if (animal_name %in% ensembl_animal_database[,1]) {
  ensembl <- useDataset(animal_name, mart = ensembl_animal)
} else if (plant_name %in% ensembl_plant_database[,1]) {
  ensembl <- useDataset(plant_name, mart = ensembl_plant)
} else {
  stop('animal not in ensembl database')
}


## go annotation
goids <- getBM(attributes=c('ensembl_gene_id','go_id'), filters='ensembl_gene_id', values=gene_id, mart=ensembl)
write.table(goids, file = paste(out_name, ".go.txt", sep = ""), quote = F, sep = ',' , row.names = F )

## go detail information
go_file <- subset(goids, go_id != "")
go_list <- unique(go_file$go_id)
go_out <- as.data.frame(go_list)
go_term <- Term(go_list)
go_ontology <- Ontology(go_list)
go_out$go_ontology <- go_ontology
go_out$go_term <- go_term
write.table(go_out, file = paste(out_name, ".go_detail.txt", sep = ""), quote = F,  sep = "\t", row.names = F)

## gene interpro description
gene_inf <- getBM(attributes=c('ensembl_gene_id','external_gene_name','interpro'), filters='ensembl_gene_id', values=gene_id,, mart=ensembl)
write.table(gene_inf, file = paste(out_name, ".des.raw.txt", sep = ""), quote = F, sep = ',' , row.names = F )
