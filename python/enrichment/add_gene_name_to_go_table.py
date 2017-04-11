'''
Usage:
    add_gene_name_to_go_table.py <go.txt> <biomart.anno.txt> <diff.gene.list> <go.enrich.result> <go.anno.table>
'''

from docopt import docopt
import pandas as pd

arguments = docopt(__doc__, version = 'v1')

go_file = arguments['<go.txt>']
biomart_anno = arguments['<biomart.anno.txt>']
diff_gene_file = arguments['<diff.gene.list>']
go_enrich_file = arguments['<go.enrich.result>']
go_anno_table = arguments['<go.anno.table>']

go_file_tb = pd.read_table(go_file, sep = ',')
go_file_tb.columns = ['gene', 'go']
go_file_tb = go_file_tb.loc[go_file_tb.go != '', :]

biomart_anno_tb = pd.read_table(biomart_anno)
biomart_anno_tb.columns = ['gene', 'gene_name', 'interpro_id', 'interpro_des']

diff_gene_file_tb = pd.read_table(diff_gene_file, header = False)
diff_gene_file_tb.columns = ['gene']

diff_gene_go_tb = pd.merge(diff_gene_file_tb, go_file_tb)
diff_gene_go_anno_tb = pd.merge(diff_gene_go_tb, biomart_anno_tb)

diff_gene_go_anno_dict = {}
diff_gene_go_anno_groupby_gene = diff_gene_go_anno_tb.groupby('go')['gene']
diff_gene_go_anno_groupby_name = diff_gene_go_anno_tb.groupby('go')['gene_name']
diff_gene_go_anno_dict['gene'] = map('|'.join, diff_gene_go_anno_groupby_gene.unique())
diff_gene_go_anno_dict['gene_name'] = map('|'.join, diff_gene_go_anno_groupby_name.unique())
