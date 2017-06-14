'''
Usage:
    get_gene_pos.py <gtf> <gene.pos>

'''

import sys
from docopt import docopt

RNASEQ_LIBPATH = '/home/public/scripts/RNAseq/python/lib/'
sys.path.insert(0, RNASEQ_LIBPATH)
import RNAseq_tools

arguments = docopt(__doc__, version='gtf to gene position v1')

gtf_file = arguments['<gtf>']
pos_file = arguments['<gene.pos>']

tr_dict = RNAseq_tools.get_transcript_info(gtf_file)
gene_dict = RNAseq_tools.get_gene_info(tr_dict)

with open(pos_file, 'w') as pos_file_inf:
    pos_file_inf.write('gene_id\tlocus\tstrand\n')
    for each_gene in gene_dict:
        each_start = min(gene_dict[each_gene]['tss'])
        each_end = max(gene_dict[each_gene]['tts'])
        each_chrom = gene_dict[each_gene]['chrom']
        each_strand = gene_dict[each_gene]['strand']
        pos_file_inf.write('{0}\t{1}:{2}-{3}\t{4}\n'.format(each_gene, each_chrom, each_start, each_end, each_strand))
