'''
Usage:
    exp_table_split_by_group.py <group.sample> <gene.exp> <out.dir>

'''

from docopt import docopt
import pandas as pd
from os import path

arguments = docopt(__doc__, version='split expression table by group v1')

group_sample_file = arguments['<group.sample>']
gene_exp_file = arguments['<gene.exp>']
out_dir = arguments['<out.dir>']

group_sample_df = pd.read_table(group_sample_file, header = None, index_col = 0)
gene_exp_df = pd.read_table(gene_exp_file, index_col = 0)
for each_group in group_sample_df.index:
    each_group_exp_out = path.join(out_dir,'{}.tpm.txt'.format(each_group))
    each_group_sample = group_sample_df.loc[each_group, 1]
    if isinstance(each_group_sample, str):
        gene_exp_df.to_csv(each_group_exp_out, columns = [each_group_sample], sep = '\t')
    else:
        gene_exp_df.to_csv(each_group_exp_out, columns = each_group_sample, sep = '\t')
