'''
Usage:
    diff_summary.py <quant_dir>

'''

from docopt import docopt
import pandas as pd
from os import listdir
from os import path

if __name__ == '__main__':
    arguments = docopt(__doc__, version='v1')
    quant_dir = arguments['<quant_dir>']
    diff_dir = path.join(quant_dir, 'differential_analysis')
    compare_list = listdir(diff_dir)
    gene_count_file =
    for each_compare in compare_list:
        re_list = each_compare.split('_vs_')
        each_compare_diff_gene_list = []
        each_compare_status_list = []
        for each_re in re_list:
            each_re_file = path.join(diff_dir, each_compare, '{0}.{1}-UP.edgeR.DE_results.txt')
            each_re_gene_list = [each.strip() for each in open(each_re_file)]
            each_re_status_list = ['{}-UP'.format(each_re)] * len(each_re_gene_list)
            each_compare_diff_gene_list.extend(each_re_gene_list)
            each_compare_status_list.extend(each_re_status_list)
        each_compare_status_series = pd.Series(each_compare_status_list, index=each_compare_diff_gene_list, name = each_compare)
