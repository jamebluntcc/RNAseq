'''
Usage:
    star_mapping_stats.py <sample_inf> <star_mapping_dir> <star_mapping_stats_prefix>

'''

from docopt import docopt
import pandas as pd
from os import path

def read_star_mapping_log(star_log_file):
    star_df = pd.read_table(star_log_file, header = None, sep = '|', index_col = 0)
    star_df = star_df.dropna()
    star_df = star_df.ix[4:]
    star_df.loc[:,1] = [each.strip() for each in star_df.loc[:,1]]
    return star_df


if __name__ == '__main__':
    arguments = docopt(__doc__, version = 'v1')
    sample_inf = arguments['<sample_inf>']
    star_mapping_dir = arguments['<star_mapping_dir>']
    star_mapping_stats_prefix = arguments['<star_mapping_stats_prefix>']

    sample_list = [each.strip().split()[1] for each in open(sample_inf)]
    star_log_file_list = [path.join(star_mapping_dir, each_sample, 'Log.final.out') for each_sample in sample_list]
    star_log_df_list = map(read_star_mapping_log, star_log_file_list)
    star_log_df = pd.concat(star_log_df_list, axis = 1)
    star_log_df.columns = sample_list
    star_log_out_df = star_log_df.T
    star_log_out_df.index.name = 'Sample'
    star_log_out_df.columns = [each.strip() for each in star_log_out_df.columns]

    ## output
    star_mapping_stats_txt = '{}.txt'.format(star_mapping_stats_prefix)
    star_log_out_df.to_csv(star_mapping_stats_txt, sep = '\t')
    star_mapping_stats_excel = '{}.xlsx'.format(star_mapping_stats_prefix)
    pd.formats.format.header_style = None
    writer = pd.ExcelWriter(star_mapping_stats_excel, engine='xlsxwriter',options={'strings_to_urls': False})
    star_log_out_df.to_excel(writer, 'star_mapping_stats')
    writer.save()

    ## output plot data
    plot_data = star_log_out_df.loc[:,['Number of input reads', 'Uniquely mapped reads number', 'Number of reads mapped to multiple loci']]
    plot_data.columns = ['total_reads', 'unique_mapped_reads', 'multiple_mapped_reads']
    star_mapping_stats_plot = '{}.plot'.format(star_mapping_stats_prefix)
    plot_data.to_csv(star_mapping_stats_plot, sep = '\t')
