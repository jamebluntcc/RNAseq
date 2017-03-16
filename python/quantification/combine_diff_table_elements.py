'''
Usage:
    combine_diff_table_elements.py <files_for_diff_table> <merged.file>

'''

from docopt import docopt
import pandas as pd

arguments = docopt(__doc__, version='combine_diff_table_elements v1')

all_file_list = arguments['<files_for_diff_table>'].split(',')
merged_file = arguments['<merged.file>']

for n, each_file in enumerate(all_file_list):
    each_file_df = pd.read_table(each_file, index_col = 0, sep = '\t')
    if n == 0:
        global merged_df
        merged_df = each_file_df
    else:
        merged_df = pd.merge(merged_df, each_file_df, left_index = True, right_index = True, how = 'outer')

filter_merged_df = merged_df.loc[merged_df.FDR.notnull(), ]
filter_merged_df = filter_merged_df.fillna('--')
filter_merged_df.index.name = 'Gene_id'
sorted_filter_merged_df = filter_merged_df.sort_values(by=['FDR', 'logFC'])
sorted_filter_merged_df.to_csv(merged_file, sep = '\t')
