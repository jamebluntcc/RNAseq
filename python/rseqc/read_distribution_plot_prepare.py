'''
Usage:
    read_distribution_plot.py <sample_inf> <read_distribution_data_dir>
'''

from docopt import docopt
import pandas as pd
from os import path

if __name__ == '__main__':
    ## read arguments
    arguments = docopt(__doc__, version = 'plot read distribution v1')
    sample_inf = arguments['<sample_inf>']
    data_dir = arguments['<read_distribution_data_dir>']

    output_list = ['Group\tTotal_bases\tTag_count\tTags/Kb\tSample']
    group_sample_df = pd.read_table(sample_inf, header = None, index_col = 1)
    read_flag = 0
    for each_sample in group_sample_df.index:
        each_data_file = path.join(data_dir, '{}.read_distribution.txt'.format(each_sample))
        if not path.exists(each_data_file):
            print 'WARNING: {} not exists!'.format(each_data_file)
            continue
        ## strip rseqc output head and tail, add sample information for R plot
        with open(each_data_file) as each_data_file_inf:
            for eachline in each_data_file_inf:
                if eachline.startswith('Group'):
                    read_flag = 1
                    continue
                if eachline.startswith('='):
                    read_flag = 0
                if read_flag == 1:
                    eachline = '\t'.join(eachline.strip().split())
                    eachline = '{0}\t{1}'.format(eachline, each_sample)
                    output_list.append(eachline)

    ## output plot data
    plot_data = path.join(data_dir, 'read_distribution.summary.txt')
    with open(plot_data, 'w') as plot_data_inf:
        for each_record in output_list:
            plot_data_inf.write('{}\n'.format(each_record))
