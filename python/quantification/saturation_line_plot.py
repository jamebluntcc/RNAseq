'''
Usage:
    saturation_line_plot.py <rseqc_saturation_output> <saturation_plot_prefix>
'''

from __future__ import division
from docopt import docopt
import pandas as pd
from bisect import bisect
import subprocess

arguments = docopt(__doc__, version='saturation line plot v1')

rpkm_cutoff_dict = {
    0:'0-1',
    1:'1-5',
    5:'5-10',
    10:'10-50',
    50:'50-100',
    100:'>100'
}

def run_cmd(cmd):
    p = subprocess.Popen(cmd, shell=False, universal_newlines=True, stdout=subprocess.PIPE)
    ret_code = p.wait()
    output = p.communicate()[0]
    return output

rpkm_cutoff_list = [0, 1, 5, 10, 50, 100]

saturation_data = arguments['<rseqc_saturation_output>']
saturation_prefix = arguments['<saturation_plot_prefix>']

## prepare plot data
saturation_data_df = pd.read_table(saturation_data, index_col = 3)
saturation_data_exp_df = saturation_data_df.loc[saturation_data_df['100%'] > 0,:]

for x in range(5, 100, 5):
    each_name_data = '{}%'.format(x)
    each_name_stat = 'p_{}'.format(x)
    saturation_data_exp_df.loc[:, each_name_stat] = saturation_data_exp_df[each_name_data]/saturation_data_exp_df['100%']
    saturation_data_exp_df.loc[abs(saturation_data_exp_df[each_name_stat] - 1) <= 0.15, each_name_stat] = 1
    saturation_data_exp_df.loc[abs(saturation_data_exp_df[each_name_stat] - 1) > 0.15, each_name_stat] = 0

for each_transcript in saturation_data_exp_df.index:
    each_transcript_exp = saturation_data_exp_df.loc[each_transcript, '100%']
    each_transcript_exp_flag = bisect(rpkm_cutoff_list, each_transcript_exp)
    saturation_data_exp_df.loc[each_transcript, 'exp_region'] = rpkm_cutoff_dict[rpkm_cutoff_list[each_transcript_exp_flag-1]]

saturation_stat_dict = {}
saturation_stat_index = [rpkm_cutoff_dict[each] for each in rpkm_cutoff_list]
saturation_stat_columns = ['p_{}'.format(each) for each in range(5, 100, 5)]
for x in range(5, 100, 5):
    each_name_stat = 'p_{}'.format(x)
    each_name_stat_group = saturation_data_exp_df.groupby('exp_region')[each_name_stat]
    for each_fpkm in rpkm_cutoff_list:
        each_cutoff = rpkm_cutoff_dict[each_fpkm]
        each_saturation_value = each_name_stat_group.value_counts()[each_cutoff][1]/float(sum(each_name_stat_group.value_counts()[each_cutoff]))
        if each_name_stat not in saturation_stat_dict:
            saturation_stat_dict[each_name_stat] = [each_saturation_value]
        else:
            saturation_stat_dict[each_name_stat].append(each_saturation_value)
saturation_stat_df = pd.DataFrame(saturation_stat_dict, index = saturation_stat_index, columns = saturation_stat_columns)
saturation_stat_df.index.name = 'rpkm'
saturation_stat_df.columns = range(5, 100, 5)
saturation_stat_df.to_csv('{}.saturation.txt'.format(saturation_prefix), sep = '\t')

## plot
run_cmd(['/home/public/software/R-3.3.2/bin/Rscript',
        '/home/public/scripts/RNAseq/R/quantification/saturation.plot.R',
        '{}.saturation.txt'.format(saturation_prefix),
        saturation_prefix])
