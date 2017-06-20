#! /usr/bin/python

'''
Usage:
    circ_pipe.py <circ_pred_files> <mbs_pred> <circ_exp_matrix>

'''

import pandas as pd
from docopt import docopt
from os import path


class circ_inf:

    def __init__(self):
        self.bed = ""
        self.type = ""
        self.gene = ""
        self.name = "--"
        self.sample = []
        self.count = []
        self.miRNA = []


class circ_parser:

    def __init__(self):
        self.circ_pred = ""
        self.circ_dict = {}
        self.junction_count_limit = 2

    def circexplorer_parser(self, each_line_array):
        bedline_elements = each_line_array[:12]
        bedline_elements_list = [str(each) for each in bedline_elements]
        circ_name = '{0}:{1}-{2},{3}'.format(
            bedline_elements_list[0], bedline_elements_list[1], bedline_elements_list[2], bedline_elements_list[5])
        bedline_elements_list[3] = circ_name
        bedline = '\t'.join(bedline_elements_list)
        junction_count = each_line_array[12]
        if junction_count >= self.junction_count_limit:
            self.circ_dict[circ_name] = circ_inf()
            self.circ_dict[circ_name].bed = bedline
            self.circ_dict[circ_name].type = each_line_array[13]
            self.circ_dict[circ_name].gene = each_line_array[14]

    def get_info(self):
        circ_pred_df = pd.read_table(self.circ_pred, header=None)
        for each in circ_pred_df.index:
            self.circexplorer_parser(circ_pred_df.loc[each, ])


if __name__ == '__main__':
    arguments = docopt(__doc__, version='1.0')
    circ_pred_files = arguments['<circ_pred_files>']
    mbs_pred = arguments['<mbs_pred>']
    circ_exp_matrix = arguments['<circ_exp_matrix>']

    all_circ_dict = {}
    circ_pred_file_list = [each.strip() for each in open(circ_pred_files)]
    for each_file in circ_pred_file_list:
        each_circ_parser = circ_parser()
        each_circ_parser.circ_pred = each_file
        each_circ_parser.get_info()
        all_circ_dict.update(each_circ_parser.circ_dict)

    # get miRNA binding info
    mbs_df = pd.read_table(mbs_pred, header=None, index_col=1)
    for each_circ in mbs_df.index:
        if isinstance(mbs_df.loc[each_circ, 0], str):
            all_circ_dict[each_circ].miRNA = [mbs_df.loc[each_circ, 0]]
        else:
            all_circ_dict[each_circ].miRNA = list(mbs_df.loc[each_circ, 0])

    # add circRNA information
    circRNA_inf_dict = {}
    circ_list = all_circ_dict.keys()
    for each_circ in circ_list:
        each_circ_type = all_circ_dict[each_circ].type
        each_circ_miRNAs = all_circ_dict[each_circ].miRNA
        each_circ_miRNA_num = len(each_circ_miRNAs)
        each_circ_miRNA_mbs = len(set(each_circ_miRNAs))
        each_circ_miRNAs = ','.join(set(each_circ_miRNAs))
        if not each_circ_miRNAs:
            each_circ_miRNAs = '--'
        circRNA_inf_dict.setdefault('circRNA_type', []).append(each_circ_type)
        circRNA_inf_dict.setdefault(
            'miRNA_binding_number', []).append(each_circ_miRNA_num)
        circRNA_inf_dict.setdefault(
            'miRNA_binding_site_number', []).append(each_circ_miRNA_mbs)
        circRNA_inf_dict.setdefault('miRNA_ids', []).append(each_circ_miRNAs)
    circRNA_inf_dict['circRNA_id'] = circ_list
    circRNA_inf_df = pd.DataFrame(circRNA_inf_dict)
    circ_exp_df = pd.read_table(circ_exp_matrix)
    circ_exp_inf_df = pd.merge(circ_exp_df, circRNA_inf_df, left_on="chrom:beg-end",
                               right_on="circRNA_id", how="left").drop('circRNA_id', axis=1)
    circ_exp_matrix_name = path.basename(circ_exp_matrix)
    circ_exp_matrix_dir = path.dirname(circ_exp_matrix)
    circ_exp_matrix_name_pre, circ_exp_matrix_name_suf = path.splitext(
        circ_exp_matrix_name)
    circ_exp_inf_out_path = path.join(circ_exp_matrix_dir, '{0}.ann.{1}'.format(
        circ_exp_matrix_name_pre, circ_exp_matrix_name_suf))
    circ_exp_inf_df.to_csv(circ_exp_inf_out_path, sep='\t', index=False)
