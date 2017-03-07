'''
Usage:
    anoroad_data_rename.py <ano_clean_data_dir> <rename_fq_dir>

'''

from os import path
from os import listdir
from os import system
import glob
import re
from docopt import docopt


if __name__ == '__main__':
    arguments = docopt(__doc__, version='anoroad data rename 1.0')
    ano_clean_data_dir = path.abspath(arguments['<ano_clean_data_dir>'])
    rename_fq_dir = path.abspath(arguments['<rename_fq_dir>'])

    if not path.exists(rename_fq_dir):
        system('mkdir -p {}'.format(rename_fq_dir))
    sample_list = listdir(ano_clean_data_dir)
    for each_sample in sample_list:
        each_sample_path = path.join(ano_clean_data_dir, each_sample)
        each_sample_fq_files = glob.glob('{}/*gz'.format(each_sample_path))
        for each_fq in each_sample_fq_files:
            each_fq_name = path.basename(each_fq)
            new_fq_name = re.sub("R(\w).fq", "\g<1>.clean.fq",each_fq_name)
            system('ln -s {0} {1}/{2}'.format(each_fq, rename_fq_dir, new_fq_name))
