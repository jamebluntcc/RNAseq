import subprocess
import luigi
from os import path
from os import system
import sys
import pandas as pd
import itertools

script_path = path.dirname(path.abspath(__file__))
RNAseq_lib_path = path.join(script_path, '..')
sys.path.insert(0, RNAseq_lib_path)
from RNAseq_lib import run_cmd

splicing_difference_cutoff = '0.0001'
analysis_type = 'U'
library_type = 'fr-unstranded'
novel_splice = '0'


class prepare(luigi.Task):

    def run(self):
        log_dir = path.join(OutDir, 'logs')
        rmats_dir = path.join(OutDir, 'rmats')
        plot_dir = path.join(OutDir, 'sashimiplot')

        tmp = run_cmd(['mkdir',
        '-p',
        log_dir,
        rmats_dir,
        plot_dir])

        with self.output().open('w') as prepare_log:
            prepare_log.write(tmp)

    def output(self):
        return luigi.LocalTarget('{}/logs/prepare.log'.format(OutDir))

class run_rmats(luigi.Task):
    '''
    transcript splicing analysis using rMATS
    '''

    compare = luigi.Parameter()

    def requires(self):
        return prepare()

    def run(self):
        group1, group2 = self.compare
        group1_bam_list = ['{0}/{1}.bam'.format(BamDir, each_sample) for each_sample in group_sample_df.loc[group1][1]]
        group2_bam_list = ['{0}/{1}.bam'.format(BamDir, each_sample) for each_sample in group_sample_df.loc[group2][1]]
        out_dir = path.join(OutDir, 'rmats', '{0}_vs_{1}'.format(group1, group2))

        tmp = run_cmd(['RNASeq-MATS.py',
        '-b1',
        '{}'.format(','.join(group1_bam_list)),
        '-b2',
        '{}'.format(','.join(group2_bam_list)),
        '-t',
        'paired',
        '-len',
        '150',
        '-gtf',
        '{}'.format(Gtf),
        '-o',
        '{}'.format(out_dir),
        '-c',
        '{}'.format(splicing_difference_cutoff),
        '-analysis',
        '{}'.format(analysis_type),
        '-libType',
        '{}'.format(library_type),
        '-novelSS',
        '{0}'.format(novel_splice)])

        with self.output().open('w') as run_rmats_log:
            run_rmats_log.write(tmp)

    def output(self):
        out_name = '_vs_'.join(self.compare)
        return luigi.LocalTarget('{0}/logs/{1}.run_rmats.log'.format(OutDir, out_name))


class rmats_collection(luigi.Task):

    BamDir = luigi.Parameter()
    SampleInf = luigi.Parameter()
    OutDir = luigi.Parameter()
    Gtf = luigi.Parameter()

    def requires(self):
        global BamDir, OutDir, Gtf, group_sample_df, compare_list
        BamDir = self.BamDir
        OutDir = self.OutDir
        Gtf = self.Gtf
        group_sample_df = pd.read_table(self.SampleInf, header = None, index_col = 0)
        compare_list = itertools.combinations(group_sample_df.index.unique(), 2)

        return [run_rmats(compare = each_compare) for each_compare in compare_list]

    def run(self):
        pass

    def output(self):
        pass

if __name__ == '__main__':
    luigi.run()
