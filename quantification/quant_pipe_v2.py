import luigi
from os import path
from os import walk
import sys
import pandas as pd
import itertools

script_path = path.dirname(path.abspath(__file__))
RNAseq_lib_path = path.join(script_path, '..')
sys.path.insert(0, RNAseq_lib_path)
from RNAseq_lib import run_cmd
from RNAseq_lib import KALLISTO_TO_TABLE
from RNAseq_lib import DIFF_ANALYSIS
from RNAseq_lib import QUANT_REPORT
from RNAseq_lib import rsync_pattern_to_file
from RNAseq_lib import txt_to_excel
from python_tools import write_obj_to_file


class prepare(luigi.Task):
    '''
    prepare directory and others
    '''
    OutDir = luigi.Parameter()

    def run(self):
        log_dir = path.join(OutDir, 'logs')
        kallisto_dir = path.join(OutDir, 'kallisto')
        diff_dir = path.join(OutDir, 'differential_analysis')
        exp_dir = path.join(OutDir, 'expression_summary')

        tmp = run_cmd(['mkdir',
                       '-p',
                       log_dir,
                       kallisto_dir,
                       diff_dir,
                       exp_dir])
        with self.output().open('w') as prepare_logs:
            prepare_logs.write(tmp)

    def output(self):
        return luigi.LocalTarget('{}/logs/prepare.log'.format(OutDir))


class run_kallisto(luigi.Task):
    '''
    run kallisto
    '''
    OutDir = luigi.Parameter()
    sample = luigi.Parameter()

    def requires(self):
        return prepare(OutDir=OutDir)

    def run(self):
        tmp = run_cmd(['kallisto',
                       'quant',
                       '-i',
                       '{}.kallisto_idx'.format(Transcript),
                       '--bootstrap-samples=100',
                       '--output-dir={0}/kallisto/{1}'.format(
                           OutDir, self.sample),
                       '{0}/{1}_1.clean.fq.gz'.format(CleanDir, self.sample),
                       '{0}/{1}_2.clean.fq.gz'.format(CleanDir, self.sample)])

        with self.output().open('w') as run_kallisto_logs:
            run_kallisto_logs.write(tmp)

    def output(self):
        return luigi.LocalTarget('{0}/logs/{1}_kallisto.log'.format(OutDir, self.sample))


class kallisto_to_matrix(luigi.Task):

    '''
    generate gene expression matrix from kallisto quantification results
    '''

    OutDir = luigi.Parameter()

    def requires(self):
        return [run_kallisto(sample=each_sample, OutDir=OutDir) for each_sample in sample_list]

    def run(self):
        gener_mat_cmd = ['Rscript',
                         KALLISTO_TO_TABLE,
                         '--kallisto_dir',
                         '{}/kallisto'.format(OutDir),
                         '--sample_inf',
                         SampleInf,
                         '--gene2tr',
                         Gene2Tr,
                         '--out_dir',
                         OutDir]

        gener_mat_inf = run_cmd(gener_mat_cmd)
        with self.output().open('w') as gener_mat_logs:
            gener_mat_logs.write(gener_mat_inf)

    def output(self):
        return luigi.LocalTarget('{0}/logs/kallisto_to_matrix.log'.format(OutDir))


class run_diff(luigi.Task):
    '''
    run diff analysis for each compare
    '''

    compare = luigi.Parameter()
    OutDir = luigi.Parameter()

    def requires(self):
        return kallisto_to_matrix(OutDir=OutDir)

    def run(self):
        diff_cmd = ['Rscript',
                    DIFF_ANALYSIS,
                    '--kallisto_dir',
                    '{}/kallisto'.format(OutDir),
                    '--tpm_table',
                    '{}/expression_summary/Gene.tpm.txt'.format(OutDir),
                    '--compare',
                    self.compare,
                    '--sample_inf',
                    SampleInf,
                    '--gene2tr',
                    Gene2Tr,
                    '--out_dir',
                    '{0}/differential_analysis/{1}'.format(OutDir, self.compare)]

        diff_inf = run_cmd(diff_cmd)
        with self.output().open('w') as diff_log:
            diff_log.write(diff_inf)

    def output(self):
        return luigi.LocalTarget('{0}/logs/diff_analysis_{1}.log'.format(OutDir, self.compare))


class get_excel_table(luigi.Task):

    '''
    generate excel format table
    '''
    OutDir = luigi.Parameter()

    def requires(self):
        return [run_diff(compare=each_compare, OutDir=OutDir) for each_compare in compare_name_list]

    def run(self):

        for dirpath, dirnames, filenames in walk(OutDir):
            for each_file in filenames:
                if each_file.endswith('.txt'):
                    each_file_path = path.join(dirpath, each_file)
                    txt_to_excel(each_file_path)

        with self.output().open('w') as get_excel_table_log:
            get_excel_table_log.write('txt to excel finished')

    def output(self):
        return luigi.LocalTarget('{0}/logs/get_excel_table.log'.format(OutDir))


class report_data(luigi.Task):
    '''
    generate table and plots for report
    '''
    OutDir = luigi.Parameter()

    def requires(self):
        return get_excel_table(OutDir=OutDir)

    def run(self):
        report_tb_cmd = ['Rscript',
                         QUANT_REPORT,
                         '--quant_dir',
                         OutDir,
                         '--sample_inf',
                         SampleInf]

        report_tb_inf = run_cmd(report_tb_cmd)
        with self.output().open('w') as report_log:
            report_log.write(report_tb_inf)

    def output(self):
        return luigi.LocalTarget('{0}/logs/report_data.log'.format(OutDir))


class quant_collection(luigi.Task):

    OutDir = luigi.Parameter()
    SampleInf = luigi.Parameter()
    CleanDir = luigi.Parameter()
    Transcript = luigi.Parameter()
    Gene2Tr = luigi.Parameter()

    def requires(self):
        global OutDir, SampleInf, CleanDir, sample_list
        global Transcript, Gene2Tr, compare_name_list
        OutDir = self.OutDir
        SampleInf = self.SampleInf
        CleanDir = self.CleanDir
        Transcript = self.Transcript
        Gene2Tr = self.Gene2Tr
        group_sample_df = pd.read_table(
            self.SampleInf, header=None, index_col=0)
        compare_list = itertools.combinations(
            group_sample_df.index.unique(), 2)
        compare_name_list = ['{0}_vs_{1}'.format(
            each_compare[0], each_compare[1]) for each_compare in compare_list]
        sample_list = [each.strip().split()[1] for each in open(SampleInf)]
        return report_data(OutDir=OutDir)

    def run(self):
        ignore_files = ['.ignore', 'logs', 'kallisto/*/run_info.json',
                        '.report_files', 'Rplots.pdf']
        report_files_pattern = ['expression_summary/*.png',
                        'differential_analysis/*/*png',
                        'expression_summary/*Gene.tpm.txt',
                        'expression_summary/*example.diff.table.txt',
                        'differential_analysis/*/*.edgeR.DE_results.txt']
        report_files = rsync_pattern_to_file(self.OutDir, report_files_pattern)
        report_ini = path.join(self.OutDir, '.report_files')
        write_obj_to_file(report_files, report_ini)
        with self.output().open('w') as ignore_inf:
            for each_file in ignore_files:
                ignore_inf.write('{}\n'.format(each_file))

    def output(self):
        return luigi.LocalTarget('{}/.ignore'.format(self.OutDir))


if __name__ == '__main__':
    luigi.run()
