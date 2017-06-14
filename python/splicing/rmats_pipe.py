import luigi
from os import path
import pandas as pd
import itertools
from RNAseq_lib import run_cmd
from RNAseq_lib import get_diff_splicing_table
from RNAseq_lib import get_diff_as_plot_cmd
from RNAseq_lib import GET_AS_SUMMARY_PLOT_DATA
from RNAseq_lib import SIG_AS_PLOT

splicing_difference_cutoff = '0.0001'
analysis_type = 'U'
library_type = 'fr-unstranded'
novel_splice = '0'
splicing_type = ['SE', 'RI', 'MXE', 'A5SS', 'A3SS']


class prepare(luigi.Task):

    def run(self):
        log_dir = path.join(OutDir, 'logs')
        rmats_dir = path.join(OutDir, 'rmats')

        tmp = run_cmd(['mkdir',
                       '-p',
                       log_dir,
                       rmats_dir])

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
        group1_bam_list = ['{0}/{1}.bam'.format(BamDir, each_sample)
                           for each_sample in group_sample_df.loc[group1][1]]
        group2_bam_list = ['{0}/{1}.bam'.format(BamDir, each_sample)
                           for each_sample in group_sample_df.loc[group2][1]]
        out_dir = path.join(
            OutDir, 'rmats', '{0}_vs_{1}'.format(group1, group2))

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


class extract_diff_splicing(luigi.Task):

    compare = luigi.Parameter()

    def requires(self):

        return run_rmats(compare=self.compare)

    def run(self):

        compare_name = '{0}_vs_{1}'.format(self.compare[0], self.compare[1])
        junction_only_results = ['{0}/rmats/{1}/MATS_output/{2}.MATS.JunctionCountOnly.txt'.format(
            OutDir, compare_name, each_as_type) for each_as_type in splicing_type]
        t_and_j_results = ['{0}/rmats/{1}/MATS_output/{2}.MATS.ReadsOnTargetAndJunctionCounts.txt'.format(
            OutDir, compare_name, each_as_type) for each_as_type in splicing_type]
        out_dir_list = [
            '{0}/{1}'.format(OutDir, compare_name)] * len(splicing_type)
        map(get_diff_splicing_table, junction_only_results, out_dir_list)
        map(get_diff_splicing_table, t_and_j_results, out_dir_list)
        summary_file = path.join(
            OutDir, 'rmats', compare_name, 'summary.txt')
        summary_file_out = path.join(OutDir, compare_name, 'summary.txt')
        summary_file_out_inf = open(summary_file_out, 'w')
        headers = splicing_type[:]
        headers.append('EventType')
        with open(summary_file) as summary_file_inf:
            for eachline in summary_file_inf:
                eachline_inf = eachline.strip().split('\t')
                if eachline_inf[0] in headers:
                    summary_file_out_inf.write(eachline)
        summary_file_out_inf.close()

        with self.output().open('w') as extract_diff_splicing_log:
            extract_diff_splicing_log.write(
                '{} extracted diff splicing finished'.format(self.compare))

    def output(self):
        return luigi.LocalTarget('{0}/logs/{1}_vs_{2}.extract_diff_splicing.log'.format(OutDir, self.compare[0], self.compare[1]))


class rmats_sashimiplot(luigi.Task):

    compare = luigi.Parameter()

    def requires(self):
        return extract_diff_splicing(compare=self.compare)

    def run(self):

        group1, group2 = self.compare
        compare_name = '{0}_vs_{1}'.format(group1, group2)
        diff_t_and_j_results = ['{0}/{1}/diff.{2}.MATS.ReadsOnTargetAndJunctionCounts.txt'.format(
            OutDir, compare_name, each_as_type) for each_as_type in splicing_type]
        diff_junction_only_results = ['{0}/{1}/diff.{2}.MATS.JunctionCountOnly.txt'.format(
            OutDir, compare_name, each_as_type) for each_as_type in splicing_type]
        t_and_j_plot_dir = [
            '{0}/{1}/plot/ReadsOnTargetAndJunctionCounts/{2}'.format(OutDir, compare_name, each_as_type) for each_as_type in splicing_type]
        junction_only_plot_dir = [
            '{0}/{1}/plot/JunctionCountOnly/{2}'.format(OutDir, compare_name, each_as_type) for each_as_type in splicing_type]
        # get group configure file for plot
        group_cf = path.join(OutDir, 'rmats', compare_name, 'group.cf')
        group1_sample_number = len(group_sample_df.loc[group1][1])
        group2_sample_number = len(group_sample_df.loc[group2][1])
        with open(group_cf, 'w') as group_cf_inf:
            group_cf_inf.write(
                '{0}: 1-{1}\n'.format(group1, group1_sample_number))
            group2_start = group1_sample_number + 1
            group2_end = group1_sample_number + group2_sample_number
            group_cf_inf.write(
                '{0}: {1}-{2}\n'.format(group2, group2_start, group2_end))
        # get plot cmd
        group1_bam_list = ['{0}/{1}.bam'.format(BamDir, each_sample)
                           for each_sample in group_sample_df.loc[group1][1]]
        group2_bam_list = ['{0}/{1}.bam'.format(BamDir, each_sample)
                           for each_sample in group_sample_df.loc[group2][1]]
        bam_file_list = [','.join(group1_bam_list), ','.join(group2_bam_list)]
        t_and_j_plot_cmd = []
        for n, each_as_type in enumerate(splicing_type):
            each_t_and_j_cmd = get_diff_as_plot_cmd(
                diff_t_and_j_results[n], self.compare, bam_file_list, each_as_type, t_and_j_plot_dir[n], group_cf)
            each_junction_only_cmd = get_diff_as_plot_cmd(
                diff_junction_only_results[n], self.compare, bam_file_list, each_as_type, junction_only_plot_dir[n], group_cf)
            t_and_j_plot_cmd.extend([each_t_and_j_cmd, each_junction_only_cmd])
        # run cmds
        # stop after plot one record? try to put script into a sh script
        sh_script = path.join(OutDir, 'rmats', compare_name, 'plot.sh')
        with open(sh_script, 'w') as sh_script_inf:
            for each_cmd in t_and_j_plot_cmd:
                sh_script_inf.write('{}\n'.format(' '.join(each_cmd)))
        run_sh_script = ['sh', sh_script]
        plot_cmd_log_inf = run_cmd(run_sh_script)
        with self.output().open('w') as plot_cmd_log:
            plot_cmd_log.write(plot_cmd_log_inf)

    def output(self):
        return luigi.LocalTarget('{0}/logs/{1}_vs_{2}.diff_splicing_plot.log'.format(OutDir, self.compare[0], self.compare[1]))


class rmats_summary_plot(luigi.Task):

    compare = luigi.Parameter()

    def requires(self):
        return extract_diff_splicing(compare=self.compare)

    def run(self):
        group1, group2 = self.compare
        compare_name = '{0}_vs_{1}'.format(group1, group2)
        summary_file = path.join(OutDir, 'rmats', compare_name, 'summary.txt')
        summary_plot_file = path.join(
            OutDir, 'rmats', compare_name, '{}.summary.plot.txt'.format(compare_name))
        plot_out_dir = path.join(OutDir, compare_name)
        get_plot_file_cmd = ['sh',
                             GET_AS_SUMMARY_PLOT_DATA,
                             summary_file,
                             summary_plot_file]

        plot_cmd = ['Rscript',
                    SIG_AS_PLOT,
                    '--as_summary',
                    summary_plot_file,
                    '--out_dir',
                    plot_out_dir]
        plot_cmd_log_inf = run_cmd([get_plot_file_cmd, plot_cmd])
        with self.output().open('w') as plot_cmd_log:
            plot_cmd_log.write(plot_cmd_log_inf)

    def output(self):
        return luigi.LocalTarget('{0}/logs/{1}_vs_{2}.summary_plot.log'.format(OutDir, self.compare[0], self.compare[1]))


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
        group_sample_df = pd.read_table(
            self.SampleInf, header=None, index_col=0)
        compare_list = itertools.combinations(
            group_sample_df.index.unique(), 2)

        return [rmats_summary_plot(compare=each_compare) for each_compare in compare_list]

    def run(self):
        ignore_files = ['.ignore', 'logs', 'rmats']
        with self.output().open('w') as ignore_files_inf:
            for each_file in ignore_files:
                ignore_files_inf.write('{}\n'.format(each_file))

    def output(self):
        return luigi.LocalTarget('{}/logs/rmats_collection.log'.format(self.OutDir))


if __name__ == '__main__':
    luigi.run()
