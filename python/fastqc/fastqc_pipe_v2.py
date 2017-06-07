#! /usr/bin/python

import luigi
from os import path
import sys

script_path = path.dirname(path.abspath(__file__))
RNAseq_lib_path = path.join(script_path, '..')
sys.path.insert(0, RNAseq_lib_path)
from RNAseq_lib import FASTQC_SUMMERY
from RNAseq_lib import FASTQC
from RNAseq_lib import GC_PLOT
from RNAseq_lib import RQ_PLOT
from RNAseq_lib import run_cmd
from python_tools import write_obj_to_json


class prepare(luigi.Task):
    '''
    prepare directory and others
    '''
    OutDir = luigi.Parameter()

    def run(self):
        log_dir = path.join(OutDir, 'logs')
        fastqc_outdir = path.join(OutDir, 'fastqc_results')
        tmp = run_cmd(['mkdir',
                       '-p',
                       log_dir,
                       fastqc_outdir])
        with self.output().open('w') as prepare_logs:
            prepare_logs.write(tmp)

    def output(self):
        return luigi.LocalTarget('{}/logs/prepare.log'.format(OutDir))

# run fastqc


class run_fastqc(luigi.Task):
    '''
    run fastqc
    '''
    sample = luigi.Parameter()
    OutDir = luigi.Parameter()

    def requires(self):
        return prepare(OutDir=OutDir)

    def output(self):
        return luigi.LocalTarget('{0}/logs/{1}.log'.format(OutDir, self.sample))

    def run(self):
        tmp = run_cmd([FASTQC,
                       '{0}/{1}_1.clean.fq.gz'.format(CleanDir, self.sample),
                       '{0}/{1}_2.clean.fq.gz'.format(CleanDir, self.sample),
                       '--extract',
                       '-o',
                       '{}/fastqc_results'.format(OutDir)])
        with self.output().open('w') as qc_logs:
            qc_logs.write(tmp)


class fastqc_summary(luigi.Task):
    '''
    generate summary table
    '''
    OutDir = luigi.Parameter()

    def requires(self):
        return [run_fastqc(sample=sample, OutDir=OutDir) for sample in sample_list]

    def run(self):
        tmp = run_cmd(['python',
                       FASTQC_SUMMERY,
                       '{}'.format(SampleInf),
                       '{}'.format(OutDir),
                       '{}/fastqc_general_stats'.format(OutDir)])
        with self.output().open('w') as qc_summary:
            qc_summary.write(tmp)

    def output(self):
        return luigi.LocalTarget('{}/logs/fastqc_summary.log'.format(OutDir))


class gc_plot(luigi.Task):
    '''
    plot gc graph
    '''
    OutDir = luigi.Parameter()

    def requires(self):
        return fastqc_summary(OutDir=OutDir)

    def run(self):
        gc_dir = path.join(OutDir, 'gc_plot')
        tmp = run_cmd(['Rscript',
                       GC_PLOT,
                       '--gc_dir',
                       gc_dir,
                       '--out_dir',
                       gc_dir,
                       '--sample_inf',
                       SampleInf])
        with self.output().open('w') as gc_plot_log:
            gc_plot_log.write(tmp)

    def output(self):
        return luigi.LocalTarget('{0}/logs/gc_plot.log'.format(OutDir))


class reads_quality_plot(luigi.Task):
    '''
    plot reads quality
    '''
    OutDir = luigi.Parameter()

    def requires(self):
        return fastqc_summary(OutDir=OutDir)

    def run(self):
        tmp = run_cmd(['Rscript',
                       RQ_PLOT,
                       '{}/reads_quality_plot/'.format(OutDir),
                       SampleInf])
        with self.output().open('w') as reads_quality_plot_log:
            reads_quality_plot_log.write(tmp)

    def output(self):
        return luigi.LocalTarget('{}/logs/reads_quality_plot.log'.format(OutDir))


class fastqc_collection(luigi.Task):

    OutDir = luigi.Parameter()
    SampleInf = luigi.Parameter()
    CleanDir = luigi.Parameter()

    def requires(self):
        global OutDir, SampleInf, CleanDir, sample_list
        OutDir = self.OutDir
        SampleInf = self.SampleInf
        CleanDir = self.CleanDir
        sample_list = [each.strip().split()[1] for each in open(SampleInf)]
        return [reads_quality_plot(OutDir=OutDir), gc_plot(OutDir=OutDir)]

    def run(self):
        ignore_files = ['.ignore', 'logs', 'fastqc_results/*zip', '.report_files']
        pdf_report_files = ['fastqc_general_stats.txt',
                            'gc_plot/*gc_distribution.line.png', 'reads_quality_plot/*reads_quality.bar.png']
        pdf_report_ini = path.join(self.OutDir, '.report_files')
        write_obj_to_json(pdf_report_files, pdf_report_ini)
        with self.output().open('w') as ignore_files_inf:
            for each_file in ignore_files:
                ignore_files_inf.write('{}\n'.format(each_file))

    def output(self):
        return luigi.LocalTarget('{}/.ignore'.format(self.OutDir))


if __name__ == '__main__':
    luigi.run()
